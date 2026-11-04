# -*- coding: utf-8 -*-
"""패널 영역 리드만 모아 작은 풀 BAM 을 만듦
검체마다 BAM 을 한 번 훑어 패널에 걸치는 리드만 뽑고 정렬, 색인

사용: python 11_pool.py --panel <패널 csv> --use <검체 목록> --out <풀.bam> [--pad 200]
"""

# 스크립트 전용. import 하면 본체가 돌아 자료를 덮어씀
if __name__ != '__main__':
    raise ImportError(__file__ + ' 은(는) 스크립트다 — import 하지 않는다')

import os as _os, sys as _sys
# 0_setup 찾기: 위로 올라가며 찾고, 없으면 METH_CONF_DIR 사용
_h = _os.path.dirname(_os.path.abspath(__file__))
for _ in range(5):
    _c = _os.path.join(_h, '0_setup')
    if _os.path.isdir(_c):
        _sys.path.insert(0, _c); break
    _h = _os.path.dirname(_h)
else:
    _e = _os.environ.get('METH_CONF_DIR')
    if _e and _os.path.isdir(_e):
        _sys.path.insert(0, _e)
try:
    import meth_config as _cfg
    _M, _D = _cfg.METH_ROOT, _cfg.DATASET_ROOT
except ImportError:
    raise SystemExit(
        '설정을 못 찾았습니다. 0_setup/meth_config.py 가 있는 자리를 찾지 못했습니다.\n'
        '  이 파일은 정리본 안에서 도는 것이라 기본값으로 넘어가지 않습니다.\n'
        '  트리 밖에서 돌리려면 METH_CONF_DIR 로 0_setup 자리를 주십시오.')


import os, re, sys, glob, subprocess
import pandas as pd, pysam

# 이 단계 경로. 폴더 구조 다르면 여기만 수정
#   NORMAL_BAM_DIR (정상 BAM), GBM_BAM_DIR (암 BAM) 환경변수로 바꿀 수 있음
NB = _os.environ.get('NORMAL_BAM_DIR') or _D + '/for_in_silico_test/normal_cfDNA_public/aligned_bam'
GB = _os.environ.get('GBM_BAM_DIR')    or _D + '/training/GBM_cell-line/aligned_bam'
BLOCK = 100
try:   # 01~03 과 같은 블록 크기 사용
    BLOCK = int(_cfg.BLOCK_SIZE)
except Exception:
    pass

a = sys.argv
PANEL = a[a.index('--panel')+1]
OUT   = a[a.index('--out')+1]
USE   = [x for x in open(a[a.index('--use')+1]).read().split() if x]
PAD   = int(a[a.index('--pad')+1]) if '--pad' in a else 200

# LLR_HELDOUT 이면 혼합군 검체만 풀에 넣음 (참조군은 05_llrref 용)
#   나누는 단위: 정상은 사람, 암은 세포주 (samples.csv 의 group)
if int(getattr(_cfg, 'LLR_HELDOUT', 0)):
    _cand = _os.environ.get('CAND_DIR') or (_M + '/results/dmr/j_candidates' + str(_cfg.GEN))
    _S = pd.read_csv(_cand + '/samples.csv')
    _g = dict(zip(_S['sample'].astype(str), _S['group'].astype(str)))
    _t = dict(zip(_S['sample'].astype(str), _S['type'].astype(str)))
    _miss = [s for s in USE if s not in _t]
    if _miss:
        sys.exit('held-out: samples.csv 에 없는 검체 — 어느 군인지 정할 수 없다: %s' % _miss[:5])
    _kinds = sorted({_t[s] for s in USE})
    if len(_kinds) != 1:
        sys.exit('held-out: 한 목록에 정상과 암이 섞여 있다: %s' % _kinds)
    _k = _kinds[0]
    _unit = (lambda s: s) if _k == 'Normal' else (lambda s: _g[s])
    _ref, _mix = _cfg.split_ref_mix(_unit(s) for s in _t if _t[s] == _k)
    _before = len(USE)
    USE = [s for s in USE if _unit(s) in set(_mix)]
    print('[held-out] %s 혼합군만 씀 — 단위 %d/%d · 검체 %d/%d'
          % (_k, len(_mix), len(_ref) + len(_mix), len(USE), _before))
    if not USE:
        sys.exit('held-out: 혼합군에 남은 검체가 없다')

# 파일 이름 → 검체 이름. 기본은 첫 _ 앞까지
#   규칙이 다르면 SID_RE 로 지정 (첫 괄호가 이름, 안 맞으면 기본 규칙)
#   예: SID_RE='(.+?)__'  H1876_0d__SRR... → H1876_0d
_SRE = _os.environ.get('SID_RE')
_srx = re.compile(_SRE) if _SRE else None
def sid(p):
    b = os.path.basename(p)
    for suf in ('.bismark.cov.gz', '.cov.gz', '.bam'):
        if b.endswith(suf): b = b[:-len(suf)]; break
    if _srx is not None:
        m = _srx.match(b)
        if m and m.groups(): return m.group(1)
    return b.split('_')[0]

P = pd.read_csv(PANEL)
keep = set(zip(P['chr'].astype(str), P['blk'].astype('int64')))
print('패널 %d칸 · 검체 %d개' % (len(keep), len(USE)))

bams = []
for s in USE:
    # 정상, 암 BAM 폴더가 같을 수 있어서 실경로로 중복 제거
    _seen, hit = set(), []
    for f in glob.glob(NB+'/*.bam') + glob.glob(GB+'/*.bam'):
        r = os.path.realpath(f)
        if r in _seen or sid(f) != s: continue
        _seen.add(r); hit.append(f)
    if len(hit) != 1:
        sys.exit('BAM 을 못 찾음: %s (%d개)' % (s, len(hit)))
    bams.append((s, hit[0]))

_cfg.check_same_refs([b for _, b in bams], '풀 입력')
os.makedirs(os.path.dirname(OUT) or '.', exist_ok=True)
tmp = OUT + '.unsorted.bam'
hdr = pysam.AlignmentFile(bams[0][1], 'rb', check_sq=False).header
out = pysam.AlignmentFile(tmp, 'wb', header=hdr)
tot = wr = 0
for i, (s, b) in enumerate(bams, 1):
    f = pysam.AlignmentFile(b, 'rb', check_sq=False)
    n = m = 0
    for r in f:
        if r.is_unmapped or r.is_secondary or r.is_supplementary: continue
        n += 1
        e = r.reference_end
        if e is None: continue
        ch = r.reference_name
        b0 = (r.reference_start + 1 - PAD) // BLOCK
        b1 = (e + PAD) // BLOCK
        if any((ch, k) in keep for k in range(b0, b1 + 1)):
            out.write(r); m += 1
    f.close(); tot += n; wr += m
    print('  [%2d/%d] %-16s 리드 %10d → %8d (%.3f%%)'
          % (i, len(bams), s, n, m, 100.0*m/max(1, n)), flush=True)
out.close()
print('\n전체 %d → 뽑은 리드 %d (%.3f%%)' % (tot, wr, 100.0*wr/max(1, tot)))

print('정렬·색인 중...')
subprocess.check_call(['samtools', 'sort', '-@', '4', '-o', OUT, tmp])
subprocess.check_call(['samtools', 'index', OUT])
os.remove(tmp)
print('저장: %s  (%.1f MB)' % (OUT, os.path.getsize(OUT) / 1e6))
