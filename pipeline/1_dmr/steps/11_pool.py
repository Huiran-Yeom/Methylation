# -*- coding: utf-8 -*-
"""S3-0 패널 영역 리드만 모아 가벼운 합본 BAM 을 만든다 (2026-08-20).

39G 를 통째로 합쳐 정렬할 필요가 없다. 패널 300칸은 유전체의 30kb 뿐이다.
검체마다 한 번 훑으며 패널에 걸치는 리드만 뽑아 모은 뒤 정렬·색인한다.

  정상을 A·B 풀로 나누면 음성 대조가 별도 작업이 아니라 같은 파이프라인의 한 조건이 된다.
"""

# 이 파일은 스크립트다. 다른 코드가 import 하면 본체가 그대로 돌아
# 실제 자료를 덮어쓴다. 실행은 `python <파일>` 로만 한다.
if __name__ != '__main__':
    raise ImportError(__file__ + ' 은(는) 스크립트다 — import 하지 않는다')

import os as _os, sys as _sys
# 0_setup 을 찾는다. 이 파일이 「복사되어」 정리본 밖에서 돌 수도 있으므로
#   ① 자기 위에서 위로 올라가며 찾고 ② 환경변수 METH_CONF_DIR ③ 그래도 없으면
#   meth_config 없이 돌 수 있게 기본값으로 간다(검증이 검체 폴더로 복사해 쓰는 경로다).
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

# ╔══ 이 단계가 쓰는 경로 ═══════════════════════════════════════════╗
#   배치가 다르면 「여기만」 고친다. 0_setup/config.conf 에는 뿌리 둘
#   (METH_ROOT=_M · DATASET_ROOT=_D)과 실험 설정만 둔다.
#   환경변수나 명령줄 인자로 그때그때 덮어쓸 수 있다(아래 참조).
#   이 단계가 읽는 자료 폴더는 환경변수로 옮길 수 있다.
#     NORMAL_BAM_DIR  정상 BAM
#     GBM_BAM_DIR     암 BAM
NB = _os.environ.get('NORMAL_BAM_DIR') or _D + '/for_in_silico_test/normal_cfDNA_public/aligned_bam'
GB = _os.environ.get('GBM_BAM_DIR')    or _D + '/training/GBM_cell-line/aligned_bam'
BLOCK = 100
try:   # 01~03 과 같은 격자. 어긋나면 풀이 패널과 다른 자리를 담는다.
    BLOCK = int(_cfg.BLOCK_SIZE)
except Exception:
    pass

a = sys.argv
PANEL = a[a.index('--panel')+1]
OUT   = a[a.index('--out')+1]
# ╚══════════════════════════════════════════════════════════════════╝
USE   = [x for x in open(a[a.index('--use')+1]).read().split() if x]
PAD   = int(a[a.index('--pad')+1]) if '--pad' in a else 200

# 파일 이름에서 검체 이름을 뽑는다. 기본은 첫 _ 앞까지 — cov 와 BAM 의 뒤쪽
#   꼬리표가 달라도 짝이 맞게 하기 위해서다.
#   자료마다 이름 규칙이 다르다. 안 맞으면 «자료를 고치지 말고» SID_RE 를 준다.
#   첫 괄호가 검체 이름이 되고, 안 걸리면 기본 규칙으로 떨어진다. 예:
#     SID_RE='(.+?)__'   H1876_0d__SRR…  -> H1876_0d   (밑줄 둘로 가르는 자료)
#                        H1876__SRR…     -> H1876
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
    # 두 폴더가 같을 수 있다(정상·암 BAM 을 한 곳에 두는 배치). 그대로 두면 같은
    #   파일이 두 번 걸려 len(hit)!=1 로 건너뛴다. 문자열이 아니라 «실경로» 로
    #   묶는다 — /data/bam 과 /data/bam/ 은 다른 문자열이지만 같은 폴더다.
    _seen, hit = set(), []
    for f in glob.glob(NB+'/*.bam') + glob.glob(GB+'/*.bam'):
        r = os.path.realpath(f)
        if r in _seen or sid(f) != s: continue
        _seen.add(r); hit.append(f)
    if len(hit) != 1:
        sys.exit('BAM 을 못 찾음: %s (%d개)' % (s, len(hit)))
    bams.append((s, hit[0]))

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
