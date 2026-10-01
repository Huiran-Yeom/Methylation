# -*- coding: utf-8 -*-
"""S0b 리드 무늬 세기 (2026-08-20): 라벨을 쓰지 않는다.

표준 CpG 자리(02_blockcpg)로 창을 만들고 BAM 을 통째로 훑는다.
창을 모두 덮은 리드의 무늬를 센다. 리드를 버리지 않고 전부 센다.

  - 창은 모든 검체에서 동일하다 (JSD 가 성립하려면 그래야 한다)
  - 창당 최소 리드 수를 여기서 걸지 않는다. 문턱은 카운트를 보고 정한다
  - 페어드엔드 두 짝은 한 조각으로 합친다 (같은 분자를 두 번 세지 않는다)
"""

# 이 파일은 스크립트다. 다른 코드가 import 하면 본체가 그대로 돌아
# 실제 자료를 덮어쓴다. 실행은 `python <파일>` 로만 한다.
if __name__ != '__main__':
    raise ImportError(__file__ + ' 은(는) 스크립트다 — import 하지 않는다')

import os as _os, sys as _sys
# 0_setup 을 찾는다. 이 파일이 「복사되어」 정리본 밖에서 돌 수도 있으므로
#   ① 자기 위에서 위로 올라가며 찾고 ② 환경변수 METH_CONF_DIR ③ 그래도 없으면
#   설정 없이 돌 수 있게 기본값으로 간다(검증이 검체 폴더로 복사해 쓰는 경로다).
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


import os, sys, glob
import numpy as np, pandas as pd, pysam

# 산출 폴더는 GEN(패널 세대)을 따른다. 비우면 세대 없는 옛 폴더를 가리킨다.
_GEN = _os.environ.get('GEN', '')
if not _GEN:
    try:
        _GEN = _cfg.GEN
    except Exception:
        raise SystemExit('GEN(패널 세대)을 못 읽었습니다 — config.conf 의 GEN 을 채우거나 '
                         '환경변수로 주십시오. 비워 두면 옛 세대 폴더를 가리킵니다.')
# ╔══ 이 단계가 쓰는 경로 ═══════════════════════════════════════════╗
#   배치가 다르면 「여기만」 고친다. 0_setup/config.conf 에는 뿌리 둘
#   (METH_ROOT=_M · DATASET_ROOT=_D)과 실험 설정만 둔다.
#   환경변수나 명령줄 인자로 그때그때 덮어쓸 수 있다(아래 참조).
#   이 단계가 읽는 자료 폴더는 환경변수로 옮길 수 있다.
#     NORMAL_SET      정상 cov (검체 목록을 여기서 뽑는다)
#     GBM_COV_DIR     암 cov (같음)
#     NORMAL_BAM_DIR  정상 BAM (실제로 읽는 곳)
#     GBM_BAM_DIR     암 BAM (같음)
CAND = _M + '/results/dmr/j_candidates' + _GEN
OUT  = _M + '/results/dmr/j_jsdcount' + _GEN
NB   = _os.environ.get('NORMAL_BAM_DIR') or _D + '/for_in_silico_test/normal_cfDNA_public/aligned_bam'
GB   = _os.environ.get('GBM_BAM_DIR')    or _D + '/training/GBM_cell-line/aligned_bam'
ND   = _os.environ.get('NORMAL_SET')     or _M + '/sample_data/_input/normal_pub15'
GD   = _os.environ.get('GBM_COV_DIR')    or _D + '/training/GBM_cell-line/cov'
# ╚══════════════════════════════════════════════════════════════════╝
BLOCK, KS = 100, (3, 4, 5)
try:   # 2026-09-29: BLOCK_SIZE 가 config.conf 에 있는데 1단계는 100 을 박고 있었다.
    BLOCK = int(_cfg.BLOCK_SIZE)
except Exception:
    pass

a = sys.argv
if '--nd'   in a: ND   = a[a.index('--nd')+1]        # 2026-09-10 학습 정상 폴더를 CLI 로
if '--cand' in a: CAND = a[a.index('--cand')+1]
if '--out'  in a: OUT  = a[a.index('--out')+1]
ONLY  = a[a.index('--sample')+1] if '--sample' in a else None
LIMIT = int(a[a.index('--limit')+1]) if '--limit' in a else 0

def sid(p):
    b = os.path.basename(p)
    for suf in ('.bismark.cov.gz', '.cov.gz', '.bam'):
        if b.endswith(suf): b = b[:-len(suf)]; break
    return b.split('_')[0]

P = pd.read_parquet(CAND + '/blocks_cpg.parquet').sort_values(['chr','blk','pos'])
gp = P.groupby(['chr','blk'], observed=True)['pos'].apply(list)
wins   = {K: [] for K in KS}
blkwin = {}
for (ch, bk), pos in gp.items():
    for K in KS:
        for s in range(len(pos) - K + 1):
            w = tuple(pos[s:s+K])
            blkwin.setdefault((ch,bk), {}).setdefault(K, []).append((len(wins[K]), w))
            wins[K].append((ch, bk, s, w))
print('칸 %d개 · 창 %s' % (len(gp), ' · '.join('K%d %d' % (K, len(wins[K])) for K in KS)))
C = {K: np.zeros((len(wins[K]), 1 << K), dtype=np.int32) for K in KS}

def flush(frag):
    for k, d in frag.items():
        for K, lst in blkwin[k].items():
            arr = C[K]
            for wid, w in lst:
                v = 0
                for p in w:
                    b = d.get(p)
                    if b is None: v = -1; break
                    v = (v << 1) | b
                if v >= 0: arr[wid, v] += 1

def scan(bam):
    for K in KS: C[K][:] = 0
    f = pysam.AlignmentFile(bam, 'rb', check_sq=False)
    last, frag, nr, nh = None, {}, 0, 0
    for r in f:
        if r.is_unmapped or r.is_secondary or r.is_supplementary: continue
        nr += 1
        if LIMIT and nr > LIMIT: break
        if r.query_name != last:
            if frag: flush(frag)
            frag, last = {}, r.query_name
        ch, e = r.reference_name, r.reference_end
        if e is None: continue
        if not any((ch, b) in blkwin for b in range((r.reference_start+1)//BLOCK, e//BLOCK + 1)):
            continue
        nh += 1
        try: xm = r.get_tag('XM')
        except KeyError: continue
        for qp, rp in r.get_aligned_pairs(matches_only=True):
            c = xm[qp]
            if c != 'z' and c != 'Z': continue
            p = rp + 1
            k = (ch, p // BLOCK)
            if k in blkwin: frag.setdefault(k, {}).setdefault(p, 1 if c == 'Z' else 0)
    if frag: flush(frag)
    f.close()
    return nr, nh

os.makedirs(OUT, exist_ok=True)
for K in KS:
    fp = OUT + '/windows_K%d.parquet' % K
    if not os.path.exists(fp):
        d = pd.DataFrame({'chr':[x[0] for x in wins[K]], 'blk':[x[1] for x in wins[K]],
                          'wi':[x[2] for x in wins[K]]})
        for j in range(K): d['pos%d' % (j+1)] = [x[3][j] for x in wins[K]]
        d.to_parquet(fp, index=False); print('창 표 저장 K%d' % K)

jobs = [(sid(f), NB, 'Normal') for f in sorted(glob.glob(ND+'/*.cov.gz'))] + \
       [(sid(f), GB, 'GBM')    for f in sorted(glob.glob(GD+'/*.cov.gz'))]
if ONLY: jobs = [j for j in jobs if j[0] == ONLY]
if not jobs: sys.exit('대상 검체 없음')

for i, (s, bd, ty) in enumerate(jobs, 1):
    hit = [f for f in glob.glob(bd+'/*.bam') if sid(f) == s]
    if len(hit) != 1:
        print('  [%2d/%d] %-16s BAM %d개: 건너뜀' % (i, len(jobs), s, len(hit))); continue
    nr, nh = scan(hit[0])
    msg = []
    for K in KS:
        tot = C[K].sum(1)
        nz  = int((tot > 0).sum())
        med = int(np.median(tot[tot > 0])) if nz else 0
        msg.append('K%d %d창(%.0f%%) 중앙%d' % (K, nz, 100.0*nz/max(1,len(tot)), med))
        cd = pd.DataFrame(C[K], columns=['p%d' % j for j in range(1 << K)])
        cd.insert(0, 'wid', np.arange(len(cd), dtype=np.int32))
        cd = cd[tot > 0]
        cd.to_parquet(OUT + '/count_%s_K%d.parquet' % (s, K), index=False)
    print('  [%2d/%d] %-16s 리드 %10d · 후보걸침 %9d · %s'
          % (i, len(jobs), s, nr, nh, ' | '.join(msg)), flush=True)
print('\n저장: ' + OUT)
