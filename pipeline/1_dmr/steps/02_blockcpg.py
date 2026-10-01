# -*- coding: utf-8 -*-
"""S0a 후보 칸의 표준 CpG 자리 (2026-08-20) — 라벨을 쓰지 않는다.

JSD 창은 모든 검체에서 같아야 한다. 검체마다 관측된 자리로 창을 만들면
서로 다른 자리의 분포를 맞대게 되어 JSD 가 성립하지 않는다.
그래서 칸마다 '충분한 검체가 함께 보는 자리'를 먼저 확정한다.
기준은 후보 칸과 같다. 정상 60% AND 암 60%(종당 파일 MINREP개 이상).
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


import os, sys, glob, math
import numpy as np, pandas as pd

# ╔══ 이 단계가 쓰는 경로 ═══════════════════════════════════════════╗
#   배치가 다르면 「여기만」 고친다. 0_setup/config.conf 에는 뿌리 둘
#   (METH_ROOT=_M · DATASET_ROOT=_D)과 실험 설정만 둔다.
#   환경변수나 명령줄 인자로 그때그때 덮어쓸 수 있다(아래 참조).
#   이 단계가 읽는 자료 폴더는 환경변수로 옮길 수 있다.
#     NORMAL_SET      정상 cov (*.cov.gz)
#     GBM_COV_DIR     암 cov
ND   = _os.environ.get('NORMAL_SET')  or _M + '/sample_data/_input/normal_pub15'
GD   = _os.environ.get('GBM_COV_DIR') or _D + '/training/GBM_cell-line/cov'
# 산출 폴더는 GEN(패널 세대)을 따른다. 비우면 세대 없는 옛 폴더를 가리킨다.
_GEN = _os.environ.get('GEN', '')
if not _GEN:
    try:
        _GEN = _cfg.GEN
    except Exception:
        raise SystemExit('GEN(패널 세대)을 못 읽었습니다 — config.conf 의 GEN 을 채우거나 '
                         '환경변수로 주십시오. 비워 두면 옛 세대 폴더를 가리킵니다.')
CAND = _M + '/results/dmr/j_candidates' + _GEN
BLOCK, MINCOV, FRAC, MINREP = 100, 9, 0.60, 2
try:   # 2026-09-29: BLOCK_SIZE 가 config.conf 에 있는데 1단계는 100 을 박고 있었다.
    BLOCK = int(_cfg.BLOCK_SIZE)
except Exception:
    pass

a = sys.argv
if '--nd'   in a: ND   = a[a.index('--nd')+1]        # 2026-09-10 학습 정상 폴더를 CLI 로
if '--cand' in a: CAND = a[a.index('--cand')+1]
if '--depth' in a: MINCOV = int(a[a.index('--depth')+1]) - 1
OUT = a[a.index('--out')+1] if '--out' in a else CAND
# ╚══════════════════════════════════════════════════════════════════╝

def sid(p):
    b = os.path.basename(p)
    for suf in ('.bismark.cov.gz', '.cov.gz'):
        if b.endswith(suf):
            b = b[:-len(suf)]; break
    return b.split('_')[0]
def grp(s):  return s.split('-')[0] if s.upper().startswith('SNU') and '-' in s else s

files = [(f,'Normal') for f in sorted(glob.glob(ND+'/*.cov.gz'))] + \
        [(f,'GBM')    for f in sorted(glob.glob(GD+'/*.cov.gz'))]
if not files: sys.exit('cov 파일 없음')
_nn = sum(1 for f, _t in files if _t == 'Normal')
if _nn == 0:          sys.exit('정상 cov 가 없다 — NORMAL_SET 또는 --nd 를 확인하라: ' + ND)
if _nn == len(files): sys.exit('암 cov 가 없다 — GBM_COV_DIR 을 확인하라: ' + GD)
if not files: sys.exit('cov 파일 없음')

u_norm = len({grp(sid(f)) for f,t in files if t=='Normal'})
u_gbm  = len({grp(sid(f)) for f,t in files if t=='GBM'})
need_n = int(math.ceil(FRAC*u_norm))
need_g = int(math.ceil(FRAC*u_gbm))
_nrep = {}
for _f,_t in files:
    if _t == 'GBM':
        _g = grp(sid(_f)); _nrep[_g] = _nrep.get(_g,0)+1
_req = {g: min(MINREP,n) for g,n in _nrep.items()}
print('정상 %d명 · 암 %d종 · 자리 기준 정상 %d AND 암 %d'
      % (u_norm, u_gbm, need_n, need_g))

B  = pd.read_csv(CAND + '/blocks_frac%03d.csv' % round(FRAC*100), dtype={'chr':str})
BM = pd.MultiIndex.from_arrays([B['chr'].values, B['blk'].values])
print('후보 칸 %d개' % len(BM))

rows = []
for i,(p,t) in enumerate(files,1):
    s = sid(p)
    d = pd.read_csv(p, sep='\t', header=None, usecols=[0,1,4,5],
                    names=['chr','start','tn','cn'], compression='gzip',
                    dtype={'chr':str,'start':np.int64,'tn':np.int32,'cn':np.int32})
    d = d[(d.tn + d.cn) > MINCOV]
    d['blk'] = d.start // BLOCK
    d = d[pd.MultiIndex.from_arrays([d['chr'].values, d['blk'].values]).isin(BM)]
    print('  [%2d/%d] %-16s %8d' % (i,len(files),s,len(d)))
    if d.empty: continue
    d = d[['chr','blk','start']].copy()
    d['grp'] = grp(s); d['ty'] = t
    rows.append(d)

A = pd.concat(rows, ignore_index=True)
del rows

N  = A[A.ty=='Normal'].groupby(['chr','blk','start'], observed=True)['grp'].nunique()
G  = A[A.ty=='GBM'   ].groupby(['chr','blk','start','grp'], observed=True).size()
Gv = np.array([_req[g] for g in G.index.get_level_values('grp')])
Go = pd.Series((G.values >= Gv).astype(np.int32),
               index=G.index).groupby(level=[0,1,2]).sum()
del A, G, Gv

ix = N.index.union(Go.index)
N  = N.reindex(ix, fill_value=0)
Go = Go.reindex(ix, fill_value=0)
keep = ix[(N.values >= need_n) & (Go.values >= need_g)]
print('\n표준 자리 %d개 / 관측 자리 %d개' % (len(keep), len(ix)))

P = pd.DataFrame({'chr': keep.get_level_values(0),
                  'blk': keep.get_level_values(1),
                  'pos': keep.get_level_values(2)}).sort_values(['chr','blk','pos'])
c = P.groupby(['chr','blk'], observed=True).size()
print('자리가 남은 칸 %d개 / 후보 %d개' % (len(c), len(BM)))
print('칸당 자리 수  최소 %d · 중앙 %d · 평균 %.1f · 최대 %d'
      % (c.min(), int(c.median()), c.mean(), c.max()))
print()
for K in (3,4,5):
    ok = c[c >= K]
    print('  K=%d   칸 %7d (%5.1f%%)   창 %9d'
          % (K, len(ok), 100.0*len(ok)/len(BM), int((ok-K+1).sum())))

os.makedirs(OUT, exist_ok=True)
P.to_parquet(OUT+'/blocks_cpg.parquet', index=False)
txt = ['표준 자리 %d / 칸 %d' % (len(P), len(c)),
       '기준 정상 %d명 AND 암 %d종 (종당 파일 %d개 이상)' % (need_n, need_g, MINREP),
       '칸당 자리 최소%d 중앙%d 평균%.1f 최대%d' % (c.min(), int(c.median()), c.mean(), c.max())]
for K in (3,4,5):
    ok = c[c >= K]
    txt.append('K=%d 칸 %d 창 %d' % (K, len(ok), int((ok-K+1).sum())))
open(OUT+'/blocks_cpg_summary.txt','w',encoding='utf-8').write('\n'.join(txt)+'\n')
print('\n저장: '+OUT+'/blocks_cpg.parquet')
