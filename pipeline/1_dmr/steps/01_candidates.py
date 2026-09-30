# -*- coding: utf-8 -*-
"""S1 후보 블록 생성: 라벨을 쓰지 않는다 (2026-08-12)"""

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
ND  = _M + '/sample_data/_input/normal_pub15'
GD  = _D + '/training/GBM_cell-line/cov'
# 산출 폴더는 GEN(패널 세대)을 따른다. 비우면 세대 없는 옛 폴더를 가리킨다.
_GEN = _os.environ.get('GEN', '')
if not _GEN:
    try:
        _GEN = _cfg.GEN
    except Exception:
        raise SystemExit('GEN(패널 세대)을 못 읽었습니다 — config.conf 의 GEN 을 채우거나 '
                         '환경변수로 주십시오. 비워 두면 옛 세대 폴더를 가리킵니다.')
OUT = _M + '/results/dmr/j_candidates' + _GEN
BLOCK, MINCOV, MINCPG, FRAC = 100, 9, 3, 0.60   # MINCOV=9 → 리드 10 이상
try:   # 02·03·11 과 같은 격자를 써야 한다. 여기만 100 을 박으면 04 가
    BLOCK = int(_cfg.BLOCK_SIZE)   #   「쓸 수 있는 창이 0」 으로 죽는다.
except Exception:
    pass
MINREP = 2   # 2026-08-20: 한 독립단위(세포주)는 파일 이만큼 관측돼야 관측으로 친다

a = sys.argv
if '--frac'  in a: FRAC   = float(a[a.index('--frac')+1])
if '--depth' in a: MINCOV = int(a[a.index('--depth')+1]) - 1
if '--nd'    in a: ND     = a[a.index('--nd')+1]        # 2026-09-10 학습 정상 폴더를 CLI 로
if '--cpg'   in a: MINCPG = int(a[a.index('--cpg')+1])
if '--out' in a: OUT = a[a.index('--out')+1]
DRY = '--dry' in a

def sid(p):
    b = os.path.basename(p)
    for suf in ('.bismark.cov.gz', '.cov.gz'):
        if b.endswith(suf):
            b = b[:-len(suf)]; break
    return b.split('_')[0]
def grp(s):  return s.split('-')[0] if s.upper().startswith('SNU') and '-' in s else s

def read_one(p):
    d = pd.read_csv(p, sep='\t', header=None, usecols=[0,1,3,4,5],
                    names=['chr','start','beta','tn','cn'], compression='gzip',
                    dtype={'chr':str,'start':np.int64,'beta':np.float32,
                           'tn':np.int32,'cn':np.int32})
    d = d[(d.tn + d.cn) > MINCOV]
    if d.empty: return None
    d['blk'] = d.start // BLOCK
    d['reads'] = d.tn + d.cn
    g = d.groupby(['chr','blk'], observed=True).agg(
        n_cpg=('beta','size'), reads=('reads','sum'), beta=('beta','mean')).reset_index()
    return g[g.n_cpg >= MINCPG]

files = [(f,'Normal') for f in sorted(glob.glob(ND+'/*.cov.gz'))] + \
        [(f,'GBM')    for f in sorted(glob.glob(GD+'/*.cov.gz'))]
if not files: sys.exit('cov 파일 없음')
if DRY: files = files[:1] + files[-1:]   # 2026-08-20: 정상 1 + 암 1 (암 경로도 지나가게)
# 커버 기준을 집단별로 나눈다.
#   한 덩어리로 세면 한쪽으로 몰린 블록이 통과해 t검정에서 NaN 이 된다.
#   검체 수가 바뀔 때 문턱이 요동치는 문제도 같이 없앤다.
# 문턴을 파일 수가 아니라 독립단위 수로 잡는다.
#   SNU 30파일 = 세포주 10종 x 3복제. 파일로 세면 암 표본이 3배로 부푼다.
n_norm = sum(1 for _, t in files if t == 'Normal')
n_gbm  = len(files) - n_norm
u_norm = len({grp(sid(f)) for f, t in files if t == 'Normal'})
u_gbm  = len({grp(sid(f)) for f, t in files if t == 'GBM'})
assert u_norm == n_norm, '정상에 복제 파일이 있다 — 커버 계산을 다시 봐야 한다'
need_n = int(math.ceil(FRAC * u_norm))
need_g = int(math.ceil(FRAC * u_gbm))
need   = need_n + need_g
print('파일 %d개 (정상 %d · 암 %d) · 독립단위 %d (정상 %d · 암 %d)'
      % (len(files), n_norm, n_gbm, u_norm + u_gbm, u_norm, u_gbm))
print('리드 %d이상 · CpG %d이상' % (MINCOV+1, MINCPG))
print('커버 %.0f%%: 정상 %d명 이상 AND 암 %d종 이상 (한 종은 파일 %d개 이상)'
      % (FRAC*100, need_n, need_g, MINREP))

parts, meta = [], []
for i,(p,t) in enumerate(files,1):
    s = sid(p); g = read_one(p)
    if g is None or g.empty: print('  [%d] %s 통과 0' % (i,s)); continue
    g['sample'] = s
    parts.append(g); meta.append({'sample':s,'group':grp(s),'type':t,'blocks':len(g),'file':p})
    print('  [%2d/%d] %-16s %8d' % (i,len(files),s,len(g)))

A = pd.concat(parts, ignore_index=True)
# 커버를 파일이 아니라 독립단위로 센다.
#   한 세포주는 파일 MINREP개 이상 관측돼야 '관측된 종'으로 친다.
_gmap = {m['sample']: m['group'] for m in meta}
_nset = {m['sample'] for m in meta if m['type'] == 'Normal'}
A['_n'] = A['sample'].isin(_nset)

# 세포주별 파일 수. 파일이 MINREP보다 적으면 있는 만큼만 요구한다.
_nrep = {}
for _f, _t in files:
    if _t == 'GBM':
        _g = grp(sid(_f)); _nrep[_g] = _nrep.get(_g, 0) + 1
_req   = {g: min(MINREP, n) for g, n in _nrep.items()}
_short = sorted(g for g, n in _nrep.items() if n < MINREP)
if _short:
    print('  [주의] 파일 %d개 미만인 세포주: %s' % (MINREP, ','.join(_short)))

n_obs = A[A._n].groupby(['chr','blk'], observed=True).size()
_gg   = A[~A._n][['chr','blk','sample']].copy()
_gg['_grp'] = _gg['sample'].map(_gmap)
_rep  = _gg.groupby(['chr','blk','_grp'], observed=True).size()
_reqv = np.array([_req[g] for g in _rep.index.get_level_values('_grp')])
g_obs = pd.Series((_rep.values >= _reqv).astype(np.int32),
                  index=_rep.index).groupby(level=[0,1]).sum()
del _gg, _rep, _reqv

idx   = n_obs.index.union(g_obs.index)
n_obs = n_obs.reindex(idx, fill_value=0)
g_obs = g_obs.reindex(idx, fill_value=0)
keep  = idx[(n_obs.values >= need_n) & (g_obs.values >= need_g)]

cnt = A.groupby(['chr','blk'], observed=True).size()
A = A.drop(columns=['_n'])
_oldneed = int(math.ceil(FRAC * len(files)))
print('\n후보 %d개 / 관측 %d개' % (len(keep), len(cnt)))
print('  (파일 단위 옛 기준: 전체 %d파일 이상 → %d개)'
      % (_oldneed, int((cnt >= _oldneed).sum())))
if len(keep) == 0: sys.exit('조건 만족 블록 없음')

A = A.set_index(['chr','blk'])
A = A.loc[A.index.isin(keep)].reset_index()
B = A.pivot_table(index=['chr','blk'], columns='sample', values='beta',  aggfunc='first')
D = A.pivot_table(index=['chr','blk'], columns='sample', values='reads', aggfunc='first')
# ╚══════════════════════════════════════════════════════════════════╝
C = A.pivot_table(index=['chr','blk'], columns='sample', values='n_cpg', aggfunc='first')

os.makedirs(OUT, exist_ok=True)
B.to_parquet(OUT+'/beta_matrix.parquet');  D.to_parquet(OUT+'/depth_matrix.parquet')
C.to_parquet(OUT+'/cpg_matrix.parquet')
M = pd.DataFrame(meta); M.to_csv(OUT+'/samples.csv', index=False)
pd.DataFrame({'chr':B.index.get_level_values(0),'blk':B.index.get_level_values(1)}
             ).to_csv(OUT+'/blocks_frac%03d.csv' % round(FRAC*100), index=False)
txt = ['검체 %d (정상 %d · 암 %d)' % (len(meta), (M.type=='Normal').sum(), (M.type=='GBM').sum()),
       '독립단위 정상 %d · 암 %d' % (M[M.type=='Normal'].group.nunique(), M[M.type=='GBM'].group.nunique()),
       '기준 리드%d · CpG%d · 커버%.0f%% (정상 %d명 AND 암 %d종 · 종당 파일 %d개 이상)'
       % (MINCOV+1, MINCPG, FRAC*100, need_n, need_g, MINREP),
       '후보 블록 %d' % len(B), '결측률 %.1f%%' % (B.isna().mean().mean()*100)]
open(OUT+'/summary.txt','w',encoding='utf-8').write('\n'.join(txt)+'\n')
print('\n'+'\n'.join(txt)+'\n저장: '+OUT)
