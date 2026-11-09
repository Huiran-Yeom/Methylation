# -*- coding: utf-8 -*-
"""후보 블록 만들기 (정상, 암 cov 에서 커버 기준 통과 블록, 라벨 안 씀)

사용: python 01_candidates.py [--nd <정상 cov 폴더>] [--out ..] [--frac F] [--depth N] [--cpg N] [--dry]
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


import os, re, sys, glob, math
import numpy as np, pandas as pd

# 이 단계 경로. 폴더 구조 다르면 여기만 수정 (환경변수, 인자로 덮어쓰기 가능)
#   NORMAL_SET      정상 cov (*.cov.gz)
#   CANCER_COV_DIR  암 cov (옛 이름 GBM_COV_DIR)
ND  = _os.environ.get('NORMAL_SET')  or _M + '/sample_data/_input/normal_pub15'
GD  = _cfg.cancer_env('CANCER_COV_DIR', 'GBM_COV_DIR') or _D + '/training/GBM_cell-line/cov'
# 산출 폴더 이름에 GEN(패널 세대) 붙임. 못 읽으면 중단
_sys.path.insert(0, _os.path.dirname(_os.path.abspath(__file__)))
import dmr_common   # noqa: E402  (GEN, 암 검사, 단위 평균)
_GEN = dmr_common.gen(_cfg)
OUT = _M + '/results/dmr/j_candidates' + _GEN
BLOCK, MINCOV, MINCPG, FRAC = 100, 9, 3, 0.60   # MINCOV=9 → 리드 10 이상
try:   # 02, 03, 11 과 같은 블록 격자 사용 (BLOCK_SIZE)
    BLOCK = int(_cfg.BLOCK_SIZE)
except Exception:
    pass
MINREP = 2   # 세포주 하나가 관측으로 치려면 필요한 파일 수

a = sys.argv
if '--frac'  in a: FRAC   = float(a[a.index('--frac')+1])
if '--depth' in a: MINCOV = int(a[a.index('--depth')+1]) - 1
if '--nd'    in a: ND     = a[a.index('--nd')+1]        # 학습 정상 폴더를 CLI 로
if '--cpg'   in a: MINCPG = int(a[a.index('--cpg')+1])
if '--out' in a: OUT = a[a.index('--out')+1]
DRY = '--dry' in a

# 검체 이름, 단위 규칙은 sample_id.py (SID_RE, GROUP_RE)
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import sample_id as _sid
sid = _sid.sid
_CANCER = getattr(_cfg, 'CANCER_LABEL', 'GBM')   # samples.csv 의 type 값


def grp(s, t=_CANCER):
    # GROUP_RE 는 암 쪽에만 사용. 정상은 한 사람당 한 파일 전제
    return _sid.grp(s, t == _CANCER)

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
        [(f,_CANCER)    for f in sorted(glob.glob(GD+'/*.cov.gz'))]
if not files: sys.exit('cov 파일 없음')
# 같은 sid 가 둘 이상이면 중단 (표 열 이름이라 겹치면 pivot_table 이 하나를 버림)
_byid = {}
for _f, _t in files: _byid.setdefault(sid(_f), []).append(_os.path.basename(_f))
_dup = {k: v for k, v in _byid.items() if len(v) > 1}
if _dup:
    _msg = ['중단: 검체 이름이 겹칩니다 (첫 _ 앞까지가 이름입니다).']
    for _k, _v in sorted(_dup.items())[:10]:
        _msg.append('  %s : %s' % (_k, ' · '.join(_v)))
    if len(_dup) > 10: _msg.append('  … 외 %d개' % (len(_dup) - 10))
    _msg.append('  겹치면 표의 열이 덮여 한 쪽이 조용히 사라집니다.')
    _msg.append('  파일 이름의 첫 _ 앞을 서로 다르게 바꾸십시오.')
    sys.exit(chr(10).join(_msg))
# 정상, 암 한 쪽만 비어도 중단 (need 가 0 이 돼서 한 쪽만으로 후보가 나옴)
_nn = sum(1 for _, _t in files if _t == 'Normal')
if _nn == 0:        sys.exit('정상 cov 가 없다 — NORMAL_SET 또는 --nd 를 확인하라: ' + ND)
if _nn == len(files): sys.exit('암 cov 가 없다 — CANCER_COV_DIR (옛 이름 GBM_COV_DIR) 을 확인하라: ' + GD)
if DRY: files = files[:1] + files[-1:]   # 정상 1 + 암 1 (암 경로도 지나가게)
# 커버 기준은 정상, 암 따로 (합쳐 세면 한쪽에 몰린 블록이 통과해 t검정 NaN)
# 문턱은 파일 수 말고 세포주 수 기준 (복제 파일로 세면 암 표본이 부풂)
n_norm = sum(1 for _, t in files if t == 'Normal')
n_gbm  = len(files) - n_norm
u_norm = len({grp(sid(f), 'Normal') for f, t in files if t == 'Normal'})
u_gbm  = len({grp(sid(f)) for f, t in files if t == _CANCER})
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
    parts.append(g); meta.append({'sample':s,'group':grp(s, t),'type':t,'blocks':len(g),'file':p})
    print('  [%2d/%d] %-16s %8d' % (i,len(files),s,len(g)))

A = pd.concat(parts, ignore_index=True)
# 커버는 세포주 단위로 셈 (파일 MINREP개 이상 관측돼야 관측으로 침)
_gmap = {m['sample']: m['group'] for m in meta}
_nset = {m['sample'] for m in meta if m['type'] == 'Normal'}
A['_n'] = A['sample'].isin(_nset)

# 세포주별 파일 수. MINREP 보다 적으면 있는 만큼만 요구
_nrep = {}
for _f, _t in files:
    if _t == _CANCER:
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
C = A.pivot_table(index=['chr','blk'], columns='sample', values='n_cpg', aggfunc='first')

os.makedirs(OUT, exist_ok=True)
B.to_parquet(OUT+'/beta_matrix.parquet');  D.to_parquet(OUT+'/depth_matrix.parquet')
C.to_parquet(OUT+'/cpg_matrix.parquet')
M = pd.DataFrame(meta); M.to_csv(OUT+'/samples.csv', index=False)
pd.DataFrame({'chr':B.index.get_level_values(0),'blk':B.index.get_level_values(1)}
             ).to_csv(OUT+'/blocks_frac%03d.csv' % round(FRAC*100), index=False)
txt = ['검체 %d (정상 %d · 암 %d)' % (len(meta), (M.type=='Normal').sum(), (M.type==_CANCER).sum()),
       '독립단위 정상 %d · 암 %d' % (M[M.type=='Normal'].group.nunique(), M[M.type==_CANCER].group.nunique()),
       '기준 리드%d · CpG%d · 커버%.0f%% (정상 %d명 AND 암 %d종 · 종당 파일 %d개 이상)'
       % (MINCOV+1, MINCPG, FRAC*100, need_n, need_g, MINREP),
       '후보 블록 %d' % len(B), '결측률 %.1f%%' % (B.isna().mean().mean()*100)]
open(OUT+'/summary.txt','w',encoding='utf-8').write('\n'.join(txt)+'\n')
print('\n'+'\n'.join(txt)+'\n저장: '+OUT)
