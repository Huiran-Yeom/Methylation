# -*- coding: utf-8 -*-
"""JSD 상위 N 블록을 패널로 뽑음
q, |dbeta| 는 보고용으로만 계산 (선정에 안 씀)

사용: python 06_select.py --jsd <jsd_K?_minr?.parquet> [--cand ..] [--out ..] [--n N] [--col ..]
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


import os, sys
import numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from moderated_t import moderated_ttest, bh

# 산출 폴더 이름에 GEN(패널 세대) 붙임. 못 읽으면 중단
_GEN = _os.environ.get('GEN', '')
if not _GEN:
    try:
        _GEN = _cfg.GEN
    except Exception:
        raise SystemExit('GEN(패널 세대)을 못 읽었습니다 — config.conf 의 GEN 을 채우거나 '
                         '환경변수로 주십시오. 비워 두면 옛 세대 폴더를 가리킵니다.')
# 이 단계 경로. 폴더 구조 다르면 여기만 수정 (인자로 덮어쓰기 가능)
# --jsd 는 필수 - 기본값 두면 다른 JSD 표로 패널을 고를 수 있음
CAND = _M + '/results/dmr/j_candidates' + _GEN
JSD  = None
OUT  = _M + '/results/dmr/j_panel' + _GEN
N, COL = 300, 'jsd_mean'
a = sys.argv
if '--cand' in a: CAND = a[a.index('--cand')+1]
if '--jsd'  in a: JSD  = a[a.index('--jsd')+1]
if '--out'  in a: OUT  = a[a.index('--out')+1]
if '--n'    in a: N    = int(a[a.index('--n')+1])
if '--col'  in a: COL  = a[a.index('--col')+1]
# 빈 문자열(--jsd "")도 막음. 파일 없으면 run_dmr.sh 가 빈 값을 넘김
if not JSD:
    sys.exit('사용: python 06_select.py --jsd <jsd_K?_minr?.parquet> [--cand ..] [--out ..] [--n N]'
             '   ·  --jsd 에 기본값을 두지 않는다. minr 이 다른 표를 읽으면 다른 패널이 나온다.'
             '   ·  09-10 에 쓴 값: --jsd %s/results/dmr/j_jsd%s/jsd_K3_minr6.parquet --n 200'
             % (_M, _GEN))

B = pd.read_parquet(CAND + '/beta_matrix.parquet')
M = pd.read_csv(CAND + '/samples.csv')
X = B.T.reindex(M['sample']).values.astype(np.float32)
_CANCER = getattr(_cfg, 'CANCER_LABEL', 'GBM')   # samples.csv 의 type 값

y = (M['type'] == _CANCER).astype(int).values
# CANCER_LABEL 과 samples.csv type 이 안 맞으면 암 0개 → 중단
if int(y.sum()) == 0:
    sys.exit('samples.csv 에 type=%s 인 검체가 없습니다. CANCER_LABEL 과 1단계가 쓴 값이'
             ' 어긋났습니다 (현재 CANCER_LABEL=%s).' % (_CANCER, _CANCER))
_c = [str(i) for i in range(X.shape[1])]
_D = pd.DataFrame(X, columns=_c); _D['_g'] = M['group'].values; _D['_y'] = y
_A = _D.groupby('_g', sort=True).mean()
y = (_A['_y'].values > 0.5).astype(int)
X = _A[_c].values.astype(np.float32)
mu = np.nanmean(X, axis=0); mu = np.where(np.isfinite(mu), mu, 0.0)
X = np.where(np.isfinite(X), X, mu)
print('단위 %d (암 %d) · 블록 %d' % (len(y), int(y.sum()), X.shape[1]))

d, t, p, dft, d0, s0 = moderated_ttest(X[y == 1], X[y == 0])
q = bh(p)
T = pd.DataFrame({'chr': B.index.get_level_values(0).astype(str),
                  'blk': B.index.get_level_values(1).astype('int64'),
                  'dbeta': d, 'q': q})

J = pd.read_parquet(JSD)
J['chr'] = J['chr'].astype(str); J['blk'] = J['blk'].astype('int64')
Z = T.merge(J, on=['chr', 'blk'], how='inner')
print('JSD 계산된 블록 %d / 후보 %d (%.1f%%)'
      % (len(Z), len(T), 100.0 * len(Z) / len(T)))

# 위치순 정렬 후 안정 정렬 - 동률 블록에서 N 경계가 매번 같게
Z = (Z.sort_values(['chr', 'blk'])
       .sort_values(COL, ascending=False, kind='mergesort').reset_index(drop=True))
P = Z.head(N).copy()

def stat(D, tag):
    nq = int((D.q < 0.05).sum()); nd = int((D.dbeta.abs() >= 30).sum())
    nb = int(((D.q < 0.05) & (D.dbeta.abs() >= 30)).sum())
    print('  %-12s q<0.05 %5d (%5.1f%%) · |dbeta|>=30 %5d (%5.1f%%) · 둘 다 %5d (%5.1f%%)'
          % (tag, nq, 100.0*nq/len(D), nd, 100.0*nd/len(D), nb, 100.0*nb/len(D)))

print('\n기존 기준으로 보면 (선정에는 안 씀)')
stat(P, '상위 %d' % N)
rs = np.random.RandomState(0)
stat(Z.iloc[rs.choice(len(Z), min(N, len(Z)), replace=False)], '무작위 %d' % N)
stat(Z, '전체')

print('\n상위 %d 의 %s   최소 %.3f · 중앙 %.3f · 최대 %.3f'
      % (N, COL, P[COL].min(), P[COL].median(), P[COL].max()))
print('저메틸(암<정상) %d · 고메틸 %d'
      % (int((P.dbeta < 0).sum()), int((P.dbeta > 0).sum())))
print('칸당 창 중앙 %d · 창 리드 중앙 %d' % (int(P.nwin.median()), int(P.rd.median())))

os.makedirs(OUT, exist_ok=True)
f = OUT + '/panel_jsd_%s_n%d.csv' % (COL, N)
P.to_csv(f, index=False)
print('\n저장: ' + f)
