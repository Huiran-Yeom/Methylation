# -*- coding: utf-8 -*-
"""S0e JSD 순위로 패널 뽑기 (2026-08-20)

지시대로 JSD 순위 상위 N 개를 그대로 패널로 쓴다.
q값과 |dbeta| 는 '기존 기준으로 보면 얼마나 겹치나' 를 보고하려고 함께 계산만 한다.
선정에는 쓰지 않는다.
"""

# 이 파일은 스크립트다. 다른 코드가 import 하면 본체가 그대로 돌아
# 실제 자료를 덮어쓴다. 실행은 `python <파일>` 로만 한다.
if __name__ != '__main__':
    raise ImportError(__file__ + ' 은(는) 스크립트다 — import 하지 않는다')

import os as _os, sys as _sys
# 0_setup 을 찾는다. 이 파일이 «복사되어» 정리본 밖에서 돌 수도 있으므로
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
        '설정을 못 찾았습니다 — 0_setup/meth_config.py 가 있는 자리를 찾지 못했습니다.\n'
        '  이 파일은 정리본 안에서 도는 것이라 기본값으로 넘어가지 않습니다.\n'
        '  트리 밖에서 돌리려면 METH_CONF_DIR 로 0_setup 자리를 주십시오.')


import os, sys
import numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from moderated_t import moderated_ttest, bh

# 2026-09-26: 기본값에 세대 접미사가 없어 «옛 세대» 를 가리켰다.
#   09-23 사고가 옛 세대 폴더를 덮은 원인이 바로 이것이다. GEN 을 따른다.
_GEN = _os.environ.get('GEN', '')
if not _GEN:
    try:
        _GEN = _cfg.GEN
    except Exception:
        raise SystemExit('GEN(패널 세대)을 못 읽었습니다 — config.conf 의 GEN 을 채우거나 '
                         '환경변수로 주십시오. 비워 두면 옛 세대 폴더를 가리킵니다.')
# ╔══ 이 단계가 쓰는 경로 ═══════════════════════════════════════════╗
#   배치가 다르면 «여기만» 고칩니다. 0_setup/config.conf 에는 뿌리 둘
#   (METH_ROOT=_M · DATASET_ROOT=_D)과 실험 설정만 있습니다.
#   환경변수나 명령줄 인자로 그때그때 덮어쓸 수도 있습니다(아래 참조).
# 2026-09-26: 기본값이 `j_candidates_d5` · `jsd_K3_minr4` 였다 — **깊이5 감도 변종**이고
#   minr 도 4 다. 실제로 쓴 것은 `j_candidates<GEN>` 과 `jsd_K3_minr6` 다.
#   그대로 두면 인자를 빠뜨린 사람이 «다른 JSD 표»로 패널을 고른다. 필수로 바꾼다.
CAND = _M + '/results/dmr/j_candidates' + _GEN
JSD  = None
OUT  = _M + '/results/dmr/j_panel' + _GEN
# ╚══════════════════════════════════════════════════════════════════╝
N, COL = 300, 'jsd_mean'
a = sys.argv
if '--cand' in a: CAND = a[a.index('--cand')+1]
if '--jsd'  in a: JSD  = a[a.index('--jsd')+1]
if '--out'  in a: OUT  = a[a.index('--out')+1]
if '--n'    in a: N    = int(a[a.index('--n')+1])
if '--col'  in a: COL  = a[a.index('--col')+1]
# 2026-09-29: `is None` 이면 드라이버가 넘기는 빈 문자열(--jsd "")이 관문을 비켜간다.
#   run_dmr.sh 는 `ls ... | head -1` 로 만들어 파일이 없으면 빈 값을 넘긴다.
if not JSD:
    sys.exit('사용: python 06_select.py --jsd <jsd_K?_minr?.parquet> [--cand ..] [--out ..] [--n N]'
             '   ·  --jsd 에 기본값을 두지 않는다 — minr 이 다른 표를 읽으면 다른 패널이 나온다.'
             '   ·  09-10 에 쓴 값: --jsd %s/results/dmr/j_jsd%s/jsd_K3_minr6.parquet --n 200'
             % (_M, _GEN))

B = pd.read_parquet(CAND + '/beta_matrix.parquet')
M = pd.read_csv(CAND + '/samples.csv')
X = B.T.reindex(M['sample']).values.astype(np.float32)
y = (M['type'] == 'GBM').astype(int).values
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

# 2026-09-29: 기본 quicksort 는 불안정해 동률 블록의 N=200 경계가 빌드에 달린다.
#   04_coverage_check 와 같은 방식으로 인덱스 먼저·안정 정렬로 결정적으로 만든다.
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
