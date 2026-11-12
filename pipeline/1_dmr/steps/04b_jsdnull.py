# -*- coding: utf-8 -*-
"""JSD 귀무 분포 (순열 검정, 문턱 정하기용)
라벨을 섞어 무작위 단위를 암으로 지정하고 JSD 다시 계산, NPERM 번 반복
섞는 단위는 사람 또는 세포주 (복제본은 함께 움직임), 집단 크기는 유지

사용: python 04b_jsdnull.py [--dir ..] [--cand ..] [--out ..] [--k K] [--minr N] [--frac F]
                            [--need N] [--needg N] [--nperm N]
"""

# 스크립트 전용. import 하면 본체가 돌아 자료를 덮어씀
if __name__ != '__main__':
    raise ImportError(__file__ + ' 은(는) 스크립트다. import 하지 않는다')

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
from jsd_resample import setup   # noqa: E402  (04b · 04c 공통: 설정 · 자료 · 관측 JSD)
S = setup(sys.argv, _cfg)

rs = np.random.RandomState(S.SEED)
nullw, nullb = [], []
for p in range(S.NPERM):
    perm = rs.permutation(S.allu)
    J, ok = S.jsd_of(perm[:S.NG], perm[S.NG:])
    nullw.append(J[ok]); nullb.append(S.by_block(J, ok).values)
    if (p + 1) % 20 == 0: print('  순열 %d/%d' % (p + 1, S.NPERM), flush=True)
nw_all = np.concatenate(nullw); nb_all = np.concatenate(nullb)

print('\n귀무 창JSD  중앙 %.3f · 99%% %.3f · 99.9%% %.3f · 최대 %.3f'
      % (np.median(nw_all), np.percentile(nw_all, 99),
         np.percentile(nw_all, 99.9), nw_all.max()))
print('귀무 칸JSD  중앙 %.3f · 99%% %.3f · 99.9%% %.3f · 최대 %.3f'
      % (np.median(nb_all), np.percentile(nb_all, 99),
         np.percentile(nb_all, 99.9), nb_all.max()))

print('\n문턱 후보 (칸 기준)')
print('  귀무백분위      t       관측 통과 칸   기대 가짜')
for q in (99.0, 99.5, 99.9, 99.95, 99.99):
    t  = float(np.percentile(nb_all, q))
    hit = int((S.Bo >= t).sum())
    exp = len(S.Bo) * (1 - q / 100.0)
    print('  %7.2f%%   %.4f   %9d   %9.1f' % (q, t, hit, exp))

os.makedirs(S.OUT, exist_ok=True)
pd.DataFrame({'null_block_jsd': nb_all}).to_parquet(
    S.OUT + '/null_K%d_minr%d.parquet' % (S.K, S.MINR), index=False)
print('\n저장: %s/null_K%d_minr%d.parquet' % (S.OUT, S.K, S.MINR))
