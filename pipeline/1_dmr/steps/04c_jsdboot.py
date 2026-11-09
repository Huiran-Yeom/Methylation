# -*- coding: utf-8 -*-
"""JSD 부트스트랩 순위 안정성
정상, 암 안에서 단위(사람 또는 세포주)를 복원추출해 JSD 다시 계산, 상위 TOPN 에 든 횟수(freq)를 셈
산출: <out>/boot_K<K>_minr<MINR>.parquet (열: chr, blk, jsd_obs, freq) → 06b 가 jsdb200 선정에 사용

사용: python 04c_jsdboot.py [--dir ..] [--cand ..] [--out ..] [--k K] [--minr N] [--frac F]
                            [--need N] [--needg N] [--nperm N] [--topn N]
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
from jsd_resample import setup   # noqa: E402  (04b · 04c 공통: 설정 · 자료 · 관측 JSD)
S = setup(sys.argv, _cfg)

# 안정성: 라벨은 그대로, 단위만 복원추출. 매 회 JSD 상위 TOPN 에 든 블록에 +1
TOPN = int(S.a[S.a.index('--topn') + 1]) if '--topn' in S.a else 1000
rs = np.random.RandomState(S.SEED)
print('안정성: 복원추출 %d회 · JSD 상위 %d 에 든 횟수' % (S.NPERM, TOPN))
freq = pd.Series(0, index=S.Bo.index, dtype=int)
gi, ni = S.allu[S.NN:], S.allu[:S.NN]
for p in range(S.NPERM):
    pg = rs.choice(gi, len(gi), replace=True)
    pn = rs.choice(ni, len(ni), replace=True)
    J, ok = S.jsd_of(pg, pn)
    hit = S.by_block(J, ok).sort_values(ascending=False).head(TOPN).index
    freq[freq.index.isin(hit)] += 1
    if (p + 1) % 20 == 0:
        print('  복원추출 %d/%d' % (p + 1, S.NPERM), flush=True)
out = pd.DataFrame({'jsd_obs': S.Bo, 'freq': freq}).reset_index()
os.makedirs(S.OUT, exist_ok=True)   # 출력 폴더 없으면 만듦
fp = '%s/boot_K%d_minr%d.parquet' % (S.OUT, S.K, S.MINR)
out.to_parquet(fp)
print('freq: 만점 %d칸 · 90회이상 %d칸 · 50회이상 %d칸 · 1회이상 %d칸'
      % (int((freq == S.NPERM).sum()), int((freq >= 0.9 * S.NPERM).sum()),
         int((freq >= 0.5 * S.NPERM).sum()), int((freq > 0).sum())))
print('저장: %s' % fp)
