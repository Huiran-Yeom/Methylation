# -*- coding: utf-8 -*-
"""귀무B 세 판 (randb1-3): bl200 의 CpG 수 분포에 맞춰 층화 추출
bl200 비교용 (창 안 서는 칸 비율까지 bl200 과 같게). jsd 계열은 귀무A
모집단 blocks_frac060.csv (Baseline 과 같은 모집단), CpG 8개 이상은 한 버킷
씨앗 20260914 + 100 + i 고정 (귀무A 와 안 겹치게)

  사용
    python 10_null_b.py                  # 실제 자리에 씀
    python 10_null_b.py --out /tmp/test   # 다른 자리에 써서 대조만
    python 10_null_b.py --cand <후보폴더>
"""

# 스크립트 전용. import 하면 본체가 돌아 자료를 덮어씀
if __name__ != '__main__':
    raise ImportError(__file__ + ' 은(는) 스크립트다 — import 하지 않는다')

import os as _os, sys as _sys
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
    _M = _cfg.METH_ROOT
except ImportError:
    raise SystemExit(
        '설정을 못 찾았습니다. 0_setup/meth_config.py 가 있는 자리를 찾지 못했습니다.\n'
        '  이 파일은 정리본 안에서 도는 것이라 기본값으로 넘어가지 않습니다.\n'
        '  트리 밖에서 돌리려면 METH_CONF_DIR 로 0_setup 자리를 주십시오.')


import collections
import os
import sys
import numpy as np
import pandas as pd

# 폴더 이름에 GEN(패널 세대) 붙임. 박아 두면 GEN 바꿔도 옛 세대 후보로 귀무를 만듦
_GEN = _os.environ.get('GEN', '')
if not _GEN:
    try:
        _GEN = _cfg.GEN
    except Exception:
        raise SystemExit('GEN(패널 세대)을 못 읽었습니다 — config.conf 의 GEN 을 채우거나 '
                         '환경변수로 주십시오. 비워 두면 옛 세대 폴더를 가리킵니다.')
# 이 단계 경로. 폴더 구조 다르면 여기만 수정 (인자로 덮어쓰기 가능)
R = _M + '/results/dmr'
OUT = R + '/j_panel_dmr' + _GEN + '_panels'
CAND = R + '/j_candidates' + _GEN
SEED0 = 20260914 + 100    # 씨앗 i = SEED0 + i (고정)
CAP = 8                   # CpG 8개 이상은 한 버킷
a = sys.argv
if '--cand' in a:
    CAND = a[a.index('--cand') + 1].rstrip('/')
if '--out' in a:
    OUT = a[a.index('--out') + 1].rstrip('/')

uni = pd.read_csv(CAND + '/blocks_frac060.csv')
U = [(str(c), int(b)) for c, b in zip(uni['chr'], uni['blk'])]
B = pd.read_parquet(CAND + '/blocks_cpg.parquet')
cnt = collections.Counter()
for (c, b), n in B.groupby(['chr', 'blk']).size().items():
    cnt[(str(c), int(b))] = int(n)

# bl200 의 CpG 수 분포를 그대로 따름. bl200 파일 없으면 중단
_bl = '%s/bl200/DMR_confirmed_bl200.csv' % OUT
if not os.path.isfile(_bl):
    raise SystemExit(
        ('입력이 없다: %s' % _bl)
        + '  |  08_export_panels.py 가 만든다 (그 앞은 07_panel.py).'
        + '  귀무B 는 bl200 의 CpG 수 분포를 맞춰 뽑으므로 bl200 이 먼저 있어야 한다.')
bl = pd.read_csv(_bl)
blk = [(x.split('_')[0], int(x.split('_')[1])) for x in bl['Feature']]
hist = collections.Counter(min(cnt.get(k, 0), CAP) for k in blk)
print('bl200 CpG 수 분포(8+ 묶음):', dict(sorted(hist.items())))
NPICK = len(blk)          # 칸 수도 bl200 에 맞춤 (보통 200, 교집합 패널이면 그 크기)

bucket = collections.defaultdict(list)
for k in U:
    bucket[min(cnt.get(k, 0), CAP)].append(k)
print('모집단 버킷 크기:', {k: len(v) for k, v in sorted(bucket.items())})

for i in (1, 2, 3):
    rs = np.random.RandomState(SEED0 + i)
    pick = []
    for c, n in sorted(hist.items()):
        pool = bucket[c]
        idx = rs.choice(len(pool), min(n, len(pool)), replace=False)
        pick += [pool[j] for j in idx]
    # 버킷이 모자라 칸 수 못 채우면 중단
    assert len(pick) == NPICK, (len(pick), NPICK)
    d = pd.DataFrame({'Feature': ['%s_%d' % (c, b) for c, b in pick],
                      'dbeta': [0.0] * NPICK, 'q': [1.0] * NPICK})
    o = '%s/randb%d200' % (OUT, i)
    os.makedirs(o, exist_ok=True)
    d.to_csv('%s/DMR_confirmed_randb%d200.csv' % (o, i), index=False)
    with open(o + '/_origin.txt', 'w') as h:
        h.write('무작위 패널 귀무B 씨앗 %d: 2026-09-14\n'
                '  모집단: blocks_frac060.csv %s칸 (Baseline 이 고르는 그 모집단)\n'
                '  bl200 의 CpG 수 분포에 맞춰 층화 추출. bl200 의 칸 수 불리함을 그대로 진다.\n'
                '  bl200 의 잣대로 쓴다. jsd/jsdb 의 잣대는 귀무A(rand*200).\n' % (i, format(len(U), ',')))
    have3 = sum(1 for k in pick if cnt.get(k, 0) >= 3)
    print('randb%d200  %d칸 · K=3창 서는 칸 %d · CpG<3 인 칸 %d'
          % (i, NPICK, have3, sum(1 for k in pick if cnt.get(k, 0) < 3)))
