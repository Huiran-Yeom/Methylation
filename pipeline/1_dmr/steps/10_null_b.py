# -*- coding: utf-8 -*-
"""귀무B 세 판 — bl200 의 CpG 수 분포에 맞춰 층화 추출 (2026-09-14 · 파일로 굳힘 2026-09-26)

귀무A 만 두면 bl200 이 못 넘겼을 때 「bl200 블록이 쓸모없다」와
「bl200 은 칸이 절반이다」를 구별할 수 없다. 귀무B 는 bl200 의 칸 수 불리함을
**그대로 진다** — 그래서 bl200 의 잣대다. jsd/jsdb 의 잣대는 귀무A. 사전등록 3-1절.

모집단은 `blocks_frac060.csv` 23,532칸 — Baseline 이 고르는 그 모집단이다.
CpG 수로 버킷을 만들고(8 이상은 한 칸으로 묶음) bl200 의 분포를 그대로 따라 뽑는다.
실측으로 세 씨앗 모두 창 서는 칸 104 · CpG<3 인 칸 96 — bl200 과 정확히 일치한다.

씨앗은 `20260914 + 100 + i` 다. 귀무A(`20260914 + i`)와 겹치지 않게 100 을 띄웠다.

원래 인라인으로 돌려 파일이 없었다. 09-26 에 세션 기록에서 찾아 그대로 옮겼다.

  사용
    python 10_null_b.py                  # 실제 자리에 쓴다
    python 10_null_b.py --out /tmp/test   # 다른 자리에 써서 대조만
"""

# 이 파일은 스크립트다. 다른 코드가 import 하면 본체가 그대로 돌아
# 실제 자료를 덮어쓴다. 실행은 `python <파일>` 로만 한다.
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
        '설정을 못 찾았습니다 — 0_setup/meth_config.py 가 있는 자리를 찾지 못했습니다.\n'
        '  이 파일은 정리본 안에서 도는 것이라 기본값으로 넘어가지 않습니다.\n'
        '  트리 밖에서 돌리려면 METH_CONF_DIR 로 0_setup 자리를 주십시오.')


import collections
import os
import sys
import numpy as np
import pandas as pd

# 2026-09-26: 세대를 박아 두면 GEN 을 바꿨을 때 «이전 세대» 후보로 귀무를 만든다.
#   그러면 실제 패널과 귀무 패널의 세대가 갈리는데, 여기가 대조의 기준선이다.
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
R = _M + '/results/dmr'
OUT = R + '/j_panel_dmr' + _GEN + '_panels'
CAND = R + '/j_candidates' + _GEN
# ╚══════════════════════════════════════════════════════════════════╝
SEED0 = 20260914 + 100    # 씨앗 i 는 SEED0 + i. 고정값이다
CAP = 8                   # CpG 수 8 이상은 한 버킷으로 묶는다
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

# bl200 의 CpG 수 분포를 그대로 베낀다 — 그래야 칸 수 불리함이 같아진다
# 2026-09-28: bl200 의 CpG 수 분포를 베끼는 것이 이 판의 핵심이다.
#   그 파일이 없으면 «누가 만드는지» 를 말하고 멈춘다.
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
    # 버킷이 모자라 200칸을 못 채우면 분포가 어긋난 것이다. 멈춘다.
    assert len(pick) == 200, len(pick)
    d = pd.DataFrame({'Feature': ['%s_%d' % (c, b) for c, b in pick],
                      'dbeta': [0.0] * 200, 'q': [1.0] * 200})
    o = '%s/randb%d200' % (OUT, i)
    os.makedirs(o, exist_ok=True)
    d.to_csv('%s/DMR_confirmed_randb%d200.csv' % (o, i), index=False)
    with open(o + '/_origin.txt', 'w') as h:
        h.write('무작위 패널 귀무B 씨앗 %d — 2026-09-14\n'
                '  모집단: blocks_frac060.csv 23,532칸 (Baseline 이 고르는 그 모집단)\n'
                '  bl200 의 CpG 수 분포에 맞춰 층화 추출. bl200 의 칸 수 불리함을 그대로 진다.\n'
                '  bl200 의 잣대로 쓴다. jsd/jsdb 의 잣대는 귀무A(rand*200).\n' % i)
    have3 = sum(1 for k in pick if cnt.get(k, 0) >= 3)
    print('randb%d200  200칸 · K=3창 서는 칸 %d · CpG<3 인 칸 %d'
          % (i, have3, sum(1 for k in pick if cnt.get(k, 0) < 3)))
