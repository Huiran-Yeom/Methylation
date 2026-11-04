# -*- coding: utf-8 -*-
"""귀무A 세 판 (rand1-3): K=3 창이 서는 블록에서 무작위 균등 추출
jsd 계열 판과 비교용. bl200 은 귀무B(10_null_b.py)와 비교
씨앗 20260914 + i 고정

  사용
    python 09_null_a.py                  # 실제 자리에 씀
    python 09_null_a.py --out /tmp/test   # 다른 자리에 써서 대조만
    python 09_null_a.py --cand <후보폴더>
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
SEED0 = 20260914          # 씨앗 i = SEED0 + i (고정)
NPICK = 200
MINCPG = 3                # K=3 창이 서려면 블록에 CpG 3개 이상
a = sys.argv
if '--cand' in a:
    CAND = a[a.index('--cand') + 1].rstrip('/')
if '--out' in a:
    OUT = a[a.index('--out') + 1].rstrip('/')

# 칸 수는 jsd200 패널에 맞춤 (보통 200, 교집합 패널이면 그 크기)
_js = '%s/jsd200/DMR_confirmed_jsd200.csv' % OUT
if os.path.isfile(_js):
    NPICK = len(pd.read_csv(_js))
# NULL_A_N 주면 그 칸 수 사용
if _os.environ.get('NULL_A_N'):
    NPICK = int(_os.environ['NULL_A_N'])
print('뽑을 칸 수 %d (jsd200 패널 %s%s)' % (NPICK, '있음' if os.path.isfile(_js) else '없음 — 기본 200', ' · NULL_A_N 지정' if _os.environ.get('NULL_A_N') else ''))

B = pd.read_parquet(CAND + '/blocks_cpg.parquet')
cnt = B.groupby(['chr', 'blk']).size()
pool = [(str(c), int(b)) for (c, b), n in cnt.items() if n >= MINCPG]
print('고를 수 있는 블록(K=3 창 서는 것) %d개' % len(pool))

used = set()
for i in (1, 2, 3):
    rs = np.random.RandomState(SEED0 + i)
    idx = rs.choice(len(pool), NPICK, replace=False)
    pick = [pool[j] for j in idx]
    d = pd.DataFrame({'Feature': ['%s_%d' % (c, b) for c, b in pick],
                      'dbeta': [0.0] * NPICK, 'q': [1.0] * NPICK})
    o = '%s/rand%d200' % (OUT, i)
    os.makedirs(o, exist_ok=True)
    d.to_csv('%s/DMR_confirmed_rand%d200.csv' % (o, i), index=False)
    with open(o + '/_origin.txt', 'w') as h:
        h.write('무작위 패널 씨앗 %d: 2026-09-14\n'
                '  K=3 창이 서는 블록 %d개에서 %d개를 비복원 추출.\n'
                '  암과 무관하게 골랐다. 실제 패널과 견줄 무작위 기준 패널로 쓴다.\n'
                '  창이 서는 블록에서 고르므로 리드 기반 피처를 전 칸 받는다 —\n'
                '  bl200(104칸)에 대해서는 보수적인(넘기 어려운) 귀무다.\n'
                % (i, len(pool), NPICK))
    used |= set(pick)
    print('rand%d200  %d칸' % (i, NPICK))
print('세 씨앗 합집합 %d칸 (겹침 %d)' % (len(used), 3 * NPICK - len(used)))
