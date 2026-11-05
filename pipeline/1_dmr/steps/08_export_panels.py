# -*- coding: utf-8 -*-
"""실제 패널 세 판을 2단계가 읽는 형식(Feature,dbeta,q)으로 내보냄
Feature = <chr>_<blk>, N = PANEL_N (기본 200)

  bl<N>   j_panel_dmr<GEN>/bl<N>_d10_cellline/DMR_confirmed_m.csv   (07_panel.py)
  jsd<N>  j_panel<GEN>/panel_jsd_jsd_mean_n<N>.csv                  (06_select.py)
  jsdb<N> j_panel<GEN>_supplement/panel_jsd_supplement_n<N>.csv     (06b_select_supplement.py)

  사용
    python 08_export_panels.py                  # 실제 자리에 씀
    python 08_export_panels.py --out /tmp/test    # 다른 자리에 써서 대조만
    python 08_export_panels.py --r <results/dmr>  # 입력 폴더 바꿈
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
    _M = _cfg.METH_ROOT
except ImportError:
    raise SystemExit(
        '설정을 못 찾았습니다. 0_setup/meth_config.py 가 있는 자리를 찾지 못했습니다.\n'
        '  이 파일은 정리본 안에서 도는 것이라 기본값으로 넘어가지 않습니다.\n'
        '  트리 밖에서 돌리려면 METH_CONF_DIR 로 0_setup 자리를 주십시오.')


import itertools
import os
import sys
import pandas as pd

# 폴더 이름에 GEN(패널 세대) 붙임. 박아 두면 GEN 바꿔도 옛 세대를 읽음
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
a = sys.argv
if '--r' in a:
    R = a[a.index('--r') + 1].rstrip('/')
    OUT = R + '/j_panel_dmr' + _GEN + '_panels'
if '--out' in a:
    OUT = a[a.index('--out') + 1].rstrip('/')

N = _cfg.PANEL_N    # 패널 블록 수 (config.conf PANEL_N)
SRC = {'bl%d' % N:   R + '/j_panel_dmr' + _GEN + '/bl%d_d10_cellline/DMR_confirmed_m.csv' % N,
       'jsd%d' % N:  R + '/j_panel' + _GEN + '/panel_jsd_jsd_mean_n%d.csv' % N,
       'jsdb%d' % N: R + '/j_panel' + _GEN + '_supplement/panel_jsd_supplement_n%d.csv' % N}

# 입력 없을 때 어느 단계를 먼저 돌릴지 알려 주는 용도
_MADE_BY = {'bl%d' % N:   '07_panel.py  (Baseline 순위)',
            'jsd%d' % N:  '06_select.py (JSD 상위 N)',
            'jsdb%d' % N: '06b_select_supplement.py 를 먼저 돌리십시오 (j_panel<GEN>_supplement 에 씁니다)'}
# 옛 트리(200 블록)는 bl200, jsdb200 입력이 한글 이름. 새 이름 없을 때만 옛 이름 읽음
_OLD = {}
if N == 200:
    _OLD = {'bl200':   SRC['bl200'].replace('bl200_d10_cellline', 'bl200_d10_세포주'),
            'jsdb200': SRC['jsdb200'].replace('_supplement/panel_jsd_supplement_n200.csv',
                                              '_보완/panel_jsd보완_n200.csv')}
for _p, _o in _OLD.items():
    if not os.path.isfile(SRC[_p]) and os.path.isfile(_o):
        print('  옛 이름으로 읽는다: %s' % _o)
        SRC[_p] = _o

# 설정 PANELS 에 없는 판은 건너뜀 (예: N 칸 못 채운 jsdb)
_skip = [p for p in SRC if p not in _cfg.PANELS]
if _skip:
    print('  PANELS(%s) 에 없어 건너뜀: %s' % (','.join(_cfg.PANELS), ','.join(_skip)))
SRC = {p: v for p, v in SRC.items() if p in _cfg.PANELS}
if not SRC:
    raise SystemExit('PANELS 에 내보낼 판이 없다: %s' % _cfg.PANELS)

for panel, src in SRC.items():
    if not os.path.isfile(src):
        raise SystemExit(
            ('입력이 없다: %s' % src)
            + ('  |  %s 판은 ' % panel) + _MADE_BY.get(panel, '앞 단계')
            + ' 가 만든다. 그 단계를 먼저 돌리십시오.')

for panel, src in SRC.items():
    d = pd.read_csv(src)
    f = pd.DataFrame({'Feature': d['chr'].astype(str) + '_' + d['blk'].astype(str),
                      'dbeta':   d['dbeta'].astype(float),
                      'q':       d['q'].astype(float)})
    # N 칸 아니거나 중복 칸 있으면 중단
    assert len(f) == N and f['Feature'].is_unique, (panel, len(f))
    os.makedirs(OUT + '/' + panel, exist_ok=True)
    f.to_csv('%s/%s/DMR_confirmed_%s.csv' % (OUT, panel, panel), index=False)
    with open('%s/%s/_origin.txt' % (OUT, panel), 'w') as h:
        h.write('원본 %s\n변환 2026-09-14 · Feature=chr_blk · %d블록\n' % (src, N))
    print('%-8s %3d블록  <- %s' % (panel, len(f), '/'.join(src.split('/')[-2:])))

# 판끼리 겹침 출력 (많이 겹치면 판 비교 의미가 약함)
F = {panel: set(pd.read_csv('%s/%s/DMR_confirmed_%s.csv' % (OUT, panel, panel))['Feature'])
     for panel in SRC}
for x, y in itertools.combinations(SRC, 2):
    n = len(F[x] & F[y])
    print('겹침 %-8s x %-8s %3d칸 (%.1f%%)' % (x, y, n, 100 * n / N))
print('합집합 %d칸' % len(set().union(*F.values())))
