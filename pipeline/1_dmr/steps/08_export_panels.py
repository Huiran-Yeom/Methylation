# -*- coding: utf-8 -*-
"""실판 셋을 검증이 읽는 형식으로 내보낸다 (2026-09-14 · 파일로 굳힘 2026-09-26)

세 판은 서로 다른 코드가 만들어서 열 이름도 자리도 다르다. 검증(`config.py` 의
`DMR_CSV`)은 `Feature,dbeta,q` 한 가지만 읽으므로 여기서 한 형식으로 맞춘다.
`Feature` 는 `<chr>_<blk>` 다.

  bl200   j_panel_dmr15/bl200_d10_cellline/DMR_confirmed_m.csv   (07_panel.py 산출)
  jsd200  j_panel15/panel_jsd_jsd_mean_n200.csv                (06_select.py 산출)
  jsdb200 j_panel15_supplement/panel_jsd_supplement_n200.csv

원래 인라인으로 돌려 파일이 없었다. 09-26 에 세션 기록에서 찾아 그대로 옮겼다.
계산은 한 글자도 바꾸지 않았다. 바꾸면 기존 패널과 달라진다.

  사용
    python 08_export_panels.py                  # 실제 자리에 쓴다
    python 08_export_panels.py --out /tmp/test    # 다른 자리에 써서 대조만
"""

# 이 파일은 스크립트다. 다른 코드가 import 하면 본체가 그대로 돌아
# 실제 자료를 덮어쓴다. 실행은 `python <파일>` 로만 한다.
if __name__ != '__main__':
    raise ImportError(__file__ + ' 은(는) 스크립트다 — import 하지 않는다')

import os as _os, sys as _sys
# 0_setup 을 찾는다. 위로 올라가며 → METH_CONF_DIR → 기본값.
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

# 세대를 박아 두면 GEN 을 바꿨을 때 「이전 세대」 후보로 패널을 내보낸다.
#   그러면 실제 패널과 귀무 패널의 세대가 갈리는데, 여기가 대조의 기준선이다.
_GEN = _os.environ.get('GEN', '')
if not _GEN:
    try:
        _GEN = _cfg.GEN
    except Exception:
        raise SystemExit('GEN(패널 세대)을 못 읽었습니다 — config.conf 의 GEN 을 채우거나 '
                         '환경변수로 주십시오. 비워 두면 옛 세대 폴더를 가리킵니다.')
# ╔══ 이 단계가 쓰는 경로 ═══════════════════════════════════════════╗
#   배치가 다르면 「여기만」 고친다. 0_setup/config.conf 에는 뿌리 둘
#   (METH_ROOT=_M · DATASET_ROOT=_D)과 실험 설정만 둔다.
#   환경변수나 명령줄 인자로 그때그때 덮어쓸 수 있다(아래 참조).
R = _M + '/results/dmr'
OUT = R + '/j_panel_dmr' + _GEN + '_panels'
# ╚══════════════════════════════════════════════════════════════════╝
a = sys.argv
if '--r' in a:
    R = a[a.index('--r') + 1].rstrip('/')
    OUT = R + '/j_panel_dmr' + _GEN + '_panels'
if '--out' in a:
    OUT = a[a.index('--out') + 1].rstrip('/')

SRC = {'bl200':   R + '/j_panel_dmr' + _GEN + '/bl200_d10_cellline/DMR_confirmed_m.csv',
       'jsd200':  R + '/j_panel' + _GEN + '/panel_jsd_jsd_mean_n200.csv',
       'jsdb200': R + '/j_panel' + _GEN + '_supplement/panel_jsd_supplement_n200.csv'}

# 없는 입력에서 pandas traceback 만 남으면 「누가 만드는지」 를 모른다.
#   어느 단계를 먼저 돌려야 하는지 같이 말한다.
_MADE_BY = {'bl200':   '07_panel.py  (Baseline 순위)',
            'jsd200':  '06_select.py (JSD 상위 N)',
            'jsdb200': '06b_select_supplement.py 를 먼저 돌리십시오 (j_panel<GEN>_supplement 에 씁니다)'}
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
    # 200칸이 아니거나 같은 칸이 두 번 들어오면 멈춘다. 조용히 넘기면
    # 학습과 검증이 서로 다른 패널을 본다.
    assert len(f) == 200 and f['Feature'].is_unique, (panel, len(f))
    os.makedirs(OUT + '/' + panel, exist_ok=True)
    f.to_csv('%s/%s/DMR_confirmed_%s.csv' % (OUT, panel, panel), index=False)
    with open('%s/%s/_origin.txt' % (OUT, panel), 'w') as h:
        h.write('원본 %s\n변환 2026-09-14 · Feature=chr_blk · 200블록\n' % src)
    print('%-8s %3d블록  <- %s' % (panel, len(f), '/'.join(src.split('/')[-2:])))

# 판끼리 얼마나 겹치는지. 겹침이 크면 「판 비교」 가 묶음 비교가 된다
F = {panel: set(pd.read_csv('%s/%s/DMR_confirmed_%s.csv' % (OUT, panel, panel))['Feature'])
     for panel in SRC}
for x, y in itertools.combinations(SRC, 2):
    n = len(F[x] & F[y])
    print('겹침 %-8s x %-8s %3d칸 (%.1f%%)' % (x, y, n, 100 * n / 200))
print('합집합 %d칸' % len(set().union(*F.values())))
