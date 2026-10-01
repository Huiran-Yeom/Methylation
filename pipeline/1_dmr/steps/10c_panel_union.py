# -*- coding: utf-8 -*-
"""11단계 풀 BAM 이 쓰는 입력 다섯을 만든다.

  panel_union3_cellline.csv   bl200 · jsd200 · jsdb200 의 합집합
  panel_union_rand3.csv       rand1-3 의 합집합
  panel_union_randb3.csv      randb1-3 의 합집합
  use<GEN>_normal.txt         풀에 넣을 정상 검체 이름
  use<GEN>_gbm.txt            풀에 넣을 암 검체 이름

전에는 이 다섯을 손으로 만들었다. 그래서 새 자료로 1단계를 처음부터 돌리면
11단계가 「패널 목록이 없다」 로 멈췄다. 판 파일과 samples.csv 에서 전부
계산되는 값이므로 여기서 만든다.

합집합은 판의 Feature(<chr>_<blk>)를 chr·blk 로 쪼개 중복을 없앤 것이다.
정렬은 chr 를 글자로, blk 를 수로 본다. 줄 끝은 LF 다 — 손으로 만든 기존
다섯 파일과 바이트 단위로 같은 것을 확인했다(2026-10-01).
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
except Exception as _x:
    raise SystemExit(
        '설정을 못 찾았습니다. 0_setup/meth_config.py 가 있는 자리를 찾지 못했습니다.\n'
        '  METH_CONF_DIR 로 알려 주십시오.  (%s)' % _x)

import csv
import os
import sys

NL = chr(10)

GEN = os.environ.get('GEN') or getattr(_cfg, 'GEN', '15')
R = os.environ.get('DMR_OUT') or (_M + '/results/dmr')

a = sys.argv
if '--gen' in a: GEN = a[a.index('--gen') + 1]
if '--dmr' in a: R = a[a.index('--dmr') + 1].rstrip('/')

PANELS = list(getattr(_cfg, 'PANELS', ['bl200', 'jsd200', 'jsdb200']))
NULLS = list(getattr(_cfg, 'NULL_PANELS', []))
_CANCER = getattr(_cfg, 'CANCER_LABEL', 'GBM')

PDIR = R + '/j_panel_dmr' + GEN + '_panels'
OUT = R + '/j_panel' + GEN
CAND = R + '/j_candidates' + GEN

# 판을 풀 종류로 가른다. 2_train/steps/config.py 가 PANEL 로 풀을 고르는 규칙과 같다.
#   어긋나면 귀무판이 실판 풀을 보고 200칸 중 몇 칸으로 조용히 돈다.
GROUPS = {'cellline': [], 'rand': [], 'randb': []}
for _p in PANELS + NULLS:
    if _p.startswith('randb'):  GROUPS['randb'].append(_p)
    elif _p.startswith('rand'): GROUPS['rand'].append(_p)
    else:                       GROUPS['cellline'].append(_p)

NAME = {'cellline': 'panel_union3_cellline.csv',
        'rand':     'panel_union_rand3.csv',
        'randb':    'panel_union_randb3.csv'}


def read_panel(p):
    f = PDIR + '/' + p + '/DMR_confirmed_' + p + '.csv'
    if not os.path.exists(f):
        sys.exit('판 파일이 없습니다 — 08·09·10 단계를 먼저 돌리십시오: ' + f)
    rows = []
    with open(f, newline='', encoding='utf-8') as fh:
        rd = csv.DictReader(fh)
        if 'Feature' not in (rd.fieldnames or []):
            sys.exit('Feature 열이 없습니다: ' + f)
        for r in rd:
            v = str(r['Feature'])
            if '_' not in v:
                sys.exit('Feature 가 <chr>_<blk> 형식이 아닙니다: %s (%s)' % (v, f))
            c, b = v.rsplit('_', 1)
            try:
                rows.append((c, int(b)))
            except ValueError:
                sys.exit('blk 가 수가 아닙니다: %s (%s)' % (v, f))
    if not rows:
        sys.exit('판이 비었습니다: ' + f)
    return rows


os.makedirs(OUT, exist_ok=True)
for g in ('cellline', 'rand', 'randb'):
    ps = GROUPS[g]
    if not ps:
        sys.exit('%s 에 해당하는 판이 없습니다 — config.conf 의 PANELS·NULL_PANELS 를 보십시오' % g)
    seen, uni = set(), []
    for p in ps:
        for ch, bk in read_panel(p):
            if (ch, bk) in seen: continue
            seen.add((ch, bk)); uni.append((ch, bk))
    uni.sort(key=lambda t: (t[0], t[1]))
    f = OUT + '/' + NAME[g]
    # csv 기본 종결자는 CRLF 다. 기존 파일은 LF 라 그냥 두면 전부 다르게 보인다.
    with open(f, 'w', newline='', encoding='utf-8') as fh:
        w = csv.writer(fh, lineterminator=NL)
        w.writerow(['chr', 'blk'])
        for ch, bk in uni:
            w.writerow([ch, bk])
    print('  %-24s %4d칸  <- %s' % (NAME[g], len(uni), ' '.join(ps)))

# ── 풀에 넣을 검체 이름 ──────────────────────────────────────────
sc = CAND + '/samples.csv'
if not os.path.exists(sc):
    sys.exit('samples.csv 가 없습니다 — 01 단계를 먼저 돌리십시오: ' + sc)
norm, canc = [], []
with open(sc, newline='', encoding='utf-8') as fh:
    rd = csv.DictReader(fh)
    for c in ('sample', 'type'):
        if c not in (rd.fieldnames or []):
            sys.exit('samples.csv 에 %s 열이 없습니다: %s' % (c, sc))
    _odd = set()
    for r in rd:
        t = r['type']
        if t == 'Normal':   norm.append(r['sample'])
        elif t == _CANCER:  canc.append(r['sample'])
        else:               _odd.add(t)
    # 'Normal' 이 아닌 것을 전부 암으로 세면 CANCER_LABEL 이 어긋나도 조용히 돈다.
    if _odd:
        sys.exit('samples.csv 에 모르는 type 이 있습니다: %s (아는 것은 Normal · %s). '
                 'CANCER_LABEL 이 1단계가 쓴 값과 맞는지 보십시오: %s'
                 % (', '.join(sorted(_odd)), _CANCER, sc))
if not norm: sys.exit('samples.csv 에 정상이 없습니다: ' + sc)
if not canc: sys.exit('samples.csv 에 암이 없습니다 (type=%s): %s' % (_CANCER, sc))

for nm, lst in (('use%s_normal.txt' % GEN, norm), ('use%s_gbm.txt' % GEN, canc)):
    f = OUT + '/' + nm
    with open(f, 'w', encoding='utf-8') as fh:
        fh.write(NL.join(lst) + NL)
    print('  %-24s %4d명' % (nm, len(lst)))

print('  저장: ' + OUT)
