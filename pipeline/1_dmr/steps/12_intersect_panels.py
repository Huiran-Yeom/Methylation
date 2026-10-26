"""두 대비의 패널을 판마다 교집합한다 (다음 판 2안).

  DMR1  암 대 정상 조직 (다른 트리의 1단계 08 산출)
  DMR2  암 대 정상 cfDNA (이 트리의 1단계 08 산출과 같은 것)
  판마다 같은 이름끼리 교집합한다: bl200 ∩ bl200 · jsd200 ∩ jsd200 · …(설정 PANELS)
  행은 DMR2 의 순서 · 값(dbeta · q)을 그대로 둔다.

  방향: dbeta 는 두 대비 모두 «암 − 정상» 이다. 판마다 부호가 같은 칸 · 다른 칸을 센다.
    --same-direction 이면 부호가 같은 칸만 남긴다(종양이 혈액과도 · 정상 조직과도 같은 쪽으로 다른 칸).
    --report-only 면 수만 찍고 아무것도 쓰지 않는다(패널 목록만 보는 것이라 성능 결과와 무관).
  교집합이 MIN(기본 50)칸 미만인 판은 «계산 불가» 다. 그 판은 쓰지 않고, 남은 판 목록을
  OUT/_intersect_panels.txt 에 적는다. 남은 판이 없으면 멈춘다.
  쓰는 법:
    python 12_intersect_panels.py --a <DMR1 의 j_panel_dmr<GEN>_panels> --b <DMR2 의 같은 폴더> --out <이 트리의 같은 폴더> [--min 50]
  무작위 패널(09 · 10)은 이 단계 뒤에 돌린다. 09 는 jsd200 칸 수, 10 은 bl200 칸 수에 맞춰 뽑는다.
"""
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
except ImportError:
    raise SystemExit('설정을 못 찾았습니다. 0_setup/meth_config.py 가 있는 자리를 찾지 못했습니다.')

import os
import sys
import pandas as pd

a = sys.argv
for k in ('--a', '--b', '--out'):
    if k not in a:
        raise SystemExit('사용: 12_intersect_panels.py --a <DMR1 판 폴더> --b <DMR2 판 폴더> --out <출력 판 폴더> [--min 50]')
A = a[a.index('--a') + 1].rstrip('/')
B = a[a.index('--b') + 1].rstrip('/')
OUT = a[a.index('--out') + 1].rstrip('/')
MIN = int(a[a.index('--min') + 1]) if '--min' in a else 50
SAME = '--same-direction' in a
ONLY = '--report-only' in a
if os.path.realpath(OUT) in (os.path.realpath(A), os.path.realpath(B)):
    raise SystemExit('출력이 입력 폴더와 같다 — 입력 패널을 덮는다: %s' % OUT)

print('DMR1 %s\nDMR2 %s\n출력 %s · 최소 %d칸 · 방향 %s%s'
      % (A, B, OUT, MIN, '같은 쪽만' if SAME else '안 봄', ' · 보고만' if ONLY else ''))
keep, report = [], []
for p in _cfg.PANELS:
    fa = '%s/%s/DMR_confirmed_%s.csv' % (A, p, p)
    fb = '%s/%s/DMR_confirmed_%s.csv' % (B, p, p)
    for f in (fa, fb):
        if not os.path.isfile(f):
            raise SystemExit('입력이 없다: %s — 그 트리의 1단계 08 을 먼저 돌리십시오' % f)
    da, db = pd.read_csv(fa), pd.read_csv(fb)
    for n, d in (('DMR1', da), ('DMR2', db)):
        if not d['Feature'].is_unique:
            raise SystemExit('%s %s 에 같은 칸이 두 번 있다' % (n, p))
    sa = dict(zip(da['Feature'], da['dbeta']))
    d = db[db['Feature'].isin(sa)].reset_index(drop=True)
    same = [(sa[f] > 0 and b > 0) or (sa[f] < 0 and b < 0) for f, b in zip(d['Feature'], d['dbeta'])]
    # NaN · 0 은 방향을 판단할 수 없다 — «다름» 과 따로 센다(--same-direction 이면 빠진다)
    undec = int(sum(not (abs(sa[f] * b) > 0) for f, b in zip(d["Feature"], d["dbeta"])))
    n_all, n_same = len(d), int(sum(same))
    if SAME:
        d = d[same].reset_index(drop=True)
    ok = len(d) >= MIN
    report.append('%s  DMR1 %d · DMR2 %d · 교집합 %d (방향 같음 %d · 다름 %d · 판단 불가 %d) · 쓰는 칸 %d%s'
                  % (p, len(da), len(db), n_all, n_same, n_all - n_same - undec, undec, len(d),
                     '' if ok else '  → «계산 불가»(%d칸 미만)' % MIN))
    if not ok or ONLY:
        continue
    o = '%s/%s' % (OUT, p)
    os.makedirs(o, exist_ok=True)
    d.to_csv('%s/DMR_confirmed_%s.csv' % (o, p), index=False)
    with open(o + '/_origin.txt', 'w', encoding='utf-8') as h:
        h.write('%s 교집합 패널 (12_intersect_panels.py)\n  DMR1 %s\n  DMR2 %s\n'
                '  DMR1 %d칸 ∩ DMR2 %d칸 = %d칸 (방향 같음 %d) · 방향 %s · 쓰는 칸 %d · 행은 DMR2 의 순서와 값\n'
                % (p, fa, fb, len(da), len(db), n_all, n_same, '같은 쪽만' if SAME else '안 봄', len(d)))
    keep.append(p)
for r in report:
    print('  ' + r)
if ONLY:
    print('보고만 — 아무것도 쓰지 않았다')
    sys.exit(0)
os.makedirs(OUT, exist_ok=True)
with open(OUT + '/_intersect_panels.txt', 'w', encoding='utf-8') as h:
    h.write(','.join(keep) + '\n')
    h.write('\n'.join(report) + '\n')
if not keep:
    raise SystemExit('모든 판이 «계산 불가» 다(교집합 %d칸 미만)' % MIN)
print('쓸 판: %s  (기록 %s/_intersect_panels.txt)' % (','.join(keep), OUT))
