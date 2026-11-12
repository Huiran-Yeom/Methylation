"""두 대비에서 모두 통과한 블록으로 패널을 만듦
  DMR1: 암 대 정상 조직 (다른 트리의 1단계 산출), DMR2: 암 대 정상 cfDNA (이 트리)
  bl200   두 트리 07 순위 목록(DMR_ranked_m.csv)에 다 있는 칸 → DMR2 순위로 상위 N
  jsd200  두 트리 다 JSD 가 귀무 PCT 백분위 넘는 칸 (06b 와 같은 기준) → DMR2 jsd_obs 순으로 상위 N
  dbeta, q 는 06 과 같은 방법으로 트리마다 계산, 부호 같은 칸 셈
  MIN 칸 미만인 판은 계산 불가로 빼고, 남은 판 목록을 OUT/_intersect_panels.txt 에 씀

  사용
    python 12b_passset_panels.py --a <DMR1 results/dmr> --b <DMR2 results/dmr> --out <이 트리 j_panel_dmr<GEN>_panels>
                                 [--n 200] [--pct 99] [--min 50] [--report-only]
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
import numpy as np
import pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from moderated_t import moderated_ttest, bh
import dmr_common   # noqa: E402

a = sys.argv
for k in ('--a', '--b', '--out'):
    if k not in a:
        raise SystemExit('사용: 12b_passset_panels.py --a <DMR1 results/dmr> --b <DMR2 results/dmr> --out <출력 판 폴더> [--n 200] [--pct 99] [--min 50] [--report-only]')
RA = a[a.index('--a') + 1].rstrip('/')
RB = a[a.index('--b') + 1].rstrip('/')
OUT = a[a.index('--out') + 1].rstrip('/')
PN = _cfg.PANEL_N    # 판 이름, 07 폴더 이름의 크기 (config.conf PANEL_N)
N = int(a[a.index('--n') + 1]) if '--n' in a else PN
PCT = float(a[a.index('--pct') + 1]) if '--pct' in a else 99.0
MIN = int(a[a.index('--min') + 1]) if '--min' in a else 50
ONLY = '--report-only' in a
GEN = str(_cfg.GEN)
K, MINR = int(_cfg.PATTERN_K), int(_cfg.MIN_READS_PER_WINDOW)
if os.path.realpath(OUT).startswith(os.path.realpath(RA) + os.sep) or os.path.realpath(OUT).startswith(os.path.realpath(RB) + os.sep):
    raise SystemExit('출력이 입력 트리 안이다. 입력을 덮을 수 있다: %s' % OUT)
_CANCER = getattr(_cfg, 'CANCER_LABEL', 'GBM')


def stats(R):
    """06 과 같은 dbeta, q (세포주 평균, moderated t, BH). 칸 이름별 표로 돌려줌"""
    C = '%s/j_candidates%s' % (R, GEN)
    B = pd.read_parquet(C + '/beta_matrix.parquet')
    M = pd.read_csv(C + '/samples.csv')
    X, y, _ = dmr_common.unit_matrix(B, M, _CANCER, '%s 에 type=%s 인 검체가 없다' % (C, _CANCER))
    X = dmr_common.fill_col_mean(X)
    d, t, p, dft, d0, s0 = moderated_ttest(X[y == 1], X[y == 0])
    q = bh(p)
    f = B.index.get_level_values(0).astype(str) + '_' + B.index.get_level_values(1).astype('int64').astype(str)
    return pd.DataFrame({'dbeta': d, 'q': q}, index=f)


def ranked(R):
    d = pd.read_csv('%s/j_panel_dmr%s/bl%d_d10_cellline/DMR_ranked_m.csv' % (R, GEN, PN))
    return list(d['chr'].astype(str) + '_' + d['blk'].astype(str))


def jsd_pass(R):
    boot = pd.read_parquet('%s/j_jsd%s/boot_K%d_minr%d.parquet' % (R, GEN, K, MINR))
    null = pd.read_parquet('%s/j_jsd%s/null_K%d_minr%d.parquet' % (R, GEN, K, MINR))['null_block_jsd'].values
    t = float(np.percentile(null, PCT))
    boot['F'] = boot['chr'].astype(str) + '_' + boot['blk'].astype('int64').astype(str)
    z = boot[boot['jsd_obs'] > t].sort_values('jsd_obs', ascending=False)
    return list(z['F']), t


print('DMR1 %s\nDMR2 %s\n출력 %s · 판당 최대 %d칸 · jsd 귀무 %.1f%% · 최소 %d칸%s'
      % (RA, RB, OUT, N, PCT, MIN, ' · 보고만' if ONLY else ''))
SA, SB = stats(RA), stats(RB)
la, lb = ranked(RA), ranked(RB)
ja, ta = jsd_pass(RA)
jb, tb = jsd_pass(RB)
BL, JS = 'bl%d' % PN, 'jsd%d' % PN
sel = {BL: [f for f in lb if f in set(la)],
       JS: [f for f in jb if f in set(ja)]}
note = {BL: 'DMR1 통과 %d · DMR2 통과 %d · 둘 다 %d' % (len(la), len(lb), len(sel[BL])),
        JS: 'DMR1 통과 %d(t %.4f) · DMR2 통과 %d(t %.4f) · 둘 다 %d' % (len(ja), ta, len(jb), tb, len(sel[JS]))}
keep, report = [], []
for p in _cfg.PANELS:
    if p not in sel:
        raise SystemExit('이 단계는 %s · %s 만 만든다: PANELS 에 %s' % (BL, JS, p))
    F = sel[p][:N]
    da, db = SA.reindex(F)['dbeta'], SB.reindex(F)['dbeta']
    same = int(((da > 0) & (db > 0)).sum() + ((da < 0) & (db < 0)).sum())
    undec = int((da.isna() | db.isna() | (da == 0) | (db == 0)).sum())
    ok = len(F) >= MIN
    report.append('%s  %s → 쓰는 칸 %d (방향 같음 %d · 다름 %d · 판단 불가 %d)%s'
                  % (p, note[p], len(F), same, len(F) - same - undec, undec, '' if ok else '  → «계산 불가»(%d칸 미만)' % MIN))
    if not ok or ONLY:
        continue
    d = pd.DataFrame({'Feature': F, 'dbeta': SB.reindex(F)['dbeta'].round(2).values, 'q': SB.reindex(F)['q'].values})
    assert d['Feature'].is_unique
    o = '%s/%s' % (OUT, p)
    os.makedirs(o, exist_ok=True)
    d.to_csv('%s/DMR_confirmed_%s.csv' % (o, p), index=False)
    with open(o + '/_origin.txt', 'w', encoding='utf-8') as h:
        h.write('%s 둘 다 통과 패널 (12b_passset_panels.py)\n  DMR1 %s\n  DMR2 %s\n  %s · 쓰는 칸 %d · 행은 DMR2 의 순서 · dbeta · q 는 DMR2\n'
                % (p, RA, RB, note[p], len(F)))
    keep.append(p)
for r in report:
    print('  ' + r)
if ONLY:
    print('보고만. 아무것도 쓰지 않았다')
    sys.exit(0)
os.makedirs(OUT, exist_ok=True)
with open(OUT + '/_intersect_panels.txt', 'w', encoding='utf-8') as h:
    h.write(','.join(keep) + '\n')
    h.write('\n'.join(report) + '\n')
if not keep:
    raise SystemExit('모든 판이 «계산 불가» 다(%d칸 미만)' % MIN)
print('쓸 판: %s  (기록 %s/_intersect_panels.txt)' % (','.join(keep), OUT))
