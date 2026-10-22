"""두 대비에서 «둘 다 통과한» 블록으로 패널을 만든다 (다음 판 2안-b).

  DMR1  암 대 정상 조직 (다른 트리의 1단계 산출 · 정상 = 정상 폐 조직)
  DMR2  암 대 정상 cfDNA (이 판의 1안 트리 1단계 산출)
  bl200   두 트리의 07 순위 목록(DMR_ranked_m.csv · q<0.05 · |dbeta|>=30)에 모두 있는 칸 → DMR2 순위 상위 N
  jsd200  두 트리 모두 JSD 가 귀무 PCT 백분위를 넘는 칸(06b 와 같은 정의 · jsd_obs > t) → DMR2 jsd_obs 순서 상위 N
  dbeta · q 는 06 과 같은 방법(세포주 단위 평균 · moderated t · BH)으로 각 트리에서 낸다. 방향(dbeta 부호)은 두 대비를 맞대어 센다.
  MIN(기본 50)칸 미만인 판은 «계산 불가» 다. 남은 판 목록을 OUT/_intersect_panels.txt 에 적는다.

  쓰는 법:
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

a = sys.argv
for k in ('--a', '--b', '--out'):
    if k not in a:
        raise SystemExit('사용: 12b_passset_panels.py --a <DMR1 results/dmr> --b <DMR2 results/dmr> --out <출력 판 폴더> [--n 200] [--pct 99] [--min 50] [--report-only]')
RA = a[a.index('--a') + 1].rstrip('/')
RB = a[a.index('--b') + 1].rstrip('/')
OUT = a[a.index('--out') + 1].rstrip('/')
N = int(a[a.index('--n') + 1]) if '--n' in a else 200
PCT = float(a[a.index('--pct') + 1]) if '--pct' in a else 99.0
MIN = int(a[a.index('--min') + 1]) if '--min' in a else 50
ONLY = '--report-only' in a
GEN = str(_cfg.GEN)
K, MINR = int(_cfg.PATTERN_K), int(_cfg.MIN_READS_PER_WINDOW)
if os.path.realpath(OUT).startswith(os.path.realpath(RA) + os.sep) or os.path.realpath(OUT).startswith(os.path.realpath(RB) + os.sep):
    raise SystemExit('출력이 입력 트리 안이다 — 입력을 덮을 수 있다: %s' % OUT)
_CANCER = getattr(_cfg, 'CANCER_LABEL', 'GBM')


def stats(R):
    """06 과 같은 dbeta · q (세포주 단위 평균 · moderated t · BH). 칸 이름 → (dbeta, q)."""
    C = '%s/j_candidates%s' % (R, GEN)
    B = pd.read_parquet(C + '/beta_matrix.parquet')
    M = pd.read_csv(C + '/samples.csv')
    X = B.T.reindex(M['sample']).values.astype(np.float32)
    y = (M['type'] == _CANCER).astype(int).values
    if int(y.sum()) == 0:
        raise SystemExit('%s 에 type=%s 인 검체가 없다' % (C, _CANCER))
    cols = [str(i) for i in range(X.shape[1])]
    D = pd.DataFrame(X, columns=cols); D['_g'] = M['group'].values; D['_y'] = y
    G = D.groupby('_g', sort=True).mean()
    y = (G['_y'].values > 0.5).astype(int)
    X = G[cols].values.astype(np.float32)
    mu = np.nanmean(X, axis=0); mu = np.where(np.isfinite(mu), mu, 0.0)
    X = np.where(np.isfinite(X), X, mu)
    d, t, p, dft, d0, s0 = moderated_ttest(X[y == 1], X[y == 0])
    q = bh(p)
    f = B.index.get_level_values(0).astype(str) + '_' + B.index.get_level_values(1).astype('int64').astype(str)
    return pd.DataFrame({'dbeta': d, 'q': q}, index=f)


def ranked(R):
    d = pd.read_csv('%s/j_panel_dmr%s/bl200_d10_cellline/DMR_ranked_m.csv' % (R, GEN))
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
sel = {'bl200': [f for f in lb if f in set(la)],
       'jsd200': [f for f in jb if f in set(ja)]}
note = {'bl200': 'DMR1 통과 %d · DMR2 통과 %d · 둘 다 %d' % (len(la), len(lb), len(sel['bl200'])),
        'jsd200': 'DMR1 통과 %d(t %.4f) · DMR2 통과 %d(t %.4f) · 둘 다 %d' % (len(ja), ta, len(jb), tb, len(sel['jsd200']))}
keep, report = [], []
for p in _cfg.PANELS:
    if p not in sel:
        raise SystemExit('이 단계는 bl200 · jsd200 만 만든다: PANELS 에 %s' % p)
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
    print('보고만 — 아무것도 쓰지 않았다')
    sys.exit(0)
os.makedirs(OUT, exist_ok=True)
with open(OUT + '/_intersect_panels.txt', 'w', encoding='utf-8') as h:
    h.write(','.join(keep) + '\n')
    h.write('\n'.join(report) + '\n')
if not keep:
    raise SystemExit('모든 판이 «계산 불가» 다(%d칸 미만)' % MIN)
print('쓸 판: %s  (기록 %s/_intersect_panels.txt)' % (','.join(keep), OUT))
