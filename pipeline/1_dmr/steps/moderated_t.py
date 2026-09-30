#!/usr/bin/env python
# coding: utf-8
"""moderated t-test (limma eBayes) + BH-FDR. R 의존성 없음.

블록마다 분산을 그 블록 표본만으로 추정하면 표본이 적을 때 불안정하다.
우연히 분산이 작게 나온 블록이 가짜 유의가 된다. 전체 블록의 분산 분포에서
사전분산 (d0, s0^2) 을 추정해 개별 추정치를 그쪽으로 수축시킨다 (Smyth 2004).

    s~^2 = (d0*s0^2 + d*s^2) / (d0 + d)
    t    = diff / sqrt(s~^2 * (1/n1 + 1/n2)),   df = d + d0

사용: python moderated_t.py --selftest
      python moderated_t.py <candidate_beta.csv> [--q 0.05] [--out DMR.csv]
"""
import os, sys, numpy as np, pandas as pd
from scipy.special import psi, polygamma
from scipy import stats


def _trigamma_inv(x):
    x = np.asarray(x, dtype=float)
    y = 0.5 + 1.0 / x
    for _ in range(50):
        tri = polygamma(1, y)
        dif = tri * (1 - tri / x) / polygamma(2, y)
        y = y + dif
        if np.max(np.abs(dif / y)) < 1e-8:
            break
    return y


def fit_prior(s2, df):
    ok = np.isfinite(s2) & (s2 > 0) & (df > 0)
    z = np.log(s2[ok])
    e = z - psi(df[ok] / 2.0) + np.log(df[ok] / 2.0)
    evar = e.var(ddof=1) - np.mean(polygamma(1, df[ok] / 2.0))
    if evar > 0:
        d0 = float(2 * _trigamma_inv(evar))
        s0_2 = float(np.exp(e.mean() + psi(d0 / 2.0) - np.log(d0 / 2.0)))
    else:
        d0, s0_2 = 1e6, float(np.exp(e.mean()))
    return min(d0, 1e6), s0_2


def moderated_ttest(G, N):
    n1 = np.sum(np.isfinite(G), axis=0).astype(float)
    n2 = np.sum(np.isfinite(N), axis=0).astype(float)
    with np.errstate(invalid='ignore', divide='ignore'):
        m1, m2 = np.nanmean(G, axis=0), np.nanmean(N, axis=0)
        v1, v2 = np.nanvar(G, axis=0, ddof=1), np.nanvar(N, axis=0, ddof=1)
        df = n1 + n2 - 2
        s2 = ((n1 - 1) * v1 + (n2 - 1) * v2) / df
    d0, s0_2 = fit_prior(s2, df)
    s2p = (d0 * s0_2 + df * s2) / (d0 + df)
    with np.errstate(invalid='ignore', divide='ignore'):
        se = np.sqrt(s2p * (1.0 / n1 + 1.0 / n2))
        t = (m1 - m2) / se
        dft = df + d0
        p = 2 * stats.t.sf(np.abs(t), dft)
    return m1 - m2, t, p, dft, d0, s0_2


def bh(p):
    p = np.asarray(p, dtype=float)
    q = np.full_like(p, np.nan)
    ok = np.flatnonzero(np.isfinite(p))
    if not len(ok):
        return q
    o = ok[np.argsort(p[ok])]
    adj = p[o] * len(o) / np.arange(1, len(o) + 1)
    q[o] = np.minimum.accumulate(adj[::-1])[::-1].clip(max=1.0)
    return q


def efp(q, sel):
    """BH-FDR 로 뽑은 sel 중 기대 위양성 수. q 합이 표준 추정치이고,
       0.05 x 개수 는 상한이다 (표 3단계의 '1,000개 중 약 50개'가 이 상한)."""
    import numpy as np
    n = len(sel)
    if n == 0: return 0.0, 0.0
    return float(np.nansum(q[sel])), 0.05 * n

def selftest():
    rng = np.random.RandomState(0)
    nb = 3000
    N = rng.normal(0.5, 0.10, (40, nb))
    G = rng.normal(0.5, 0.10, (30, nb))
    G[:, 0] += 0.30
    N[:, 1] = 0.5 + rng.normal(0, 1e-4, 40)
    G[:, 1] = 0.5 + rng.normal(0, 1e-4, 30) + 0.002
    _, _, pm, _, d0, _ = moderated_ttest(G, N)
    _, po = stats.ttest_ind(G, N, axis=0, equal_var=True)
    assert np.isfinite(d0) and d0 > 0, 'd0 추정 실패'
    assert pm[0] < 1e-5, '진짜 차이를 못 잡음 %g' % pm[0]
    assert pm[1] > po[1] * 100, '분산 과소추정 블록을 못 눌렀음 %g vs %g' % (pm[1], po[1])
    q = bh(pm)
    assert np.nanmax(q) <= 1.0 and q[0] < 0.05, 'BH 이상'
    print('selftest OK  d0=%.1f  진짜 p=%.2e  가짜: 보통 %.2e -> moderated %.2e'
          % (d0, pm[0], po[1], pm[1]))


def main():
    if '--selftest' in sys.argv:
        return selftest()
    if len(sys.argv) < 2:
        sys.exit(__doc__)
    src = sys.argv[1]
    Q = float(sys.argv[sys.argv.index('--q') + 1]) if '--q' in sys.argv else 0.05
    X = pd.read_csv(src, index_col=0)
    y = X.pop('Type').values
    A = X.values.astype(float)
    diff, t, p, dft, d0, s0_2 = moderated_ttest(A[y == 1], A[y == 0])
    q = bh(p)
    sel = np.isfinite(q) & (q < Q)
    d = pd.DataFrame({'Feature': X.columns, 'diff': diff, 't': t,
                      'p_value': p, 'q_value': q})[sel].sort_values('q_value')
    n = A.shape[1]
    print('후보 %d · d0=%.1f · s0^2=%.2e' % (n, d0, s0_2))
    print('BH-FDR q<%.3g  통과 %d개 (%.1f%%)' % (Q, len(d), 100.0 * len(d) / n))
    print('참고 Bonferroni p<%.2e  %d개' % (0.05 / n, int(np.nansum(p < 0.05 / n))))
    out = sys.argv[sys.argv.index('--out') + 1] if '--out' in sys.argv else \
        os.path.join(os.path.dirname(src), 'DMR_modt_' + os.path.basename(src))
    d.to_csv(out, index=False)
    print('저장 %s' % out)


if __name__ == '__main__':
    sys.exit(main())
