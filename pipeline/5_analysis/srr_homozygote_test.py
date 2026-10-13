#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
<정상A>/607 이 같은 사람인지 — 유전형에 특이한 통계로 다시.

왜 상관으로는 안 되나
  앞 검사에서 상위 다섯 쌍이 전부 사내 4명(<사내1>·<사내2>·<사내3>·<사내4>)이었다.
  이들은 「서로 다른 사람」이고 같은 실험실·같은 프로토콜이다. 즉 베타 상관은
  배치가 지배한다. 배치가 지배하는 축에서 「상관이 높다」는 같은 사람의
  증거가 되지 못한다.

더 특이한 통계
  배치 효과는 베타를 「매끄럽게」 민다. 유전형은 0 / 0.5 / 1 세 덩어리를 만든다.
  그래서 「한쪽은 0 쪽, 다른 쪽은 1 쪽」인 반대 동형접합 자리를 센다.
  같은 사람이면 이 비율이 측정오차 수준이어야 하고, 다른 사람이면 자리가 쌓인다.
  배치가 같아도 사람이 다르면 쌓인다 — 사내 4명이 그 음성 대조다.

읽는 법
  사내 4명 6쌍의 반대동형 비율이 「같은 배치·다른 사람」의 바닥이다.
  두 SRR 이 그 바닥보다 뚜렷이 낮으면 같은 사람이 의심된다.
  바닥 안에 있으면 「같은 사람이라는 증거가 없다」까지만 말한다.
  이 자료에 「같은 사람을 두 번 잰」 양성 대조가 없어서 상한은 못 세운다 —
  그 한계를 같이 적는다.
"""
import itertools
import sys

import numpy as np
import pandas as pd

D = '<METH_ROOT>/results/dmr/j_candidates15/'
PAIR = ('<정상A>', '<정상B>')
INHOUSE = ['<사내1>', '<사내2>', '<사내3>', '<사내4>']
MIN_DEPTH = 10            # 동형접합 판정이니 깊이를 올린다
HOM_LO, HOM_HI = 0.10, 0.90
MIN_PAIR_SITES = 200


def main():
    s = pd.read_csv(D + 'samples.csv')
    normals = sorted(s.loc[s.type == 'Normal', 'sample'])
    beta = pd.read_parquet(D + 'beta_matrix.parquet')[normals]
    dep = pd.read_parquet(D + 'depth_matrix.parquet')[normals]

    cov = dep >= MIN_DEPTH
    bm = beta.where(cov)
    # 자리 후보: 덮인 검체 중 누군가 0 쪽, 누군가 1 쪽 — 유전형이 갈릴 자리
    ncov = cov.sum(axis=1)
    cand = (ncov >= 8) & (bm.min(axis=1) < HOM_LO) & (bm.max(axis=1) > HOM_HI)
    g, gc = bm.loc[cand], cov.loc[cand]
    print('깊이 %d · 8명 이상 덮임 · 양끝이 갈리는 자리: %d' % (MIN_DEPTH, len(g)))
    if len(g) < 200:
        raise SystemExit('자리가 적어 판정하지 않는다: %d' % len(g))

    rows = []
    for x, y in itertools.combinations(normals, 2):
        both = gc[x] & gc[y]
        n = int(both.sum())
        if n < MIN_PAIR_SITES:
            continue
        bx, by = g.loc[both, x].values, g.loc[both, y].values
        opp = int(np.sum(((bx < HOM_LO) & (by > HOM_HI)) |
                         ((bx > HOM_HI) & (by < HOM_LO))))
        rows.append((opp / n, x, y, n, opp))
    rows.sort()

    tgt = [t for t in rows if {t[1], t[2]} == set(PAIR)]
    if not tgt:
        raise SystemExit('두 SRR 쌍이 자리 부족으로 빠졌다 — 판정하지 않는다.')
    tgt = tgt[0]
    inh = [t for t in rows if t[1] in INHOUSE and t[2] in INHOUSE]
    rank = [i for i, t in enumerate(rows, 1) if {t[1], t[2]} == set(PAIR)][0]

    print()
    print('=' * 72)
    print('반대 동형접합 비율 — 낮을수록 같은 사람 쪽')
    print('=' * 72)
    print('%-32s %9s %8s %8s' % ('쌍', '비율', '자리', '반대'))
    print('%-32s %8.4f%% %8d %8d   <<< 두 SRR'
          % ('%s / %s' % PAIR, tgt[0] * 100, tgt[3], tgt[4]))
    print()
    print('음성 대조 — 사내 4명 (같은 배치 · 서로 다른 사람) %d쌍' % len(inh))
    for p, x, y, n, o in inh:
        print('%-32s %8.4f%% %8d %8d' % ('%s / %s' % (x, y), p * 100, n, o))
    if inh:
        floor = min(p for p, *_ in inh)
        print('  바닥(가장 낮은 사내 쌍)  %.4f%%' % (floor * 100))
    print()
    print('전체 %d쌍 중 낮은 쪽 다섯' % len(rows))
    for p, x, y, n, o in rows[:5]:
        mark = '  <<< 두 SRR' if {x, y} == set(PAIR) else ''
        print('  %8.4f%%  %-14s %-14s n=%-6d%s' % (p * 100, x, y, n, mark))
    print()
    print('두 SRR 순위: %d쌍 중 %d위 (낮은 쪽부터)' % (len(rows), rank))

    print()
    if inh and tgt[0] < floor * 0.5:
        print('판정: 사내 바닥의 절반 아래다. 같은 사람으로 의심된다 —')
        print('      그렇다면 held-out 분할이 그 한 사람에 대해 사람 단위가 아니다.')
        return 2
    print('판정: 같은 사람이라는 증거가 없다.')
    print('      이 자료에 「같은 사람을 두 번 잰」 양성 대조가 없어서')
    print('      「같은 사람이 아님」을 증명한 것은 아니다. 거기까지가 한계다.')
    return 0


if __name__ == '__main__':
    sys.exit(main())
