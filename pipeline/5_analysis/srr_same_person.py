#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
학습 정상 15명 중 <정상A> 과 <정상B> 이 «같은 사람»인지 잰다.

왜 재나
  held-out 분할은 사람 단위여야 뜻이 있다. 두 SRR 은 번호가 연속이라
  같은 사람의 두 런일 수 있고, 그렇다면 하나는 참조군(607) 하나는 혼합군(606)
  으로 갈려 있으므로 그 사람에 대해서는 누수가 «안 고쳐진다».

방법 (EGAD 공여자 가르기와 같은 원리)
  C>T SNP 자리는 메틸화 0% 처럼 보인다. 그래서 「어떤 사람은 0 에 가깝고
  어떤 사람은 안 그런」 CpG 는 유전형 자리일 가능성이 높다. 그 자리들만 골라
  쌍별 베타 상관을 낸다. 같은 사람의 두 런이면 상관이 다른 105쌍 중
  뚜렷한 바깓값이어야 한다.

이 판정은 「상관이 높다」가 아니라 «다른 모든 쌍과 비교해 바깓값인가» 로 한다.
바깓값이 아니면 「같은 사람이라는 증거가 없다」로 적는다 — 「다른 사람이다」가
아니다. 두 문장은 다르다.
"""
import itertools
import sys

import numpy as np
import pandas as pd

NL = chr(10)
D = '<METH_ROOT>/results/dmr/j_candidates15/'
PAIR = ('<정상A>', '<정상B>')
MIN_DEPTH = 5
MIN_COVERED = 8          # 자리를 「유전형처럼 갈린다」고 부를 최소 덮인 검체 수
MIN_PAIR_SITES = 200     # 이보다 적으면 그 쌍은 상관을 내지 않는다
LOW, HIGH = 0.10, 0.60


def main():
    s = pd.read_csv(D + 'samples.csv')
    normals = sorted(s.loc[s.type == 'Normal', 'sample'])
    if len(normals) != 15:
        raise SystemExit('정상이 15명이 아니다: %d' % len(normals))

    beta = pd.read_parquet(D + 'beta_matrix.parquet')[normals]
    dep = pd.read_parquet(D + 'depth_matrix.parquet')[normals]

    # RRBS 커버리지가 검체마다 크게 달라서 「15명 전원」을 걸면 자리가 32개로
    # 말라붙는다. 그래서 자리 정의는 «덮인 검체만» 보고, 쌍별 상관은 그 쌍이
    # 둘 다 덮은 자리에서 낸다. n 을 같이 찍어 쌍마다 다른 것을 숨기지 않는다.
    cov = dep >= MIN_DEPTH
    bm = beta.where(cov)
    ncov = cov.sum(axis=1)
    lo = bm.min(axis=1)
    hi = bm.max(axis=1)
    snpish = (ncov >= MIN_COVERED) & (lo < LOW) & (hi > HIGH)
    g, gc = bm.loc[snpish], cov.loc[snpish]
    print('깊이 %d 기준 · %d명 이상 덮인 자리 중 유전형처럼 갈리는 것: %d / %d'
          % (MIN_DEPTH, MIN_COVERED, len(g), len(beta)))
    if len(g) < 200:
        raise SystemExit('자리가 너무 적어 판정하지 않는다: %d' % len(g))

    rows = []
    for x, y in itertools.combinations(normals, 2):
        both = gc[x] & gc[y]
        n = int(both.sum())
        if n < MIN_PAIR_SITES:
            rows.append((float('nan'), x, y, n))
            continue
        r = float(np.corrcoef(g.loc[both, x].values, g.loc[both, y].values)[0, 1])
        rows.append((r, x, y, n))
    thin = [t for t in rows if np.isnan(t[0])]
    rows = [t for t in rows if not np.isnan(t[0])]
    rows.sort(reverse=True)
    print('자리가 %d개 미만이라 뺀 쌍: %d' % (MIN_PAIR_SITES, len(thin)))
    for _r, x, y, n in thin:
        print('   %-14s %-14s n=%d' % (x, y, n))

    hit = [t for t in rows if {t[1], t[2]} == set(PAIR)]
    if not hit:
        raise SystemExit('두 SRR 쌍이 자리 부족으로 빠졌다 — 판정하지 않는다.')
    target = hit[0]
    others = [t[0] for t in rows if {t[1], t[2]} != set(PAIR)]
    mu, sd = float(np.mean(others)), float(np.std(others, ddof=1))
    rank = [i for i, t in enumerate(rows, 1) if {t[1], t[2]} == set(PAIR)][0]

    print()
    print('=' * 68)
    print('%s 대 %s' % PAIR)
    print('=' * 68)
    print('  상관            %.4f  (자리 %d개)' % (target[0], target[3]))
    print('  나머지 %d쌍    중앙 %.4f · 평균 %.4f · 표준편차 %.4f · 폭 %.4f~%.4f'
          % (len(others), float(np.median(others)), mu, sd,
             min(others), max(others)))
    print('  z               %+.2f' % ((target[0] - mu) / sd))
    print('  순위            %d쌍 중 %d위' % (len(rows), rank))
    print()
    print('가장 높은 다섯 쌍')
    for r, x, y, n in rows[:5]:
        mark = '  <<< 두 SRR' if {x, y} == set(PAIR) else ''
        print('  %.4f  %-14s %-14s n=%-6d%s' % (r, x, y, n, mark))

    print()
    if rank == 1 and (target[0] - mu) / sd > 4:
        print('판정: 같은 사람으로 보인다. 분할이 사람 단위가 아니다 — 고쳐야 한다.')
        return 2
    if rank <= 3 and (target[0] - mu) / sd > 3:
        print('판정: 의심스럽다. 추가 근거 없이 넘기지 않는다.')
        return 2
    print('판정: 바깓값이 아니다 — 같은 사람이라는 증거가 없다.')
    print('      (「다른 사람임을 증명했다」가 아니다. 이 자료로는 거기까지다.)')
    return 0


if __name__ == '__main__':
    sys.exit(main())
