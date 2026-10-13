#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
검증 AUC 의 「재어진」 우연 기준선. 사전등록 추기 38항.

읽는 것은 검증 「음성」 행의 y_score 뿐이다. 양성 행은 읽는 즉시 버린다.
설계·읽는 규칙·쓰는 규칙은 추기 38항에 돌리기 전에 적었다.

파일 이름 세대
  판마다 두 형식이 있다 — `y_prob_all_<판>.csv` 와 `y_prob_all_<판>_모든조합.csv`.
  bl200 만 기본형이 164개 있고(부분 실행), 세 판 전부 `_모든조합` 이 178개다.
  **`_모든조합` 을 쓴다.** 어느 쪽을 읽었는지 찍는다 — make_tables.py 의
  이름세대 경고와 같은 함정이다.

귀무 두 가지 (38-4)
  (가) 비율 섞음   검체의 음성 70행 중 무작위 10행을 위장 양성
  (나) 비율 고정   검체의 음성 중 「비율 0.005 인 10행」을 위장 양성, 나머지를 음성
                   → 실제 통계와 모양이 같다(양성은 한 비율, 음성은 전 비율)
  (나)를 주로 본다. (가)는 민감도다.
"""
import glob
import sys

import numpy as np
import pandas as pd

VAL = '<METH_ROOT>/patient178_j51_val'
PANELS = ['j15bl200', 'j15jsd200', 'j15jsdb200']
DEPTH = '5k_uni'
COMBO = 'mean|entropy|jsd|llr|pdr'
RATIO = 0.005
NPOS_PER = 10
NREP = 200
SEED = 42
REPORTED = {'j15bl200': 0.5227, 'j15jsd200': 0.5420, 'j15jsdb200': 0.5254}


def auc(y, s):
    r = pd.Series(s).rank().to_numpy()
    p, n = int((y == 1).sum()), int((y == 0).sum())
    if p == 0 or n == 0:
        return np.nan
    return (r[y == 1].sum() - p * (p + 1) / 2.0) / (p * n)


def load_negatives(panel):
    """검체 -> (전체 음성 점수, 비율 0.005 인 음성 마스크). 양성은 즉시 버린다."""
    d = '%s_%s' % (panel, DEPTH)
    base = '%s/*/results/%s/step04_ML_classifier/' % (VAL, d)
    cands = [('모든조합', base + 'y_prob_all_%s_모든조합.csv' % d),
             ('기본', base + 'y_prob_all_%s.csv' % d)]
    used, files = None, []
    for name, pat in cands:
        f = sorted(glob.glob(pat))
        if len(f) > len(files):
            used, files = name, f
    if not files:
        raise SystemExit('점수 파일이 없다: ' + panel)
    per = {}
    for f in files:
        t = pd.read_csv(f, usecols=['sample', 'feature_set', 'ratio',
                                    'y_true', 'y_score'])
        t = t[(t.feature_set == COMBO) & (t.y_true == 0)]   # 주의: 양성 버림
        if not len(t):
            continue
        smp = str(t['sample'].iloc[0])
        per[smp] = (t['y_score'].to_numpy(dtype=float),
                    np.isclose(t['ratio'].to_numpy(dtype=float), RATIO))
    return per, used, len(files)


def run_null(per, samples, rng, fixed):
    fp, fn = [], []
    for s in samples:
        v, at_ratio = per[s]
        if fixed:
            m = at_ratio.copy()
            if m.sum() == 0 or m.sum() >= len(v):
                continue
        else:
            k = min(NPOS_PER, max(1, len(v) // 2))
            idx = rng.choice(len(v), k, replace=False)
            m = np.zeros(len(v), bool)
            m[idx] = True
        fp.append(v[m])
        fn.append(v[~m])
    if not fp:
        return np.nan
    fp, fn = np.concatenate(fp), np.concatenate(fn)
    y = np.r_[np.zeros(len(fn)), np.ones(len(fp))]
    return auc(y, np.r_[fn, fp])


def main():
    print('검증 음성만 읽는다 · 조합 %s · 깊이 %s · 비율 %.3f'
          % (COMBO, DEPTH, RATIO))
    print('양성 행은 읽는 즉시 버린다 (추기 38-1 · 38-5)')
    print()
    data = {}
    for panel in PANELS:
        per, used, nf = load_negatives(panel)
        data[panel] = per
        print('%-12s 파일 형식 「%s」 · %d개 · 검체 %d'
              % (panel, used, nf, len(per)))
    common = sorted(set.intersection(*[set(data[p]) for p in PANELS]))
    print('세 판 공통 검체 %d' % len(common))
    print()

    print('%-12s %9s %9s %9s %9s %9s %-6s'
          % ('판', '음성행', '귀무중앙', '귀무95', '귀무99', '보고값', '판정'))
    print('— (나) 비율 고정 · 주 —')
    for panel in PANELS:
        per = data[panel]
        rng = np.random.default_rng(SEED)
        v = run_null(per, common, rng, True)
        tot = sum(len(per[s][0]) for s in common)
        rep = REPORTED[panel]
        # 비율 고정은 값이 하나다. 반복은 의미가 없다 — 그대로 적는다.
        print('%-12s %9d %9.4f %9s %9s %9.4f %-6s'
              % (panel, tot, v, '—', '—', rep, '밖' if rep > v else '안'))
    print('   비율 고정은 추출이 없어 값이 하나다. 분위수가 없으므로')
    print('   「보고값 > 귀무」로만 읽는다.')
    print()
    print('— (가) 비율 섞음 · 민감도 · 반복 %d회 —' % NREP)
    for panel in PANELS:
        per = data[panel]
        rng = np.random.default_rng(SEED)
        a = np.array([run_null(per, common, rng, False) for _ in range(NREP)])
        a = a[np.isfinite(a)]
        tot = sum(len(per[s][0]) for s in common)
        rep = REPORTED[panel]
        q95, q99 = np.percentile(a, 95), np.percentile(a, 99)
        print('%-12s %9d %9.4f %9.4f %9.4f %9.4f %-6s'
              % (panel, tot, np.median(a), q95, q99, rep,
                 '밖' if rep > q95 else '안'))
    print()
    print('보고값은 4_report 가 낸 것을 그대로 인용했다 (여기서 재계산 안 함).')
    print('판정이 「밖」이어도 순위는 매기지 않는다 — 0.70 미달과 사전등록 금지가')
    print('그대로 걸린다 (추기 38-1 · 38-4).')
    return 0


if __name__ == '__main__':
    sys.exit(main())
