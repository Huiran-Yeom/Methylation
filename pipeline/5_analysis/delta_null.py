#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Δ 의 귀무 분포 — 수정판 귀무와 정본 귀무를 짝지어 뺀다. 추기 42-7. 2026-10-07

왜
  42-2 는 Δ 를 「한 번 뽑은 값의 단측 폭(0.081)」으로 판정했다. 틀렸다.
  Δ 는 두 번 뽑은 것의 차이고 양측 질문이다. σ√2 를 추정할 수도 있지만
  그럴 필요가 없다 — 두 귀무를 반복 i 끼리 빼면 Δ 의 귀무가 «직접» 나온다.

내는 것
  중앙(Δ_귀무)  = n_eff 35→18 의 실측 효과   (42-3 이 재려던 것)
  2.5~97.5 퍼센타일 = Δ 를 판정할 올바른 문턱  (42-7 에 박은 규칙)

읽는 규칙은 추기 42-3 과 42-7 에 숫자 보기 전에 적었다. 여기서 다시 정하지 않는다.
"""
import sys

import numpy as np
import pandas as pd
from sklearn.ensemble import RandomForestClassifier
from sklearn.impute import SimpleImputer
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import StratifiedKFold, cross_val_predict
from sklearn.pipeline import Pipeline
from sklearn.svm import SVC
from xgboost import XGBClassifier

SEED = 42
NREP = 120
NF = 5
HALF = 50
MR = '25'

CANON = '<METH_ROOT>/results/판'
HELD = '<METH_ROOT_HELDOUT>/results/판'
PANELS = ['j15bl200', 'j15jsd200', 'j15jsdb200']
SETS = {
    'entropy': ['entropy'],
    'llr 뺀 조합': ['mean', 'entropy', 'jsd', 'pdr'],
    '1차 조합': ['mean', 'entropy', 'jsd', 'llr', 'pdr'],
}
# 실측 Δ (추기 40-1·40-3). 여기서 다시 계산하지 않고 인용한다.
OBS = {
    ('j15bl200', 'entropy'): None,
    ('j15bl200', 'llr 뺀 조합'): +0.1102,
    ('j15bl200', '1차 조합'): +0.0117,
    ('j15jsd200', 'llr 뺀 조합'): -0.1020,
    ('j15jsd200', '1차 조합'): -0.1439,
    ('j15jsdb200', 'llr 뺀 조합'): +0.0608,
    ('j15jsdb200', '1차 조합'): +0.0074,
}


def models():
    return {
        'LR': LogisticRegression(random_state=SEED),
        'RF': RandomForestClassifier(random_state=SEED, n_jobs=1),
        'XGB': XGBClassifier(random_state=SEED, eval_metric='logloss', n_jobs=1),
        'SVM': SVC(random_state=SEED, probability=False),
    }


def auc(y, s):
    r = pd.Series(s).rank().to_numpy()
    p, n = int((y == 1).sum()), int((y == 0).sum())
    return (r[y == 1].sum() - p * (p + 1) / 2.0) / (p * n)


def cv_max(X, y):
    cv = StratifiedKFold(NF, shuffle=True, random_state=SEED)
    best = -1.0
    for mk in models().values():
        pipe = Pipeline([('imp', SimpleImputer(strategy='mean',
                                               keep_empty_features=True)),
                         ('m', mk)])
        meth = ('predict_proba' if hasattr(mk, 'predict_proba')
                else 'decision_function')
        try:
            o = np.asarray(cross_val_predict(pipe, X, y, method=meth, cv=cv,
                                             n_jobs=1), dtype=float)
        except Exception:
            continue
        a = auc(y, o[:, 1] if o.ndim == 2 else o)
        if np.isfinite(a):
            best = max(best, a)
    return best


def load_neg(base, panel, groups):
    d = panel + '_5k_uni'
    cols = []
    for g in groups:
        f = '%s/%s/%s/step04_ML_classifier/%s_%s_%s.csv' % (base, panel, d, MR, g, d)
        t = pd.read_csv(f, index_col=0)
        if 'Type' not in t.columns:
            raise SystemExit('Type 열이 없다: ' + f)
        t = t.drop(columns=['Type'])
        t = t.loc[t.index.to_series().str.startswith('Normal')]
        t.columns = ['%s::%s' % (g, c) for c in t.columns]
        cols.append(t.sort_index())
    X = cols[0]
    for e in cols[1:]:
        X = X.join(e, how='inner')
    if len(X) != 100:
        raise SystemExit('음성이 100행이 아니다: %d' % len(X))
    return X.to_numpy(dtype=float)


def null_series(X, rng):
    """같은 라벨 순서를 두 트리에 쓰려면 라벨을 밖에서 만들어 넘겨야 한다."""
    out = []
    for y in rng:
        out.append(cv_max(X, y))
    return np.array(out, dtype=float)


def main():
    print('Δ 의 귀무 — 수정판 귀무 − 정본 귀무 · 반복 %d회 · 50/50' % NREP)
    print('라벨과 행 추출을 두 트리에 «같게» 주어 짝을 맞춘다')
    print()
    print('%-12s %-13s %9s %9s %9s %9s %9s %-6s'
          % ('판', '특징', 'Δ귀무중앙', 'Δ귀무2.5%', 'Δ귀무97.5%', '실측Δ', '정본귀무중앙', '판정'))
    for panel in PANELS:
        for name, groups in SETS.items():
            try:
                Xc = load_neg(CANON, panel, groups)
                Xh = load_neg(HELD, panel, groups)
            except Exception as e:
                print('%-12s %-13s 못 읽음 %s' % (panel, name, e))
                continue
            # ★ 라벨을 미리 만들어 두 트리에 같은 순서로 쓴다
            rg = np.random.default_rng(SEED)
            labels = []
            for _ in range(NREP):
                y = np.zeros(100, int)
                y[rg.choice(100, HALF, replace=False)] = 1
                labels.append(y)
            nc = null_series(Xc, labels)
            nh = null_series(Xh, labels)
            d = nh - nc
            d = d[np.isfinite(d)]
            obs = OBS.get((panel, name))
            lo, hi = np.percentile(d, [2.5, 97.5])
            if obs is None:
                verdict = '실측없음'
            else:
                verdict = '밖' if (obs < lo or obs > hi) else '안'
            print('%-12s %-13s %+9.4f %+9.4f %+9.4f %9s %9.4f %-6s'
                  % (panel, name, np.median(d), lo, hi,
                     ('%+.4f' % obs) if obs is not None else '—',
                     np.median(nc), verdict))
    print()
    print('Δ귀무중앙 = n_eff 35→18 의 실측 효과 (추기 42-3)')
    print('2.5~97.5 퍼센타일 = Δ 판정 문턱 (추기 42-7 에 숫자 보기 전에 박았다)')
    print('한 번 뽑은 값의 폭(0.081)은 쓰지 않는다.')
    return 0


if __name__ == '__main__':
    sys.exit(main())
