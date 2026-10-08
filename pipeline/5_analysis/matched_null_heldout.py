#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
관측 CV 와 귀무 바닥을 «같은 클래스 크기»에서 비교한다. 2026-10-07

왜 필요한가
  33-4 의 귀무는 음성 100행을 50/50 으로 가른다. 실제 보고 CV 는 100/100 이다.
  AUC 의 귀무 폭은 n_pos·n_neg 에 달려 있어 2500 대 10000 이면 귀무가 약 두 배
  넓다. 그래서 귀무 95%(0.65~0.68) 를 관측(0.686 등)과 바로 못 비교한다.

  해결: «관측 쪽을 50/50 으로 내린다». 음성 50 + 양성 50 을 뽑아 같은 파이프라인을
  돌려 관측 분포를 내고, 같은 크기의 귀무 분포(음성 50/50)와 나란히 놓는다.
  두 분포의 중앙 차이가 「신호」이고, 겹침이 「구분 못 함」이다.

읽는 규칙 (숫자 보기 전에 적는다)
  관측 중앙이 귀무 95% 밖    그 CV 는 파이프라인 잡음으로 설명되지 않는다
  관측 중앙이 귀무 95% 안    구분 못 한다. 「신호 없음」이 아니라 「이 설계로 못 가림」
  관측 중앙 < 귀무 중앙      관측이 귀무보다 «나쁘다». 모델이 실제 구조를
                            못 쓰고 있다는 뜻이고, 과적합(34-3)과 같은 방향이다
  분포 겹침 비율도 같이 적는다 — 중앙만으로는 폭을 못 본다
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

SETS = {
    'entropy': ['entropy'],
    'llr 뺀 조합': ['mean', 'entropy', 'jsd', 'pdr'],
    '1차 조합': ['mean', 'entropy', 'jsd', 'llr', 'pdr'],
}
# 추기 32-3(다) 의 설계: 정본과 수정판에서 각각 귀무를 내 차이를 n_eff 효과로 읽는다.
JOBS = [
    ('정본', 'j15bl200', '<METH_ROOT>/results/판', '25'),
    ('수정판', 'j15bl200', '<METH_ROOT_HELDOUT>/results/판', '25'),
    ('정본', 'j15jsd200', '<METH_ROOT>/results/판', '25'),
    ('수정판', 'j15jsd200', '<METH_ROOT_HELDOUT>/results/판', '25'),
    ('정본', 'j15jsdb200', '<METH_ROOT>/results/판', '25'),
    ('수정판', 'j15jsdb200', '<METH_ROOT_HELDOUT>/results/판', '25'),
]


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


def load(base, panel, mr, groups):
    d = panel + '_5k_uni'
    cols = []
    for g in groups:
        f = '%s/%s/%s/step04_ML_classifier/%s_%s_%s.csv' % (base, panel, d, mr, g, d)
        df = pd.read_csv(f, index_col=0)
        if 'Type' not in df.columns:
            raise SystemExit('Type 열이 없다: ' + f)
        df = df.drop(columns=['Type'])
        df.columns = ['%s::%s' % (g, c) for c in df.columns]
        cols.append(df)
    X = cols[0]
    for e in cols[1:]:
        X = X.join(e, how='inner')
    neg = X.index.to_series().str.startswith('Normal').to_numpy()
    return X.to_numpy(dtype=float), neg


def main():
    print('관측 대 귀무 — 둘 다 %d/%d · 반복 %d회' % (HALF, HALF, NREP))
    print('%-4s %-11s %-5s %-13s %7s %7s %7s %7s %7s %6s'
          % ('자료', '판', '팔', '특징', '관측중앙', '관측95',
             '귀무중앙', '귀무95', '차이', '칸수'))
    for tag, panel, base, mr in JOBS:
        arm = '0.5%' if mr == '25' else '0.1%'
        for name, groups in SETS.items():
            try:
                X, neg = load(base, panel, mr, groups)
            except Exception as e:
                print('%-4s %-11s %-5s %-13s 못 읽음 %s' % (tag, panel, arm, name, e))
                continue
            ni = np.where(neg)[0]
            pi = np.where(~neg)[0]
            rng = np.random.default_rng(SEED)
            obs, nul = [], []
            for _ in range(NREP):
                a = rng.choice(ni, HALF, replace=False)
                b = rng.choice(pi, HALF, replace=False)
                idx = np.concatenate([a, b])
                y = np.concatenate([np.zeros(HALF, int), np.ones(HALF, int)])
                obs.append(cv_max(X[idx], y))
                # 귀무: 음성만 100행 중 50/50 가짜 라벨
                c = rng.choice(ni, 2 * HALF, replace=False) if ni.size >= 2 * HALF \
                    else ni
                yn = np.zeros(c.size, int)
                yn[rng.choice(c.size, c.size // 2, replace=False)] = 1
                nul.append(cv_max(X[c], yn))
            o = np.array(obs)
            u = np.array(nul)
            print('%-4s %-11s %-5s %-13s %7.4f %7.4f %7.4f %7.4f %+7.4f %6d'
                  % (tag, panel, arm, name, np.median(o), np.percentile(o, 95),
                     np.median(u), np.percentile(u, 95),
                     np.median(o) - np.median(u), X.shape[1]))
    print()
    print('차이 = 관측 중앙 − 귀무 중앙. 관측 중앙이 귀무 95% 안이면 구분 못 한다.')
    return 0


if __name__ == '__main__':
    sys.exit(main())
