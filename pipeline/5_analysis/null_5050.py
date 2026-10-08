#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
이 파이프라인이 «신호가 전혀 없을 때» 낼 수 있는 AUC — 귀무 바닥. 2026-10-07

왜
  entropy 단독 다변량 CV 가 0.1% 와 0.5% 에서 0.6011 / 0.6012 였다. 비율 무관이라
  종양 기원이 아니라고 적었는데, 원인이 「기술 차이」가 아니라 **파이프라인 자체의
  귀무 바닥**일 수 있다. 197칸 × 학습 160행에 모델 4개 중 최댓값을 쓰면 신호가
  0 이어도 0.5 보다 높게 나온다.

설계 (검토 세션 제안)
  음성 100행(mut_0_reads)만 쓴다. 생물은 전부 같다 — 같은 풀에서 뽑은 복제본이다.
  거기에 무작위 50/50 가짜 라벨을 붙이고 «같은 CV 파이프라인»을 돌린다.
  반복해서 귀무 분포를 낸다.

  ★ 라벨 섞기(음성·양성을 뒤섞는 것)를 쓰지 않는다. 그러면 같은 클래스 안의
    리드 공유가 라벨과 어긋나 귀무가 0.5 로 주저앉는다. 재려는 팽창은
    «한 풀에서 온 행끼리의 공유가 라벨과 정렬될 때» 생긴다.

맞춘 것
  모델 4개 (LR · RF · XGB · SVM) · random_state=SEED · 결측은 평균 채움
  5접기 StratifiedKFold(shuffle, random_state=SEED) · cross_val_predict
  조합당 «모델 4개 중 최댓값» 까지 같다

한계
  50/50 이라 실제 100/100 보다 클래스가 작다. 자릿수 추정으로 읽는다.
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
NREP = 200
NF = 5
MR = '25'                       # 어느 팔의 CSV 를 읽어도 음성 100행은 같다

SETS = {
    'entropy': ['entropy'],
    'mean|entropy|jsd|pdr': ['mean', 'entropy', 'jsd', 'pdr'],
}
TARGETS = [
    ('GBM 정본', 'j15bl200', '<METH_ROOT>/results/판'),
    ('GBM 정본', 'j15jsd200', '<METH_ROOT>/results/판'),
    ('GBM 정본', 'j15jsdb200', '<METH_ROOT>/results/판'),
    ('폐', 'j16jsdb200', '<METH_ROOT>/lung/results/panels'),
]


def models(njobs=1):
    return {
        'LR': LogisticRegression(random_state=SEED),
        'RF': RandomForestClassifier(random_state=SEED, n_jobs=njobs),
        'XGB': XGBClassifier(random_state=SEED, eval_metric='logloss', n_jobs=njobs),
        # 08 은 SVM_PROBABILITY=False 다 (2_train/steps/config.py:66).
        #   probability=True 로 두면 Platt 스케일링이 들어가 점수가 달라지고
        #   AUC 가 어긋난다 — 폐 jsdb200 0.5% 에서 0.568 대 0.601 로 갈렸다.
        'SVM': SVC(random_state=SEED, probability=False),
    }


def auc(y, s):
    r = pd.Series(s).rank().to_numpy()
    p, n = int((y == 1).sum()), int((y == 0).sum())
    if p == 0 or n == 0:
        return np.nan
    return (r[y == 1].sum() - p * (p + 1) / 2.0) / (p * n)


def load_negatives(base, panel, groups):
    d = panel + '_5k_uni'
    cols = []
    for g in groups:
        f = '%s/%s/%s/step04_ML_classifier/%s_%s_%s.csv' % (base, panel, d, MR, g, d)
        df = pd.read_csv(f, index_col=0)
        # ★ Type 은 라벨이다. 08_ML_classifier 도 drop 한다.
        if 'Type' not in df.columns:
            raise SystemExit('Type 열이 없다 — 파일 형식이 바뀌었다: ' + f)
        df = df.drop(columns=['Type'])
        df = df.loc[df.index.to_series().str.startswith('Normal')]
        df.columns = ['%s::%s' % (g, c) for c in df.columns]
        cols.append(df.sort_index())
    X = pd.concat(cols, axis=1)
    if len(X) != 100:
        raise SystemExit('음성이 100행이 아니다: %d' % len(X))
    return X


def one_rep(X, rng):
    y = np.zeros(len(X), dtype=int)
    y[rng.choice(len(X), len(X) // 2, replace=False)] = 1
    cv = StratifiedKFold(NF, shuffle=True, random_state=SEED)
    best = -1.0
    for mk in models().values():
        pipe = Pipeline([('imp', SimpleImputer(strategy='mean',
                                               keep_empty_features=True)),
                         ('m', mk)])
        try:
            meth = ('predict_proba' if hasattr(mk, 'predict_proba')
                    else 'decision_function')
            o = cross_val_predict(pipe, X.to_numpy(dtype=float), y,
                                  method=meth, cv=cv, n_jobs=1)
            o = np.asarray(o, dtype=float)
            s = o[:, 1] if o.ndim == 2 else o
        except Exception:
            continue
        a = auc(y, s)
        if np.isfinite(a):
            # 08_ML_classifier 는 roc_auc_score 를 그대로 쓴다 — 0.5 아래를
            # 뒤집지 않는다. 뒤집으면 귀무가 크게 부풀어 비교가 안 된다.
            best = max(best, a)
    return best


def main():
    print('귀무 바닥 — 음성 100행에 무작위 50/50 라벨 · 반복 %d회' % NREP)
    print('모델 4개 중 최댓값까지 08_ML_classifier 와 같게 맞췄다')
    print()
    print('%-10s %-12s %-22s %8s %8s %8s %8s'
          % ('자료', '판', '특징', '중앙', '95%', '최대', '>0.60'))
    for tag, panel, base in TARGETS:
        for name, groups in SETS.items():
            try:
                X = load_negatives(base, panel, groups)
            except Exception as e:
                print('%-10s %-12s %-22s 못 읽음: %s' % (tag, panel, name, e))
                continue
            rng = np.random.default_rng(SEED)
            vals = np.array([one_rep(X, rng) for _ in range(NREP)])
            v = vals[np.isfinite(vals)]
            print('%-10s %-12s %-22s %8.4f %8.4f %8.4f %7.0f%%'
                  % (tag, panel, '%s (%d칸)' % (name, X.shape[1]),
                     np.median(v), np.percentile(v, 95), v.max(),
                     100.0 * (v > 0.60).mean()))
    print()
    print('이 값이 실제 CV 와 비슷하면 그 CV 는 신호가 아니라 바닥이다.')
    return 0


if __name__ == '__main__':
    sys.exit(main())
