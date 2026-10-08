#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
칸 평균 «한 칸» 을 같은 CV 파이프라인에 넣으면 197칸보다 잘하나. 2026-10-07

왜
  33-3 에서 entropy 칸 평균 한 축의 AUC 가 0.5% 에서 0.669(GBM bl200)인데
  같은 특징의 197칸 다변량 CV 는 0.52 였다. 두 수의 단위가 다르다 —
  앞은 적합 없는 통계량이고 뒤는 교차검증된 모델이다. 그래서 그냥 못 비교한다.

  칸 평균을 «한 칸 특징»으로 만들어 같은 4모델·5접기·최댓값 파이프라인에 넣으면
  단위가 맞는다. 그 값이 197칸보다 높으면 197칸은 과적합이다.

읽는 규칙 (숫자 보기 전에)
  1칸 CV > 197칸 CV   197칸이 과적합이다. 「비-LLR 특징에 신호가 없다」가 아니라
                      「이 파이프라인이 그 신호를 못 쓴다」가 맞는 문장이다
  1칸 CV ≈ 197칸 CV   과적합이 아니다. 칸 평균 축의 0.669 는 적합 없는 통계량의
                      낙관이었던 것이다
  1칸 CV < 197칸 CV   197칸이 칸 평균보다 많은 것을 쓰고 있다. 과적합 가설 기각
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
NF = 5
ARMS = [('0.1%', '5'), ('0.5%', '25')]
TARGETS = [
    ('GBM', 'j15bl200', '<METH_ROOT>/results/판'),
    ('GBM', 'j15jsdb200', '<METH_ROOT>/results/판'),
    ('폐', 'j16jsdb200', '<METH_ROOT>/lung/results/panels'),
    ('폐', 'j16bl200', '<METH_ROOT>/lung/results/panels'),
]


def models():
    return {
        'LR': LogisticRegression(random_state=SEED),
        'RF': RandomForestClassifier(random_state=SEED, n_jobs=1),
        'XGB': XGBClassifier(random_state=SEED, eval_metric='logloss', n_jobs=1),
        # 08 은 SVM_PROBABILITY=False 다 (2_train/steps/config.py:66).
        #   probability=True 로 두면 Platt 스케일링이 들어가 점수가 달라지고
        #   AUC 가 어긋난다 — 폐 jsdb200 0.5% 에서 0.568 대 0.601 로 갈렸다.
        'SVM': SVC(random_state=SEED, probability=False),
    }


def auc(y, s):
    r = pd.Series(s).rank().to_numpy()
    p, n = int((y == 1).sum()), int((y == 0).sum())
    return (r[y == 1].sum() - p * (p + 1) / 2.0) / (p * n)


def cv_max(X, y):
    """4모델 중 최댓값. 08_ML_classifier 와 같게 roc 를 뒤집지 않는다."""
    cv = StratifiedKFold(NF, shuffle=True, random_state=SEED)
    best, who = -1.0, None
    for nm, mk in models().items():
        pipe = Pipeline([('imp', SimpleImputer(strategy='mean',
                                               keep_empty_features=True)),
                         ('m', mk)])
        meth = ('predict_proba' if hasattr(mk, 'predict_proba')
                else 'decision_function')
        o = np.asarray(cross_val_predict(pipe, X, y, method=meth, cv=cv,
                                         n_jobs=1), dtype=float)
        a = auc(y, o[:, 1] if o.ndim == 2 else o)
        if a > best:
            best, who = a, nm
    return best, who


def main():
    print('%-4s %-12s %-6s %8s %8s %8s %8s %-6s'
          % ('자료', '판', '팔', '1칸CV', '다칸CV', '적합없음', '칸수', '1칸모델'))
    for tag, panel, base in TARGETS:
        d = panel + '_5k_uni'
        for lab, mr in ARMS:
            f = '%s/%s/%s/step04_ML_classifier/%s_entropy_%s.csv' % (
                base, panel, d, mr, d)
            try:
                df = pd.read_csv(f, index_col=0)
            except Exception as e:
                print('%-4s %-12s %-6s 못 읽음 %s' % (tag, panel, lab, e))
                continue
            # ★ Type 은 라벨이다. 08_ML_classifier 도 drop 한다.
            if 'Type' not in df.columns:
                raise SystemExit('Type 열이 없다 — 파일 형식이 바뀌었다: ' + f)
            df = df.drop(columns=['Type'])
            y = (~df.index.to_series().str.startswith('Normal')).astype(int).to_numpy()
            Xm = df.to_numpy(dtype=float)
            one = np.nanmean(Xm, axis=1).reshape(-1, 1)
            a1, w1 = cv_max(one, y)
            am, _ = cv_max(Xm, y)
            raw = auc(y, one.ravel())
            print('%-4s %-12s %-6s %8.4f %8.4f %8.4f %8d %-6s'
                  % (tag, panel, lab, a1, am, raw, Xm.shape[1], w1))
    print()
    print('1칸CV   칸 평균 한 칸을 같은 4모델·5접기·최댓값 파이프라인에 넣은 값')
    print('다칸CV  같은 특징 전체 칸을 넣은 값')
    print('적합없음 칸 평균 축으로 바로 낸 AUC (교차검증 아님 · 참고용)')
    return 0


if __name__ == '__main__':
    sys.exit(main())
