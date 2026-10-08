#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
바닥 0.55 를 두 조각으로 가른다 · 36-2 의 「안 쟀다」를 잰다. 2026-10-07

왜
  35-1 의 귀무 바닥이 중앙 0.545~0.561 이다. 그 높이가 어디서 오는지 모르면
  다음 판에서 어디를 고칠지 모른다.
      바닥 = 0.5 + (상관된 복제본)  + (모델 4개 중 최댓값)
                   ↑ 단일 모델 귀무     ↑ 최댓값 귀무 − 단일 모델 귀무
  대부분이 최댓값 쪽이면 «모델을 하나로 고정»하면 바닥이 내려간다.
  대부분이 복제본 쪽이면 «n_eff 를 키워야» 한다. 고칠 곳이 다르다.

  그리고 36-2 에 적은 논증(「신호가 있으면 4모델이 모여 팽창이 작다」)을
  같이 잰다. 모델별 AUC 를 남겨 관측과 귀무에서 퍼짐을 비교하면 된다.
  36-2 에 「안 쟀다」고 적었으니 여기서 잰다.

읽는 규칙 (숫자 보기 전에)
  단일 모델 귀무 중앙이 0.5 근처      바닥은 거의 전부 최댓값 선택 탓이다
                                      → 다음 판은 모델을 미리 하나로 고정한다
  단일 모델 귀무 중앙이 0.54 근처      바닥은 복제본 상관 탓이다
                                      → n_eff 를 키워야 한다. 모델 고정은 소용없다
  관측 퍼짐 < 귀무 퍼짐                36-2 의 논증이 성립한다. 판정이 보수적이다
  관측 퍼짐 ≥ 귀무 퍼짐                논증이 성립하지 않는다. 36-2 를 정정한다
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

CELLS = [
    ('GBM', 'j15bl200', '<METH_ROOT>/results/판', '25', 'entropy',
     ['entropy']),
    ('GBM', 'j15bl200', '<METH_ROOT>/results/판', '25', '1차 조합',
     ['mean', 'entropy', 'jsd', 'llr', 'pdr']),
    ('폐', 'j16jsdb200', '<METH_ROOT>/lung/results/panels', '25',
     'entropy', ['entropy']),
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


def per_model(X, y):
    """모델별 AUC 를 전부 돌려준다. 최댓값으로 뭉개지 않는다."""
    cv = StratifiedKFold(NF, shuffle=True, random_state=SEED)
    out = {}
    for nm, mk in models().items():
        pipe = Pipeline([('imp', SimpleImputer(strategy='mean',
                                               keep_empty_features=True)),
                         ('m', mk)])
        meth = ('predict_proba' if hasattr(mk, 'predict_proba')
                else 'decision_function')
        try:
            o = np.asarray(cross_val_predict(pipe, X, y, method=meth, cv=cv,
                                             n_jobs=1), dtype=float)
        except Exception:
            out[nm] = np.nan
            continue
        out[nm] = auc(y, o[:, 1] if o.ndim == 2 else o)
    return out


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


def summarize(rows, label):
    """rows = 반복별 {모델: auc} 목록."""
    names = list(models().keys())
    M = np.array([[r.get(n, np.nan) for n in names] for r in rows], dtype=float)
    mx = np.nanmax(M, axis=1)
    spread = mx - np.nanmedian(M, axis=1)
    out = ['%-6s 최댓값 중앙 %.4f · 95%% %.4f' % (label, np.nanmedian(mx),
                                                 np.nanpercentile(mx, 95))]
    for j, n in enumerate(names):
        out.append('       %-4s 중앙 %.4f · 95%% %.4f'
                   % (n, np.nanmedian(M[:, j]), np.nanpercentile(M[:, j], 95)))
    out.append('       퍼짐(최댓값 − 모델중앙) 중앙 %.4f' % np.nanmedian(spread))
    return out, float(np.nanmedian(mx)), {n: float(np.nanmedian(M[:, j]))
                                          for j, n in enumerate(names)}, \
        float(np.nanmedian(spread))


def main():
    print('바닥 조각내기 · 모델별 AUC 를 남긴다 · 50/50 · 반복 %d회' % NREP)
    for tag, panel, base, mr, sname, groups in CELLS:
        X, neg = load(base, panel, mr, groups)
        ni, pi = np.where(neg)[0], np.where(~neg)[0]
        rng = np.random.default_rng(SEED)
        obs, nul = [], []
        for _ in range(NREP):
            a = rng.choice(ni, HALF, replace=False)
            b = rng.choice(pi, HALF, replace=False)
            y = np.concatenate([np.zeros(HALF, int), np.ones(HALF, int)])
            obs.append(per_model(X[np.concatenate([a, b])], y))
            c = rng.choice(ni, 2 * HALF, replace=False)
            yn = np.zeros(c.size, int)
            yn[rng.choice(c.size, c.size // 2, replace=False)] = 1
            nul.append(per_model(X[c], yn))
        print()
        print('=' * 70)
        print('%s %s · %s · %s · %d칸' % (tag, panel, '0.5%' if mr == '25'
                                          else '0.1%', sname, X.shape[1]))
        print('=' * 70)
        lo, _, nmed, nspr = summarize(nul, '귀무')
        lp, _, omed, ospr = summarize(obs, '관측')
        for ln in lo:
            print(ln)
        for ln in lp:
            print(ln)
        single = min(nmed.values())
        best_single = max(nmed.values())
        print('  바닥 조각   단일 모델 중 가장 낮은 귀무 중앙 %.4f' % single)
        print('              단일 모델 중 가장 높은 귀무 중앙 %.4f' % best_single)
        print('              최댓값 선택이 더한 양 %+.4f'
              % (np.nanmedian([max(r.values()) for r in nul]) - best_single))
        print('  퍼짐 비교   관측 %.4f 대 귀무 %.4f → %s'
              % (ospr, nspr, '관측이 작다 (36-2 논증 성립)' if ospr < nspr
                 else '관측이 크거나 같다 (36-2 정정 필요)'))
    return 0


if __name__ == '__main__':
    sys.exit(main())
