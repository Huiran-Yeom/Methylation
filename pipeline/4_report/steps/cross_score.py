# -*- coding: utf-8 -*-
"""교차 채점 — A 트리의 모델로 B 트리의 학습 특징을, B 의 모델로 A 를 채점한다

  split_ab.py 로 나눈 두 무리를 따로 학습한 트리 둘을 받는다. 학습 CV(W)는 같은 세포주 ·
  정상인이 접기 양쪽에 들어가지만, 교차(G)는 분류기가 처음 보는 세포주 · 정상인이다.
  주의: 패널 선정과 jsd · pdr 의 정상 기준 분포는 A · B 를 모두 보고 만들어졌으므로, G 는
        «파이프라인 전체가 처음 보는» 경우가 아니다.
  주의: G 는 모든 비율을 합쳐 학습한 pooled BEST 의 점수이고, W 는 --mut 비율만으로 학습한
        교차검증 값이다. 그래서 G − W 에는 «처음 보는 무리» 효과 말고도 pooled · 비율별 학습의
        차이가 섞인다. W 는 BEST 와 같은 모델의 값(W)과 모델 넷 중 최고(W4)를 함께 낸다.

  채점은 3_validate/steps/score_one_sample.py 의 run_model 과 같다:
    특징 = <피처>_<열> 접두 · 피처끼리 안쪽 결합 · 모델 features 순서로 맞춤 ·
    빈 칸은 모델의 train_mean → 0 으로 채움 · predict_proba 폭 < 1e-6 이면 decision_function.
    채점 묶음(모델 하나 · 자료 하나)마다 단위가 하나라 점수 단위가 섞이지 않는다.

사용:  python cross_score.py --a <트리 A 의 METH_ROOT> --b <트리 B 의 METH_ROOT>
           [--gen 16] [--depth 5k] [--mut 25] [--panels bl200,jsd200,jsdb200]
           [--combos 'mean|entropy|jsd|llr|pdr,mean|entropy|jsd|pdr']
       --mut 은 섞은 종양 리드 수다(5k 의 0.5% = 25).
"""
# 이 파일은 스크립트다. 다른 코드가 import 하면 본체가 그대로 돈다.
if __name__ != '__main__':
    raise ImportError(__file__ + ' 은(는) 스크립트다 — import 하지 않는다')

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
import meth_config as _cfg

import argparse
import os
import sys

import joblib   # 이 파이프라인의 09 단계가 만든 BEST.joblib 만 연다(외부 파일 아님)
import numpy as np
import pandas as pd
from sklearn.metrics import roc_auc_score

ap = argparse.ArgumentParser()
ap.add_argument('--a', required=True)
ap.add_argument('--b', required=True)
ap.add_argument('--gen', default=str(_cfg.GEN))
ap.add_argument('--depth', default='5k')
ap.add_argument('--mut', type=int, default=25)
ap.add_argument('--panels', default='bl200,jsd200,jsdb200')
ap.add_argument('--combos', default='mean|entropy|jsd|llr|pdr,mean|entropy|jsd|pdr')
a = ap.parse_args()
T = {'A': a.a.rstrip('/'), 'B': a.b.rstrip('/')}
if os.path.realpath(T['A']) == os.path.realpath(T['B']):
    sys.exit('--a 와 --b 가 같은 트리다')
PANELS = a.panels.split(',')
COMBOS = a.combos.split(',')


def ml(g, p):
    v = 'j%s%s_%s_uni' % (a.gen, p, a.depth)
    return '%s/results/panels/j%s%s/%s/step04_ML_classifier' % (T[g], a.gen, p, v), v


def score(model_g, data_g, p, combo):
    md, _ = ml(model_g, p)
    dd, v = ml(data_g, p)
    bf = '%s/models/%s/pooled/BEST.joblib' % (md, combo)
    if not os.path.exists(bf):
        sys.exit('모델이 없다: %s — 그 트리에서 2단계(09)를 먼저 돌린다' % bf)
    B = joblib.load(bf)
    keys = B['feature_set'].split('|')
    Xs, y = [], None
    for k in keys:
        ff = '%s/%d_%s_%s.csv' % (dd, a.mut, k, v)
        if not os.path.exists(ff):
            sys.exit('특징 표가 없다: %s — --gen · --depth · --mut 을 확인한다' % ff)
        d = pd.read_csv(ff, index_col=0)
        y = d['Type'].astype(int) if y is None else y
        Xs.append(d.drop(columns=['Type']).add_prefix(k + '_'))
    X0 = Xs[0]
    for e in Xs[1:]:
        X0 = X0.join(e, how='inner')
    y = y.loc[X0.index]
    feats = list(B['features'])
    tm = pd.Series(B.get('train_mean', {})).reindex(feats)
    have = sum(1 for f in feats if f in X0.columns)
    X = X0.reindex(columns=feats).fillna(tm).fillna(0.0)
    m = B['model']
    sc = m.predict_proba(X)[:, 1] if hasattr(m, 'predict_proba') else m.decision_function(X)
    sat = ''
    if float(np.max(sc) - np.min(sc)) < 1e-6 and hasattr(m, 'decision_function'):
        sc = m.decision_function(X); sat = ' (포화 → df)'
    return roc_auc_score(y, sc), len(y), int(y.sum()), have, len(feats), B.get('model_name', type(m).__name__) + sat


def within_cv(g, p, combo):
    """(W, W4): BEST 와 같은 모델의 교차검증 값, 모델 넷 중 최고."""
    md, v = ml(g, p)
    c = pd.read_csv('%s/cv_%d_%s.csv' % (md, a.mut, v))
    c = c[c.feature_set == combo]
    best = joblib.load('%s/models/%s/pooled/BEST.joblib' % (md, combo)).get('model_name')
    same = c.loc[c.model == best, 'cv_auc']
    return (float(same.max()) if len(same) else float('nan')), float(c['cv_auc'].max())


def main():
    print('A = %s\nB = %s' % (T['A'], T['B']))
    for combo in COMBOS:
        print('== %s · %s · 섞은 종양 리드 %d' % (combo, a.depth, a.mut))
        for p in PANELS:
            r = {}
            for mg, dg in (('A', 'B'), ('B', 'A')):
                auc, n, npos, have, nf, mn = score(mg, dg, p, combo)
                r[mg + dg] = auc
                print('  %-8s %s→%s  AUC %.4f  (행 %d · 양성 %d · 특징 %d/%d · %s)' % (p, mg, dg, auc, n, npos, have, nf, mn))
            (wa, wa4), (wb, wb4) = within_cv('A', p, combo), within_cv('B', p, combo)
            g, w, w4 = (r['AB'] + r['BA']) / 2, (wa + wb) / 2, (wa4 + wb4) / 2
            print('  %-8s G=%.4f · W=%.4f (A %.4f · B %.4f) · G−W %+.4f · W4=%.4f (A %.4f · B %.4f) · G−W4 %+.4f'
                  % (p, g, w, wa, wb, g - w, w4, wa4, wb4, g - w4))


if __name__ == '__main__':
    sys.exit(main())
