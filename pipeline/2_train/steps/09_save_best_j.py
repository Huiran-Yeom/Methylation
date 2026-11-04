# -*- coding: utf-8 -*-
"""검증에 넘길 모델과 문턱 저장 (피처 조합마다 BEST.joblib)

조합마다 모델 4개 교차검증해서 1등 저장, 문턱은 정상 OOF 점수 95백분위
  models/<조합>/pooled/BEST.joblib          비율 합친 것 (검증은 이것만 씀)
  models/<조합>/mut_<p>_reads/BEST.joblib   비율마다 따로 (--pooled 면 안 만듦)
합칠 때 정상은 한 번만 넣음 (비율마다 같은 0% 행이라), 암만 비율마다 모음

사용법: 09_save_best_j.py [판이름 ...] [--feat mean,entropy ...] [--allcombos] [--pooled]
"""

# 스크립트 전용 - import 하면 본체가 돌아 자료를 덮어씀
if __name__ != '__main__':
    raise ImportError(__file__ + ' 은(는) 스크립트다 — import 하지 않는다')

import os as _os, sys as _sys
# 0_setup 찾기: 위 폴더로 올라가며 → METH_CONF_DIR (못 찾으면 중단)
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
try:
    import meth_config as _cfg
    _M, _D = _cfg.METH_ROOT, _cfg.DATASET_ROOT
except ImportError:
    raise SystemExit(
        '설정을 못 찾았습니다. 0_setup/meth_config.py 가 있는 자리를 찾지 못했습니다.\n'
        '  이 파일은 정리본 안에서 도는 것이라 기본값으로 넘어가지 않습니다.\n'
        '  트리 밖에서 돌리려면 METH_CONF_DIR 로 0_setup 자리를 주십시오.')


import os
import sys
import glob
import numpy as np
import pandas as pd
import joblib
from sklearn.pipeline import Pipeline as _Pipe
from sklearn.impute import SimpleImputer as _Imputer
from sklearn.ensemble import RandomForestClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.svm import SVC
from xgboost import XGBClassifier
from sklearn.model_selection import (cross_val_predict, StratifiedKFold,
                                     StratifiedGroupKFold)
from sklearn.metrics import roc_auc_score
from sklearn.base import clone

R = _M + '/results'
SEED, NF = 42, 5
try:   # SEED 는 config.conf 값 사용
    SEED = int(_cfg.SEED)
except Exception:
    pass
FIXED = str(getattr(_cfg, 'FIXED_MODEL', '') or '').strip()   # 비면 모델 4개 비교

_a = sys.argv[1:]
FEATS = []
while '--feat' in _a:
    i = _a.index('--feat')
    FEATS.append([x.strip() for x in _a[i + 1].split(',') if x.strip()])
    del _a[i:i + 2]
# --allcombos: 있는 피처의 모든 조합 (08 과 같은 규칙)
# --pooled: 비율별 모델 건너뜀 (검증은 pooled 만 씀)
ALLC = '--allcombos' in _a
POOLED_ONLY = '--pooled' in _a
_a = [x for x in _a if x not in ('--allcombos', '--pooled')]
if not FEATS and not ALLC:
    FEATS = [['mean', 'entropy'], ['mean', 'entropy', 'jsd']]
VERS = [v for v in _a if not v.startswith('-')]
if not VERS:
    VERS = sorted(os.path.basename(d) for d in
                  glob.glob(R + '/panels/*/j*_uni') + glob.glob(R + '/panels/*/j*_cov'))


def models():
    ms = {'Logistic Regression': LogisticRegression(random_state=SEED, max_iter=2000),
            'Random Forest': RandomForestClassifier(random_state=SEED, n_jobs=1),
            'XGBoost': XGBClassifier(random_state=SEED, eval_metric='logloss', n_jobs=1),
            'SVM': SVC(random_state=SEED)}
    # FIXED_MODEL 로 모델 하나 고정
    if FIXED:
        if FIXED not in ms:
            raise SystemExit('FIXED_MODEL=%r 은 모델 이름이 아니다: %s' % (FIXED, ' · '.join(ms)))
        return {FIXED: ms[FIXED]}
    return ms


def load(D, V, p, feats):
    """한 비율의 피처들을 옆으로 붙임 (열 이름 앞에 피처 이름)"""
    Xs = None
    y = None
    for k in feats:
        f = '%s/%d_%s_%s.csv' % (D, p, k, V)
        if not os.path.exists(f):
            return None, None
        d = pd.read_csv(f, index_col=0)
        if 'Type' not in d.columns:
            return None, None
        if y is None:
            y = d['Type'].astype(int)
        x = d.drop(columns=['Type']).add_prefix(k + '_')
        Xs = x if Xs is None else Xs.join(x, how='inner')
    return Xs, y


def fit_save(X, y, md, tag):
    """모델 4개 교차검증해서 1등 저장"""
    mu = X.mean()
    # 빈칸은 fold 안에서 학습 fold 평균으로 채움 (이 문턱이 검증 판정에 쓰여서 누출 막음)
    # mu 는 train_mean 으로 저장
    # pooled 는 같은 샘플이 여러 비율로 들어있어 "<샘플>@r<n>" 의 @ 앞으로 그룹 나눔
    # 그룹이 다 다르면(비율별) 그냥 StratifiedKFold
    _grp = [str(ix).split("@")[0] for ix in X.index]
    if len(set(_grp)) < len(_grp):
        cv = StratifiedGroupKFold(NF, shuffle=True, random_state=SEED)
    else:
        cv = StratifiedKFold(NF, shuffle=True, random_state=SEED)
        _grp = None
    best = None
    for mn, mk in models().items():
        meth = 'predict_proba' if hasattr(mk, 'predict_proba') else 'decision_function'
        try:
            _mk = _Pipe([('imp', _Imputer(strategy='mean', keep_empty_features=True)),
                         ('m', clone(mk))])
            o = cross_val_predict(_mk, X, y, method=meth, n_jobs=1, cv=cv,
                                  groups=_grp)
            o = np.asarray(o, dtype=float)
            o = o[:, 1] if o.ndim > 1 else o
            auc = float(roc_auc_score(y, o))
        except Exception as e:
            print('      %-20s 실패 (%s)' % (mn, type(e).__name__))
            failures.append('모델 %s 실패 (%s): %s' % (mn, type(e).__name__, tag))
            continue
        if best is None or auc > best[1]:
            best = (mn, auc, o, meth)
    if best is None:
        print('      %s: 네 모델 다 실패' % tag)
        failures.append('네 모델 다 실패: %s' % tag)
        return None
    mn, auc, o, meth = best
    m = models()[mn]
    m = _Pipe([('imp', _Imputer(strategy='mean', keep_empty_features=True)), ('m', m)])
    m.fit(X, y)
    yv = np.asarray(y)
    thr = float(np.percentile(o[yv == 0], 95))
    # OOF 정상, 암 점수 분포도 저장 (정상 수가 적어 문턱이 거칠어서 해석용)
    _neg = np.asarray(o[yv == 0], dtype=float)
    _pos = np.asarray(o[yv == 1], dtype=float)
    def _dist(a):
        if not len(a):
            return {}
        return {'n': int(len(a)),
                'min': float(np.min(a)), 'q05': float(np.percentile(a, 5)),
                'q25': float(np.percentile(a, 25)), 'median': float(np.median(a)),
                'q75': float(np.percentile(a, 75)), 'q95': float(np.percentile(a, 95)),
                'max': float(np.max(a)), 'mean': float(np.mean(a)),
                'sd': float(np.std(a, ddof=1)) if len(a) > 1 else 0.0}
    negdist, posdist = _dist(_neg), _dist(_pos)
    # 정상 점수 상위 20개 저장 (문턱 신뢰구간 계산용)
    negtop = sorted(float(x) for x in _neg)[-20:]
    os.makedirs(md, exist_ok=True)
    joblib.dump({'model': m, 'features': list(X.columns),
                 'feature_set': tag.split(' ')[0], 'model_name': mn,
                 'score_type': meth, 'cv_auc': round(auc, 4),
                 'threshold_spec95': thr,
                 'oof_neg_dist': negdist, 'oof_pos_dist': posdist,
                 'oof_neg_top20': negtop,
                 'cv_auc_주의': '모든 비율을 합친 값이다. 특정 농도의 값이 아니다 — 5% 이상이 1.000 이라 평균이 올라간다. 보고에 그대로 쓰지 않는다.',
                 'sens_at_spec95': float((o[yv == 1] > thr).mean()),
                 'spec_on_cv': float((o[yv == 0] <= thr).mean()),
                 'train_mean': {k: float(v) for k, v in mu.items()},
                 'n_train': int(len(yv)), 'n_pos': int((yv == 1).sum()),
                 'seed': SEED, 'cv_folds': NF},
                md + '/BEST.joblib')
    print('      %-22s %-20s AUC %.4f · 문턱 %.4f · 민감도 %.3f'
          % (tag, mn, auc, thr, float((o[yv == 1] > thr).mean())))
    print('        OOF 음성 n=%d · 중앙 %.4f · 95분위 %.4f · 최대 %.4f  |  양성 중앙 %.4f'
          % (negdist.get('n', 0), negdist.get('median', float('nan')),
             negdist.get('q95', float('nan')), negdist.get('max', float('nan')),
             posdist.get('median', float('nan'))))
    return auc


made = 0
made_by_V = {}      # 판별 저장 수
expect_by_V = {}    # 판별 기대 수
failures = []   # 실패 기록, 끝에서 있으면 중단
_conv = _cfg.count_convergence()
for V in VERS:
    D = '%s/panels/%s/%s/step04_ML_classifier' % (R, V.split('_')[0], V)
    if not os.path.isdir(D):
        continue
    ps = sorted(int(os.path.basename(f).split('_')[0])
                for f in glob.glob('%s/*_mean_%s.csv' % (D, V)))
    if not ps:
        print('%s 건너뜀: mean 행렬 없음' % V)
        continue
    print('')
    print('== %s   비율 %d개' % (V, len(ps)))
    _F = FEATS
    if ALLC:
        import itertools as _it
        _av = [k for k in ('mean', 'entropy', 'jsd', 'llr', 'pdr')
               if os.path.exists('%s/%d_%s_%s.csv' % (D, ps[0], k, V))]
        _F = [list(c) for n in range(1, len(_av) + 1) for c in _it.combinations(_av, n)]
        print('   전수 조합 %d개 (피처 %s)' % (len(_F), ','.join(_av)))
    expect_by_V[V] = len(_F) * (1 if POOLED_ONLY else len(ps) + 1)
    for feats in _F:
        tag = '|'.join(feats)
        parts = []
        for p in ps:
            X, y = load(D, V, p, feats)
            if X is None:
                continue
            if not POOLED_ONLY:
                fit_save(X, y, '%s/models/%s/mut_%d_reads' % (D, tag, p),
                         '%s mut_%d' % (tag, p))
                made += 1
                made_by_V[V] = made_by_V.get(V, 0) + 1
            parts.append((p, X, y))
        if not parts:
            continue
        cols = parts[0][1].columns
        Xa = [parts[0][1][parts[0][2] == 0]]
        ya = [parts[0][2][parts[0][2] == 0]]
        for p, X, y in parts:
            Xa.append(X.reindex(columns=cols)[y == 1])
            ya.append(y[y == 1])
        Xp = pd.concat(Xa, axis=0)
        yp = pd.concat(ya, axis=0)
        # 행 이름은 <샘플>@r<번호> - 샘플 이름을 남겨야 fit_save 가 그룹으로 나눌 수 있음
        Xp.index = ['%s@r%d' % (str(ix), k)
                    for k, ix in enumerate(Xp.index)]
        yp.index = Xp.index
        fit_save(Xp, yp, '%s/models/%s/pooled' % (D, tag), '%s pooled' % tag)
        made += 1
        made_by_V[V] = made_by_V.get(V, 0) + 1

print('')
print('끝. 저장한 모델 %d개' % made)
# 판마다 기대 수와 저장 수 비교, 다르면 실패 처리
for _V in sorted(made_by_V):
    _n = made_by_V[_V]
    _e = expect_by_V.get(_V)
    if _e is None:
        continue
    print('   %-24s 저장 %d / 기대 %d%s'
          % (_V, _n, _e, '' if _n == _e else '   <-- 어긋남'))
    if _n != _e:
        failures.append('%s: 저장 %d / 기대 %d — 어긋남' % (_V, _n, _e))
_cfg.convergence_check(_conv, failures, '09')
if failures:
    print('')
    print('실패 %d건:' % len(failures))
    for f in failures[:20]:
        print('  - %s' % f)
    sys.exit(1)
