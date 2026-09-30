# -*- coding: utf-8 -*-
"""검증에 넘길 모델과 문턱을 굳힌다. j 판용 (2026-08-21)

없어진 save_best_pool.py 자리다. 전신 ~/tools/bin/save_best_mean.py (08-14) 와
같은 방식을 쓰되 셋을 바꿨다.

  ① 특징 조합을 인자로 받는다      mean|entropy · mean|entropy|jsd
  ② 비율을 합친 판(pooled)도 만든다  m 판에서 비율마다 AUC 가 전부 1.0 이라
                                    비율별로는 모델을 못 골랐다
  ③ 판을 여러 개 한 번에 돈다

저장 위치
  models/<조합>/pooled/BEST.joblib          비율을 합쳐 하나로
  models/<조합>/mut_<p>_reads/BEST.joblib   비율마다 따로

문턱은 정상(y=0)의 교차검증 점수 95백분위. 전신과 같다.

합칠 때 정상을 한 번만 쓴다
  파일마다 Normal 100판이 들어 있는데 전부 같은 0% 판이다. 그대로 쌓으면
  같은 정상이 아홉 번 들어가 정상 쪽으로 치우친다. GBM 만 비율마다 모은다.

명령줄:  09_save_best_j.py [판이름 ...] [--feat mean,entropy] [--feat mean,entropy,jsd]
"""

# 이 파일은 스크립트다. 다른 코드가 import 하면 본체가 그대로 돌아
# 실제 자료를 덮어쓴다. 실행은 `python <파일>` 로만 한다.
if __name__ != '__main__':
    raise ImportError(__file__ + ' 은(는) 스크립트다 — import 하지 않는다')

import os as _os, sys as _sys
# 0_setup 을 찾는다. 이 파일이 「복사되어」 정리본 밖에서 돌 수도 있으므로
#   ① 자기 위에서 위로 올라가며 찾고 ② 환경변수 METH_CONF_DIR ③ 그래도 없으면
#   설정 없이 돌 수 있게 기본값으로 간다(검증이 검체 폴더로 복사해 쓰는 경로다).
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
try:   # 2026-09-29: SEED 가 config.conf 에 있는데 여기서만 박혀 있었다.
    SEED = int(_cfg.SEED)
except Exception:
    pass

_a = sys.argv[1:]
FEATS = []
while '--feat' in _a:
    i = _a.index('--feat')
    FEATS.append([x.strip() for x in _a[i + 1].split(',') if x.strip()])
    del _a[i:i + 2]
# --allcombos 면 있는 피처의 공집합 아닌 부분집합을 전부 굳힌다.
#   2_train/steps/08 의 조합 목록과 같은 규칙이다. 손으로 적지 않는다.
#   --pooled 면 비율별 적합을 건너뛴다. 검증은 pooled 모델만 쓴다
#   (score_one_sample.py:46  models/<조합>/pooled/BEST.joblib).
#   127조합 x 10비율 x 4모델 = 5,080 적합이 127 x 4 = 508 로 줄어든다.
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
    return {'Logistic Regression': LogisticRegression(random_state=SEED, max_iter=2000),
            'Random Forest': RandomForestClassifier(random_state=SEED, n_jobs=1),
            'XGBoost': XGBClassifier(random_state=SEED, eval_metric='logloss', n_jobs=1),
            'SVM': SVC(random_state=SEED)}


def load(D, V, p, feats):
    """한 비율의 여러 특징을 옆으로 이어붙인다. 열 이름에 특징 이름을 붙인다."""
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
    """모델 4종을 교차검증으로 겨뤄 1등을 저장한다."""
    mu = X.mean()
    # [실행차단 교정] 여기서 전체 자료 평균으로 채우고 있었다.
    #   아래 교차검증이 만드는 thr 이 pooled/BEST.joblib 에 실리고,
    #   score_one_sample.py 가 그 thr_spec95 로 외부 검증 특이도를 판정한다.
    #   즉 실제로 출하되는 문턱이 누출된 값이었다.
    #   채우기를 접기 안으로 넣어 접기마다 학습 접기 평균으로만 채운다.
    #   (mu 는 train_mean 으로 저장되므로 계산은 남긴다.)
    # [교정] pooled 는 같은 검체의 여러 비율판이 쌓여 있다.
    #   StratifiedKFold 는 그것을 학습 접기와 검증 접기에 갈라 넣는다 — 누출이다.
    #   인덱스가 "<검체>@r<n>" 이므로 @ 앞을 그룹으로 쓴다.
    #   비율별 판(중복 없음)에서는 그룹이 전부 달라 기존과 같게 동작한다.
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
    # [진단] 문턱이 어디에 서 있는지 남긴다.
    #   문턱은 OOF 음성의 95분위인데, 음성이 100개뿐이라 알갱이가 1% 다.
    #   게다가 실측 문턱이 5k 0.9871 · 15k 0.8610 으로 [0,1] 의 맨 끝이고
    #   깊이마다 0.126 이나 다르다. 분포가 없으면 나중에 해석할 수 없다.
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
    # [문턱 오차] 분위수 하나로는 "그 문턱이 얼마나 흔들리나" 를 못 낸다.
    #   음성 n=100 · p=0.95 일 때 95% 분위수의 95% 신뢰구간은 순서통계량으로 나온다
    #: 정렬한 음성 점수의 90~100번째. Bin(100,0.95) 의 평균 95 · SD 2.18 에서 온다.
    #   여유 있게 상위 20개를 남긴다. 이게 있어야 특이도 주장에 구간이 붙는다.
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
made_by_V = {}      # 2026-09-15 판마다 센다 (made 는 전체 누적이라 판별 부족을 못 잡는다)
expect_by_V = {}    # 판마다 기대치
failures = []   # 2026-09-14 print 만 하면 아무도 못 읽는다. 기록하고 끝에서 멈춘다.
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
        # 여기서 검체 이름을 r0,r1,... 로 덮어쓰고 있었다.
        #   그러면 같은 검체의 다른 비율판이 학습 접기와 검증 접기에 갈라 들어가도
        #   막을 수가 없다(GroupKFold 가 구조적으로 불가능해진다).
        #   이름을 남긴다. 중복은 비율을 붙여 구분하되 검체 부분을 보존한다.
        Xp.index = ['%s@r%d' % (str(ix), k)
                    for k, ix in enumerate(Xp.index)]
        yp.index = Xp.index
        fit_save(Xp, yp, '%s/models/%s/pooled' % (D, tag), '%s pooled' % tag)
        made += 1
        made_by_V[V] = made_by_V.get(V, 0) + 1

print('')
print('끝. 저장한 모델 %d개' % made)
# 만들어야 할 수와 실제 만든 수를 대조한다.
#   한 조합이 BEST 를 못 만들면 val04_2 의 pick() 이 조용히 다음 후보로 내려간다.
# [고침] 판마다 센다. 그리고 방향을 단정하지 않는다.
#   앞서는 made(판 루프 누적)를 _expect(판 하나 기준)와 비교했다. 범위가 달랐다.
#   5k+15k 를 굳히니 254 vs 127, 네 판이면 508 vs 127 로 정상인데 실패가 떴다.
#   그리고 문구가 "조용히 적게" 였는데 실제로는 많았다. 반대로 알려준 것이다.
#   [배운 것] 관문이 거짓으로 울면 그 관문은 죽는다. 254 때 이미 울었는데 지나쳤다.
#   멈추는 관문이어도 한 번 거짓이면 두 번째부터 안 읽힌다.
for _V in sorted(made_by_V):
    _n = made_by_V[_V]
    _e = expect_by_V.get(_V)
    if _e is None:
        continue
    print('   %-24s 저장 %d / 기대 %d%s'
          % (_V, _n, _e, '' if _n == _e else '   <-- 어긋남'))
    if _n != _e:
        failures.append('%s: 저장 %d / 기대 %d — 어긋남' % (_V, _n, _e))
if failures:
    print('')
    print('실패 %d건:' % len(failures))
    for f in failures[:20]:
        print('  - %s' % f)
    sys.exit(1)
