# -*- coding: utf-8 -*-
"""학습에서 굳힌 모델로 검체 채점 (j 판용, 재학습 없음)

  모델: models/<조합>/pooled/BEST.joblib (09_save_best_j.py 가 만듦)
  검체에 없는 행렬이 든 조합은 건너뛰고 더 짧은 조합 사용 (pick())
  없는 블록은 train_mean 으로 채우고 채운 비율 기록

  사용: python score_one_sample.py                (검체 폴더 안에서)
        python score_one_sample.py --model 'mean|entropy|jsd'
        python score_one_sample.py --allmodels    (굳힌 조합 전부)
"""

# import 금지 - 스크립트 본체가 그대로 돎
if __name__ != '__main__':
    raise ImportError(__file__ + ' 은(는) 스크립트다. import 하지 않는다')

import os as _os, sys as _sys
# 0_setup 찾기: 위 폴더로 올라가며 찾고, 없으면 METH_CONF_DIR, 그래도 없으면 중단
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


import os, sys, glob
import numpy as np, pandas as pd, joblib
from sklearn.metrics import roc_auc_score

sys.path.insert(0, os.getcwd())
import config as C

a = sys.argv
WANT = a[a.index('--model') + 1] if '--model' in a else None
# 굳힌 조합 전부 한 번에 채점 (검체 피처는 한 번만 읽음)
ALL = '--allmodels' in a
import glob as _glob

def _md5f(p):
    """모델 파일 md5 앞 8자리 (채점 행에 남김)"""
    import hashlib
    h = hashlib.md5()
    with open(p, 'rb') as _fh:
        for _c in iter(lambda: _fh.read(1 << 20), b''):
            h.update(_c)
    return h.hexdigest()[:8]
V, ML = C.VERSION, C.ML_OUT
TRAIN = _M + '/results/panels/%s/%s/step04_ML_classifier/models' % (V.split('_')[0], V)
# MODELS_DIR 주면 모델 자리만 바꿈 (검체 자료, 산출 자리는 그대로)
TRAIN = os.environ.get('MODELS_DIR') or TRAIN
SM = os.path.basename(os.getcwd())
RATIOS = [r for r in C.MUT_RATIOS if r > 0]
if not RATIOS:
    # 리드 부족으로 모든 비율이 건너뛰어진 경우 (채점할 칸 없음)
    sys.exit('%s: 풀 부족으로 채점할 비율이 없다. 이 깊이에서는 채점 안 함' % SM)
_SKIP = []        # (조합, 비율), 특징 못 읽어 건너뛴 칸

# 이미 채점된 검체 건너뛰기 - 표의 모델 지문(md5)이 현재 모델과 같을 때만
#   모델 다시 굳히면 지문이 달라져 다시 채점
def _already_done():
    if not ALL:
        return False                      # 조합 하나만 채점할 때는 안 건너뜀
    p = '%s/y_prob_all_%s_all_combos.csv' % (ML, V)
    if not os.path.exists(p) or os.path.getsize(p) == 0:
        return False
    try:
        _t = pd.read_csv(p, usecols=['model_file', 'model_md5'])
    except Exception:
        return False                      # 지문 열 없는 표 → 다시 채점
    if _cfg.SCORE_UNIT == 'consistent' and 'y_df' not in pd.read_csv(p, nrows=0).columns:
        return False                      # y_df 없는 표 → 다시 채점
    # 특징이 점수 파일보다 새것이면 다시 채점
    _feat = [f for f in _glob.glob('%s/[0-9]*_*_%s.csv' % (ML, V))
             if not os.path.basename(f).startswith('y_prob')]
    if any(os.path.getmtime(f) > os.path.getmtime(p) for f in _feat):
        return False
    # 조합마다 모든 비율이 있어야 완료
    _want = {round(x, 6) for x in RATIOS}
    _r = pd.read_csv(p, usecols=['model_file', 'ratio'])
    for _mf, _g in _r.groupby('model_file'):
        if not _want <= {round(x, 6) for x in _g['ratio'].unique()}:
            return False
    have = set(map(tuple, _t.drop_duplicates().values))
    now = set()
    for _f in _glob.glob(TRAIN + '/*/pooled/BEST.joblib'):
        now.add((os.path.basename(os.path.dirname(os.path.dirname(_f))), _md5f(_f)))
    if have == now and now:
        print('%s: 이미 채점됨. 조합 %d개 지문이 현재 모델과 같습니다. 건너뜁니다.'
              % (SM, len(now)))
        return True
    return False


if _already_done():
    sys.exit(0)

def have(key, p):
    return os.path.exists('%s/%d_%s_%s.csv' % (ML, p, key, V))

def pick():
    """검체에 있는 행렬만으로 쓸 수 있는 모델 고름"""
    # 기본 후보는 1차 조합(mean|entropy|jsd|llr|pdr)부터 긴 순서로
    cands = [WANT] if WANT else ['mean|entropy|jsd|llr|pdr',
                                 'mean|entropy|jsd|llr',
                                 'mean|entropy|jsd',
                                 'mean|entropy']
    p0 = int(RATIOS[0] * C.TARGET_DEPTH_READS)
    for name in cands:
        f = '%s/%s/pooled/BEST.joblib' % (TRAIN, name)
        if not os.path.exists(f): continue
        B = joblib.load(f)
        B['_md5'] = _md5f(f)
        keys = B['feature_set'].split('|')
        miss = [k for k in keys if not have(k, p0)]
        if miss:
            print('  %s 못 씀: 검체에 %s 행렬 없음' % (name, ','.join(miss)))
            continue
        return name, B, keys
    return None, None, None

name, B, keys = (None, None, None) if ALL else pick()
if B is None and not ALL:
    # 모델이 없는 건지 검체 특징이 없는 건지 구분 (후보 목록은 pick() 과 같음)
    _cands = [WANT] if WANT else ['mean|entropy|jsd|llr|pdr', 'mean|entropy|jsd|llr',
                                  'mean|entropy|jsd', 'mean|entropy']
    _has = any(os.path.exists('%s/%s/pooled/BEST.joblib' % (TRAIN, n)) for n in _cands)
    if _has:
        print('%s: 검체 특징이 없습니다. 아직 안 만들어졌을 수 있습니다 (판이 끝난 뒤 다시)' % SM)
    else:
        print('%s: 학습 모델이 없습니다 (%s). 09_save_best_j.py 를 먼저 돌리세요' % (SM, TRAIN))
    sys.exit(1)
_cache = {}


def _load(p, k):
    """(비율, 피처) CSV 캐시해서 읽음"""
    if (p, k) not in _cache:
        _cache[(p, k)] = pd.read_csv('%s/%d_%s_%s.csv' % (ML, p, k, V), index_col=0)
    return _cache[(p, k)]


def run_model(name, B, keys, rows):
    feats = list(B['features'])
    tm = pd.Series(B.get('train_mean', {})).reindex(feats)
    thr = float(B.get('threshold_spec95', np.nan))
    st = B.get('score_type', 'predict_proba')
    m = B['model']
    _tag = 'pooled'
    for r in RATIOS:
        p = int(r * C.TARGET_DEPTH_READS)
        try:
            Xs, y = [], None
            for k in keys:
                d = _load(p, k)
                y = d['Type'].astype(int) if y is None else y
                Xs.append(d.drop(columns=['Type']).add_prefix(k + '_'))
            X0 = Xs[0]
            for e in Xs[1:]:
                X0 = X0.join(e, how='inner')
            y = y.loc[X0.index]
        except Exception as e:
            print('  비율 %g%% 건너뜀 (%s)' % (r * 100, type(e).__name__))
            _SKIP.append((name, r))
            continue

        have_n = sum(1 for f in feats if f in X0.columns)
        X = X0.reindex(columns=feats).fillna(tm).fillna(0.0)
        miss = 100.0 * (len(feats) - have_n) / len(feats)

        # 포화 우회: 깊은 판에서 predict_proba 가 전부 1.0 에 붙으면 순위가 사라져 AUC 0.500
        #   그때만 decision_function 값으로 바꿈 (모델은 그대로, 같은 z 를 다르게 읽을 뿐)
        _sat = ''
        _thr_z = None      # 포화로 df 쓸 때의 문턱 (z 공간)
        # consistent: 전환과 상관없이 decision_function 도 같이 기록
        _ydf = None
        if _cfg.SCORE_UNIT == 'consistent':
            _ydf = m.decision_function(X) if hasattr(m, 'decision_function') else np.full(len(X), np.nan)
        if st == 'predict_proba' and hasattr(m, 'predict_proba'):
            sc = m.predict_proba(X)[:, 1]
            # 폭(max-min)이 1e-6 미만이면 뭉친 것으로 봄 (고유값 개수로는 부족)
            _rng = float(np.max(sc) - np.min(sc)) if len(sc) else 0.0
            if len(X) > 1 and _rng < 1e-6 and hasattr(m, 'decision_function'):
                sc = m.decision_function(X)
                _sat = '+df(포화)'
                # 문턱도 z 공간으로 옮김: 로지스틱이면 logit(thr) 이 같은 지점
                #   LogisticRegression 일 때만, 다른 모형은 y_pred 를 -1 로 둠
                _e = m
                if hasattr(_e, 'named_steps'):
                    _e = list(_e.named_steps.values())[-1]
                if type(_e).__name__ == 'LogisticRegression' and thr == thr and 0.0 < thr < 1.0:
                    _thr_z = float(np.log(thr / (1.0 - thr)))
        elif hasattr(m, 'decision_function'):
            sc = m.decision_function(X)
        else:
            sc = np.full(len(X), np.nan)
        for i, sid in enumerate(X.index):
            rows.append({'sample': SM, 'version': V, 'model_file': name,
                         'model_depth': _tag,
                         'model_md5': B.get('_md5', ''),
                         'model': B['model_name'], 'feature_set': B['feature_set'],
                         'ratio': r, 'fraction_pct': r * 100, 'mut_reads': p,
                         'rep': sid, 'y_true': int(y.loc[sid]),
                         'y_score': float(sc[i]), 'score_type': st + _sat,
                         'thr_spec95': (thr if _thr_z is None else _thr_z),
                         'y_pred': (int(sc[i] > thr) if (thr == thr and not _sat)
                                    else (int(sc[i] > _thr_z) if _thr_z is not None else -1)),
                         'feat_missing_pct': round(miss, 2),
                         **({} if _ydf is None else {'y_df': float(_ydf[i])})})
    return rows

# 채점할 모델 고르기. --allmodels 면 굳힌 조합 전부 (피처 CSV 는 _load 캐시)
_models = []
if ALL:
    _p0 = int(RATIOS[0] * C.TARGET_DEPTH_READS)
    for _f in sorted(_glob.glob(TRAIN + '/*/pooled/BEST.joblib')):
        _nm = os.path.basename(os.path.dirname(os.path.dirname(_f)))
        _B = joblib.load(_f)
        _B['_md5'] = _md5f(_f)
        _ks = _B['feature_set'].split('|')
        if any(not have(_k, _p0) for _k in _ks):
            continue
        _models.append((_nm, _B, _ks))
    if not _models:
        print('%s: 쓸 수 있는 굳힌 모델이 없습니다 (%s)' % (SM, TRAIN))
        sys.exit(1)
    print('%s · 조합 %d개 채점' % (SM, len(_models)))
else:
    print('%s · %s · %s · %s (%d열)'
          % (SM, name, B['model_name'], B['feature_set'], len(B['features'])))
    _models = [(name, B, keys)]

rows = []
for _nm, _B, _ks in _models:
    run_model(_nm, _B, _ks, rows)
thr = float(_models[-1][1].get('threshold_spec95', float('nan')))
name = _models[-1][0] if not ALL else 'all_combos'

if not rows:
    print('%s: 채점할 것이 없습니다. val04_1 이 행렬을 만들었는지 확인' % SM); sys.exit(1)
T = pd.DataFrame(rows)
os.makedirs(ML, exist_ok=True)
# --model 주면 조합 이름을 파일명에 넣음 (다른 결과 덮어쓰기 방지)
_tag = '_all_combos' if ALL else (('_' + name.replace('|', '-')) if WANT else '')
out = '%s/y_prob_all_%s%s.csv' % (ML, V, _tag)
# 임시 파일에 쓰고 이름 바꿈 (중간에 죽어도 반쪽 표가 안 남음)
_tmp = out + '.writing'   # 산출 이름은 ASCII 만
T.to_csv(_tmp, index=False)
os.replace(_tmp, out)
pos = T[T.y_true == 1]
print('  %d행 · 결측 특징 %.1f%% · 비율 %d개' % (len(T), T.feat_missing_pct.mean(), T.ratio.nunique()))
neg = T[T.y_true == 0]
def _auc(r):
    q = pos[pos.ratio == r]
    if not len(q) or not len(neg): return float('nan')
    return roc_auc_score(np.r_[np.zeros(len(neg)), np.ones(len(q))],
                         np.r_[neg.y_score.values, q.y_score.values])
print('  검출률 %s' % ' · '.join('%g%%:%.2f' % (r * 100, pos[pos.ratio == r].y_pred.mean())
                                 for r in sorted(pos.ratio.unique())))
print('  AUC    %s' % ' · '.join('%g%%:%.3f' % (r * 100, _auc(r))
                                 for r in sorted(pos.ratio.unique())))
print('  정상 오탐 %.2f (문턱 %.4f)' % (neg.y_pred.mean(), thr))
if _SKIP:
    # 표는 남기고 실패로 끝냄 (다음 실행에서 빠진 비율 보고 다시 채점)
    print('!! 특징을 못 읽어 건너뛴 칸 %d개: %s' % (len(_SKIP), ' · '.join('%s %g%%' % (n, r * 100) for n, r in _SKIP[:10])))
    sys.exit(1)
print('  저장 %s' % out)
