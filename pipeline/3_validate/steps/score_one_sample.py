# -*- coding: utf-8 -*-
"""학습에서 저장한 모델을 환자 검체에 적용: j 판용 (2026-08-24)

  ~/tools/bin/val04_2_apply.py 를 j 판 모델 구조에 맞춘 것.
  바꾼 것은 모델을 찾는 자리뿐이다.
      m 판   models/BEST_POOL.joblib
      j 판   models/<조합>/pooled/BEST.joblib   (09_save_best_j.py 가 만든다)
  조합은 mean|entropy(Baseline·JSD-1) 와 mean|entropy|jsd(JSD-2) 둘.
  jsd 행렬이 없는 검체는 앞의 것으로 내려간다.


  다시 학습하지 않는다. joblib 의 모델·특징순서·문턱·학습셋평균을 그대로 쓴다.
  모델은 비율을 합쳐 학습한 것이라 하나뿐이다. 검체의 암 비율을 몰라도 쓸 수 있다.

  모델 고르기
    BEST_POOL      특징 전부. 검체에 그 행렬이 다 있어야 쓴다
    BEST_POOLMEAN  mean 단독. 위가 안 되면 이것
  없는 블록은 joblib 의 train_mean 으로 채우고, 채운 비율을 기록한다.

  사용: python val04_2_apply.py            (검체 폴더 안에서)
        python val04_2_apply.py --model BEST_POOLMEAN
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


import os, sys, glob
import numpy as np, pandas as pd, joblib
from sklearn.metrics import roc_auc_score

sys.path.insert(0, os.getcwd())
import config as C

a = sys.argv
WANT = a[a.index('--model') + 1] if '--model' in a else None
# 모든 조합을 한 번에 채점한다. 검체 피처는 한 번만 읽는다.
ALL = '--allmodels' in a
import glob as _glob

def _md5f(p):
    """2026-09-15 채점 행에 모델 파일 지문을 남긴다. 어느 모델이 이 점수를 냈는지."""
    import hashlib
    h = hashlib.md5()
    with open(p, 'rb') as _fh:
        for _c in iter(lambda: _fh.read(1 << 20), b''):
            h.update(_c)
    return h.hexdigest()[:8]
V, ML = C.VERSION, C.ML_OUT
TRAIN = _M + '/results/panels/%s/%s/step04_ML_classifier/models' % (V.split('_')[0], V)
# 모델만 다른 트리에서 읽어야 할 때가 있다 — 같은 검체 피처에 새로 굳힌 모델을
#   적용해 보는 경우다. 검체 자료와 산출 자리는 그대로 두고 모델 자리만 옮긴다.
TRAIN = os.environ.get('MODELS_DIR') or TRAIN
SM = os.path.basename(os.getcwd())
RATIOS = [r for r in C.MUT_RATIOS if r > 0]

# ─── [이미 된 일 관문] 잎에 한 번 두어 런처마다 기억하지 않는다 ───
#   오늘 이미 끝난 5k 178검체를 다시 채점시켜 3~6시간을 버릴 뻔했다.
#   런처들에는 건너뛰기가 있는데
#   잎에는 하나도 없었다. 호출하는 쪽을 전부 고치는 방식은 계속 실패한다.
#
#   "파일 있으면 건너뛰기" 가 아니라 지문이 현재 모델과 같으면 건너뛴다.
#   모델을 다시 굳히면 지문이 달라져 스스로 다시 채점한다.
def _already_done():
    if not ALL:
        return False                      # 조합 하나만 채점할 때는 건너뛰지 않는다
    p = '%s/y_prob_all_%s_all_combos.csv' % (ML, V)
    if not os.path.exists(p) or os.path.getsize(p) == 0:
        return False
    try:
        _t = pd.read_csv(p, usecols=['model_file', 'model_md5'])
    except Exception:
        return False                      # 지문 열이 없는 옛 표 → 다시 한다
    have = set(map(tuple, _t.drop_duplicates().values))
    now = set()
    for _f in _glob.glob(TRAIN + '/*/pooled/BEST.joblib'):
        now.add((os.path.basename(os.path.dirname(os.path.dirname(_f))), _md5f(_f)))
    if have == now and now:
        print('%s: 이미 채점됨 — 조합 %d개 지문이 현재 모델과 같습니다. 건너뜁니다.'
              % (SM, len(now)))
        return True
    return False


if _already_done():
    sys.exit(0)
# ──────────────────────────────────────────────────────────────────────

def have(key, p):
    return os.path.exists('%s/%d_%s_%s.csv' % (ML, p, key, V))

def pick():
    """검체에 있는 행렬만으로 쓸 수 있는 모델을 고른다."""
    # j 판은 조합마다 폴더가 따로다. jsd 가 든 것을 먼저 본다.
    # [실행차단 교정] 기본 후보에 1차 조합이 없었다.
    #   --model 을 안 주면 영원히 'mean|entropy|jsd' 로 검증한다.
    #   그러면 사전등록의 1차 지표(mean|entropy|jsd|llr|pdr)가 통째로 안 나온다.
    #   검증이 다 돌아간 뒤에야 드러날 자리였다. 긴 것부터 본다.
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
    # [교정] cands 는 pick() 지역변수라 여기서 NameError 가 났다.
    #   진짜 원인 메시지가 그 예외에 덮였다. 목록을 밖에 둔다.
    _cands = [WANT] if WANT else ['mean|entropy|jsd|llr|pdr', 'mean|entropy|jsd|llr',
                                  'mean|entropy|jsd', 'mean|entropy']
    _has = any(os.path.exists('%s/%s/pooled/BEST.joblib' % (TRAIN, n)) for n in _cands)
    if _has:
        print('%s: 검체 특징이 없습니다 — 아직 안 만들어졌을 수 있습니다 (판이 끝난 뒤 다시)' % SM)
    else:
        print('%s: 학습 모델이 없습니다 (%s) — 09_save_best_j.py 를 먼저 돌리세요' % (SM, TRAIN))
    sys.exit(1)
_cache = {}


def _load(p, k):
    """같은 (비율, 피처) CSV 를 조합마다 다시 읽지 않는다."""
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
            continue

        have_n = sum(1 for f in feats if f in X0.columns)
        X = X0.reindex(columns=feats).fillna(tm).fillna(0.0)
        miss = 100.0 * (len(feats) - have_n) / len(feats)

        # [포화 우회] predict_proba 가 float64 한계로 1.0 에 붙으면
        #   검체마다 값이 같아져 AUC 가 정확히 0.500 이 된다(순위 소멸).
        #   원인: LLR 원값이 리드 수에 비례해 깊은 판에서 decision_function 이 60~106 까지 간다.
        #   z>36 이면 1/(1+e^-z) 는 1.0 과 구별되지 않는다.
        #   실측: bl200 50k z=106 · jsd200 100k z=60 → 전 검체 proba 1.0 → AUC 0.500
        #              bl200 5k z=8 → 0.9997 (정상) · jsd200 50k z=29 → 0.992 (아슬아슬)
        #   AUC 는 순위만 보므로 뭉친 것이 확인되면 decision_function 으로 바꾼다.
        #   모델은 그대로다. 같은 z 를 다른 방식으로 읽을 뿐이라 판정 규칙을 안 바꾼다.
        _sat = ''
        _thr_z = None      # 포화로 df 를 쓸 때의 문턱 (z 공간)
        if st == 'predict_proba' and hasattr(m, 'predict_proba'):
            sc = m.predict_proba(X)[:, 1]
            # 「고유값 2개 미만」 만으로는 부족하다.
            #   1차 조합에서 고유값 76개인데 범위가 1.000~1.000 인 경우를 봤다 —
            #   순위가 부동소수점 끝자리에만 남아 위태롭다. 폭으로 판정한다.
            #   폭(max-min)이 1e-6 미만이면 사실상 뭉친 것으로 본다.
            _rng = float(np.max(sc) - np.min(sc)) if len(sc) else 0.0
            if len(X) > 1 and _rng < 1e-6 and hasattr(m, 'decision_function'):
                sc = m.decision_function(X)
                _sat = '+df(포화)'
                # [문턱도 옮긴다] proba 문턱을 그냥 버리면 y_pred 가 공백이 되고,
                #   사전등록이 요구하는 특이도·위양성(⑪·②)을 그 칸에서 못 낸다.
                #   로지스틱은 predict_proba = sigmoid(decision_function) 이므로
                #   logit(thr) 이 정확히 같은 지점이다. 재학습도 재보정도 필요 없다.
                #   (2026-09 시점에 다섯 칸 전부 LogisticRegression 이었다.
                #    모델이 바뀌면 이 말은 안 맞을 수 있다 — 그래서 아래에서
                #    런타임에 type 을 확인한다. 주석을 믿지 않는다.
                #    다른 모형은 이 등식이 안 서므로 그대로 공백으로 둔다.)
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
                         'feat_missing_pct': round(miss, 2)})
    return rows

# ---- 어느 모델들을 채점할지 ----
# --allmodels 면 굳혀 둔 조합을 전부 채점한다.
#   검체 피처 CSV 는 _load 가 캐시하므로 조합이 늘어도 파일을 다시 읽지 않는다.
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
    print('%s: 채점할 것이 없습니다 — val04_1 이 행렬을 만들었는지 확인' % SM); sys.exit(1)
T = pd.DataFrame(rows)
os.makedirs(ML, exist_ok=True)
# --model 을 주면 파일명에 넣어 나눈다. 안 그러면 +JSD 결과를 덮어쓴다.
_tag = '_all_combos' if ALL else (('_' + name.replace('|', '-')) if WANT else '')
out = '%s/y_prob_all_%s%s.csv' % (ML, V, _tag)
# 원자적 쓰기. 중간에 죽으면 머리줄만 남은 표가 생기고,
#   score_parallel.sh 의 "머리에 model_md5 있으면 건너뛴다"가 그걸 완료로 본다.
_tmp = out + '.writing'   # 2026-09-29: 산출 이름은 ASCII 로 (한글 파일명 금지)
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
print('  저장 %s' % out)
