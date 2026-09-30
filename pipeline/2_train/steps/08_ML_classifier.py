#!/usr/bin/env python
# coding: utf-8
"""
2_train/steps/08: ML 분류 + 분류점수(y_prob) 저장.  step04_1 과 2_train/steps/08 의 feature 를 통합한다.

원본 : Methylation/scripts/step04_ML_classifier.py 의 `3) 모델링`

회의록 250219 의 feature 조합 계획을 그대로 구현한다.
   1) Mean 으로만                    (step04_1  {p}_mean_v3.csv)
   2) Entropy 로만 (site 수준)       (step04_1  {p}_entropy_v3.csv)
   3) Site별 엔트로피 → read 합치기  (2_train/steps/08  {p}_ent_site_read_v3.csv)
   4) Read별 엔트로피 → site 합치기  (2_train/steps/08  {p}_ent_read_site_v3.csv)
   5) Site+Read 한번에 (패턴)        (2_train/steps/08  {p}_ent_pattern_v3.csv)
   조합: 1+2, 1+3, 1+4, 1+5

변경점
  1. **분류점수(y_prob)를 샘플마다 저장한다.** 원본은 Accuracy 하나만 남겨서 논문
     피규어(비율별 점그래프·혼동행렬·ROC)를 그릴 수 없었다.
  2. 2_train/steps/08 의 read-level feature 를 조합에 포함한다.
  3. feature 를 **샘플 id(index) 기준으로 병합**한다. 원본은 위치 기반
     `pd.concat(axis=1)` 이라 파일 순서에 의존했다. id 병합이 더 안전하고,
     순서가 같으면 결과도 같다.
  4. 예외를 삼키지 않고 traceback 을 찍고 실패 목록을 요약한다.

바뀌지 않은 것
  - train_test_split(test_size=0.2, random_state=SEED), 모델 4종과 하이퍼파라미터,
    LogisticRegression 기본 max_iter(원본 동작 — 수렴 경고가 나는 것도 원본과 같음)
  - stratify 는 기본 꺼짐(원본과 동일). config.STRATIFY=True 로 켤 수 있다.

SVM 점수
  원본 SVC 는 probability=False 라 predict_proba 가 없다. 기본값에서는
  decision_function 을 점수로 쓴다. 원본과 학습·예측이 완전히 동일하다.
  단 0~1 범위가 아니므로 `score_type` 열을 보고 **모델별로 따로** 문턱값을 잡아야 한다.
"""

# 2026-09-26: 09-23 사고에서 가드 없는 이 부류가 importlib 로 「확인」 되다 실제로 돌았다.
#   82개 파일이 덮였다(내용은 결정적이라 같았지만 그건 운이었다).
if __name__ != '__main__':
    raise ImportError(__file__ + ' 은(는) 스크립트다 — import 하지 않는다')

import itertools
import os
import sys
import traceback

import joblib
import numpy as np
import pandas as pd
from sklearn.ensemble import RandomForestClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import accuracy_score, roc_auc_score
from sklearn.base import clone
from sklearn.pipeline import Pipeline as _Pipe
from sklearn.impute import SimpleImputer as _Imputer
from sklearn.model_selection import (train_test_split,
                                     cross_val_predict, StratifiedKFold)
from sklearn.svm import SVC
from xgboost import XGBClassifier

import config as C

if '--ratios' in sys.argv:
    C.MUT_READ_COUNTS = [int(float(x) * C.TARGET_DEPTH_READS) for x in
                         sys.argv[sys.argv.index('--ratios') + 1].split(',')]
    print(f'[CLI] 대상 비율 -> {C.MUT_READ_COUNTS}')

# feature 이름 -> (폴더, 파일 접미사)
SOURCES = {
    # 최종 스펙([논문] 슬라이드 1 "block + mean/entropy")
    'mean':    (C.ML_OUT, 'mean'),
    'entropy': (C.ML_OUT, 'entropy'),
}

COMBOS = []

# 모든 조합: 있는 특징의 공집합 아닌 부분집합 전부 (2026-08-13)
# 손으로 적지 않는다. 특징이 늘어도 이 줄은 그대로 둔다.
import glob as _g, itertools as _it
_all = ['mean', 'entropy']
# 2026-09-16 (마) readent·mhl 을 뺀다. 사용자 결정.
#   채움률이 1차 지표 깊이(5k)에서 47.9% · 48.7% 다(음성 행만 보면 83.3% · 83.4%).
#   칸 절반이 그 검체에서 잰 값이 아니라 다른 검체들의 평균이다.
#   근거는 **채움률**이지 성능이 아니다. 성능을 근거로 삼지 않는다.
#   조합 2^7-1=127 -> 2^5-1=31. 옛 127 결과는 지우지 않는다.
for _k in ('jsd', 'llr', 'pdr'):
    if _g.glob('%s/*_%s_%s.csv' % (C.ML_OUT, _k, C.VERSION)):
        SOURCES[_k] = (C.ML_OUT, _k)
        _all.append(_k)
# 2026-09-14 [관문] _g.glob 이 조용히 비면 조합 수가 줄어든 채 **에러 없이** 끝난다.
#   실제로 2_train/steps/08 가 jsd·pdr 이 만들어지기 16분 전에 돌아 2^5-1=31 조합만 학습했다.
#   R109 에서 잡은 llr 비대칭(본판 63 · 대조 31)과 똑같은 구조가 판 안에서 재현됐다.
#   세는 곳이 없으면 또 놓친다. 여기서 멈춘다.
# 2026-09-16 (마) 이후 **5종**. 관문은 끄지 않는다 — 4종이나 6종이면 여전히 잡는다.
_want = ['mean', 'entropy', 'jsd', 'llr', 'pdr']
if len(_all) != len(_want):
    _miss = [x for x in _want if x not in _all]
    sys.exit('중단: 특징 %d종만 보입니다 (%s) — 없는 것 %s. 앞 단계가 안 끝났습니다.'
             % (len(_all), ','.join(_all), ','.join(_miss)))
COMBOS = [tuple(c) for n in range(1, len(_all) + 1) for c in _it.combinations(_all, n)]
print('[특징] %d조합 · %s' % (len(COMBOS), ' · '.join('|'.join(c) for c in COMBOS)))


# 선택 실험: optional/read_entropy_experiment.py 를 돌렸다면 아래 주석을 해제
# SOURCES.update({'rd_site_read': (C.READ_ENT_OUT,'ent_site_read'),
#                 'rd_read_site': (C.READ_ENT_OUT,'ent_read_site'),
#                 'rd_pattern':   (C.READ_ENT_OUT,'ent_pattern')})
# COMBOS += [('rd_site_read',), ('rd_read_site',), ('rd_pattern',), ('mean','rd_pattern')]


def make_models():
    """원본과 동일한 모델 4종.

    XGBClassifier 의 use_label_encoder=False 는 넣지 않았다. 라벨이 이미 0/1 정수라
    학습 결과에 영향이 없고(xgboost 1.x 에서 경고만 없애는 역할), 2.x 에서는 에러다.
    """
    return {
        'Logistic Regression': LogisticRegression(random_state=C.SEED),
        'Random Forest': RandomForestClassifier(random_state=C.SEED, n_jobs=C.N_JOBS),
        'XGBoost': XGBClassifier(
            random_state=C.SEED, eval_metric='logloss', n_jobs=C.N_JOBS
        ),
        'SVM': SVC(random_state=C.SEED, probability=C.SVM_PROBABILITY),
    }


def get_scores(model, X):
    if hasattr(model, 'predict_proba'):
        return model.predict_proba(X)[:, 1], 'predict_proba'
    if hasattr(model, 'decision_function'):
        return model.decision_function(X), 'decision_function'
    return np.full(len(X), np.nan), 'none'


def load_one(p, key):
    folder, suffix = SOURCES[key]
    path = f'{folder}/{p}_{suffix}_{C.VERSION}.csv'
    if not os.path.exists(path):
        raise FileNotFoundError(path)
    df = pd.read_csv(path, index_col=0)  # index = sample id (Normal1.. / GBM1..)
    y = df['Type'].astype(int)
    X = df.drop(columns=['Type']).add_prefix(f'{key}_')
    return X, y


def load_combo(p, keys):
    Xs, y0 = [], None
    for k in keys:
        X, y = load_one(p, k)
        Xs.append(X)
        if y0 is None:
            y0 = y
        else:
            common = y0.index.intersection(y.index)
            if not y0.loc[common].equals(y.loc[common]):
                raise ValueError(f'{k}: 라벨이 다릅니다')
    # 샘플 id 기준 병합 (원본의 위치 기반 concat 보다 안전)
    X = Xs[0]
    for extra in Xs[1:]:
        X = X.join(extra, how='inner')
    y = y0.loc[X.index]
    return X, y


def main():
    """피처 조합마다 모델 넷을 학습·교차검증하고 y_prob 를 남긴다. 242줄.

    받는 것
      ML_OUT/<비율>_<피처>_<판>.csv   피처 행렬 (mean·entropy·jsd·pdr·llr)
      SEED · NF                      StratifiedKFold(NF, shuffle, random_state=SEED)

    내는 것
      ML_OUT/y_prob_<비율>_<판>.csv · y_prob_all_<판>.csv
      ML_OUT/cv_<비율>_<판>.csv · Acc_<비율>_<판>.xlsx
      ML_OUT/models/<조합>/{pooled,mut_<n>_reads}/BEST.joblib

    구역
      1) 피처 원본 확인   SOURCES 에 실제 파일이 있는 조합만 남긴다
      2) 조합 만들기      _want 로 종 수를 고정한다 (5종 아니면 멈춘다)
      3) 비율 루프        비율마다 행렬을 붙여 X·y 를 만든다
      4) 모델 루프        LR·RF·XGB·SVC 를 교차검증 -> cv 표
      5) 저장             y_prob · BEST.joblib · 문턱(spec95)
    """
    C.guard_not_original()
    C.banner('2_train/steps/08  ML 분류 + y_prob 저장')
    C.ensure_dirs()

    print(f'SVM probability={C.SVM_PROBABILITY} · stratify={C.STRATIFY} · '
          f'train 점수 저장={C.SAVE_TRAIN_SCORES}\n')

    all_scores, failures = [], []
    _cvcount = {}   # 2026-09-14 비율별 cv 행수

    for p in C.MUT_READ_COUNTS:
        if p == 0:
            continue
        frac = p / C.TARGET_DEPTH_READS * 100
        print(f'=== mut_{p}_reads  (종양 비율 {frac:g}%) ===')
        sheet, rows_p, cv_rows = {}, [], []

        for keys in COMBOS:
            name = '|'.join(keys)
            try:
                X, y = load_combo(p, keys)
                strat = y if C.STRATIFY else None
                Xtr, Xva, ytr, yva = train_test_split(
                    X, y, test_size=C.TEST_SIZE, random_state=C.SEED, stratify=strat
                )
                # 결측은 **학습셋 평균으로만** 채운다. 검증셋 값을 섞으면 누출이다.
                # 행렬이 이미 채워져 온 경우에는 아무 일도 일어나지 않는다.
                trmean = Xtr.mean()
                if Xtr.isna().any().any() or Xva.isna().any().any():
                    Xtr = Xtr.fillna(trmean).fillna(0.0)
                    Xva = Xva.fillna(trmean).fillna(0.0)
                print('   [%s] 학습 %d(암 %d) · 검증 %d(암 %d)'
                      % (name, len(ytr), int(ytr.sum()), len(yva), int(yva.sum())))

                splits = [('val', Xva, yva)]
                if C.SAVE_TRAIN_SCORES:
                    splits.append(('train', Xtr, ytr))

                res = []
                for mname, model in make_models().items():
                    model.fit(Xtr, ytr)
                    acc = accuracy_score(yva, model.predict(Xva))
                    res.append({'Model': mname, 'Feature Set': name, 'Accuracy': acc})

                    # --- 학습된 모델 저장 (나중에 검증에 그대로 쓰기 위해) ---
                    # 특징 순서와 특이도 95% 문턱값을 함께 남긴다. 둘 중 하나만 없어도
                    # 같은 판정을 재현할 수 없다.
                    # 문턱은 두 가지를 남긴다.
                    #  thr_train : 학습셋 정상군의 95백분위. **새 검체에 적용할 때 쓴다.**
                    #  thr_val   : 검증셋 기준. 보고용이지만 같은 셋으로 평가하므로 낙관적이다.
                    sva, stype0 = get_scores(model, Xva)
                    sva = np.asarray(sva, dtype=float); yv = np.asarray(yva)
                    yt = np.asarray(ytr)
                    neg, pos = sva[yv == 0], sva[yv == 1]
                    thr = float(np.percentile(neg, 95)) if len(neg) else float('nan')
                    sens95 = float((pos > thr).mean()) if len(neg) and len(pos) else float('nan')

                    # 재사용 문턱은 **교차검증 예측(out-of-fold)** 으로 구한다.
                    # 학습셋 in-sample 점수를 쓰면 트리 모델이 학습셋을 외우는 탓에
                    # 문턱이 지나치게 낮아져 새 데이터에서 특이도가 무너진다.
                    # 📊 2026-08-06 실측(v4 0.1%): RF 특이도 0.005 · XGB 0.106
                    #    (로지스틱·SVM 은 0.89~0.91 로 정상: 외우지 못하기 때문)
                    try:
                        cvm = 'predict_proba' if hasattr(model, 'predict_proba') else 'decision_function'
                        oof = cross_val_predict(
                            clone(model), Xtr, ytr, method=cvm, n_jobs=1,
                            cv=StratifiedKFold(5, shuffle=True, random_state=C.SEED))
                        oof = np.asarray(oof, dtype=float)
                        oof = oof[:, 1] if oof.ndim > 1 else oof
                        negc = oof[yt == 0]
                        thr_cv = float(np.percentile(negc, 95)) if len(negc) else float('nan')
                    except Exception as _e:
                        print(f'      ! OOF 문턱 실패 ({type(_e).__name__}) — NaN 으로 둔다')
                        thr_cv = float('nan')
                    sens_cv = float((pos > thr_cv).mean()) if len(pos) and thr_cv == thr_cv else float('nan')
                    spec_cv = float((neg <= thr_cv).mean()) if len(neg) and thr_cv == thr_cv else float('nan')
                    mdir = f'{C.ML_OUT}/models/mut_{p}_reads'
                    os.makedirs(mdir, exist_ok=True)
                    joblib.dump({
                        'model': model, 'features': list(X.columns),
                        'feature_set': name, 'model_name': mname, 'score_type': stype0,
                        'threshold_spec95': thr_cv,          # 교차검증 문턱 — **새 검체에는 이것**
                        'sens_at_spec95': sens_cv,           # 그 문턱을 검증셋에 적용한 민감도
                        'spec_on_val': spec_cv,              # 실제로 지켜진 특이도 (0.95 근처여야 정상)
                        'threshold_spec95_val': thr,         # 검증셋 기준(낙관적, 보고서 수치)
                        'sens_at_spec95_val': sens95,
                        'accuracy_val': float(acc), 'version': C.VERSION,
                        'mut_reads': p, 'fraction_pct': frac,
                        'depth_reads': C.TARGET_DEPTH_READS, 'seed': C.SEED,
                        'n_train': int(len(Xtr)), 'n_val': int(len(Xva)),
                        # 새 검체에 없는 블록을 채울 값. 없으면 적용 단계에서
                        # 검체 자기 평균을 쓰게 되어 같은 판정을 재현할 수 없다.
                        'train_mean': {k: float(v) for k, v in trmean.items()},
                        'stratify': bool(C.STRATIFY),
                        'n_pos_val': int((yva == 1).sum()), 'n_neg_val': int((yva == 0).sum()),
                    }, f'{mdir}/{name.replace("|", "_")}__{mname.replace(" ", "")}.joblib')
                    # ── 5-fold 교차검증(전체 200개): 한 번 나눈 것에 결과가 달리지 않게.
                    #    최고 모델을 고르는 기준은 이 AUC 다 (2026-08-14).
                    try:
                        cm2 = 'predict_proba' if hasattr(model, 'predict_proba') else 'decision_function'
                        # 2026-09-14 [실행차단 교정] 여기서 **안 채운 X** 를 그대로 넣고 있었다.
                        #   결측이 있는 피처(jsd·pdr·llr)가 든 조합은 sklearn 이 ValueError 를
                        #   내고 except 로 빠져 **cv 표에 아예 안 들어갔다.** 실측: 0.1% 에서
                        #   127조합 중 31조합만 남았다. 사라진 96개가 전부 jsd·pdr 조합이고
                        #   거기에 1차 조합(mean|entropy|jsd|llr|pdr)이 들어 있다.
                        #   채우기를 파이프라인 안에 넣어 **접기마다 학습 접기 평균으로만** 채운다.
                        #   밖에서 미리 채우면 검증 접기 값이 섞여 누출이다.
                        _mdl = _Pipe([('imp', _Imputer(strategy='mean', keep_empty_features=True)),
                                      ('m', clone(model))])
                        oa = cross_val_predict(_mdl, X, y, method=cm2, n_jobs=1,
                                               cv=StratifiedKFold(5, shuffle=True, random_state=C.SEED))
                        oa = np.asarray(oa, dtype=float)
                        oa = oa[:, 1] if oa.ndim > 1 else oa
                        ya = np.asarray(y)
                        cv_rows.append({'mut_reads': p, 'fraction_pct': frac,
                                        'feature_set': name, 'model': mname,
                                        'cv_auc': round(float(roc_auc_score(ya, oa)), 4),
                                        'n': int(len(ya))})
                    except Exception as _e:
                        # 2026-09-14: print 만 하면 아무도 못 읽는다. 기록한다.
                        failures.append('교차검증 %s %s 비율%d - %s'
                                        % (name, mname, p, type(_e).__name__))

                    for split, Xs, ys in splits:
                        score, stype = get_scores(model, Xs)
                        pred = model.predict(Xs)
                        for j, sid in enumerate(Xs.index):
                            rows_p.append(
                                {
                                    'mut_reads': p,
                                    'fraction_pct': frac,
                                    'feature_set': name,
                                    'model': mname,
                                    'split': split,
                                    'sample_id': sid,
                                    'y_true': int(ys.iloc[j]),
                                    'y_score': float(score[j]),
                                    'score_type': stype,
                                    'y_pred': int(pred[j]),
                                }
                            )
                    print(f'   {name:22s} {mname:20s} acc={acc:.4f}')
                sheet[name] = pd.DataFrame(res)

            except FileNotFoundError as e:
                msg = f'p={p}, {name}: 입력 없음 → {os.path.basename(str(e))}'
                print(f'   X {msg}')
                failures.append(msg)
            except Exception as e:
                msg = f'p={p}, {name}: {type(e).__name__}: {e}'
                print(f'   X {msg}')
                traceback.print_exc()
                failures.append(msg)

        if rows_p:
            df = pd.DataFrame(rows_p)
            out = f'{C.ML_OUT}/y_prob_{p}_{C.VERSION}.csv'
            df.to_csv(out, index=False)
            print(f'   -> 점수 {len(df):,}행 저장 {out}')
            all_scores.extend(rows_p)

        if cv_rows:
            CV = pd.DataFrame(cv_rows).sort_values('cv_auc', ascending=False)
            CV.to_csv(f'{C.ML_OUT}/cv_{p}_{C.VERSION}.csv', index=False)
            _cvcount[p] = len(CV)
            b = CV.iloc[0]
            print(f"   ★ 최고 {b['feature_set']} + {b['model']}  교차검증 AUC {b['cv_auc']:.4f}")
            try:
                Xb, yb = load_combo(p, tuple(b['feature_set'].split('|')))
                # 2026-09-14 [교정] 전체 자료 평균으로 채우고 있었다.
                #   이 CV 가 만드는 thr 이 BEST.joblib 에 실려
                #   외부 검증 특이도를 판정한다. 출하되는 문턱이 누출이었다.
                mb = Xb.mean()
                best = _Pipe([('imp', _Imputer(strategy='mean', keep_empty_features=True)),
                              ('m', make_models()[b['model']])])
                best.fit(Xb, yb)
                cmb = 'predict_proba' if hasattr(best, 'predict_proba') else 'decision_function'
                ob = cross_val_predict(clone(best), Xb, yb, method=cmb, n_jobs=1,
                                       cv=StratifiedKFold(5, shuffle=True, random_state=C.SEED))
                ob = np.asarray(ob, dtype=float); ob = ob[:, 1] if ob.ndim > 1 else ob
                yb2 = np.asarray(yb)
                thr = float(np.percentile(ob[yb2 == 0], 95))
                mdir = f'{C.ML_OUT}/models/mut_{p}_reads'
                os.makedirs(mdir, exist_ok=True)
                joblib.dump({
                    'model': best, 'features': list(Xb.columns),
                    'feature_set': b['feature_set'], 'model_name': b['model'],
                    'score_type': cmb, 'cv_auc': float(b['cv_auc']),
                    'threshold_spec95': thr,
                    'sens_at_spec95': float((ob[yb2 == 1] > thr).mean()),
                    'spec_on_cv': float((ob[yb2 == 0] <= thr).mean()),
                    'train_mean': {k: float(v) for k, v in mb.items()},
                    'n_train': int(len(yb2)), 'stratify': bool(C.STRATIFY),
                    'version': C.VERSION, 'mut_reads': p, 'fraction_pct': frac,
                    'depth_reads': C.TARGET_DEPTH_READS, 'seed': C.SEED,
                }, f'{mdir}/BEST.joblib')
                print(f'   ★ 저장 {mdir}/BEST.joblib  (문턱 {thr:.4f})')
            except Exception as _e:
                print(f'   ! 최고 모델 저장 실패 {type(_e).__name__}: {_e}')

        if sheet:
            out = f'{C.ML_OUT}/Acc_{p}_{C.VERSION}.xlsx'
            with pd.ExcelWriter(out, engine='openpyxl') as w:
                for s, d in sheet.items():
                    d.to_excel(w, sheet_name=s.replace('|', '_')[:31], index=False)
            print(f'   -> Accuracy 저장 {out}')
        print()

    if all_scores:
        df = pd.DataFrame(all_scores)
        out = f'{C.ML_OUT}/y_prob_all_{C.VERSION}.csv'
        df.to_csv(out, index=False)
        print(f'통합 점수 {len(df):,}행 → {out}')
        print('\n비율별 val 샘플 수:')
        print(
            df[df.split == 'val']
            .groupby(['fraction_pct', 'feature_set'])
            .size()
            .unstack(fill_value=0)
            .to_string()
        )
    else:
        print('저장된 점수가 없습니다. step04_1 / 2_train/steps/08 출력을 확인하세요.')

    # 2026-09-14: 기대 행수를 코드에서 유도한다. 숫자를 박으면 낡는다.
    _expect = len(COMBOS) * len(make_models())
    for _p, _n in sorted(_cvcount.items()):
        if _n != _expect:
            failures.append('cv 행 %d / 기대 %d (비율 %d) - 조용히 적게 나왔다' % (_n, _expect, _p))
    if len(set(_cvcount.values())) > 1:
        failures.append('비율마다 cv 행수가 다르다: %s' % dict(sorted(_cvcount.items())))
    if failures:
        print(f'\n실패 {len(failures)}건:')
        for f in failures:
            print(f'  - {f}')
    else:
        print('\n실패 없음.')
    print('2_train/steps/08 완료.')
    if failures:
        # 2026-09-14: 실패가 있으면 종료 코드를 0 이 아니게 한다.
        sys.exit(1)


if __name__ == '__main__':
    sys.exit(main())
