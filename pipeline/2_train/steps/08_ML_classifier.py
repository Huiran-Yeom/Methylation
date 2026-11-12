#!/usr/bin/env python
# coding: utf-8
"""
2_train/steps/08: 피처 조합마다 모델 4개 학습, 교차검증, 점수(y_prob) 저장

피처 mean, entropy, jsd, llr, pdr 의 모든 조합 (31개), 샘플 id 기준으로 합침
SVM 은 probability=False 라 decision_function 을 점수로 씀 (score_type 열 참고)
옵션: --ratios 0.001,0.01
"""

# 스크립트 전용 - import 하면 본체가 돌아 자료를 덮어씀
if __name__ != '__main__':
    raise ImportError(__file__ + ' 은(는) 스크립트다. import 하지 않는다')

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
    # 기본 피처 (jsd, llr, pdr 은 아래에서 파일 있으면 추가)
    'mean':    (C.ML_OUT, 'mean'),
    'entropy': (C.ML_OUT, 'entropy'),
}

COMBOS = []

# 조합: 있는 피처의 공집합 아닌 부분집합 전부
import glob as _g, itertools as _it
_all = ['mean', 'entropy']
# readent, mhl 은 안 씀 (낮은 깊이에서 빈칸이 절반쯤이라)
for _k in ('jsd', 'llr', 'pdr'):
    if _g.glob('%s/*_%s_%s.csv' % (C.ML_OUT, _k, C.VERSION)):
        SOURCES[_k] = (C.ML_OUT, _k)
        _all.append(_k)
# 피처 5개가 다 없으면 중단 (앞 단계 덜 끝났을 때 일부 조합만 학습하지 않게)
_want = ['mean', 'entropy', 'jsd', 'llr', 'pdr']
if len(_all) != len(_want):
    _miss = [x for x in _want if x not in _all]
    sys.exit('중단: 특징 %d종만 보입니다 (%s). 없는 것 %s. 앞 단계가 안 끝났습니다.'
             % (len(_all), ','.join(_all), ','.join(_miss)))
COMBOS = [tuple(c) for n in range(1, len(_all) + 1) for c in _it.combinations(_all, n)]
print('[특징] %d조합 · %s' % (len(COMBOS), ' · '.join('|'.join(c) for c in COMBOS)))


# 선택 실험: optional/read_entropy_experiment.py 돌렸으면 아래 주석 해제
# SOURCES.update({'rd_site_read': (C.READ_ENT_OUT,'ent_site_read'),
#                 'rd_read_site': (C.READ_ENT_OUT,'ent_read_site'),
#                 'rd_pattern':   (C.READ_ENT_OUT,'ent_pattern')})
# COMBOS += [('rd_site_read',), ('rd_read_site',), ('rd_pattern',), ('mean','rd_pattern')]


def make_models():
    """모델 4개 (FIXED_MODEL 있으면 그 하나만)"""
    ms = {
        # 모델 설정은 09 와 같아야 함. n_jobs 만 다름 (09 는 1 고정, N_JOBS 기본 1 이면 같음)
        'Logistic Regression': LogisticRegression(random_state=C.SEED, max_iter=2000),
        'Random Forest': RandomForestClassifier(random_state=C.SEED, n_jobs=C.N_JOBS),
        'XGBoost': XGBClassifier(
            random_state=C.SEED, eval_metric='logloss', n_jobs=C.N_JOBS
        ),
        'SVM': SVC(random_state=C.SEED, probability=C.SVM_PROBABILITY),
    }
    # FIXED_MODEL 로 모델 하나 고정 (4개 중 최고를 고르면 그것만으로 값이 오름)
    fixed = (getattr(C, 'FIXED_MODEL', '') or '').strip()
    if fixed:
        if fixed not in ms:
            raise SystemExit('FIXED_MODEL=%r 은 모델 이름이 아니다: %s' % (fixed, ' · '.join(ms)))
        return {fixed: ms[fixed]}
    return ms


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
    df = pd.read_csv(path, index_col=0)  # index = 샘플 id (Normal1.. / GBM1..)
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
    # 샘플 id 기준으로 합침
    X = Xs[0]
    for extra in Xs[1:]:
        X = X.join(extra, how='inner')
    y = y0.loc[X.index]
    return X, y


def main():
    """피처 조합마다 모델 4개 학습, 교차검증, y_prob 저장

    입력
      ML_OUT/<n>_<피처>_<판>.csv   피처 행렬
      SEED                        StratifiedKFold(5, shuffle, random_state=SEED)

    출력
      ML_OUT/y_prob_<n>_<판>.csv, y_prob_all_<판>.csv
      ML_OUT/cv_<n>_<판>.csv, Acc_<n>_<판>.xlsx
      ML_OUT/models/mut_<n>_reads/<조합>__<모델>.joblib, BEST.joblib
    """
    C.guard_not_original()
    C.banner('2_train/steps/08  ML 분류 + y_prob 저장')
    C.ensure_dirs()

    print(f'SVM probability={C.SVM_PROBABILITY} · stratify={C.STRATIFY} · '
          f'train 점수 저장={C.SAVE_TRAIN_SCORES}\n')

    all_scores, failures = [], []
    _conv = C.count_convergence()
    _cvcount = {}   # 비율별 cv 행 수

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
                # 빈칸은 학습셋 평균으로만 채움 (검증셋 값 섞으면 누출)
                trmean = Xtr.mean()
                Xtr_raw = Xtr        # 채우기 전 값 (문턱용 OOF 예측은 fold 안에서 채움)
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

                    # 모델 저장 - 피처 순서, 특이도 95% 문턱, 학습 평균을 같이 저장
                    # 문턱 둘: OOF 기준(새 샘플용), 검증셋 기준(보고용, 낙관적)
                    sva, stype0 = get_scores(model, Xva)
                    sva = np.asarray(sva, dtype=float); yv = np.asarray(yva)
                    yt = np.asarray(ytr)
                    neg, pos = sva[yv == 0], sva[yv == 1]
                    thr = float(np.percentile(neg, 95)) if len(neg) else float('nan')
                    sens95 = float((pos > thr).mean()) if len(neg) and len(pos) else float('nan')

                    # 새 샘플용 문턱은 학습셋 OOF 예측으로 구함
                    # (in-sample 점수 쓰면 트리 모델이 외워서 문턱이 너무 낮아짐)
                    try:
                        cvm = 'predict_proba' if hasattr(model, 'predict_proba') else 'decision_function'
                        # 빈칸은 fold 마다 학습 fold 평균으로 채움 (아래 cv_auc 와 같게)
                        oof = cross_val_predict(
                            _Pipe([('imp', _Imputer(strategy='mean', keep_empty_features=True)), ('m', clone(model))]),
                            Xtr_raw, ytr, method=cvm, n_jobs=1,
                            cv=StratifiedKFold(5, shuffle=True, random_state=C.SEED))
                        oof = np.asarray(oof, dtype=float)
                        oof = oof[:, 1] if oof.ndim > 1 else oof
                        negc = oof[yt == 0]
                        thr_cv = float(np.percentile(negc, 95)) if len(negc) else float('nan')
                    except Exception as _e:
                        print(f'      ! OOF 문턱 실패 ({type(_e).__name__}). NaN 으로 둔다')
                        thr_cv = float('nan')
                    sens_cv = float((pos > thr_cv).mean()) if len(pos) and thr_cv == thr_cv else float('nan')
                    spec_cv = float((neg <= thr_cv).mean()) if len(neg) and thr_cv == thr_cv else float('nan')
                    mdir = f'{C.ML_OUT}/models/mut_{p}_reads'
                    os.makedirs(mdir, exist_ok=True)
                    joblib.dump({
                        'model': model, 'features': list(X.columns),
                        'feature_set': name, 'model_name': mname, 'score_type': stype0,
                        'threshold_spec95': thr_cv,          # OOF 문턱 (새 샘플에 사용)
                        'sens_at_spec95': sens_cv,           # 그 문턱의 검증셋 민감도
                        'spec_on_val': spec_cv,              # 그 문턱의 검증셋 특이도 (0.95 근처가 정상)
                        'threshold_spec95_val': thr,         # 검증셋 기준 문턱 (낙관적)
                        'sens_at_spec95_val': sens95,
                        'accuracy_val': float(acc), 'version': C.VERSION,
                        'mut_reads': p, 'fraction_pct': frac,
                        'depth_reads': C.TARGET_DEPTH_READS, 'seed': C.SEED,
                        'n_train': int(len(Xtr)), 'n_val': int(len(Xva)),
                        # 새 샘플에 없는 블록 채울 값
                        'train_mean': {k: float(v) for k, v in trmean.items()},
                        'stratify': bool(C.STRATIFY),
                        'n_pos_val': int((yva == 1).sum()), 'n_neg_val': int((yva == 0).sum()),
                    }, f'{mdir}/{name.replace("|", "_")}__{mname.replace(" ", "")}.joblib')
                    # 전체 자료 5-fold 교차검증. 최고 모델은 이 AUC 로 고름
                    try:
                        cm2 = 'predict_proba' if hasattr(model, 'predict_proba') else 'decision_function'
                        # 빈칸 채우기를 파이프라인 안에 넣음 - fold 마다 학습 fold 평균으로만 채움
                        # (안 채우면 sklearn 오류로 그 조합이 cv 표에서 빠짐)
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
                        # 실패 목록에 기록
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
            print(f"   최고 {b['feature_set']} + {b['model']}  교차검증 AUC {b['cv_auc']:.4f}")
            try:
                Xb, yb = load_combo(p, tuple(b['feature_set'].split('|')))
                # 빈칸 채우기는 파이프라인 안에서 (이 문턱이 BEST.joblib 에 들어가서 누출 막음)
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
                print(f'   저장 {mdir}/BEST.joblib  (문턱 {thr:.4f})')
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

    # 기대 cv 행 수 = 조합 수 x 모델 수
    _expect = len(COMBOS) * len(make_models())
    for _p, _n in sorted(_cvcount.items()):
        if _n != _expect:
            failures.append('cv 행 %d / 기대 %d (비율 %d) - 조용히 적게 나왔다' % (_n, _expect, _p))
    if len(set(_cvcount.values())) > 1:
        failures.append('비율마다 cv 행수가 다르다: %s' % dict(sorted(_cvcount.items())))
    C.convergence_check(_conv, failures, '08')
    if failures:
        print(f'\n실패 {len(failures)}건:')
        for f in failures:
            print(f'  - {f}')
    else:
        print('\n실패 없음.')
    print('2_train/steps/08 완료.')
    if failures:
        # 실패 있으면 종료 코드 1
        sys.exit(1)


if __name__ == '__main__':
    sys.exit(main())
