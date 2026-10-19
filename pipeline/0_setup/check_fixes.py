# -*- coding: utf-8 -*-
"""결함 수정이 실제로 그렇게 도는지 태워 본다. 자료·서버 없이 돈다.

  1. LLR 누수   split_ref_mix 가 사전등록 26-2 의 명단을 그대로 내고,
                11_pool 이 혼합군만 남기는가
  2. 모델 고정  08 · 09 의 모델 목록이 FIXED_MODEL 하나로 줄고, 틀린 이름이면 멈추는가
  3. 채우기     05 의 to_frame 이 묶음 평균으로 채우지 않고, 검증에서는 블록도 안 버리는가

단계 파일은 스크립트라 import 하면 본체가 돈다(가드가 막는다). 그래서 필요한 함수·구간만
소스에서 잘라 실행한다. 쓰는 법:  python 0_setup/check_fixes.py   (rc=0 이면 통과)
"""
import ast
import io
import os
import sys
import tempfile

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
TREE = os.path.dirname(HERE)
sys.path.insert(0, HERE)
os.environ.setdefault('METH_CONF_OPTIONAL', '1')
import meth_config as M  # noqa: E402

bad = []


def ok(cond, msg):
    print(('  통과  ' if cond else '  실패  ') + msg)
    if not cond:
        bad.append(msg)


def func_src(path, name):
    src = io.open(path, encoding='utf-8').read()
    for node in ast.parse(src).body:
        if isinstance(node, ast.FunctionDef) and node.name == name:
            return ast.get_source_segment(src, node)
    raise SystemExit('%s 에 %s 가 없다' % (path, name))


def block_src(path, start, end):
    src = io.open(path, encoding='utf-8').read()
    i = src.index(start)
    j = src.index(end, i)
    return src[i:j]


# ── 1. LLR 누수 ──
print('1. LLR 참조군 · 혼합군')
# 실제 명단은 사전등록 26-2 에 있다(공여자 이름이라 저장소에 넣지 않는다). 여기는 같은 정렬
#   성질을 가진 가짜 이름이다 — 대문자가 소문자보다 앞, 'A-3' 이 'A-24' 뒤(문자열 순).
#   기본 로케일 정렬(대소문자 무시)을 쓰면 배정이 달라져 이 검사가 실패한다.
normals = ['A-1', 'A-16', 'A-19', 'A-20', 'A-21', 'A-23', 'A-24',
           'A-3', 'A-30', 'S06', 'S07', 'ha', 'hy', 'na', 'un']
REF8 = ['A-1', 'A-19', 'A-21', 'A-24', 'A-30', 'S07', 'hy', 'un']
MIX7 = ['A-16', 'A-20', 'A-23', 'A-3', 'S06', 'ha', 'na']
ref, mix = M.split_ref_mix(reversed(normals))          # 순서를 뒤집어 넣어도 같아야 한다
ok(ref == REF8, '정상 참조군 8명 — 26-2 와 같은 규칙')
ok(mix == MIX7, '정상 혼합군 7명 — 26-2 와 같은 규칙')
lines = ['SNU%d' % i for i in range(1, 11)]
ref, mix = M.split_ref_mix(lines)
ok(ref == ['SNU1', 'SNU2', 'SNU4', 'SNU6', 'SNU8'] and mix == ['SNU10', 'SNU3', 'SNU5', 'SNU7', 'SNU9'],
   '세포주 참조군·혼합군 = 사전등록 26-2')

pool = block_src(os.path.join(TREE, '1_dmr/steps/11_pool.py'),
                 "if int(getattr(_cfg, 'LLR_HELDOUT', 0)):", "\nP = pd.read_csv(PANEL)")
tmp = tempfile.mkdtemp()
rows = [(n, 'Normal', n) for n in normals] + \
       [('SNU%d-%d' % (i, r), 'GBM', 'SNU%d' % i) for i in range(1, 11) for r in (1, 2)]
pd.DataFrame(rows, columns=['sample', 'type', 'group']).to_csv(os.path.join(tmp, 'samples.csv'), index=False)


class _C:
    LLR_HELDOUT = 1
    GEN = '15'
    split_ref_mix = staticmethod(M.split_ref_mix)


def run_pool(use):
    ns = {'_cfg': _C, '_os': os, 'pd': pd, 'sys': sys, 'USE': list(use), '_M': '/없음'}
    os.environ['CAND_DIR'] = tmp
    exec(pool, ns)
    return ns['USE']


got = run_pool(normals)
ok(got == MIX7, '11_pool 정상 목록 → 혼합군 7명만')
got = run_pool(['SNU%d-%d' % (i, r) for i in range(1, 11) for r in (1, 2)])
ok(sorted({g.split('-')[0] for g in got}) == ['SNU10', 'SNU3', 'SNU5', 'SNU7', 'SNU9'] and len(got) == 10,
   '11_pool 세포주 목록 → 혼합군 5주의 반복 10개만 (세포주 단위로 갈림)')
try:
    run_pool(['없는검체'])
    ok(False, '11_pool 모르는 검체에서 멈춘다')
except SystemExit:
    ok(True, '11_pool 모르는 검체에서 멈춘다')
import shutil  # noqa: E402
shutil.rmtree(tmp, ignore_errors=True)

# ── 2. 모델 고정 ──
print('2. 모델 고정')
from sklearn.ensemble import RandomForestClassifier  # noqa: E402
from sklearn.linear_model import LogisticRegression  # noqa: E402
from sklearn.svm import SVC  # noqa: E402
from xgboost import XGBClassifier  # noqa: E402
base = {'LogisticRegression': LogisticRegression, 'RandomForestClassifier': RandomForestClassifier,
        'SVC': SVC, 'XGBClassifier': XGBClassifier}


class _C8:
    SEED, N_JOBS, SVM_PROBABILITY = 42, 1, False
    FIXED_MODEL = ''


f08 = func_src(os.path.join(TREE, '2_train/steps/08_ML_classifier.py'), 'make_models')
f09 = func_src(os.path.join(TREE, '2_train/steps/09_save_best_j.py'), 'models')
for fixed, want in (('', 4), ('Logistic Regression', 1)):
    ns = dict(base, C=_C8, SEED=42, FIXED=fixed)
    _C8.FIXED_MODEL = fixed
    exec(f08, ns); exec(f09, ns)
    m8, m9 = ns['make_models'](), ns['models']()
    ok(len(m8) == want and len(m9) == want and (not fixed or list(m8) == list(m9) == [fixed]),
       '08·09 FIXED_MODEL=%r → 모델 %d개' % (fixed, want))
for name, fn in (('08', 'make_models'), ('09', 'models')):
    ns = dict(base, C=_C8, SEED=42, FIXED='Ridge')
    _C8.FIXED_MODEL = 'Ridge'
    exec(f08 if name == '08' else f09, ns)
    try:
        ns[fn]()
        ok(False, '%s 틀린 이름에서 멈춘다' % name)
    except SystemExit:
        ok(True, '%s 틀린 이름에서 멈춘다' % name)

# ── 3. 묶음 평균 채우기 ──
print('3. 채우기 · 블록 버리기')
ft = func_src(os.path.join(TREE, '2_train/steps/05_Feature_Matrix.py'), 'to_frame')
rows_n = {'Normal%d' % i: {'1_10': 0.1 * i, '1_20': np.nan if i == 1 else 0.5, '1_30': np.nan}
          for i in range(1, 6)}
rows_g = {'GBM%d' % i: {'1_10': 0.9, '1_20': 0.6, '1_30': np.nan if i > 1 else 0.7}
          for i in range(1, 6)}
for val in (False, True):
    ns = {'pd': pd, 'np': np, 'MIN_SAMPLE_COVERAGE': 0.8, '_VAL': val, 'COLSET_FROM': 'all'}
    exec(ft, ns)
    df = ns['to_frame'](rows_n, rows_g, 'test')
    X = df.drop(columns=['Type'])
    if val:
        ok(list(X.columns) == ['1_10', '1_20', '1_30'], '검증: 블록을 버리지 않는다 (3칸 그대로)')
        ok(X.isna().sum().sum() == 1 + 9, '검증: 빈칸을 채우지 않는다 (빈칸 10개 그대로)')
    else:
        ok(list(X.columns) == ['1_10', '1_20'], '학습: 커버율 80% 미만 블록(1_30)은 버린다')
        ok(bool(np.isnan(X.loc['Normal1', '1_20'])), '학습: 빈칸을 묶음 평균으로 채우지 않는다')
    ok(df['Type'].tolist() == [0] * 5 + [1] * 5, '라벨 열은 그대로')

# 3-2 칸을 정상 행만으로 고른다 (COLSET_FROM=normal · 사전등록 71-2 · 86)
#   블록 1_40 은 정상 다섯 행에서는 다 덮이고 암 다섯 행에서는 다 빈다.
rn = {k: {**v, '1_40': 0.3} for k, v in rows_n.items()}
rg = {k: {**v, '1_40': np.nan} for k, v in rows_g.items()}
for mode, want in (('all', ['1_10', '1_20']), ('normal', ['1_10', '1_20', '1_40'])):
    ns = {'pd': pd, 'np': np, 'MIN_SAMPLE_COVERAGE': 0.8, '_VAL': False, 'COLSET_FROM': mode}
    exec(ft, ns)
    X = ns['to_frame'](rn, rg, 'test').drop(columns=['Type'])
    ok(list(X.columns) == want, '%s: 남는 블록 %s' % (mode, want))
ns = {'pd': pd, 'np': np, 'MIN_SAMPLE_COVERAGE': 0.8, '_VAL': False, 'COLSET_FROM': 'normal'}
exec(ft, ns)
try:
    ns['to_frame']({}, rg, 'test'); ok(False, 'normal: 정상 행 0개면 멈춘다')
except SystemExit:
    ok(True, 'normal: 정상 행 0개면 멈춘다')

# ── 4. 다음 판 준비 (사전등록 65) ──
print('4. 다음 판 준비')

# 4-1 검증 판별을 환경변수로도 (65 라)
src05 = io.open(os.path.join(TREE, '2_train/steps/05_Feature_Matrix.py'), encoding='utf-8').read()
line = [l for l in src05.splitlines() if l.startswith('_VAL = ')][0]
for env, want in ((None, False), ('1', True)):
    if env is None:
        os.environ.pop('METH_VALIDATE', None)
    else:
        os.environ['METH_VALIDATE'] = env
    ns = {'os': os, '__file__': '/x/05_Feature_Matrix.py'}
    exec(line, ns)
    ok(ns['_VAL'] is want, '05: 파일 이름이 val 이 아니어도 METH_VALIDATE=%s 이면 _VAL=%s' % (env, want))
os.environ.pop('METH_VALIDATE', None)

# 4-2 05_llrref 의 held-out 분할은 --use 와 무관하게 samples.csv 전체로 (65 다)
blk = block_src(os.path.join(TREE, '1_dmr/steps/05_llrref.py'), "M = pd.read_csv(CAND + '/samples.csv')", 'import math')
with tempfile.TemporaryDirectory() as td:
    rows = [(n, n, 'Normal') for n in normals] + [('C%d-%d' % (g, r), 'C%d' % g, 'SCLC') for g in range(1, 7) for r in (1, 2)]
    pd.DataFrame(rows, columns=['sample', 'group', 'type']).to_csv(td + '/samples.csv', index=False)
    # 부분집합을 다시 가르는 방식과 전체로 가르는 방식이 다른 답을 내는 명단이어야 한다 — 우연히 같으면 시험이 못 잡는다
    use = set(normals[1:7]) | {'C2-1', 'C3-1', 'C4-1', 'C5-1'}
    class _Cf:
        LLR_HELDOUT = 1; CANCER_LABEL = 'SCLC'
        split_ref_mix = staticmethod(M.split_ref_mix)
    ns = {'pd': pd, 'sys': sys, 'CAND': td, 'USE': use, '_cfg': _Cf}
    exec(blk, ns)
    rn_full, _ = M.split_ref_mix(normals)
    rg_full, _ = M.split_ref_mix(['C%d' % g for g in range(1, 7)])
    ok(sorted(ns['norm']) == sorted(set(use) & set(rn_full)), '05: --use 를 줘도 정상 참조군은 전체 분할의 참조군 ∩ use')
    ok(sorted(ns['lines']) == sorted({'C2', 'C3', 'C4', 'C5'} & set(rg_full)), '05: 암 참조군도 전체 분할 기준')
    _old_n, _ = M.split_ref_mix(sorted(set(normals[1:7]))); _old_g, _ = M.split_ref_mix(['C2', 'C3', 'C4', 'C5'])
    ok(sorted(_old_n) != sorted(set(use) & set(rn_full)) and sorted(_old_g) != sorted({'C2', 'C3', 'C4', 'C5'} & set(rg_full)),
       '05: (시험 자체) 옛 분할과 새 분할이 이 명단에서 다르다')

# 4-3 07 의 llr 블록 집합은 LLR_KEEP_FROM(학습) 에서, 비율마다 읽는 표는 검체 것 (54-2)
blk7 = block_src(os.path.join(TREE, '2_train/steps/07_llrfeat.py'),
                 "    means = sorted(glob.glob(out_dir + '/*_mean_%s.csv' % C.VERSION),", '    keep = set()')
blk7 = '\n'.join(l[4:] if l.startswith('    ') else l for l in blk7.splitlines())
import glob as _glob, re as _re  # noqa: E402
with tempfile.TemporaryDirectory() as td:
    sd, tdir = os.path.join(td, 'sample'), os.path.join(td, 'train')
    os.makedirs(sd); os.makedirs(tdir)
    pd.DataFrame({'1_10': [0.1], '1_20': [0.2], '1_30': [0.3]}).to_csv(sd + '/5_mean_V.csv')
    pd.DataFrame({'1_10': [0.1], '1_20': [0.2]}).to_csv(sd + '/25_mean_V.csv')
    pd.DataFrame({'1_10': [0.1], '1_20': [0.2]}).to_csv(tdir + '/5_mean_V.csv')
    class _C7:
        VERSION = 'V'
    for kf, want_cols in ((None, ['1_10', '1_20', '1_30']), (tdir, ['1_10', '1_20'])):
        if kf is None:
            os.environ.pop('LLR_KEEP_FROM', None)
        else:
            os.environ['LLR_KEEP_FROM'] = kf
        ns = {'glob': _glob, 're': _re, 'os': os, 'sys': sys, 'pd': pd, 'C': _C7, 'out_dir': sd}
        exec(blk7, ns)
        got = list(ns['ref'].columns)
        ok(got == want_cols, '07: LLR_KEEP_FROM=%s → 블록 %s' % ('학습' if kf else '없음', want_cols))
        ok([os.path.basename(m) for m in ns['means']] == ['5_mean_V.csv', '25_mean_V.csv'],
           '07: 비율마다 읽는 means 는 검체 것 그대로 (%s)' % ('학습' if kf else '없음'))
    os.environ.pop('LLR_KEEP_FROM', None)

# 4-4 08 의 특이도 문턱용 접기 예측은 접기 안에서 채운다 (65 가) — 학습 반복 안이라 소스로 본다
src08 = io.open(os.path.join(TREE, '2_train/steps/08_ML_classifier.py'), encoding='utf-8').read()
seg = src08[src08.index('oof = cross_val_predict('):src08.index('oof = np.asarray(oof')]
ok('_Pipe(' in seg and '_Imputer(' in seg and 'Xtr_raw' in seg, '08: 문턱용 접기 예측이 Pipeline(Imputer) · 채우기 전 Xtr 를 쓴다')

# 4-5 검증 실행기가 두 변수를 내보낸다
sh = io.open(os.path.join(TREE, '3_validate/steps/run_one_sample.sh'), encoding='utf-8').read()
ok('export LLR_KEEP_FROM=' in sh and 'export METH_VALIDATE=1' in sh, 'run_one_sample.sh: LLR_KEEP_FROM · METH_VALIDATE 를 내보낸다')

# 4-6 08 · 09 의 로지스틱 회귀 설정이 같다
_C8.FIXED_MODEL = ''
ns8 = dict(base, C=_C8); exec(f08, ns8)
ns9 = dict(base, C=_C8, SEED=42, FIXED=''); exec(f09, ns9)
lr8, lr9 = ns8['make_models']()['Logistic Regression'], ns9['models']()['Logistic Regression']
ok(lr8.get_params()['max_iter'] == lr9.get_params()['max_iter'], '08 · 09 LogisticRegression max_iter 같음 (%s · %s)'
   % (lr8.get_params()['max_iter'], lr9.get_params()['max_iter']))

print()
# ── 5. 점수 단위 맞추기 (SCORE_UNIT · meth_config.unify_score_unit) ──
print('5. 점수 단위 맞추기')
_D = pd.DataFrame({'feature_set': ['a', 'a', 'a', 'b', 'b', 'c', 'd', 'd'],
                   'score_type': ['predict_proba', 'predict_proba+df(포화)', 'predict_proba', 'predict_proba', 'predict_proba',
                                  'decision_function', 'predict_proba', 'predict_proba'],
                   'y_score': [0.9, 7.0, 0.8, 1.0, 1.0, -1.2, 0.6, 0.7], 'y_df': [2.0, 7.0, 1.5, 40.0, 41.0, -1.2, np.nan, np.nan],
                   'thr_spec95': [0.7, 0.85, 0.7, 0.6, 0.6, 0.1, 0.5, 0.5]})
_keep = _D.copy()
_old = M.SCORE_UNIT
try:
    M.SCORE_UNIT = 'as_is'
    _o, _n = M.unify_score_unit(_D)
    ok(_o.equals(_keep) and _n == 0, 'as_is: 아무것도 바꾸지 않는다')
    M.SCORE_UNIT = 'consistent'
    _o, _n = M.unify_score_unit(_D)
    ok(_D.equals(_keep), 'consistent: 들어온 표를 제자리에서 바꾸지 않는다')
    ok(_n == 2 and _o.y_score.tolist()[:6] == [2.0, 7.0, 1.5, 40.0, 41.0, -1.2],
       'consistent: 전환 있는 조합(a)과 부분 포화 동률 조합(b)을 y_df 로 · 동률이 풀린다')
    ok(_o.y_score.tolist()[6:] == [0.6, 0.7] and _o.score_type.tolist()[6:] == ['predict_proba'] * 2,
       'consistent: decision_function 이 없는 모형(d)은 확률 그대로')
    ok(_o.score_type.tolist()[:5] == ['decision_function(통일)'] * 5 and _o.score_type.tolist()[5] == 'decision_function',
       'consistent: 바꾼 행만 표기가 바뀐다 · SVM 의 decision_function 은 전환으로 보지 않는다')
    ok(bool(np.isnan(_o.thr_spec95[0])) and _o.thr_spec95[1] == 0.85 and bool(np.isnan(_o.thr_spec95[3]))
       and _o.thr_spec95[5] == 0.1 and _o.thr_spec95[6] == 0.5,
       'consistent: 확률 단위였던 행의 문턱만 비운다')
    try:
        M.unify_score_unit(_D.drop(columns=['y_df'])); ok(False, 'consistent: y_df 열이 없으면 멈춘다')
    except SystemExit:
        ok(True, 'consistent: y_df 열이 없으면 멈춘다')
    _bad = _D.copy(); _bad.loc[0, 'y_df'] = np.nan
    try:
        M.unify_score_unit(_bad); ok(False, 'consistent: 맞출 묶음에 빈 y_df 가 있으면 멈춘다')
    except SystemExit:
        ok(True, 'consistent: 맞출 묶음에 빈 y_df 가 있으면 멈춘다')
finally:
    M.SCORE_UNIT = _old

print()
# ── 6. 08 · 09 의 수렴 경고 세기 (CONVERGENCE_GATE · meth_config.count_convergence) ──
print('6. 수렴 경고 세기')
import contextlib, warnings  # noqa: E401,E402
from sklearn.exceptions import ConvergenceWarning  # noqa: E402
_err = io.StringIO()
with warnings.catch_warnings(), contextlib.redirect_stderr(_err):
    _n = M.count_convergence()
    for _ in range(3):
        warnings.warn('같은 자리', ConvergenceWarning)   # 한 줄에서 세 번 — always 가 없으면 한 번만 온다
    warnings.warn('다른 경고', UserWarning)
ok(_n[0] == 3, '같은 자리 수렴 경고도 다 센다 (%d)' % _n[0])
ok(_err.getvalue().count('ConvergenceWarning: ') == 1, '화면에는 수렴 경고를 처음 하나만 낸다')
ok('다른 경고' in _err.getvalue(), '다른 경고는 그대로 지나간다')
_old = M.CONVERGENCE_GATE
try:
    for _g, _want in (('off', 0), ('on', 1)):
        M.CONVERGENCE_GATE = _g; _f = []
        with contextlib.redirect_stdout(io.StringIO()):
            M.convergence_check([2], _f, 't'); M.convergence_check([0], _f, 't')
        ok(len(_f) == _want, 'CONVERGENCE_GATE=%s: 경고 2건이면 실패 %d건 · 0건이면 그대로' % (_g, _want))
finally:
    M.CONVERGENCE_GATE = _old
for _st, _call in (('08_ML_classifier.py', 'C.convergence_check(_conv, failures'), ('09_save_best_j.py', '_cfg.convergence_check(_conv, failures')):
    _s = io.open(os.path.join(TREE, '2_train/steps', _st), encoding='utf-8').read()
    ok(_call in _s and _s.index(_call) < _s.rindex('if failures:'), '%s: 실패 판정 앞에서 수렴 경고를 넣는다' % _st[:2])

print('결과: %s' % ('전부 통과' if not bad else '실패 %d개' % len(bad)))
sys.exit(1 if bad else 0)
