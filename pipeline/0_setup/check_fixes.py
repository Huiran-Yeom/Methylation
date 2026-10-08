# -*- coding: utf-8 -*-
"""2026-10-08 결함 수정 세 가지가 실제로 그렇게 도는지 태워 본다. 자료·서버 없이 돈다.

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


# ── 1. LLR 누수 ────────────────────────────────────────────
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
tmp = os.path.join(HERE, '_check_cand')
os.makedirs(tmp, exist_ok=True)
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

# ── 2. 모델 고정 ───────────────────────────────────────────
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

# ── 3. 묶음 평균 채우기 ────────────────────────────────────
print('3. 채우기 · 블록 버리기')
ft = func_src(os.path.join(TREE, '2_train/steps/05_Feature_Matrix.py'), 'to_frame')
rows_n = {'Normal%d' % i: {'1_10': 0.1 * i, '1_20': np.nan if i == 1 else 0.5, '1_30': np.nan}
          for i in range(1, 6)}
rows_g = {'GBM%d' % i: {'1_10': 0.9, '1_20': 0.6, '1_30': np.nan if i > 1 else 0.7}
          for i in range(1, 6)}
for val in (False, True):
    ns = {'pd': pd, 'np': np, 'MIN_SAMPLE_COVERAGE': 0.8, '_VAL': val}
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

print()
print('결과: %s' % ('전부 통과' if not bad else '실패 %d개' % len(bad)))
sys.exit(1 if bad else 0)
