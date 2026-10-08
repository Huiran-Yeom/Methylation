# -*- coding: utf-8 -*-
"""폐 Patient 7판 대 10판 비교. 사전등록 추기 44-3 · 44-5 · 44-7. 2026-10-07

내는 것
  1차 지표 (primary_metric.py 계산법: 양성 0.5% · 음성 전 비율 합침 · 검체 단위)
    7판 · 10판 전체 · 10판 반복 1~7
  (나) 비율별 검체 AUC 567칸의 Δ(10판 − 7판) 부호
쓰는 법
  python lung10_compare.py            → 세 갈래 다
  python lung10_compare.py --only7    → 7판만 (파이프라인 값 재현 검사)
"""
import glob, sys
import numpy as np, pandas as pd
R = '<METH_ROOT>/lung'
PANELS = ['j16bl200', 'j16jsd200', 'j16jsdb200']
DEPTHS = ['5k', '50k', '100k']
PRIMARY = 'mean|entropy|jsd|llr|pdr'
RATIOS = [0.001, 0.005, 0.01, 0.02, 0.025, 0.03, 0.05]

def auc(y, s):
    r = pd.Series(s).rank().to_numpy(); p = int((y == 1).sum()); n = int((y == 0).sum())
    return np.nan if p == 0 or n == 0 else (r[y == 1].sum() - p * (p + 1) / 2.0) / (p * n)

def load(tree, ver, maxrep=None):
    fs = sorted(glob.glob('%s/%s/*/results/%s/step04_ML_classifier/y_prob_all_%s_all_combos.csv' % (R, tree, ver, ver)))
    if not fs:
        return None, 0
    d = pd.concat([pd.read_csv(f, usecols=['sample', 'feature_set', 'model', 'ratio', 'rep', 'y_true', 'y_score']) for f in fs], ignore_index=True)
    if maxrep:
        d = d[d.rep.str.replace(r'^\D+', '', regex=True).astype(int) <= maxrep]
    return d, len(fs)

def primary(d):
    x = d[d.feature_set == PRIMARY]
    pos = x[(x.y_true == 1) & np.isclose(x.ratio, 0.005)]
    neg = x[x.y_true == 0]
    return auc(np.r_[np.ones(len(pos)), np.zeros(len(neg))], np.r_[pos.y_score, neg.y_score]), len(pos), len(neg)

def per_ratio(d):
    out = {}
    for (fs, r), g in d.groupby(['feature_set', 'ratio']):
        out[(fs, round(r, 4))] = auc(g.y_true.to_numpy(), g.y_score.to_numpy())
    return out

ARMS = [('7판', 'Patient_Test_val', None)]
if '--only7' not in sys.argv:
    ARMS += [('10판 전체', 'Patient_Test10_val', None), ('10판 반복1~7', 'Patient_Test10_val', 7)]

print('== 1차 지표 (5k · 0.5% · 검체 · 음성 전 비율 합침)')
cache = {}
for name, tree, mr in ARMS:
    for p in PANELS:
        d, nf = load(tree, p + '_5k_uni', mr)
        if d is None:
            print('  %-12s %-11s 파일 없음' % (name, p)); continue
        cache[(name, p, '5k')] = d
        a, np_, nn = primary(d)
        print('  %-12s %-11s AUC %.4f  (양성 %d행 · 음성 %d행 · 검체 %d)' % (name, p, a, np_, nn, nf))

if '--only7' in sys.argv:
    sys.exit(0)

print()
print('== (나) 비율별 검체 AUC 의 Δ 부호 — 실제 3판 × 깊이 3 × 조합 × 비율 7')
for left, right in (('10판 전체', '7판'), ('10판 반복1~7', '7판'), ('10판 전체', '10판 반복1~7')):
    pos = neg = zero = 0; ds = []
    for p in PANELS:
        for dp in DEPTHS:
            v = p + '_' + dp + '_uni'
            arms = {}
            for name, tree, mr in ARMS:
                if (name, p, dp) not in cache:
                    cache[(name, p, dp)] = load(tree, v, mr)[0]
                arms[name] = cache[(name, p, dp)]
            if arms[left] is None or arms[right] is None:
                continue
            a, b = per_ratio(arms[left]), per_ratio(arms[right])
            for k in a:
                if k in b and np.isfinite(a[k]) and np.isfinite(b[k]):
                    x = a[k] - b[k]; ds.append(x)
                    pos += x > 0; neg += x < 0; zero += x == 0
    n = pos + neg + zero
    if n:
        print('  %-12s − %-12s 칸 %d · 양 %d · 음 %d · 0 %d · 같은 부호 최대 %.1f%% · |Δ| 중앙 %.4f · 최대 %.4f'
              % (left, right, n, pos, neg, zero, 100.0 * max(pos, neg) / n, np.median(np.abs(ds)), np.max(np.abs(ds))))
print()
print('귀결 (추기 44-3 나): 같은 부호 75%% 초과(567칸이면 425칸 이상)면 결과를 쓰기 전에 섞기의 반복 의존부터 본다.')
