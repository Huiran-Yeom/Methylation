# -*- coding: utf-8 -*-
"""size_cv 04: 접기 AUC 로 방식(jsd, bl)마다 패널 크기 N 고르기 (1-SE)
  접기 AUC = 1차 조합 행, 비율 RATIO 양성 행 대 음성 행 전부, 동률 0.5 (primary_metric a1 과 같은 정의)
  칸 검사: 채점 검체 = 빼 둔 단위, 양성 행 = 음성 행 = 단위 수 × NREP, model_depth 와 score_type 하나씩
  결정: 접기 평균이 가장 높은 N* (같으면 작은 N) → 평균(N*) - SE(N*) 이상인 N 중 가장 작은 N
       SE = 접기 SD / √NFOLD (ddof 1). 접기가 다 없는 N 은 고르지 않음
사용: python 04_pick.py --root <CV_ROOT> --gen <GEN> --sizes "200 300 …" --nfold 5 --ratio 0.005 --nrep 10 --depth 5000
      [--fold-root <k 를 받는 형식, 기본 {root}/fold{k}/root>] [--cohort <형식, 기본 cv_f{k}>] [--dry]
"""
if __name__ != '__main__':
    raise ImportError(__file__ + ' 은(는) 스크립트다. import 하지 않는다')
import glob
import os
import sys

import numpy as np
import pandas as pd
from sklearn.metrics import roc_auc_score

a = sys.argv


def arg(k, d=None):
    return a[a.index(k) + 1] if k in a else d


ROOT = arg('--root')
GEN = arg('--gen')
SIZES = [int(x) for x in arg('--sizes').split()]
NFOLD = int(arg('--nfold', 5))
RATIO = float(arg('--ratio', 0.005))
NREP = int(arg('--nrep', 10))
DEPTH = int(arg('--depth', 5000))
FROOT = arg('--fold-root', '{root}/fold{k}/root')
COH = arg('--cohort', 'cv_f{k}')
DRY = '--dry' in a
COMBO = 'mean|entropy|jsd|llr|pdr'
DG = '%dk' % (DEPTH // 1000) if DEPTH % 1000 == 0 else str(DEPTH)


def fold_auc(k, p):
    v = 'j%s%s_%s_uni' % (GEN, p, DG)
    root = FROOT.format(root=ROOT, k=k)
    units = [x.strip() for x in open('%s/fold%d/heldout_cancer_units.txt' % (ROOT, k)) if x.strip()]
    vd = '%s/%s_val' % (root, COH.format(k=k))
    fs = sorted(glob.glob('%s/*/results/%s/step04_ML_classifier/y_prob_all_%s_all_combos.csv' % (vd, v, v)))
    got = sorted(os.path.relpath(f, vd).replace('\\', '/').split('/')[0] for f in fs)
    if got != sorted(units):
        return None, '채점 검체 %d/%d' % (len(got), len(units))
    d = pd.concat([pd.read_csv(f) for f in fs], ignore_index=True)
    d = d[d.feature_set == COMBO]
    if d.model_depth.nunique() != 1 or d.score_type.nunique() != 1:
        return None, '점수 단위 섞임'
    pos = d[(d.y_true == 1) & np.isclose(d.ratio.astype(float), RATIO)]
    neg = d[d.y_true == 0]
    want = len(units) * NREP
    if len(pos) != want or len(neg) != want:
        return None, '행 양성 %d · 음성 %d (기대 %d)' % (len(pos), len(neg), want)
    if DRY:
        return np.nan, '행 %d+%d' % (len(pos), len(neg))
    y = np.r_[np.zeros(len(neg)), np.ones(len(pos))]
    return float(roc_auc_score(y, np.r_[neg.y_score.values, pos.y_score.values])), '행 %d+%d' % (len(pos), len(neg))


def blocks(p):
    out = []
    for k in range(1, NFOLD + 1):
        f = '%s/results/dmr/j_panel_dmr%s_panels/%s/DMR_confirmed_%s.csv' % (FROOT.format(root=ROOT, k=k), GEN, p, p)
        out.append(len(pd.read_csv(f)) if os.path.exists(f) else -1)
    return out


def pick(ok):
    """ok = [(N, 평균, SE)] → (N*, 평균*, SE*, 문턱, 고른 N)"""
    nstar, mstar, sestar = max(ok, key=lambda t: (t[1], -t[0]))
    thr = mstar - sestar
    return nstar, mstar, sestar, thr, min(n for n, m, _ in ok if m >= thr)


print('== 패널 크기 고르기 (1-SE) · 깊이 %s · 비율 %g · 접기 %d%s' % (DG, RATIO, NFOLD, ' · DRY' if DRY else ''))
for meth in ('jsd', 'bl'):
    print('\n-- %s' % meth)
    print('%-6s %-*s %7s %7s %7s  %s' % ('N', 8 * NFOLD, '접기 AUC', '평균', 'SD', 'SE', '비고'))
    ok = []
    for N in SIZES:
        p = '%s%d' % (meth, N)
        res = [fold_auc(k, p) for k in range(1, NFOLD + 1)]
        aucs = [r[0] for r in res]
        b = blocks(p)
        full = all(x is not None for x in aucs)
        cell = ' '.join('%.4f ' % x if (x is not None and not np.isnan(x)) else ('  -    ' if x is None else ' dry   ') for x in aucs)
        note = ('N 미달 %s ' % b if any(0 <= x < N for x in b) else '')
        if not full:
            note += '계산 불가: ' + '; '.join('f%d %s' % (i + 1, r[1]) for i, r in enumerate(res) if r[0] is None)
        if full and not DRY:
            m, sd = float(np.mean(aucs)), float(np.std(aucs, ddof=1))
            se = sd / np.sqrt(NFOLD)
            ok.append((N, m, se))
            print('%-6d %s %7.4f %7.4f %7.4f  %s' % (N, cell, m, sd, se, note))
        else:
            print('%-6d %s %7s %7s %7s  %s' % (N, cell, '', '', '', note))
    if DRY:
        continue
    if not ok:
        print('결정: 고를 N 없음 (접기가 다 있는 N 이 없음)')
        continue
    nstar, mstar, sestar, thr, n = pick(ok)
    print('N* = %d (평균 %.4f · SE %.4f) · 문턱 %.4f → 고른 N = %d' % (nstar, mstar, sestar, thr, n))
    if n == min(SIZES):
        print('  고른 N 이 후보 중 가장 작음. 더 작은 N 은 안 봄')
