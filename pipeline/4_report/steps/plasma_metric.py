# -*- coding: utf-8 -*-
"""환자 혈장 코호트를 희석 없이 채점한 결과의 사람 단위 판정 (탐색 · 1차 지표 아님)

  혈장은 이미 몸 안에서 희석된 검체라, 다시 정상 배경에 0.1~5% 로 섞으면 종양 리드가
  검출 한계 아래로 내려간다. 그래서 혈장 코호트는 RATIOS=0.0, 1.0 으로 돌려 사람마다
  자기 리드만으로 채점하고, 같은 플랫폼 정상 공여자도 같은 방식으로 채점해 음성으로 쓴다.
  In silico(조직·세포주) 코호트는 지금처럼 섞는다. 설정은
  3_validate/steps/cohorts/patient_plasma.conf.example 을 보십시오.

  판정 순서 (사전등록 추기 80 · 하나라도 아니면 「말하지 않는다」)
    ① 위장 대조   음성 행(배경만) 점수로 같은 AUC 를 재서 순열 귀무 95% 안이어야 한다
    ② AUC ≥ 0.70
    ③ 라벨 순열 p ≤ .05 (한쪽)
    ④ 기준선 둘(사람별 총 리드 수 · 전역 평균 메틸화)의 양쪽 AUC 순열 p ≤ .05 이면 교란.
       교란이 아니면 패널 − 기준선 짝 차이 p ≤ .05 여야 한다
  관문 (판정 전에 멈추거나 그 판을 읽지 않음)
    사람 수가 코호트 설정(SAMPLES · NEG_SAMPLES)과 다름 · 사람마다 반복 수가 NREP 와 다름
    · 점수 단위 혼합(포화 전환 decision_function 이 일부 검체에만)

  보고만 하는 것 (판정에 쓰지 않음): llr 뺀 조합 AUC · 병기별 AUC · 음성 한 명 빼기 범위

  읽는 것
    <METH_ROOT>/<COHORT>_val/<사람>/results/<판>_<깊이>_uni/step04_ML_classifier/
        y_prob_all_<판>_<깊이>_uni_all_combos.csv       (ratio 1.0 · rep GBM* 행)
    --baseline  csv: sample, reads, gmeth   (plasma_baseline.py 가 만든다)
    --stage     csv: sample, stage          (선택 · 'Stage I'/'Stage II' 를 조기로 본다)

사용:  COHORT=<혈장코호트> python plasma_metric.py --baseline b.csv [--stage s.csv]
            [--depth 5k] [--nperm 10000] [--seed 80] j16bl200 j16jsd200 j16jsdb200
       판 순서가 순열 난수 순서다. 같은 수를 다시 내려면 같은 순서로 준다.
"""
# 이 파일은 스크립트다. 다른 코드가 import 하면 본체가 그대로 돈다.
if __name__ != '__main__':
    raise ImportError(__file__ + ' 은(는) 스크립트다 — import 하지 않는다')

import os as _os, sys as _sys
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
import meth_config as _cfg

import argparse
import glob
import io
import os
import sys

import numpy as np
import pandas as pd

PR = 'mean|entropy|jsd|llr|pdr'
NOLLR = 'mean|entropy|jsd|pdr'

ap = argparse.ArgumentParser()
ap.add_argument('panels', nargs='+')
ap.add_argument('--baseline', required=True)
ap.add_argument('--stage')
ap.add_argument('--conf')
ap.add_argument('--depth', default='5k')
ap.add_argument('--nperm', type=int, default=10000)
ap.add_argument('--seed', type=int, default=80)
a = ap.parse_args()

COHORT = os.environ.get('COHORT') or sys.exit('COHORT 를 주십시오')
R = '%s/%s_val' % (_cfg.METH_ROOT, COHORT)


def load_conf():
    cands = [a.conf] if a.conf else [os.path.expanduser('~/cohorts/%s.conf' % COHORT),
                                     os.path.join(_h, '3_validate', 'steps', 'cohorts', '%s.conf' % COHORT)]
    p = next((c for c in cands if c and os.path.exists(c)), None) or sys.exit('코호트 정의 없음: ' + ' 또는 '.join(cands))
    print('코호트 정의  %s' % os.path.abspath(p))
    d = {}
    for ln in io.open(p, encoding='utf-8'):
        ln = ln.strip()
        if ln and not ln.startswith('#') and '=' in ln:
            k, v = ln.split('=', 1)
            d[k.strip()] = v.strip()
    return d


conf = load_conf()
SAMPLES = conf.get('SAMPLES', '').split()
NEG = set(conf.get('NEG_SAMPLES', '').split())
if not NEG or not NEG <= set(SAMPLES):
    sys.exit('NEG_SAMPLES 가 비었거나 SAMPLES 에 없는 이름이 있다 — 혈장 코호트 설정이 아니다')
dk = a.depth.upper()
NREP = int(conf.get('NREP_' + dk, conf.get('NREP', 10)))
if [float(x) for x in conf.get('RATIOS', '').replace(' ', '').split(',') if x] != [0.0, 1.0]:
    sys.exit('RATIOS 가 0.0, 1.0 이 아니다 — 희석 없음 채점이 아니다')
base = pd.read_csv(a.baseline).set_index('sample')
stage = pd.read_csv(a.stage).set_index('sample')['stage'] if a.stage else None
rng = np.random.default_rng(a.seed)


def auc(s, y):
    p, n = s[y == 1], s[y == 0]
    return ((p[:, None] > n[None, :]).sum() + 0.5 * (p[:, None] == n[None, :]).sum()) / (p.size * n.size)


def bounds(s, y):
    p, n = s[y == 1], s[y == 0]; t = p.size * n.size
    return (p[:, None] > n[None, :]).sum() / t, (p[:, None] >= n[None, :]).sum() / t


def scores(v, combo, rows):
    # 단위 맞추기(SCORE_UNIT=consistent)는 코호트 전체에서 정해야 해서 검체 파일을 한꺼번에 읽는다.
    fs = glob.glob('%s/*/results/%s/step04_ML_classifier/y_prob_all_%s_all_combos.csv' % (R, v, v))
    if not fs:
        return {}
    names = [f.split('/')[len(R.split('/'))] for f in fs]
    D = pd.concat([pd.read_csv(f).assign(_s=n) for f, n in zip(fs, names)], ignore_index=True)
    D, nu = _cfg.unify_score_unit(D)
    if nu and rows == 'GBM':
        print('  [점수 단위] %s · 확률과 decision_function 이 섞인 조합 %d개를 decision_function 으로 맞췄다'
              ' (아래 «df 전환 검체» 는 맞춘 뒤라 0 이고, ≥.999 줄은 이 단위에서 뜻이 없다)' % (v, nu))
    out = {}
    for s in names:                       # 파일마다 한 칸 — 머리줄만 있는 파일도 빈 칸으로 남겨 반복 수 관문에 걸리게 한다
        d = D[D._s == s]
        d = d[(d.feature_set == combo) & (d.ratio.astype(float) == 1.0) & d.rep.astype(str).str.startswith(rows)]
        out[s] = (d.y_score.mean(), len(d), d.score_type.astype(str).str.contains('df').any())
    return out


res = {}
for pan in a.panels:
    v = '%s_%s_uni' % (pan, a.depth)
    sc = scores(v, PR, 'GBM'); sh = scores(v, PR, 'Normal')
    miss = sorted(set(SAMPLES) - set(sc))
    if miss or set(sc) - set(SAMPLES):
        # 리드가 모자란 사람은 섞기 단계에서 ratio 1.0 이 조용히 빠진다. 빠짐이 깊이와 엮이므로 멈춘다.
        sys.exit('%s: 채점 파일 없는 사람 %d명 %s · 설정에 없는 사람 %s' % (v, len(miss), miss[:10], sorted(set(sc) - set(SAMPLES))[:10]))
    S = sorted(sc)
    y = np.array([0 if s in NEG else 1 for s in S])
    x = np.array([sc[s][0] for s in S]); xs = np.array([sh[s][0] for s in S])
    nrep = [sc[s][1] for s in S]; ndf = sum(sc[s][2] for s in S)
    print('== %s · 사람 %d (환자 %d · 정상 %d) · 반복 %d~%d · df 전환 검체 %d'
          % (pan[3:] if pan[:1] == 'j' else pan, len(S), y.sum(), (1 - y).sum(), min(nrep), max(nrep), ndf))
    bad = [s for s in S if sc[s][1] != NREP]
    if bad:
        res[pan] = '판정 안 함 (반복 수가 %d 가 아닌 사람 %d명: %s)' % (NREP, len(bad), bad[:10])
        print('  ' + res[pan]); continue
    if np.isnan(x).any() or np.isnan(xs).any():
        res[pan] = '판정 안 함 (점수 빠짐)'; print('  ' + res[pan]); continue
    if 0 < ndf < len(S):
        res[pan] = '말하지 않는다 (단위 혼합 — df 전환 %d/%d)' % (ndf, len(S))
        print('  ' + res[pan]); continue
    perm = [rng.permutation(y) for _ in range(a.nperm)]
    o_s = auc(xs, y); nul_s = np.array([auc(xs, yy) for yy in perm]); lo, hi = np.percentile(nul_s, [2.5, 97.5])
    tied = np.unique(np.round(xs, 12)).size == 1
    sham = tied or (lo <= o_s <= hi)
    o = auc(x, y); nul = np.array([auc(x, yy) for yy in perm]); pv = (1 + (nul >= o).sum()) / (a.nperm + 1)
    blo, bhi = bounds(x, y)
    print('  위장 대조 AUC %.3f · 귀무 [%.3f, %.3f] → %s' % (o_s, lo, hi, '성립 안 함(동률)' if tied else ('통과' if sham else '실패')))
    print('  AUC %.4f · 순열 p %.4f · 귀무 [%.3f, %.3f] · 동률 [하한 %.3f, 상한 %.3f]' % (o, pv, *np.percentile(nul, [2.5, 97.5]), blo, bhi))
    print('  사람 점수 환자 [%.4f, %.4f] · 정상 [%.4f, %.4f] · ≥.999 환자 %d · 정상 %d · 서로 다른 값 %d/%d'
          % (x[y == 1].min(), x[y == 1].max(), x[y == 0].min(), x[y == 0].max(), (x[y == 1] >= .999).sum(), (x[y == 0] >= .999).sum(),
             np.unique(np.round(x, 9)).size, len(x)))
    c4 = True; conf_ = False
    for col, lab in (('reads', '총 리드 수'), ('gmeth', '전역 평균 메틸화')):
        xb = base.loc[S, col].values.astype(float)
        two = lambda yy: max(auc(xb, yy), 1 - auc(xb, yy))
        ob = two(y); nb = np.array([two(yy) for yy in perm]); pb = (1 + (nb >= ob).sum()) / (a.nperm + 1)
        dobs = o - ob; nd = np.array([auc(x, yy) - two(yy) for yy in perm]); pdiff = (1 + (nd >= dobs).sum()) / (a.nperm + 1)
        conf_ |= pb <= 0.05; c4 &= pdiff <= 0.05
        print('  기준선 %-10s 양쪽 AUC %.4f · p %.4f · 패널 − 기준선 %+.4f · 짝 차이 p %.4f' % (lab, ob, pb, dobs, pdiff))
    why = []
    if not sham: why.append('위장 대조 실패')
    if conf_: why.append('교란')
    if o < 0.70 or pv > 0.05: why.append('기준 미달')
    if not conf_ and not c4: why.append('기준선을 못 이김')
    res[pan] = '통과' if not why else '말하지 않는다 (%s)' % ' · '.join(why)
    print('  판정 %s: %s' % (pan, res[pan]))
    xn = scores(v, NOLLR, 'GBM'); xn = np.array([xn[s][0] for s in S])
    msg = '  보고만 · llr 뺀 AUC %.4f' % auc(xn, y)
    if stage is not None:
        early = np.array([str(stage.get(s, '')) in ('Stage I', 'Stage II') for s in S])
        msg += (' · 조기(환자 %d) AUC %.4f · 진행(환자 %d) AUC %.4f'
                % ((early & (y == 1)).sum(), auc(x[early | (y == 0)], y[early | (y == 0)]),
                   ((~early) & (y == 1)).sum(), auc(x[(~early) | (y == 0)], y[(~early) | (y == 0)])))
    print(msg)
    loo = [auc(x[np.arange(len(S)) != i], y[np.arange(len(S)) != i]) for i in np.where(y == 0)[0]]
    print('  보고만 · 음성 한 명 빼기 AUC [%.3f, %.3f]' % (min(loo), max(loo)))
print('\n판정 요약:', res)
