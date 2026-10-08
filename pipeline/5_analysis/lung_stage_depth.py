# -*- coding: utf-8 -*-
"""폐 Patient 병기별 보조 분석. 사전등록 추기 45 (20:51 박음). 2026-10-07
쓰는 법: python lung_stage.py <트리> [최대반복]
  예) Patient_Test_val · Patient_Test10_val · Patient_Test10_val 7
"""
import glob, os, sys
import numpy as np, pandas as pd
R = '<METH_ROOT>/lung'
PANELS = ['j16bl200', 'j16jsd200', 'j16jsdb200']
PRIMARY = 'mean|entropy|jsd|llr|pdr'
NB, SEED = 2000, 42
tree = sys.argv[1]; DEPTH = sys.argv[2]; maxrep = int(sys.argv[3]) if len(sys.argv) > 3 else None

st = pd.read_csv(os.path.expanduser('<WORK>/tmp/gse301266_stage.tsv'), sep='\t', header=None, names=['gsm', 'title', 'stage'])
STAGE = {'G301266_' + g: s for g, s in zip(st.gsm, st.stage)}
GROUP = {'조기 I·II': ('Stage I', 'Stage II'), '진행 III·IV': ('Stage III', 'Stage IV')}

def contrib(pos, neg):
    neg = np.sort(neg)
    lo = np.searchsorted(neg, pos, 'left'); hi = np.searchsorted(neg, pos, 'right')
    return (lo + 0.5 * (hi - lo)) / len(neg)          # 양성 한 줄의 AUC 기여

print('트리 %s · 반복 %s · 1차 지표 계산법 (%s · 0.5%% · 검체 · %s)' % (tree, maxrep or '전부', DEPTH, PRIMARY))
for p in PANELS:
    v = p + '_' + DEPTH + '_uni'
    fs = sorted(glob.glob('%s/%s/*/results/%s/step04_ML_classifier/y_prob_all_%s_all_combos.csv' % (R, tree, v, v)))
    d = pd.concat([pd.read_csv(f, usecols=['sample', 'feature_set', 'ratio', 'rep', 'y_true', 'y_score']) for f in fs], ignore_index=True)
    d = d[d.feature_set == PRIMARY]
    if maxrep:
        d = d[d.rep.str.replace(r'^\D+', '', regex=True).astype(int) <= maxrep]
    miss = sorted(set(d['sample']) - set(STAGE))
    if miss:
        sys.exit('병기를 모르는 검체: %s' % miss)
    d['stage'] = d['sample'].map(STAGE)
    neg_all = d[d.y_true == 0].y_score.to_numpy()
    pos = d[(d.y_true == 1) & np.isclose(d.ratio, 0.005)]
    print()
    print('== %s  (검체 %d · 음성 공통 %d행)' % (p, d['sample'].nunique(), len(neg_all)))
    rng = np.random.default_rng(SEED)
    for gname, stages in [('전체 48', None)] + list(GROUP.items()):
        g = pos if stages is None else pos[pos.stage.isin(stages)]
        c = contrib(g.y_score.to_numpy(), neg_all)
        a = c.mean()
        # 사람 단위 재표본 (양성만 · 음성 고정)
        people = g['sample'].to_numpy(); u = np.unique(people)
        per = np.array([c[people == s].sum() for s in u]); cnt = np.array([(people == s).sum() for s in u])
        bs = []
        for _ in range(NB):
            k = rng.integers(0, len(u), len(u))
            bs.append(per[k].sum() / cnt[k].sum())
        lo, hi = np.percentile(bs, [2.5, 97.5])
        # 민감도: 음성도 그 무리 폴더에서만
        if stages is None:
            sens = '—'
        else:
            ng = d[(d.y_true == 0) & d.stage.isin(stages)].y_score.to_numpy()
            sens = '%.4f (음성 %d행)' % (contrib(g.y_score.to_numpy(), ng).mean(), len(ng))
        print('  %-11s 사람 %2d · 양성 %3d행   AUC %.4f   95%% 구간 [%.3f, %.3f]   민감도 %s'
              % (gname, len(u), len(g), a, lo, hi, sens))
