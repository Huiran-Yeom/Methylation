# -*- coding: utf-8 -*-
"""폐 Patient 병기별 전체 격자 (탐색 · 오너 요청 · 추기 45-5). 2026-10-07
5장 표와 같은 규약: 비율마다 그 비율의 양성·음성. 검체 단위 = 행 그대로, 사람 단위 = 검체폴더별 평균.
검산: 전체 48 의 값이 <WORK>/표/04_폐_Patient/patient.txt 의 1차 조합 줄과 칸마다 같아야 한다."""
import glob, re, sys, io
import numpy as np, pandas as pd
R = '<METH_ROOT>/lung'; PR = 'mean|entropy|jsd|llr|pdr'; tree = sys.argv[1]
PANELS = ['j16bl200', 'j16jsd200', 'j16jsdb200']; DEPTHS = ['5k', '50k', '100k']
RAT = [0.001, 0.005, 0.01, 0.02, 0.025, 0.03, 0.05]
st = pd.read_csv('<WORK>/tmp/gse301266_stage.tsv', sep='\t', header=None, names=['g', 't', 's'])
S = {'G301266_' + g: s for g, s in zip(st.g, st.s)}
GR = [('전체', None), ('조기', ('Stage I', 'Stage II')), ('진행', ('Stage III', 'Stage IV'))]
def auc(y, s):
    r = pd.Series(s).rank().to_numpy(); p = int((y == 1).sum()); n = int((y == 0).sum())
    return (r[y == 1].sum() - p * (p + 1) / 2.0) / (p * n)
# 검산용: 표의 1차 조합 줄
ref = {}; cur = None
for ln in io.open('<WORK>/표/04_폐_Patient/patient.txt', encoding='utf-8'):
    m = re.match(r'==\s+(\S+)\s+리드\s+·\s+(\S+)\s+·', ln)
    if m: cur = (m.group(2), m.group(1)); continue
    m = re.match(r'^ME\+JSD\+LLR\+PDR\s+\S+\s+((?:[\d.]+\s+){7})검체', ln)
    if m and cur: ref[cur + ('검체',)] = [float(x) for x in m.group(1).split()]; continue
    m = re.match(r'^\s+((?:[\d.]+\s+){7})사람', ln)
    if m and cur and cur + ('검체',) in ref and cur + ('사람',) not in ref: ref[cur + ('사람',)] = [float(x) for x in m.group(1).split()]
rows = []; bad = 0; chk = 0
for p in PANELS:
    for dp in DEPTHS:
        v = p + '_' + dp + '_uni'
        fs = sorted(glob.glob('%s/%s/*/results/%s/step04_ML_classifier/y_prob_all_%s_all_combos.csv' % (R, tree, v, v)))
        d = pd.concat([pd.read_csv(f, usecols=['sample', 'feature_set', 'ratio', 'y_true', 'y_score']) for f in fs])
        d = d[d.feature_set == PR]; d['s'] = d['sample'].map(S)
        for unit in ('검체', '사람'):
            for gname, keep in GR:
                vals = []
                for r in RAT:
                    x = d[np.isclose(d.ratio, r)]
                    if unit == '사람':
                        x = x.groupby(['sample', 'y_true', 's']).y_score.mean().reset_index()
                    if keep is not None:
                        x = x[(x.y_true == 0) | x.s.isin(keep)]   # 음성은 전체 공통
                    vals.append(auc(x.y_true.to_numpy(), x.y_score.to_numpy()))
                rows.append((p, dp, unit, gname, vals))
                if gname == '전체' and (p, dp, unit) in ref:
                    for a, b in zip(vals, ref[(p, dp, unit)]):
                        chk += 1; bad += abs(round(a, 3) - b) > 1e-9
print('검산: 전체 48 칸 %d 중 표와 다른 칸 %d' % (chk, bad))
for p, dp, unit, g, vals in rows:
    print('ROW\t%s\t%s\t%s\t%s\t%s' % (p, dp, unit, g, '\t'.join('%.3f' % v for v in vals)))
