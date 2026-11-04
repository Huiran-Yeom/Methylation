# -*- coding: utf-8 -*-
"""4단계 학습 표: 모든 판, 깊이의 학습 CV 를 CSV 하나로 모음
사용법: python cv_all.py <METH_ROOT> <GEN> <출력 CSV>"""
import glob, re, sys
import pandas as pd

r, g, out = sys.argv[1], sys.argv[2], sys.argv[3]
rows = []
for f in glob.glob(r + '/results/panels/*/*/step04_ML_classifier/cv_*_*.csv'):
    m = re.match(r'j%s(\w+?)_(\d+k)_uni$' % g, f.split('/')[-3])
    if m:
        d = pd.read_csv(f); d.insert(0, 'panel', m.group(1)); d.insert(1, 'depth', m.group(2)); rows.append(d)
pd.concat(rows).sort_values(['panel', 'depth', 'mut_reads', 'feature_set', 'model']).to_csv(out, index=False)
