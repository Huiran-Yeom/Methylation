# -*- coding: utf-8 -*-
"""4단계 학습 표: 트리의 모든 판 · 깊이의 학습 CV 를 한 CSV 로 모은다(옛 드라이버 셋의 heredoc 파이썬을 그대로 옮김).
쓰는 법: python cv_all.py <METH_ROOT> <GEN> <출력 CSV>"""
import glob, re, sys
import pandas as pd

r, g, out = sys.argv[1], sys.argv[2], sys.argv[3]
rows = []
for f in glob.glob(r + '/results/panels/*/*/step04_ML_classifier/cv_*_*.csv'):
    m = re.match(r'j%s(\w+?)_(\d+k)_uni$' % g, f.split('/')[-3])
    if m:
        d = pd.read_csv(f); d.insert(0, 'panel', m.group(1)); d.insert(1, 'depth', m.group(2)); rows.append(d)
pd.concat(rows).sort_values(['panel', 'depth', 'mut_reads', 'feature_set', 'model']).to_csv(out, index=False)
