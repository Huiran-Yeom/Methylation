# -*- coding: utf-8 -*-
"""추기 44-3(가) — 반복 1~7 의 y_score 가 옛 트리와 같은가. 하나라도 다르면 rc=1."""
import sys, glob
import pandas as pd
S, V = '<검체>', 'j16bl200_5k_uni'
def one(tree):
    f = glob.glob('<METH_ROOT>/lung/%s/%s/results/%s/step04_ML_classifier/y_prob_all_%s_all_combos.csv' % (tree, S, V, V))
    if len(f) != 1:
        sys.exit('파일이 하나가 아니다: %s %d' % (tree, len(f)))
    print('읽음', f[0])
    return pd.read_csv(f[0])
o, n = one('Patient_Test_val'), one('Patient_Test10chk_val')
K = ['feature_set', 'model', 'ratio', 'rep', 'y_true']
print('옛 반복', sorted(o.rep.unique(), key=lambda x: (x[:3], int(x.lstrip('NormalGBM')))))
print('새 반복', sorted(n.rep.unique(), key=lambda x: (x[:3], int(x.lstrip('NormalGBM')))))
n7 = n[n.rep.str.replace(r'^\D+', '', regex=True).astype(int) <= 7]
m = o.merge(n7, on=K, how='outer', suffixes=('_옛', '_새'), indicator=True)
only = m[m._merge != 'both']
both = m[m._merge == 'both']
d = (both.y_score_옛 - both.y_score_새).abs()
print('짝지은 줄 %d · 한쪽만 %d · 최대 차 %.3g · 다른 줄 %d' % (len(both), len(only), d.max(), int((d > 0).sum())))
n8 = n[n.rep.str.replace(r'^\D+', '', regex=True).astype(int) >= 8]
print('새 반복 8~10 줄 %d (기대: 반복 1~7 줄수 x 3/7)' % len(n8))
ok = len(only) == 0 and (d > 0).sum() == 0 and len(both) > 0
print('판정', '같다 — 통과' if ok else '다르다 — 멈춘다')
sys.exit(0 if ok else 1)
