# -*- coding: utf-8 -*-
"""size_cv 02b: 크기별 jsd / bl 패널을 검증 형식(Feature, dbeta, q)으로 내보내기
  jsd<N> ← j_panel<GEN>/panel_jsd_jsd_mean_n<N>.csv (06)
  bl<N>  ← j_panel_dmr<GEN>/bl<N>_cap<T>/DMR_confirmed_m.csv (07)
  모든 패널 블록의 합집합을 union.csv 로 (풀 BAM 영역)
사용: python 02b_export_sizes.py --dmr <results/dmr> --gen <GEN> --sizes "200 300 …" [--topcap N|숫자]
"""
if __name__ != '__main__':
    raise ImportError(__file__ + ' 은(는) 스크립트다 — import 하지 않는다')
import os
import sys

import pandas as pd

a = sys.argv
D = a[a.index('--dmr') + 1]
G = a[a.index('--gen') + 1]
SIZES = [int(x) for x in a[a.index('--sizes') + 1].split()]
TC = a[a.index('--topcap') + 1] if '--topcap' in a else 'N'
out = '%s/j_panel_dmr%s_panels' % (D, G)
U = set()
for N in SIZES:
    cap = N if TC == 'N' else int(TC)
    for name, src in (('jsd%d' % N, '%s/j_panel%s/panel_jsd_jsd_mean_n%d.csv' % (D, G, N)),
                      ('bl%d' % N, '%s/j_panel_dmr%s/bl%d_cap%d/DMR_confirmed_m.csv' % (D, G, N, cap))):
        d = pd.read_csv(src)
        f = pd.DataFrame({'Feature': d['chr'].astype(str) + '_' + d['blk'].astype(str),
                          'dbeta': d['dbeta'].astype(float), 'q': d['q'].astype(float)})
        if not f['Feature'].is_unique:
            sys.exit('멈춤: %s 블록 중복' % name)
        os.makedirs('%s/%s' % (out, name), exist_ok=True)
        f.to_csv('%s/%s/DMR_confirmed_%s.csv' % (out, name, name), index=False)
        short = '' if len(f) == N else ' · N 미달'
        open('%s/%s/_origin.txt' % (out, name), 'w').write('원본 %s\n%d블록%s\n' % (src, len(f), short))
        print('%-7s %4d블록%s' % (name, len(f), short))
        U |= set(zip(d['chr'].astype(str), d['blk'].astype(int)))
pd.DataFrame(sorted(U), columns=['chr', 'blk']).to_csv(D + '/union.csv', index=False)
print('합집합 %d블록' % len(U))
