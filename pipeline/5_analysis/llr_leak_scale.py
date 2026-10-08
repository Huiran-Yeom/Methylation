#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
누수 크기 +154.21 을 «척도로 나눈다». 사전등록 추기 29항 예정. 2026-10-07

왜 필요한가 (검토 세션 지적 C)
  +154.21 은 척도 없이는 뜻이 없다. 학습의 양성–음성 격차가 1,500 이면 누수가
  10% 를 설명하고, 150 이면 누수가 신호의 전부다. 전혀 다른 주장이다.

단위를 맞춰야 한다 — 그냥 나누면 안 된다
  LOO 측정(llr_loo_paired.py)은 그 검체의 «전체 리드»로 합산한다. 수십만 리드다.
  학습 특징행렬의 LLR 은 «5,000 리드»로 낸 값이다. 두 수는 단위가 다르다.
  그래서 리드당으로 환산한 뒤 5,000 을 곱해 특징행렬 단위로 옮긴다.
  옮긴 값이 맞는지는 (a) 조건을 특징행렬 음성 중앙과 대조해 검산한다.

낼 수 있는 결론의 형태
  「누수를 없애면 LLR 이 특징행렬 단위로 약 X 움직인다. 같은 조건에서
   학습 양성–음성 격차는 Y 다. X/Y = Z.」
  Z 가 1 을 넘으면 누수가 격차보다 크다는 뜻이고, 그러면 학습에서 LLR 이
  보던 분리는 누수로 설명된다.
"""
import os
import sys

import numpy as np
import pandas as pd

M = '<METH_ROOT>'
CAND = M + '/results/dmr/j_candidates15'
COUNT = M + '/results/dmr/j_jsdcount15'
REF_FULL = M + '/results/dmr/j_llrref15/llrref_K3_minr6.parquet'
W = os.path.expanduser('<WORK>/tmp/llr_loo')
NPAT = 8
TH = 0.005
DEPTH = 5000           # 1차 지표의 깊이. 특징행렬과 같은 단위로 옮기기 위해서다

GAP = {'j15bl200': 0.654, 'j15jsd200': 2.555, 'j15jsdb200': 1.752}


def ref(path):
    r = pd.read_parquet(path)
    return (r['wid'].to_numpy(),
            r[['n%d' % i for i in range(NPAT)]].to_numpy(dtype=np.float64),
            r[['g%d' % i for i in range(NPAT)]].to_numpy(dtype=np.float64))


def llr_and_reads(sample, path, keep):
    """(LLR, 실제로 쓴 리드, 공통 창의 총 리드).

    LLR_EPS=0 이라 P_Normal=0 인 무늬는 버린다. 그 집합이 참조마다 달라지므로
    «쓴 리드»는 조건마다 다르다 — 자기를 빼면 자기만 갖고 있던 무늬가 0 이 된다.
    그래서 분모는 조건과 무관한 «총 리드»로 고정하고, 쓴 리드는 따로 찍는다.
    """
    wid, PN, PG = ref(path)
    m = np.isin(wid, keep)
    wid, PN, PG = wid[m], PN[m], PG[m]
    c = (pd.read_parquet('%s/count_%s_K3.parquet' % (COUNT, sample))
         .set_index('wid').reindex(wid).fillna(0.0))
    X = c[['p%d' % i for i in range(NPAT)]].to_numpy(dtype=np.float64)
    ok = PN > 0
    pm = (1.0 - TH) * PN + TH * PG
    with np.errstate(divide='ignore', invalid='ignore'):
        v = float(np.nansum(np.where(ok, X * np.log(pm / np.where(ok, PN, 1.0)), 0.0)))
    return v, float(X[ok].sum()), float(X.sum())


def main():
    meta = pd.read_csv(CAND + '/samples.csv')
    norm = sorted(meta.loc[meta['type'] == 'Normal', 'sample'].astype(str))

    full_wid = ref(REF_FULL)[0]
    rows = []
    for s in norm:
        # j 규칙은 llr_loo_paired.py 와 같다 — 고정 norm[0], 자기가 norm[0] 이면 norm[1].
        y = norm[1] if s == norm[0] else norm[0]
        pb = '%s/ref_wo_%s/llrref_K3_minr6.parquet' % (W, s)   # 자기 제외 14명
        pc = '%s/ref_wo_%s/llrref_K3_minr6.parquet' % (W, y)   # 남 하나 제외 14명
        if not (os.path.exists(pb) and os.path.exists(pc)):
            raise SystemExit('참조가 없다: %s 또는 %s' % (pb, pc))
        keep = np.intersect1d(np.intersect1d(full_wid, ref(pb)[0]), ref(pc)[0])
        va, ua, ta = llr_and_reads(s, REF_FULL, keep)
        vb, ub, tb = llr_and_reads(s, pb, keep)
        vc, uc, tc = llr_and_reads(s, pc, keep)
        if not (ta == tb == tc):
            raise SystemExit('공통 창 총 리드가 어긋난다 — 코드 결함: %s' % s)
        rows.append((s, ta, ua, ub, uc, va, vb, vc))

    print('공통 창에서만 재고, 세 조건의 «총» 리드가 같은 것을 검산했다.')
    print()
    print('%-14s %9s %9s %9s %11s %11s %11s'
          % ('검체', '총리드', '쓴(a)', '쓴(b)', '(a)15명', '(c)14남', '(b)14자기'))
    for s, t, ua, ub, uc, va, vb, vc in rows:
        print('%-14s %9.0f %9.0f %9.0f %11.2f %11.2f %11.2f'
              % (s, t, ua, ub, va, vc, vb))

    drop = np.array([(r[2] - r[3]) for r in rows])      # 쓴(a) − 쓴(b)
    tot = np.array([r[1] for r in rows])
    print()
    print('자기를 빼서 «못 쓰게 된» 리드: 중앙 %.0f (%.3f%%) · 폭 %.0f~%.0f'
          % (np.median(drop), 100 * np.median(drop / tot), drop.min(), drop.max()))
    print('  = 그 사람만 갖고 있던 무늬의 리드. 누수의 직접 증거다.')

    nreads = tot
    dbc = np.array([r[6] - r[7] for r in rows])        # (b) − (c)
    dca = np.array([r[7] - r[5] for r in rows])        # (c) − (a)
    per = dbc / nreads

    print()
    print('=' * 72)
    print('리드당 환산 — 리드 수가 검체마다 다르므로 꼭 필요하다')
    print('=' * 72)
    print('  리드 수            중앙 %.0f · 폭 %.0f~%.0f'
          % (np.median(nreads), nreads.min(), nreads.max()))
    print('  Δ(b−c) 전체합      중앙 %+.2f · 폭 %+.2f~%+.2f'
          % (np.median(dbc), dbc.min(), dbc.max()))
    print('  Δ(b−c) 리드당      중앙 %+.6f · 폭 %+.6f~%+.6f'
          % (np.median(per), per.min(), per.max()))
    print('  Δ(c−a) 리드당      중앙 %+.6f  (참조 크기 효과 · 교란)'
          % np.median(dca / nreads))

    moved = float(np.median(per)) * DEPTH
    print()
    print('=' * 72
          )
    print('특징행렬 단위(%d리드)로 옮긴 누수 크기: %+.3f' % (DEPTH, moved))
    print('=' * 72)
    print('%-14s %10s %10s %8s' % ('판', '학습격차', '누수크기', '누수/격차'))
    for p in sorted(GAP):
        print('%-14s %10.3f %10.3f %8.2f' % (p, GAP[p], moved, moved / GAP[p]))
    print()
    print('학습격차 = 그 판 5k·0.5%·θ=0.005 에서 양성 100 중앙 − 음성 100 중앙.')
    print('누수크기 = 위 리드당 Δ(b−c) 중앙 × %d.' % DEPTH)
    print()
    print('주의 — 이 환산이 가정하는 것')
    print('  리드당 누수가 깊이에 선형이라고 본다. 참조가 바뀌면 무늬별 가중이')
    print('  바뀌므로 완전히 선형은 아니다. 그래서 이 수는 «자릿수»로 읽는다.')
    print('  그리고 학습 음성은 15명 풀에서 섞은 것이라 한 사람의 리드당 값을')
    print('  그대로 쓰는 것은 근사다. 정확한 값은 수정판 재학습이 낸다.')
    return 0


if __name__ == '__main__':
    sys.exit(main())
