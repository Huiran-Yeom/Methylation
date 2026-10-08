#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
참조 크기 15 → 8 의 «순수» 효과를 잰다. 검토 세션 지적 B. 2026-10-07

왜 필요한가
  25·27-5항의 대조 (c) 는 참조를 15→14 로만 줄였다(Δ 리드당 −0.000011).
  그런데 수정판은 15→8 로 줄인다. 참조 분포 추정은 표본 수에 비선형이라
  −0.000011 을 일곱 번 더해 쓸 수 없다. 그러면 수정판 숫자가 내려갔을 때
  「누수 제거」와 「참조 축소」를 못 가른다.

어떻게 재나 — 새로 만들 것이 없다
  held-out 참조(8명)가 이미 있다. 참조군 8명은 «자기가 그 참조 안에 있다».
  그래서 그 8명에 대해 (a) 15명-자기포함 대 (h) 8명-자기포함 을 비교하면
  누수는 양쪽에 똑같이 남고 «크기»만 달라진다 — 순수 크기 효과다.

  혼합군 7명은 자기가 참조에 없으므로 (b)-꼴이다. 같이 찍되 크기 효과로
  읽지 않는다. 그 7명은 「누수 제거 + 크기 축소」가 섞여 있다.

읽는 법
  참조군 8명의 Δ(h−a) 가 누수 크기(리드당 +0.000401)보다 훨씬 작으면
  수정판의 변화를 누수로 설명해도 된다. 비슷하면 두 원인이 섞여 있어
  수정판 숫자만으로 「누수였다」고 말할 수 없다.
"""
import os
import sys

import numpy as np
import pandas as pd

M = '<METH_ROOT>'
H = '<METH_ROOT_HELDOUT>'
CAND = M + '/results/dmr/j_candidates15'
COUNT = M + '/results/dmr/j_jsdcount15'
REF_A = M + '/results/dmr/j_llrref15/llrref_K3_minr6.parquet'
REF_H = H + '/results/dmr/j_llrref15_heldout/llrref_K3_minr6.parquet'
SPLIT = os.path.expanduser('<WORK>/tmp/heldout')
NPAT = 8
TH = 0.005
DEPTH = 5000
LEAK_PER_READ = 0.000401        # 27-5 · llr_leak_scale.py 의 Δ(b−c) 리드당 중앙


def ref(path):
    r = pd.read_parquet(path)
    return (r['wid'].to_numpy(),
            r[['n%d' % i for i in range(NPAT)]].to_numpy(dtype=np.float64),
            r[['g%d' % i for i in range(NPAT)]].to_numpy(dtype=np.float64))


def llr(sample, path, keep):
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
    return v, float(X.sum())


def main():
    meta = pd.read_csv(CAND + '/samples.csv')
    norm = sorted(meta.loc[meta['type'] == 'Normal', 'sample'].astype(str))
    with open(SPLIT + '/ref_normal.txt') as fh:
        in_ref = {ln.strip() for ln in fh if ln.strip()}
    with open(SPLIT + '/mix_normal.txt') as fh:
        in_mix = {ln.strip() for ln in fh if ln.strip()}
    if in_ref & in_mix or len(in_ref | in_mix) != 15:
        raise SystemExit('분할이 어긋난다 — 겹침 %d · 합 %d'
                         % (len(in_ref & in_mix), len(in_ref | in_mix)))
    print('참조군 %d명 · 혼합군 %d명 · 겹침 0 · 합 15 검산' % (len(in_ref), len(in_mix)))

    keep = np.intersect1d(ref(REF_A)[0], ref(REF_H)[0])
    print('공통 창 %d (정본 %d · held-out %d)'
          % (len(keep), len(ref(REF_A)[0]), len(ref(REF_H)[0])))
    print()

    out = []
    for s in norm:
        va, t = llr(s, REF_A, keep)
        vh, _ = llr(s, REF_H, keep)
        out.append((s, s in in_ref, t, va, vh))

    print('%-14s %-7s %9s %11s %11s %11s %12s'
          % ('검체', '어느쪽', '총리드', '(a)15명', '(h)8명', 'Δ(h−a)', 'Δ 리드당'))
    for s, isref, t, va, vh in out:
        print('%-14s %-7s %9.0f %11.2f %11.2f %+11.2f %+12.6f'
              % (s, '참조군' if isref else '혼합군', t, va, vh, vh - va, (vh - va) / t))

    grp = {}
    for s, isref, t, va, vh in out:
        grp.setdefault(isref, []).append((vh - va) / t)

    print()
    print('=' * 72)
    print('순수 참조 크기 효과 15 → 8 (참조군 8명 · 자기가 양쪽 참조에 있음)')
    print('=' * 72)
    a = np.array(grp.get(True, []))
    b = np.array(grp.get(False, []))
    for lab, arr in (('참조군 8명 — 크기 효과만', a),
                     ('혼합군 7명 — 누수제거+크기 섞임', b)):
        if arr.size == 0:
            continue
        print('  %-28s 리드당 중앙 %+.6f · 폭 %+.6f~%+.6f · 양수 %d/%d'
              % (lab, np.median(arr), arr.min(), arr.max(),
                 int((arr > 0).sum()), arr.size))
    print()
    if a.size:
        size_eff = float(np.median(a))
        print('  크기 효과 (리드당)   %+.6f  → %d리드에서 %+.3f'
              % (size_eff, DEPTH, size_eff * DEPTH))
        print('  누수 (리드당)        %+.6f  → %d리드에서 %+.3f'
              % (LEAK_PER_READ, DEPTH, LEAK_PER_READ * DEPTH))
        ratio = abs(size_eff) / LEAK_PER_READ
        print('  |크기| / 누수        %.3f' % ratio)
        print()
        if ratio < 0.2:
            print('판정: 크기 효과가 누수의 20% 아래다. 수정판의 변화를 누수로')
            print('      설명해도 된다 — 다만 「주로」이고 「전부」가 아니다.')
        else:
            print('판정: 크기 효과가 작지 않다. 수정판 숫자만으로 「누수였다」고')
            print('      말할 수 없다. 두 원인을 같이 적어야 한다.')
    return 0


if __name__ == '__main__':
    sys.exit(main())
