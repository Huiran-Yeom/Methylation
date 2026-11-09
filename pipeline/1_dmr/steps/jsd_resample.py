# -*- coding: utf-8 -*-
"""04b(JSD 귀무), 04c(JSD 부트스트랩)가 같이 쓰는 앞부분
setup(argv, cfg): 설정 · 인자 → samples · 창 · 단위 적재 → 관측 JSD 를 계산해 이름공간으로 돌려줌
  돌려주는 것: a NN NG NPERM SEED OUT K MINR allu Bo jsd_of by_block
  난수는 부르는 쪽이 만듦 (자기 순열 · 복원추출 순서 그대로)
"""
import math
import os
import sys
from types import SimpleNamespace

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import dmr_common   # noqa: E402


def setup(argv, _cfg):
    _M = _cfg.METH_ROOT

    # 산출 폴더 이름에 GEN(패널 세대) 붙임. 못 읽으면 중단
    _GEN = dmr_common.gen(_cfg)
    # 이 단계 경로. 폴더 구조 다르면 여기만 수정 (인자로 덮어쓰기 가능)
    D    = _M + '/results/dmr/j_jsdcount' + _GEN
    CAND = _M + '/results/dmr/j_candidates' + _GEN
    OUT  = _M + '/results/dmr/j_jsd' + _GEN
    # K, MINR 은 설정에서 읽음 (04 와 같은 값, 06b 가 이 값으로 파일 찾음). 인자가 우선
    K, MINR, NEED_N, NEED_G, MINREP, NPERM, SEED = 3, 6, 12, 6, 2, 100, 42
    try:   # SEED 는 config.conf 에서 읽음
        SEED = int(_cfg.SEED)
    except Exception:
        pass
    try:
        MINR = int(_cfg.MIN_READS_PER_WINDOW)
        K    = int(_cfg.PATTERN_K)
    except Exception:
        pass

    FRAC = 0.60          # 후보 풀과 같은 커버 기준
    a = argv
    if '--k'     in a: K     = int(a[a.index('--k')+1])
    if '--frac' in a: FRAC = float(a[a.index('--frac')+1])
    if '--need'  in a: NEED_N = int(a[a.index('--need')+1])  # 정상 커버 기준 직접 지정
    # 암 커버 기준 직접 지정
    if '--needg' in a: NEED_G = int(a[a.index('--needg')+1])
    if '--minr'  in a: MINR  = int(a[a.index('--minr')+1])
    if '--nperm' in a: NPERM = int(a[a.index('--nperm')+1])
    if '--dir'   in a: D     = a[a.index('--dir')+1]
    if '--cand'  in a: CAND  = a[a.index('--cand')+1]
    if '--out'   in a: OUT   = a[a.index('--out')+1]

    _CANCER = getattr(_cfg, 'CANCER_LABEL', 'GBM')   # samples.csv 의 type 값

    M   = pd.read_csv(CAND + '/samples.csv')
    grp = dict(zip(M['sample'], M['group'])); typ = dict(zip(M['sample'], M['type']))
    norm  = sorted(s for s in typ if typ[s] == 'Normal')
    lines = {}
    for s in sorted(t for t in typ if typ[t] == _CANCER): lines.setdefault(grp[s], []).append(s)
    dmr_common.need_cancer(len(lines), _CANCER)

    W  = pd.read_parquet(D + '/windows_K%d.parquet' % K)
    nw, m = len(W), 1 << K
    pc = ['p%d' % i for i in range(m)]
    NN, NG = len(norm), len(lines)
    # NEED 는 표본 수로 자동 계산 (01 과 같은 규칙)
    if '--need' not in a: NEED_N = int(math.ceil(FRAC * NN))
    if '--needg' not in a: NEED_G = int(math.ceil(FRAC * NG))
    print('커버 %.0f%% -> 정상 %d명 이상 · 암 %d종 이상' % (FRAC*100, NEED_N, NEED_G))
    NN, NG = len(norm), len(lines)
    print('K=%d · MINR=%d · 정상 %d명 · 암 %d종 · 순열 %d회' % (K, MINR, NN, NG, NPERM))

    def load(s):
        C = np.zeros((nw, m), dtype=np.float32)
        f = D + '/count_%s_K%d.parquet' % (s, K)
        if os.path.exists(f):
            d = pd.read_parquet(f); C[d['wid'].values] = d[pc].values
        t = C.sum(1); ok = t >= MINR
        C[ok] /= t[ok, None]; C[~ok] = 0
        return C, ok

    U = np.zeros((NN + NG, nw, m), dtype=np.float32)
    O = np.zeros((NN + NG, nw), dtype=bool)
    for i, s in enumerate(norm):
        U[i], O[i] = load(s)
    for j, (L, mem) in enumerate(sorted(lines.items())):
        Ls = np.zeros((nw, m), dtype=np.float32); nL = np.zeros(nw, np.int16)
        for s in mem:
            C, ok = load(s); Ls += C; nL += ok
        v = nL >= min(MINREP, len(mem))
        Ls[v] /= nL[v, None]; Ls[~v] = 0
        U[NN + j], O[NN + j] = Ls, v
    print('단위 %d개 적재 완료 (%.0f MB)' % (len(U), U.nbytes / 1e6))

    def jsd_of(cidx, nidx):
        nP = O[cidx].sum(0); nQ = O[nidx].sum(0)
        ok = (nQ >= NEED_N) & (nP >= NEED_G)
        J = np.full(nw, np.nan, np.float32)
        if not ok.any(): return J, ok
        P = U[cidx].sum(0)[ok] / nP[ok, None]
        Q = U[nidx].sum(0)[ok] / nQ[ok, None]
        Mm = 0.5 * (P + Q)
        with np.errstate(divide='ignore', invalid='ignore'):
            x = np.where(P > 0, P * np.log2(P / Mm), 0.0)
            y = np.where(Q > 0, Q * np.log2(Q / Mm), 0.0)
        J[ok] = 0.5 * x.sum(1) + 0.5 * y.sum(1)
        return J, ok

    def by_block(J, ok):
        T = W[['chr','blk']].copy(); T['j'] = J
        T = T[ok]
        return T.groupby(['chr','blk'], observed=True)['j'].mean()

    allu = np.arange(NN + NG)
    Jo, oko = jsd_of(allu[NN:], allu[:NN])
    Bo = by_block(Jo, oko)
    print('\n관측  창 %d개 · 칸 %d개' % (int(oko.sum()), len(Bo)))
    print('관측 칸JSD  중앙 %.3f · 95%% %.3f · 99%% %.3f · 최대 %.3f'
          % tuple(np.percentile(Bo, [50, 95, 99]).tolist() + [Bo.max()]))
    return SimpleNamespace(a=a, NN=NN, NG=NG, NPERM=NPERM, SEED=SEED, OUT=OUT, K=K, MINR=MINR, allu=allu, Bo=Bo, jsd_of=jsd_of, by_block=by_block)
