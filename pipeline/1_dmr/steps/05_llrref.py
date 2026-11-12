# -*- coding: utf-8 -*-
"""LLR 참조표: 창마다 정상, 암의 무늬 확률 저장
P_Normal 은 검체별 합 1 로 맞춘 뒤 평균, P_GBM 은 세포주 안에서 먼저 평균 → 세포주끼리 평균
04_jsd.py 와 같은 규칙 (MINR, MINREP, NEED_N, NEED_G). LLR_HELDOUT=1 이면 참조군만 사용

의사카운트는 안 넣음 (날 확률 저장, 쓰는 쪽에서 매끄럽게)
  P_smooth = (1-eps) * P + eps / 2^K
쓰는 쪽 계산 (파일로 안 만듦)
  P_mix(theta) = (1-theta) * P_Normal + theta * P_GBM
  LLR_theta    = sum_pattern  count * log( P_mix / P_Normal )

사용: python 05_llrref.py [--k 3] [--minr 6] [--dir <jsdcount>] [--cand ..] [--out <폴더>]
                          [--frac F] [--need N] [--needg N] [--use <검체 목록 파일>]
"""

# 스크립트 전용. import 하면 본체가 돌아 자료를 덮어씀
if __name__ != '__main__':
    raise ImportError(__file__ + ' 은(는) 스크립트다. import 하지 않는다')

import os as _os, sys as _sys
# 0_setup 찾기: 위로 올라가며 찾고, 없으면 METH_CONF_DIR 사용
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
try:
    import meth_config as _cfg
    _M, _D = _cfg.METH_ROOT, _cfg.DATASET_ROOT
except ImportError:
    raise SystemExit(
        '설정을 못 찾았습니다. 0_setup/meth_config.py 가 있는 자리를 찾지 못했습니다.\n'
        '  이 파일은 정리본 안에서 도는 것이라 기본값으로 넘어가지 않습니다.\n'
        '  트리 밖에서 돌리려면 METH_CONF_DIR 로 0_setup 자리를 주십시오.')


import os, sys
import numpy as np, pandas as pd

# 산출 폴더 이름에 GEN(패널 세대) 붙임. 못 읽으면 중단
_sys.path.insert(0, _os.path.dirname(_os.path.abspath(__file__)))
import dmr_common   # noqa: E402  (GEN, 암 검사, 단위 평균)
_GEN = dmr_common.gen(_cfg)
# 이 단계 경로. 폴더 구조 다르면 여기만 수정 (인자로 덮어쓰기 가능)
D    = _M + '/results/dmr/j_jsdcount' + _GEN
CAND = _M + '/results/dmr/j_candidates' + _GEN
OUT  = _M + '/results/dmr/j_llrref' + _GEN
K, MINR, NEED_N, NEED_G, MINREP = 3, 6, None, None, 2
try:   # 04, 04b, 04c 와 같은 값 쓰게 설정에서 읽음 (인자가 우선)
    MINR = int(_cfg.MIN_READS_PER_WINDOW)
    K    = int(_cfg.PATTERN_K)
except Exception:
    pass
# 창 채택 기준은 인원수 말고 비율 (01, 04 와 같은 60%)
FRAC = 0.60

a = sys.argv
if '--k'    in a: K    = int(a[a.index('--k')+1])
if '--minr' in a: MINR = int(a[a.index('--minr')+1])
if '--dir'  in a: D    = a[a.index('--dir')+1]
if '--cand' in a: CAND = a[a.index('--cand')+1]
if '--out'  in a: OUT  = a[a.index('--out')+1]
if '--frac' in a: FRAC = float(a[a.index('--frac')+1])
if '--need' in a: NEED_N = int(a[a.index('--need')+1])
if '--needg'in a: NEED_G = int(a[a.index('--needg')+1])
USE = set(open(a[a.index('--use')+1]).read().split()) if '--use' in a else None

M = pd.read_csv(CAND + '/samples.csv')
M_ALL = M          # held-out 분할은 전체 기준 (11_pool 과 같게)
if USE: M = M[M['sample'].isin(USE)]
grp = dict(zip(M['sample'], M['group'])); typ = dict(zip(M['sample'], M['type']))
norm = sorted(s for s in typ if typ[s] == 'Normal')
_CANCER = getattr(_cfg, 'CANCER_LABEL', 'GBM')   # samples.csv 의 type 값

gbm  = sorted(s for s in typ if typ[s] == _CANCER)
# CANCER_LABEL 과 samples.csv type 이 안 맞으면 암 0개 → 중단
if int(len(gbm)) == 0:
    sys.exit('samples.csv 에 type=%s 인 검체가 없습니다. CANCER_LABEL 과 1단계가 쓴 값이'
             ' 어긋났습니다 (현재 CANCER_LABEL=%s).' % (_CANCER, _CANCER))
lines = {}
for s in gbm: lines.setdefault(grp[s], []).append(s)
# LLR_HELDOUT 이면 참조표는 참조군만, 학습 풀(11)은 혼합군만 (겹치면 학습 CV 가 부풂)
# 가르는 규칙은 meth_config.split_ref_mix
if int(getattr(_cfg, 'LLR_HELDOUT', 0)):
    _ta = dict(zip(M_ALL['sample'], M_ALL['type'])); _ga = dict(zip(M_ALL['sample'], M_ALL['group']))
    _rn, _mn = _cfg.split_ref_mix([x for x in _ta if _ta[x] == 'Normal'])
    _rg, _mg = _cfg.split_ref_mix({_ga[x] for x in _ta if _ta[x] == _CANCER})
    print('[held-out] 참조군만 씀. 정상 %d/%d · 암 %d/%d종' % (len(_rn), len(norm), len(_rg), len(lines)))
    print('  정상 참조군 ' + ' '.join(_rn))
    print('  암   참조군 ' + ' '.join(_rg))
    norm = [s for s in norm if s in set(_rn)]
    lines = {g: v for g, v in lines.items() if g in set(_rg)}
    gbm = [s for g in sorted(lines) for s in lines[g]]
import math
if NEED_N is None: NEED_N = int(math.ceil(FRAC * len(norm)))
if NEED_G is None: NEED_G = int(math.ceil(FRAC * len(lines)))
print('K=%d · MINR=%d · 정상 %d명 · 암 %d종 (파일 %d개)' % (K, MINR, len(norm), len(lines), len(gbm)))
print('창 채택 문턱: 정상 %d명 이상 · 암 %d종 이상  (비율 %.0f%%)' % (NEED_N, NEED_G, 100*FRAC))

W  = pd.read_parquet(D + '/windows_K%d.parquet' % K)
nw, m = len(W), 1 << K
pc = ['p%d' % i for i in range(m)]

def load(s):
    C = np.zeros((nw, m), dtype=np.float64)
    f = D + '/count_%s_K%d.parquet' % (s, K)
    if os.path.exists(f):
        d = pd.read_parquet(f); C[d['wid'].values] = d[pc].values
    return C, C.sum(1)

def prob(C, t):
    ok = t >= MINR
    C[ok] /= t[ok, None]; C[~ok] = 0
    return ok

# 정상
Qs = np.zeros((nw, m)); nQ = np.zeros(nw, np.int16); rdN = np.zeros(nw, np.int64)
for s in norm:
    C, t = load(s); ok = prob(C, t)
    Qs += C; nQ += ok; rdN += np.where(ok, t, 0).astype(np.int64)

# 암: 세포주 안에서 먼저 평균
Ps = np.zeros((nw, m)); nP = np.zeros(nw, np.int16); rdG = np.zeros(nw, np.int64)
for L, mem in sorted(lines.items()):
    Ls = np.zeros((nw, m)); nL = np.zeros(nw, np.int16); rl = np.zeros(nw, np.int64)
    for s in mem:
        C, t = load(s); ok = prob(C, t)
        Ls += C; nL += ok; rl += np.where(ok, t, 0).astype(np.int64)
    v = nL >= min(MINREP, len(mem))
    Ls[v] /= nL[v, None]; Ls[~v] = 0
    Ps += Ls; nP += v; rdG += np.where(v, rl, 0)

ok = (nQ >= NEED_N) & (nP >= NEED_G)
PN = np.zeros((nw, m)); PG = np.zeros((nw, m))
PN[ok] = Qs[ok] / nQ[ok, None]
PG[ok] = Ps[ok] / nP[ok, None]

# 암 참조가 비면 중단 (PG 만든 뒤에 검사해야 함)
if len(lines) == 0 or not np.isfinite(PG).any():
    sys.exit('중단: 암 %d종 · P_GBM 에 유효값이 없습니다. 참조를 만들지 않습니다' % len(lines))
print('쓸 수 있는 창 %s / %s (%.1f%%)' % (format(int(ok.sum()), ','), format(nw, ','), 100.0 * ok.sum() / nw))

R = W.copy()
R['wid'] = np.arange(nw, dtype=np.int32)   # 창 번호 (count 파일, 2단계 llr 특징과 맞추는 키)
for i in range(m):
    R['n%d' % i] = PN[:, i].astype(np.float32)
    R['g%d' % i] = PG[:, i].astype(np.float32)
R['nn'], R['ng'] = nQ, nP
R['rd_n'], R['rd_g'] = rdN, rdG
R['ok'] = ok
R = R[R.ok].drop(columns=['ok']).reset_index(drop=True)

# 진단: 희귀 무늬가 LLR 을 좌우해서 미리 확인
z = (PN[ok] == 0).sum(1)
print('정상에서 확률 0 인 무늬 수: 창당 중앙 %d · 0개인 창 %.1f%%' % (int(np.median(z)), 100.0 * (z == 0).mean()))
lo = PN[ok][(PN[ok] > 0)]
print('정상 확률 분포: 1%% %.5f · 5%% %.5f · 중앙 %.4f' % tuple(np.percentile(lo, [1, 5, 50])))
d = PG[ok] - PN[ok]
print('암−정상 최대 차이: 중앙 %.3f · 90%% %.3f · 최대 %.3f' % tuple(np.percentile(d.max(1), [50, 90, 100])))

os.makedirs(OUT, exist_ok=True)
fp = OUT + '/llrref_K%d_minr%d.parquet' % (K, MINR)
R.to_parquet(fp, index=False)
print('저장: %s  (창 %s개 · 열 %d)' % (fp, format(len(R), ','), len(R.columns)))
