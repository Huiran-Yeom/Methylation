# -*- coding: utf-8 -*-
"""S0c JSD 계산 (2026-08-20) — 여기서 처음으로 암/정상 표지를 쓴다.

저장해둔 무늬 카운트로 창마다 JSD 를 잰다. 칸의 JSD 는 창들의 평균(기본).
최종 선정에서는 조각마다 다시 계산해야 한다. 이 실행은 문턱을 정하기 위한 진단이다.

  깊이 보정   검체별로 합이 1이 되게 나눈 뒤 평균  (한 검체 = 한 표)
  복제 보정   세포주 안에서 먼저 평균 → 세포주끼리 평균
"""

# 이 파일은 스크립트다. 다른 코드가 import 하면 본체가 그대로 돌아
# 실제 자료를 덮어쓴다. 실행은 `python <파일>` 로만 한다.
if __name__ != '__main__':
    raise ImportError(__file__ + ' 은(는) 스크립트다 — import 하지 않는다')

import os as _os, sys as _sys
# 0_setup 을 찾는다. 이 파일이 «복사되어» 정리본 밖에서 돌 수도 있으므로
#   ① 자기 위에서 위로 올라가며 찾고 ② 환경변수 METH_CONF_DIR ③ 그래도 없으면
#   설정 없이 돌 수 있게 기본값으로 간다(검증이 검체 폴더로 복사해 쓰는 경로다).
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
        '설정을 못 찾았습니다 — 0_setup/meth_config.py 가 있는 자리를 찾지 못했습니다.\n'
        '  이 파일은 정리본 안에서 도는 것이라 기본값으로 넘어가지 않습니다.\n'
        '  트리 밖에서 돌리려면 METH_CONF_DIR 로 0_setup 자리를 주십시오.')


import os, sys, math
import numpy as np, pandas as pd

# 2026-09-26: 기본값에 세대 접미사가 없어 «옛 세대» 를 가리켰다.
#   09-23 사고가 옛 세대 폴더를 덮은 원인이 바로 이것이다. GEN 을 따른다.
_GEN = _os.environ.get('GEN', '')
if not _GEN:
    try:
        _GEN = _cfg.GEN
    except Exception:
        raise SystemExit('GEN(패널 세대)을 못 읽었습니다 — config.conf 의 GEN 을 채우거나 '
                         '환경변수로 주십시오. 비워 두면 옛 세대 폴더를 가리킵니다.')
# ╔══ 이 단계가 쓰는 경로 ═══════════════════════════════════════════╗
#   배치가 다르면 «여기만» 고칩니다. 0_setup/config.conf 에는 뿌리 둘
#   (METH_ROOT=_M · DATASET_ROOT=_D)과 실험 설정만 있습니다.
#   환경변수나 명령줄 인자로 그때그때 덮어쓸 수도 있습니다(아래 참조).
D    = _M + '/results/dmr/j_jsdcount' + _GEN
CAND = _M + '/results/dmr/j_candidates' + _GEN
OUT  = _M + '/results/dmr/j_jsd' + _GEN
# ╚══════════════════════════════════════════════════════════════════╝
K, MINR, NEED_N, NEED_G, MINREP = 3, 6, 12, 6, 2
try:   # 2026-09-29: 04b·04c 와 같은 값을 쓰게 설정에서 읽는다 (--minr 이 이김)
    MINR = int(_cfg.MIN_READS_PER_WINDOW)
    K    = int(_cfg.PATTERN_K)
except Exception:
    pass

FRAC = 0.60          # 후보 풀과 같은 커버 기준
a = sys.argv
if '--k'    in a: K    = int(a[a.index('--k')+1])
if '--frac' in a: FRAC = float(a[a.index('--frac')+1])
if '--need' in a: NEED_N = int(a[a.index('--need')+1])   # 2026-09-10 정상 커버 기준을 CLI 로 (60%%)
# 2026-09-28: '--needg' 는 «있는지» 만 보고 값을 안 읽어서, 주면
#   자동계산만 꺼지고 기본값이 남았다. 연기시험에서 --needg 1 을 줬는데
#   «암 6종 이상» 으로 돌았다. 조용히 무시되는 인자였다.
if '--needg' in a: NEED_G = int(a[a.index('--needg')+1])
if '--minr' in a: MINR = int(a[a.index('--minr')+1])
if '--dir'  in a: D    = a[a.index('--dir')+1]
if '--cand' in a: CAND = a[a.index('--cand')+1]
if '--out'  in a: OUT  = a[a.index('--out')+1]
USE = set(open(a[a.index('--use')+1]).read().split()) if '--use' in a else None

M = pd.read_csv(CAND + '/samples.csv')
if USE: M = M[M['sample'].isin(USE)]
grp = dict(zip(M['sample'], M['group'])); typ = dict(zip(M['sample'], M['type']))
norm  = sorted(s for s in typ if typ[s] == 'Normal')
gbm   = sorted(s for s in typ if typ[s] == 'GBM')
lines = {}
for s in gbm: lines.setdefault(grp[s], []).append(s)
# 2026-09-10 NEED 를 표본 수에서 자동 계산한다. 상수로 두면 정상 인원이 바뀔 때
#   기준이 저 혼자 엄격해진다 (20명의 12 = 60%, 15명의 12 = 80%).
#   01_candidates.py:56-57 과 같은 규칙 — ceil(FRAC x 독립단위).
if '--need' not in a: NEED_N = int(math.ceil(FRAC * len(norm)))
if '--needg' not in a: NEED_G = int(math.ceil(FRAC * len(lines)))
print('커버 %.0f%% -> 정상 %d명 이상 · 암 %d종 이상' % (FRAC*100, NEED_N, NEED_G))
print('K=%d · MINR=%d · 정상 %d명 · 암 %d종' % (K, MINR, len(norm), len(lines)))

W  = pd.read_parquet(D + '/windows_K%d.parquet' % K)
nw, m = len(W), 1 << K
pc = ['p%d' % i for i in range(m)]

def load(s):
    C = np.zeros((nw, m), dtype=np.float32)
    f = D + '/count_%s_K%d.parquet' % (s, K)
    if os.path.exists(f):
        d = pd.read_parquet(f); C[d['wid'].values] = d[pc].values
    return C, C.sum(1)

def prob(C, t):
    ok = t >= MINR
    C[ok] /= t[ok, None]; C[~ok] = 0
    return ok

Qs = np.zeros((nw, m), dtype=np.float32); nQ = np.zeros(nw, np.int16); rdN = np.zeros(nw, np.int64)
for s in norm:
    C, t = load(s); ok = prob(C, t)
    Qs += C; nQ += ok; rdN += np.where(ok, t, 0).astype(np.int64)

Ps = np.zeros((nw, m), dtype=np.float32); nP = np.zeros(nw, np.int16); rdG = np.zeros(nw, np.int64)
for L, mem in sorted(lines.items()):
    Ls = np.zeros((nw, m), dtype=np.float32); nL = np.zeros(nw, np.int16); rl = np.zeros(nw, np.int64)
    for s in mem:
        C, t = load(s); ok = prob(C, t)
        Ls += C; nL += ok; rl += np.where(ok, t, 0).astype(np.int64)
    v = nL >= min(MINREP, len(mem))
    Ls[v] /= nL[v, None]; Ls[~v] = 0
    Ps += Ls; nP += v; rdG += np.where(v, rl, 0)

ok = (nQ >= NEED_N) & (nP >= NEED_G)
P = np.zeros((nw, m), np.float32); Q = np.zeros((nw, m), np.float32)
P[ok] = Ps[ok] / nP[ok, None]; Q[ok] = Qs[ok] / nQ[ok, None]

def jsd(P, Q):
    Mm = 0.5 * (P + Q)
    with np.errstate(divide='ignore', invalid='ignore'):
        x = np.where(P > 0, P * np.log2(P / Mm), 0.0)
        y = np.where(Q > 0, Q * np.log2(Q / Mm), 0.0)
    return 0.5 * x.sum(1) + 0.5 * y.sum(1)

J = np.full(nw, np.nan, np.float32); J[ok] = jsd(P[ok], Q[ok])
W['jsd'], W['rd'], W['nn'], W['ng'] = J, rdN + rdG, nQ, nP
V = W[ok]
print('쓸 수 있는 창 %d / %d (%.1f%%)' % (len(V), nw, 100.0*len(V)/nw))
# 2026-09-28: 0 이면 여기서 멈춘다. 그대로 가면 np.percentile([]) 안에서
#   터져 numpy traceback 만 남고 «왜 0인지» 를 알 수 없다.
if len(V) == 0:
    raise SystemExit(
        '쓸 수 있는 창이 0 입니다 — 기준을 채운 창이 없습니다.'
        '  정상 %d명 / 암 %d종 이 있고 기준은 정상 %d명 이상 · 암 %d종 이상입니다.'
        % (len(norm), len(lines), NEED_N, NEED_G)
        + '  --need · --needg 로 낮추거나 입력을 늘리십시오.')
print('JSD  중앙 %.3f · 75%% %.3f · 90%% %.3f · 95%% %.3f · 99%% %.3f · 최대 %.3f'
      % tuple(np.percentile(V.jsd, [50,75,90,95,99]).tolist() + [V.jsd.max()]))
print('리드(창 합) 중앙 %d · JSD 와 log리드 상관 %+.3f   ← 0에 가까울수록 깊이 편향 적음'
      % (int(V.rd.median()), np.corrcoef(V.jsd, np.log10(V.rd + 1))[0, 1]))

G = V.groupby(['chr','blk'], observed=True).agg(
        nwin=('jsd','size'), jsd_mean=('jsd','mean'), jsd_max=('jsd','max'),
        rd=('rd','median'), nn=('nn','min'), ng=('ng','min')).reset_index()
print('칸 %d개 · 칸당 창 중앙 %d' % (len(G), int(G.nwin.median())))
print('칸 JSD(평균)  중앙 %.3f · 90%% %.3f · 95%% %.3f · 99%% %.3f'
      % tuple(np.percentile(G.jsd_mean, [50,90,95,99])))
os.makedirs(OUT, exist_ok=True)
G.to_parquet(OUT + '/jsd_K%d_minr%d.parquet' % (K, MINR), index=False)
print('저장: %s/jsd_K%d_minr%d.parquet' % (OUT, K, MINR))
