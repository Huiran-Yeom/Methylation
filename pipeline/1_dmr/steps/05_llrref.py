# -*- coding: utf-8 -*-
"""S0f LLR 참조표 (2026-09-02) — 창마다 정상·암의 무늬 확률을 만들어 저장한다.

  P_Normal  정상 20명의 무늬 분포        검체별로 합=1 로 만든 뒤 평균 (한 검체 = 한 표)
  P_GBM     암 10종의 무늬 분포          세포주 안에서 먼저 평균 → 세포주끼리 평균

04_jsd.py 와 완전히 같은 규칙이다. 어긋나면 JSD 와 LLR 이 다른 것을 보게 된다.

  깊이 보정   검체별 정규화 (MINR 미만인 창은 그 검체에서 기권)
  복제 보정   세포주 안에서 먼저 평균 (MINREP 이상 반복이 살아야 그 세포주를 센다)
  창 채택     정상 NEED_N 명 이상 · 암 NEED_G 종 이상이 값을 냈을 때

의사카운트를 여기서 넣지 않는다. 날 확률을 그대로 저장하고, 매끄럽게 하는 정도는
쓰는 쪽에서 정한다 — 참조표를 다시 만들지 않고 바꿀 수 있게.

  P_smooth = (1-eps) * P + eps / 2^K

쓰는 쪽 계산 (파일로 만들지 않는다)
  P_mix(theta) = (1-theta) * P_Normal + theta * P_GBM
  LLR_theta    = sum_pattern  count * log( P_mix / P_Normal )

명령줄:  make_llr_ref.py [--k 3] [--minr 6] [--dir <jsdcount>] [--out <폴더>]
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


import os, sys
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
OUT  = _M + '/results/dmr/j_llrref' + _GEN
K, MINR, NEED_N, NEED_G, MINREP = 3, 6, None, None, 2
try:   # 2026-09-29: 04b·04c 와 같은 값을 쓰게 설정에서 읽는다 (--minr 이 이김)
    MINR = int(_cfg.MIN_READS_PER_WINDOW)
    K    = int(_cfg.PATTERN_K)
except Exception:
    pass
# 2026-09-14: NEED_N=12 는 정상 20 기준 60% 였다. 정상 15 에서는 80% 가 되어 버린다.
#   상수 대신 비율로 잡고, 후보 풀 기준(60%)과 같게 맞춘다. 04_jsd/0d 와 동일.
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
if USE: M = M[M['sample'].isin(USE)]
grp = dict(zip(M['sample'], M['group'])); typ = dict(zip(M['sample'], M['type']))
norm = sorted(s for s in typ if typ[s] == 'Normal')
gbm  = sorted(s for s in typ if typ[s] == 'GBM')
lines = {}
for s in gbm: lines.setdefault(grp[s], []).append(s)
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

# 암 — 세포주 안에서 먼저 평균
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

# 2026-09-28 [가드 자리 교정] 이 가드가 PG 를 «만들기 전» 에 있었다.
#   암이 하나라도 있으면 np.isfinite(PG) 가 평가돼 NameError 로 죽는다 —
#   가드를 넣은 09-14 08:52 이후로 이 코드는 돌 수 없었다
#   (참조표는 08:38 산출, 즉 가드 전이다). 검사 내용은 그대로 두고 자리만 옮긴다.
if len(lines) == 0 or not np.isfinite(PG).any():
    sys.exit('중단: 암 %d종 · P_GBM 에 유효값이 없습니다 — 참조를 만들지 않습니다' % len(lines))
print('쓸 수 있는 창 %s / %s (%.1f%%)' % (format(int(ok.sum()), ','), format(nw, ','), 100.0 * ok.sum() / nw))

R = W.copy()
R['wid'] = np.arange(nw, dtype=np.int32)   # 창 번호 — count 파일·04_1d 와 맞추는 열쇠
for i in range(m):
    R['n%d' % i] = PN[:, i].astype(np.float32)
    R['g%d' % i] = PG[:, i].astype(np.float32)
R['nn'], R['ng'] = nQ, nP
R['rd_n'], R['rd_g'] = rdN, rdG
R['ok'] = ok
R = R[R.ok].drop(columns=['ok']).reset_index(drop=True)
# ╚══════════════════════════════════════════════════════════════════╝

# 진단 — 희귀 무늬가 LLR 을 지배하므로 그 안정성을 미리 본다
z = (PN[ok] == 0).sum(1)
print('정상에서 확률 0 인 무늬 수 — 창당 중앙 %d · 0개인 창 %.1f%%' % (int(np.median(z)), 100.0 * (z == 0).mean()))
lo = PN[ok][(PN[ok] > 0)]
print('정상 확률 분포 — 1%% %.5f · 5%% %.5f · 중앙 %.4f' % tuple(np.percentile(lo, [1, 5, 50])))
d = PG[ok] - PN[ok]
print('암−정상 최대 차이 — 중앙 %.3f · 90%% %.3f · 최대 %.3f' % tuple(np.percentile(d.max(1), [50, 90, 100])))

os.makedirs(OUT, exist_ok=True)
fp = OUT + '/llrref_K%d_minr%d.parquet' % (K, MINR)
R.to_parquet(fp, index=False)
print('저장: %s  (창 %s개 · 열 %d)' % (fp, format(len(R), ','), len(R.columns)))
