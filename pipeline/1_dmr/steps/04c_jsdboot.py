# -*- coding: utf-8 -*-
"""S0db JSD 부트스트랩 순위 안정성 (2026-08-20)

  독립단위를 복원추출로 다시 뽑아 JSD 를 다시 재고, 「JSD 상위 TOPN 에 몇 번
  들었나」를 센다(freq). 100회 반복한다.

    섞는 단위는 독립단위다 (세포주 복제본 3개는 함께 움직인다)
    집단 크기를 유지한다

  산출  <out>/boot_K<K>_minr<MINR>.parquet   열: chr · blk · jsd_obs · freq
  쓰는 곳  06b_select_supplement.py 가 freq 로 정렬해 jsdb200 판을 고른다.

  ※ 순열 귀무분포(문턱 t)는 옆 파일 04b_jsdnull.py 가 만든다. 2026-09-29 에
    이 독스트링이 04b 의 것과 「똑같아」 부트스트랩 스크립트가 순열 검정이라고
    적혀 있었다 - 복원할 때 함께 복사된 것이다.
"""

# 이 파일은 스크립트다. 다른 코드가 import 하면 본체가 그대로 돌아
# 실제 자료를 덮어쓴다. 실행은 `python <파일>` 로만 한다.
if __name__ != '__main__':
    raise ImportError(__file__ + ' 은(는) 스크립트다 — import 하지 않는다')

import os as _os, sys as _sys
# 0_setup 을 찾는다. 이 파일이 「복사되어」 정리본 밖에서 돌 수도 있으므로
#   ① 자기 위에서 위로 올라가며 찾고 ② 환경변수 METH_CONF_DIR ③ 그래도 없으면
#   meth_config 없이 돌 수 있게 기본값으로 간다(검증이 검체 폴더로 복사해 쓰는 경로다).
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


import os, sys, math
import numpy as np, pandas as pd

# 산출 폴더는 GEN(패널 세대)을 따른다. 비우면 세대 없는 옛 폴더를 가리킨다.
_GEN = _os.environ.get('GEN', '')
if not _GEN:
    try:
        _GEN = _cfg.GEN
    except Exception:
        raise SystemExit('GEN(패널 세대)을 못 읽었습니다 — config.conf 의 GEN 을 채우거나 '
                         '환경변수로 주십시오. 비워 두면 옛 세대 폴더를 가리킵니다.')
# ╔══ 이 단계가 쓰는 경로 ═══════════════════════════════════════════╗
#   배치가 다르면 「여기만」 고친다. 0_setup/config.conf 에는 뿌리 둘
#   (METH_ROOT=_M · DATASET_ROOT=_D)과 실험 설정만 둔다.
#   환경변수나 명령줄 인자로 그때그때 덮어쓸 수 있다(아래 참조).
D    = _M + '/results/dmr/j_jsdcount' + _GEN
CAND = _M + '/results/dmr/j_candidates' + _GEN
OUT  = _M + '/results/dmr/j_jsd' + _GEN
# ╚══════════════════════════════════════════════════════════════════╝
# 기본값이 4 였다. 실제 실행은 6 이고(서버의 boot_K3_minr6.parquet),
#   06b 가 설정의 MIN_READS_PER_WINDOW(=6)로 파일을 찾는다. 기본값을 설정에서
#   읽어 세 자리(04·04b·04c)가 같은 값을 쓰게 한다. --minr 로 덮을 수 있다.
K, MINR, NEED_N, NEED_G, MINREP, NPERM, SEED = 3, 6, 12, 6, 2, 100, 42
try:   # 2026-09-29: SEED 가 config.conf 에 있는데 여기서만 박혀 있었다.
    SEED = int(_cfg.SEED)
except Exception:
    pass
try:
    MINR = int(_cfg.MIN_READS_PER_WINDOW)
    K    = int(_cfg.PATTERN_K)
except Exception:
    pass

FRAC = 0.60          # 후보 풀과 같은 커버 기준
a = sys.argv
if '--k'     in a: K     = int(a[a.index('--k')+1])
if '--frac' in a: FRAC = float(a[a.index('--frac')+1])
if '--need'  in a: NEED_N = int(a[a.index('--need')+1])  # 2026-09-10 정상 커버 기준을 CLI 로
# '--needg' 가 「있는지」 만 보고 값을 안 읽었다. 주면 자동계산만
#   꺼지고 기본값 6 이 남았다. 연기시험에서 --needg 1 을 줬는데 「암 6종」 으로 돌았다.
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
# CANCER_LABEL 과 samples.csv 의 type 이 어긋나면 양성이 0개가 된다.
#   그대로 두면 AUC 가 nan 이 되거나 전부 음성으로 학습한다. 여기서 멈춘다.
if int(len(lines)) == 0:
    sys.exit('samples.csv 에 type=%s 인 검체가 없습니다. CANCER_LABEL 과 1단계가 쓴 값이'
             ' 어긋났습니다 (현재 CANCER_LABEL=%s).' % (_CANCER, _CANCER))

W  = pd.read_parquet(D + '/windows_K%d.parquet' % K)
nw, m = len(W), 1 << K
pc = ['p%d' % i for i in range(m)]
NN, NG = len(norm), len(lines)
# NEED 를 표본 수에서 자동 계산 (01_candidates.py:56-57 과 같은 규칙)
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


# ── 안정성. 0d 와 같은 적재를 쓰되 딱지는 두고 단위를 복원추출한다.
#   Baseline stability() 와 같은 형태: 매 회 JSD 상위 TOPN 에 든 칸에 +1.
TOPN = int(a[a.index('--topn') + 1]) if '--topn' in a else 1000
rs = np.random.RandomState(SEED)
print('안정성: 복원추출 %d회 · JSD 상위 %d 에 든 횟수' % (NPERM, TOPN))
freq = pd.Series(0, index=Bo.index, dtype=int)
gi, ni = allu[NN:], allu[:NN]
for p in range(NPERM):
    pg = rs.choice(gi, len(gi), replace=True)
    pn = rs.choice(ni, len(ni), replace=True)
    J, ok = jsd_of(pg, pn)
    hit = by_block(J, ok).sort_values(ascending=False).head(TOPN).index
    freq[freq.index.isin(hit)] += 1
    if (p + 1) % 20 == 0:
        print('  복원추출 %d/%d' % (p + 1, NPERM), flush=True)
out = pd.DataFrame({'jsd_obs': Bo, 'freq': freq}).reset_index()
os.makedirs(OUT, exist_ok=True)   # 2026-09-29: 04b 에는 있고 여기만 빠져 있었다
fp = '%s/boot_K%d_minr%d.parquet' % (OUT, K, MINR)
out.to_parquet(fp)
print('freq: 만점 %d칸 · 90회이상 %d칸 · 50회이상 %d칸 · 1회이상 %d칸'
      % (int((freq == NPERM).sum()), int((freq >= 0.9 * NPERM).sum()),
         int((freq >= 0.5 * NPERM).sum()), int((freq > 0).sum())))
print('저장: %s' % fp)
