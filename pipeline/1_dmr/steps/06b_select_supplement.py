# -*- coding: utf-8 -*-
"""보완판 JSD 패널: 2026-09-10
  필터  JSD > t_null   (귀무 백분위 · 규칙: 99.9% 기본 · 200 미달이면 99.5% · 그래도 미달이면 99%)
  정렬  freq 내림  ->  JSD 내림          (freq = 복원추출 100회 중 JSD 상위 TOPN 에 든 횟수)
  컷    상위 N
  규칙은 결과를 보기 전에 정했다. q·dbeta 는 보고용으로만 계산한다.
"""

if __name__ != '__main__':
    raise ImportError(__file__ + ' 은(는) 스크립트다 — import 하지 않는다')

import os as _os, sys as _sys
# 0_setup 을 찾는다. 이 파일이 「복사되어」 정리본 밖에서 돌 수도 있으므로
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
# 표준 머리글을 끼울 때 04_jsd 의 D·OUT·K,MINR 대입이 함께 붙어왔다.
#   아래에서 다시 대입하므로 동작에는 영향이 없었지만, OUT 이 한때 j_jsd 를
#   가리키는 것은 읽는 사람을 속인다. 이 파일은 j_panel<GEN>_supplement 에 쓴다.
#   경로·상수는 아래 한 곳에서만 정한다.

import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))  # moderated_t 는 형제다
from moderated_t import moderated_ttest, bh
# ╔══ 이 단계가 쓰는 경로 ═══════════════════════════════════════════╗
#   배치가 다르면 「여기만」 고친다. 0_setup/config.conf 에는 뿌리 둘
#   (METH_ROOT=_M · DATASET_ROOT=_D)과 실험 설정만 둔다.
#   환경변수나 명령줄 인자로 그때그때 덮어쓸 수 있다(아래 참조).
R = _M + "/results/dmr"
CAND = R + "/j_candidates" + _GEN
K    = int(_os.environ.get("PATTERN_K", 0) or _cfg.PATTERN_K)
MINR = int(_os.environ.get("MIN_READS_PER_WINDOW", 0) or _cfg.MIN_READS_PER_WINDOW)
N    = int(_os.environ.get("PANEL_N", 0) or 200)
OUT = R + "/j_panel" + _GEN + "_supplement"
# ╚══════════════════════════════════════════════════════════════════╝
boot = pd.read_parquet("%s/j_jsd%s/boot_K%d_minr%d.parquet" % (R, _GEN, K, MINR))
null = pd.read_parquet("%s/j_jsd%s/null_K%d_minr%d.parquet" % (R, _GEN, K, MINR))["null_block_jsd"].values
boot["chr"] = boot["chr"].astype(str); boot["blk"] = boot["blk"].astype("int64")
print("칸 %d · freq 만점 %d" % (len(boot), int((boot.freq == boot.freq.max()).sum())))
for pct in (99.9, 99.5, 99.0):
    t = float(np.percentile(null, pct))
    n = int((boot.jsd_obs > t).sum())
    print("  귀무 %.1f%% -> t %.4f · 통과 %d칸" % (pct, t, n))
    if n >= N:
        used_pct, used_t = pct, t; break
else:
    raise SystemExit("어느 문턱으로도 200 을 못 채움")
print("=> 채택 문턱 귀무 %.1f%% · t=%.4f" % (used_pct, used_t))
Z = boot[boot.jsd_obs > used_t].copy()
Z = Z.sort_values(["freq", "jsd_obs"], ascending=[False, False]).reset_index(drop=True)
P = Z.head(N).copy()
# 보고용 q · dbeta (선정에 안 씀)
B = pd.read_parquet(CAND + "/beta_matrix.parquet"); M = pd.read_csv(CAND + "/samples.csv")
X = B.T.reindex(M["sample"]).values.astype(np.float32)
y = (M["type"] == "GBM").astype(int).values
_c = [str(i) for i in range(X.shape[1])]
_D = pd.DataFrame(X, columns=_c); _D["_g"] = M["group"].values; _D["_y"] = y
_A = _D.groupby("_g", sort=True).mean()
y2 = (_A["_y"].values > 0.5).astype(int); X2 = _A[_c].values.astype(np.float32)
mu = np.nanmean(X2, 0); mu = np.where(np.isfinite(mu), mu, 0.0)
X2 = np.where(np.isfinite(X2), X2, mu)
d, t_, p_, _, _, _ = moderated_ttest(X2[y2 == 1], X2[y2 == 0]); q = bh(p_)
T = pd.DataFrame({"chr": B.index.get_level_values(0).astype(str),
                  "blk": B.index.get_level_values(1).astype("int64"), "dbeta": d, "q": q})
P = P.merge(T, on=["chr", "blk"], how="left")
os.makedirs(OUT, exist_ok=True)
P.to_csv("%s/panel_jsd_supplement_n%d.csv" % (OUT, N), index=False)
print("\n상위 %d: freq %d~%d · JSD %.3f~%.3f" % (N, P.freq.min(), P.freq.max(), P.jsd_obs.min(), P.jsd_obs.max()))
print("  기존 기준(선정에 안 씀): q<0.05 %d (%.1f%%) · |dbeta|>=30 %d (%.1f%%)"
      % (int((P.q < 0.05).sum()), 100 * (P.q < 0.05).mean(),
         int((P.dbeta.abs() >= 30).sum()), 100 * (P.dbeta.abs() >= 30).mean()))
print("  염색체 %d개 · 상위: %s" % (P["chr"].nunique(),
      " · ".join("chr%s %d" % (c, n) for c, n in P["chr"].value_counts().head(5).items())))
w = (P["blk"] * 100 // 10_000_000 * 10)
g = P.assign(w=w).groupby(["chr", "w"]).size()
print("  10Mb 구간 %d개 · 최다 chr%s %d~%dMb %d개(%.0f%%)"
      % (len(g), g.idxmax()[0], g.idxmax()[1], g.idxmax()[1] + 10, g.max(), 100 * g.max() / N))
print("저장: %s/panel_jsd_supplement_n%d.csv" % (OUT, N))
