# -*- coding: utf-8 -*-
"""보완판 JSD 패널 (j_panel<GEN>_supplement)
  필터: JSD > 귀무 백분위 문턱 (99.9%, N 못 채우면 99.5%, 그래도 안 되면 99%)
  정렬: freq 내림 → JSD 내림 (freq = 부트스트랩 100회 중 JSD 상위에 든 횟수), 상위 N
  q, dbeta 는 보고용으로만 계산
"""

if __name__ != '__main__':
    raise ImportError(__file__ + ' 은(는) 스크립트다 — import 하지 않는다')

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

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))  # moderated_t 는 같은 폴더
from moderated_t import moderated_ttest, bh
# 이 단계 경로. 폴더 구조 다르면 여기만 수정 (인자로 덮어쓰기 가능)
R = _M + "/results/dmr"
CAND = R + "/j_candidates" + _GEN
K    = int(_os.environ.get("PATTERN_K", 0) or _cfg.PATTERN_K)
MINR = int(_os.environ.get("MIN_READS_PER_WINDOW", 0) or _cfg.MIN_READS_PER_WINDOW)
N    = _cfg.PANEL_N    # 패널 블록 수 (config.conf PANEL_N)
OUT = R + "/j_panel" + _GEN + "_supplement"
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
    raise SystemExit("어느 문턱으로도 %d 을 못 채움" % N)
print("=> 채택 문턱 귀무 %.1f%% · t=%.4f" % (used_pct, used_t))
Z = boot[boot.jsd_obs > used_t].copy()
Z = Z.sort_values(["freq", "jsd_obs"], ascending=[False, False]).reset_index(drop=True)
P = Z.head(N).copy()
# 보고용 q, dbeta (선정에 안 씀)
B = pd.read_parquet(CAND + "/beta_matrix.parquet"); M = pd.read_csv(CAND + "/samples.csv")
_CANCER = getattr(_cfg, 'CANCER_LABEL', 'GBM')   # samples.csv 의 type 값
X2, y2, _ = dmr_common.unit_matrix(B, M, _CANCER)   # 세포주 단위 평균
X2 = dmr_common.fill_col_mean(X2)
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
w = (P["blk"] * int(getattr(_cfg, "BLOCK_SIZE", 100)) // 10_000_000 * 10)
g = P.assign(w=w).groupby(["chr", "w"]).size()
print("  10Mb 구간 %d개 · 최다 chr%s %d~%dMb %d개(%.0f%%)"
      % (len(g), g.idxmax()[0], g.idxmax()[1], g.idxmax()[1] + 10, g.max(), 100 * g.max() / N))
print("저장: %s/panel_jsd_supplement_n%d.csv" % (OUT, N))
