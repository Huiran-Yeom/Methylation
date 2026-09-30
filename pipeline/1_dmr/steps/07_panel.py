# -*- coding: utf-8 -*-
"""최종 DMR 패널: 74검체 전체로 1회 (2026-08-13)

  조각을 나누지 않는다. dmr_nested_cv.py 의 fold 안 절차와 같되
  학습 자료가 74검체 전부다.

  산출
    DMR_ranked_m.csv     q<0.05 & |dbeta|>=DBETA 인 블록 전부 · 순위 매김
    DMR_confirmed_m.csv  위에서 TOPN 개              ← 이것이 패널
    panel_summary.txt

  사용
    python 07_panel.py                     # DBETA 30 · N 2000
    python 07_panel.py --topn 1000         # 자르는 위치만 변경 (재선정 안 함)
"""

# 이 파일은 스크립트다. 다른 코드가 import 하면 본체가 그대로 돌아
# 실제 자료를 덮어쓴다. 실행은 `python <파일>` 로만 한다.
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


import os, sys, time
import numpy as np
import pandas as pd
from sklearn.ensemble import RandomForestClassifier
import xgboost as xgb
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from moderated_t import moderated_ttest, bh, efp

# ╔══ 이 단계가 쓰는 경로 ═══════════════════════════════════════════╗
#   배치가 다르면 「여기만」 고친다. 0_setup/config.conf 에는 뿌리 둘
#   (METH_ROOT=_M · DATASET_ROOT=_D)과 실험 설정만 둔다.
#   환경변수나 명령줄 인자로 그때그때 덮어쓸 수 있다(아래 참조).
# 2026-09-26: 기본값이 `m_candidates`·`m_panel` 이었다 — **옛 m 판**이다.
#   `m_candidates` 는 서버에 아직 있고 `blocks_frac060.csv` 도 들어 있다. 그래서
#   인자 없이 돌리면 죽지 않고 **옛 세대 후보로 새 패널을 만든다.** 조용히.
#   이 파일은 원래 손으로 인자를 줘서 돌렸으므로(드라이버에 없다) 필수로 바꾼다.
CAND  = None
OUT   = None
# ╚══════════════════════════════════════════════════════════════════╝
DBETA = 30.0
NBOOT = 100
TOPN  = 2000
NJOB  = 8
a = sys.argv
if '--dbeta' in a: DBETA = float(a[a.index('--dbeta') + 1])
if '--boot'  in a: NBOOT = int(a[a.index('--boot')  + 1])
if '--topn'  in a: TOPN  = int(a[a.index('--topn')  + 1])
if '--cand'  in a: CAND  = a[a.index('--cand')      + 1]   # 2026-09-10 후보 풀을 CLI 로
if '--out'   in a: OUT   = a[a.index('--out')       + 1]
if CAND is None or OUT is None:
    sys.exit('사용: python 07_panel.py --cand <후보폴더> --out <산출폴더> [--topn N] [--dbeta X]\n'
             '  기본값을 두지 않는다. 옛 세대 후보로 조용히 새 패널을 만든 적이 있다.\n'
             '  09-14 에 쓴 값: --cand %s/results/dmr/j_candidates15\n'
             '                 --out  %s/results/dmr/j_panel_dmr15/bl200_d10_cellline' % (_M, _M))
# 세대가 어긋나면 멈춘다. 07_llrfeat 이 쓰는 것과 같은 관문이다.
# 2026-09-29: 06_jsdfeat·07_llrfeat 과 같은 꼴로 맞춘다. 두 구멍이 있었다 —
#   ① 「포함」(in)이라 j_candidates150 이 GEN=15 를 통과했다
#   ② `if _gen and` 라서 환경에 GEN 이 없으면 관문이 통째로 꺼졌다(손으로 돌릴 때).
_gen = _os.environ.get('GEN', '') or getattr(_cfg, 'GEN', '')
_want_cand = 'j_candidates' + _gen
if _os.path.basename(CAND.rstrip('/')) != _want_cand:
    sys.exit('중단: 세대 불일치 — GEN=%s 이면 --cand 는 %s 여야 한다 (받은 값 %s)'
             % (_gen, _want_cand, CAND))

B = pd.read_parquet(CAND + '/beta_matrix.parquet')
M = pd.read_csv(CAND + '/samples.csv')
X0 = B.T.reindex(M['sample']).values.astype(np.float32)
y0 = (M['type'] == 'GBM').astype(int).values
# 2026-09-14: **세포주 단위**로 묶는다. 문서와 코드를 맞추는 교정이다.
#   j판_DMR선정_20260820.md:441 "② 조정 t검정  정상 20 vs 암 **10**"
#                          :442 "④ |Δβ|  **세포주 평균끼리의 차이**"
#                          :449 "⑦ 부트스트랩 … (**세포주 단위 재추출**)"
#   그런데 여기는 45파일을 독립으로 세고 있었다. SNU3-1·-2·-3 은 같은 세포주다.
#   :1345 가 m 판이 부풀려진 원인으로 "파일 단위로 세서" 를 이미 지목했는데
#   같은 구조가 j 판 Baseline 에 남아 있었다.
#   06_select.py:28-31 (JSD 쪽)이 쓰는 것과 같은 묶음이다.
_c = [str(i) for i in range(X0.shape[1])]
_D = pd.DataFrame(X0, columns=_c); _D['_g'] = M['group'].values; _D['_y'] = y0
_A = _D.groupby('_g', sort=True).mean()
y = (_A['_y'].values > 0.5).astype(int)
X = _A[_c].values.astype(np.float32)
g = np.array(_A.index)                       # 세포주 단위 (정상은 각자 한 단위)
print('독립단위 %d (암 %d · 정상 %d) · 블록 %d · 그룹 %d · 부트스트랩 %d · DBETA %.0f · N %d'
      % (len(y), int(y.sum()), int((y == 0).sum()), X.shape[1], len(set(g)), NBOOT, DBETA, TOPN))


def imp(tr, Z):
    mu = np.nanmean(tr, axis=0); mu = np.where(np.isfinite(mu), mu, 0.0)
    return np.where(np.isfinite(Z), Z, mu)


def stability(Xt, yt, gt, cols, nboot):
    t0 = time.time(); rng = np.random.RandomState(42); ug = np.unique(gt)
    freq = np.zeros(len(cols)); top = min(1000, len(cols))
    for b in range(nboot):
        pick = rng.choice(ug, len(ug), replace=True)
        idx = np.concatenate([np.flatnonzero(gt == u) for u in pick])
        if len(np.unique(yt[idx])) < 2:
            continue
        Z = imp(Xt[:, cols], Xt[idx][:, cols])
        for m in (RandomForestClassifier(100, max_features='sqrt', n_jobs=NJOB, random_state=b),
                  xgb.XGBClassifier(n_estimators=100, max_depth=3, tree_method='hist',
                                    verbosity=0, n_jobs=NJOB, random_state=b)):
            m.fit(Z, yt[idx]); freq[np.argsort(-m.feature_importances_)[:top]] += 1
        if (b + 1) % 10 == 0:
            print('   부트스트랩 %d/%d  %.0f초' % (b + 1, nboot, time.time() - t0), flush=True)
    return freq


t0 = time.time()
d, t, p, dft, d0, s0 = moderated_ttest(X[y == 1], X[y == 0])
q = bh(p)
s1  = np.flatnonzero(q < 0.05)
sig = np.flatnonzero((q < 0.05) & (np.abs(d) >= DBETA))
print('q<0.05 %d  →  |dbeta|>=%.0f 까지 %d / %d' % (len(s1), DBETA, len(sig), X.shape[1]))
e1, u1 = efp(q, s1); e2, u2 = efp(q, sig)
print('  다중검정  q<0.05 %d개 중 기대 위양성 %.0f개(%.1f%%) · 상한 %.0f개'
      % (len(s1), e1, 100.0 * e1 / max(1, len(s1)), u1))
print('  효과크기 뒤 %d개 중 기대 위양성 %.0f개(%.1f%%) · 상한 %.0f개'
      % (len(sig), e2, 100.0 * e2 / max(1, len(sig)), u2))
if len(sig) < TOPN:
    print('! 유의 블록 %d 개가 N=%d 보다 적다. 전부 사용' % (len(sig), TOPN))
# 2026-09-28: 0 이면 여기서 멈춘다. 그대로 가면 stability() 안의 sklearn 이
#   'Found array with 0 feature(s)' 로 터져 「왜 0인지」 를 알 수 없다.
#   연기시험(작은 표본)에서 실제로 그 traceback 만 남았다.
if len(sig) == 0:
    raise SystemExit(
        '유의 블록이 0 입니다. q<0.05 이고 |dbeta|>=%.0f 인 블록이 없습니다. '
        '(q<0.05 는 %d개였습니다.)  DBETA 를 낮추거나(--dbeta) 입력을 늘리십시오.'
        % (DBETA, len(s1)))

fr = stability(X, y, g, sig, NBOOT)
op = np.lexsort((-np.abs(d[sig]), q[sig], -fr))       # 빈도내림 -> q오름 -> |dbeta|내림
order = sig[op]; frs = fr[op]

T = pd.DataFrame({'rank':  np.arange(1, len(order) + 1),
                  'chr':   [B.index[i][0] for i in order],
                  'blk':   [B.index[i][1] for i in order],
                  'start': [int(B.index[i][1]) * 100 for i in order],
                  'q':     np.round(q[order], 10),
                  'dbeta': np.round(d[order], 2),
                  'freq':  frs.astype(int)})
os.makedirs(OUT, exist_ok=True)
T.to_csv('%s/DMR_ranked_m.csv' % OUT, index=False)
P = T.head(min(TOPN, len(T)))
P.to_csv('%s/DMR_confirmed_m.csv' % OUT, index=False)

P2 = P.copy(); P2['chr'] = P2['chr'].astype(str)
P2['win'] = (P2['blk'].astype('int64') * 100) // 10000000
w = P2.groupby(['chr', 'win']).size().sort_values(ascending=False)
ch = P2['chr'].value_counts()
txt = [
    '최종 DMR 패널 — %s' % pd.Timestamp.now().strftime('%Y-%m-%d %H:%M'),
    '재료        %d검체 (정상 %d · 암 %d) · 후보 %d 블록'
    % (len(y), int((y == 0).sum()), int(y.sum()), X.shape[1]),
    '조건        q<0.05 · |dbeta|>=%.0f · 부트스트랩 %d회 x 2모델' % (DBETA, NBOOT),
    '정렬        빈도 내림 -> q 오름 -> |dbeta| 내림',
    '',
    'q<0.05             %d' % len(s1),
    '  기대 위양성       %.0f개 (%.1f%%)' % (e1, 100.0 * e1 / max(1, len(s1))),
    '효과크기까지        %d' % len(sig),
    '패널(N)            %d' % len(P),
    '',
    '저메틸(암<정상)     %d' % int((P['dbeta'] < 0).sum()),
    '고메틸(암>정상)     %d' % int((P['dbeta'] > 0).sum()),
    '빈도 범위          %d ~ %d' % (int(P['freq'].min()), int(P['freq'].max())),
    'q 최대             %.3g' % P['q'].max(),
    '|dbeta| 최소       %.1f' % P['dbeta'].abs().min(),
    '',
    '염색체 상위5       %s' % ', '.join('chr%s %d' % (c, n) for c, n in ch.head(5).items()),
    '10Mb 구간 수       %d' % len(w),
    '상위1구간          chr%s %d~%dMb  %d개 (%.1f%%)'
    % (w.index[0][0], w.index[0][1] * 10, w.index[0][1] * 10 + 10,
       int(w.iloc[0]), 100.0 * w.iloc[0] / len(P)),
    '상위5구간 합       %.1f%%' % (100.0 * w.head(5).sum() / len(P)),
    '',
    '소요               %.0f초' % (time.time() - t0),
]
open('%s/panel_summary.txt' % OUT, 'w', encoding='utf-8').write('\n'.join(txt) + '\n')
print('\n' + '\n'.join(txt))
print('\n저장: %s' % OUT)
