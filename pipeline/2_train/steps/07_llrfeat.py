#!/usr/bin/env python
# 2026-09-26: 09-23 사고에서 가드 없는 이 부류가 importlib 로 「확인」 되다 실제로 돌았다.
#   82개 파일이 덮였다(내용은 결정적이라 같았지만 그건 운이었다).
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
    # 검증이 이 파일을 검체 폴더로 「복사」 해 돌린다. 그때는 0_setup 이 없다.
    #   조용히 넘어가면 어느 뿌리를 썼는지 로그에 안 남는다. 알리고 넘어간다.
    _M = _os.environ.get('METH_ROOT', '/ssd_data/Methylation').rstrip('/')
    _D = _os.environ.get('DATASET_ROOT', '/ssd_data/dataset').rstrip('/')
    _sys.stderr.write('  [설정 없음] 환경변수/기본값으로 진행 — METH_ROOT=%s\n' % _M)


# coding: utf-8
"""2_train/steps/07 - Multi-θ LLR 피처.  2026-09-03 미팅에서 정한 것.

무엇을 만드나
  샘플 하나당 한 줄. 여러 종양 분획 θ 에 대한 로그우도비를 모아 놓은 것.

    P_mix(θ)(x) = (1 − θ) · P_Normal(x)  +  θ · P_GBM(x)
    LLR(θ)      = Σ_reads  ln [ P_mix(θ)(x) / P_Normal(x) ]

  x 는 창(K자리)의 메틸 무늬. 합은 판의 **모든 창**에 걸쳐 돈다.
  θ 별 reference 파일은 만들지 않는다. P_Normal · P_GBM 만 두고 즉석 계산한다.

  나오는 열
    LLR_<θ> …        θ 마다 하나
    LLR_max          θ=0 을 포함한 최댓값. LLR(0)≡0 이므로 항상 0 이상 → GLRT 통계량
    best_theta       LLR 을 최대로 만든 θ. 정상이면 0 이 된다
    Type             mean CSV 에서 그대로 복사 (2_train/steps/08 가 라벨로 쓴다)

참조표는 새로 만들지 않는다
  `_dmr/05_llrref.py` 가 2026-09-02 에 만들어 둔 것을 그대로 읽는다.
  그쪽이 더 정확하다. P_GBM 을 세포주 안에서 먼저 평균 낸 뒤 세포주끼리 평균하므로
  반복 수가 많은 세포주에 치우치지 않는다.

    <METH_ROOT>/results/dmr/j_llrref/llrref_K3_minr6.parquet
    열  chr blk wi pos1 pos2 pos3 wid  n0..n7(P_Normal)  g0..g7(P_GBM)  nn ng rd_n rd_g

평활: 기본은 **넣지 않는다**. 이름을 정확히 적는다.

    LLR_EPS = 0.0   (기본)   평활 없음.  P_Normal(x)=0 인 무늬의 리드는 항에서 뺀다
    LLR_EPS > 0              P_smooth = (1 − eps) · P + eps / 2^K

  `~/tmp/results/llr_stability_normalLOO.csv` 의 `eps=0.0` 행이 바로 이 조건이고,
  leave-one-donor-out 정상이 −2.5 ~ −180.8 (음수 = 정상 판정)로 유지된다.
  같은 검사에서 바닥값 floor=0.0005 는 +77 ~ +863 (양성)으로 부호가 뒤집힌다.
  → **채택한 것은 "무평활 + P_N=0 무늬 제외" 이고, 그것이 LOO 로 검증된 조건이다.**
     "eps 혼합을 채택했다" 고 적었던 것은 잘못된 이름이었다 (2026-09-03 정정).

  주의: P_Normal(x)=0 은 "정상에서 한 번도 안 나온 무늬" 라 우도비로는 증거력이 가장 크다.
  그것을 빼는 것은 보수적인 선택이다. 뺀 리드 수를 진단 파일에 남기므로 크기를 확인할 것.
  실측(jjsd200_50k_uni)에서는 검체당 4 / 12,606 리드 수준이었다.

무늬 규칙은 03_jsdcount.py · 2_train/steps/06 와 완전히 동일하다.
  XM 태그 'Z'=1 · 'z'=0 · 자리는 rp+1 · 첫 자리가 최상위 비트
  query_name 이 같은 연속 리드를 한 조각으로 (PE 두 짝 병합)
  창의 K 자리를 다 덮은 조각만 센다

창 번호(wid) 정합성을 런타임에 검사한다
  2_train/steps/06 는 blocks_cpg.parquet 로 wid 를 재구성하면서 검사가 없었다.
  02_blockcpg 를 다른 조건으로 다시 돌리면 조용히 어긋난다.
  여기서는 참조표의 pos1~pos3 과 재구성한 창 좌표를 대조해 다르면 즉시 멈춘다.

환경변수
  LLR_REF     참조 parquet 경로 (기본 위 경로)
  LLR_THETAS  쉼표로 구분한 θ 목록
  LLR_EPS     평활 상수 (기본 0.0: 안정성 검사에서 무한 0건이었다)
  LLR_OUT     결과만 다른 폴더에 쓴다 (시험용)

명령줄:  --limit N (비율마다 앞 N개만, 시험용)
"""
import os
import sys
import glob
import re
import numpy as np
import pandas as pd
import pysam
import config as C

# ╔══ 이 단계가 쓰는 경로 ═══════════════════════════════════════════╗
#   배치가 다르면 「여기만」 고친다. 0_setup/config.conf 에는 뿌리 둘
#   (METH_ROOT=_M · DATASET_ROOT=_D)과 실험 설정만 둔다.
#   환경변수나 명령줄 인자로 그때그때 덮어쓸 수 있다(아래 참조).
# 2026-08-27 폴더 구조 변경 반영. 2_train/steps/06 는 아직 옛 경로라 고장나 있다.
CAND = _os.environ.get('METH_ROOT', _M) + '/results/dmr/j_candidates'
if '--cand' in sys.argv:            # 2026-09-14 판마다 후보 풀이 다르다
    CAND = sys.argv[sys.argv.index('--cand') + 1]
# 2026-09-14: 검증(run_one_sample.sh)은 인자 없이 부른다. 환경변수로도 받는다.
CAND = os.environ.get('LLR_CAND', CAND)
# 2026-09-14 [세대 검사] 인자를 빠뜨리면 **조용히 옛 후보 풀**로 돈다.
#   실제로 복구 스크립트가 --cand 를 빠뜨려 정상 20명 기준으로 돌기 시작했다.
#   호출하는 쪽을 전부 맞추는 것은 계속 실패한다. 지나가는 잎에서 한 번 막는다.
# 2026-09-26: 'j15' 를 못 박으면 다음 세대에 관문이 스스로 꺼진다. GEN 을 따른다.
_G = getattr(C, 'GEN')
_want = 'j_candidates' + _G
# 2026-09-29: 06_jsdfeat 과 같은 이유로 VERSION 조건을 없앤다 — 판 이름이
#   j<GEN> 으로 시작하지 않으면 관문이 스스로 꺼지고 세대 없는 풀로 조용히 돌았다.
if os.path.basename(CAND.rstrip('/')) != _want:
    sys.exit('중단: 세대 불일치 — VERSION=%s · GEN=%s 인데 CAND=%s (기대 %s)'
             % (C.VERSION, _G, CAND, _want))
REF_DEFAULT = _os.environ.get('METH_ROOT', _M) + '/results/dmr/j_llrref/llrref_K3_minr6.parquet'

# 2026-09-14 [자기 증명] 실행 중인 코드가 자기 소스를 해시해 찍는다.
#   이 줄은 **이 코드를 실제로 돌려야만** 로그에 생긴다. 사후에 못 만든다.
#   해석된 입력 경로도 같이 찍는다. 참조표 두 개가 파일 이름이 같아 로그로
#   구분이 안 됐던 일(R114)이 경로를 안 찍어서 생겼다.
def _self_proof():
    import hashlib as _h
    try:
        _d = _h.md5(open(__file__, 'rb').read()).hexdigest()[:8]
    except Exception:
        _d = '????????'
    print('  [코드] %s %s md5 %s'
          % (getattr(C, 'VERSION', '?'), os.path.basename(__file__), _d))
    for _k in ('CAND', 'JCOUNT', 'LLR_REF', 'REF', 'REF_PATH'):
        _v = globals().get(_k)
        if isinstance(_v, str) and _v:
            print('[입력] %-8s %s' % (_k, _v))
    for _k in ('K', 'MINR', 'MAXNA'):
        _v = globals().get(_k)
        if _v is not None and not callable(_v):
            print('[설정] %-8s %s' % (_k, _v))
    sys.stdout.flush()



# 2026-09-26: getattr 의 셋째 인자를 지웠다. 기본값이 곧 조용한 오답이다 —
#   MIN_READS_PER_WINDOW 의 기본값 4 는 09-14 에 6 으로 「철회된」 값이었다.
#   복사본이 config 를 못 읽으면 06 은 4 로 세고 07 은 minr6 참조표를 읽어
#   두 피처가 다른 문턱으로 계산된다. 파일이 따로라 아무 에러도 안 난다.
K = getattr(C, 'PATTERN_K')
BLOCK = C.BLOCK_SIZE
M_PAT = 1 << K

# 미팅 화면의 θ 다섯(0.0005~0.01)에 0.02 · 0.05 를 더한 것.
#   실제 혼합비가 5% 까지 있어 그대로 두면 2~5% 검체가 전부 best_theta=0.01 에 뭉친다.
#   θ=0 은 목록에 넣지 않는다. LLR(0)≡0 이라 열로 두면 쓸모가 없다.
#   대신 LLR_max · best_theta 계산에는 0 을 후보로 참여시킨다.
THETAS_DEFAULT = [0.0005, 0.001, 0.002, 0.005, 0.01, 0.02, 0.05]

_env = os.environ.get('LLR_THETAS', '').strip()
THETAS = ([float(x) for x in _env.replace(' ', '').split(',') if x]
          if _env else list(THETAS_DEFAULT))
EPS = float(os.environ.get('LLR_EPS', 0.0))
REF_PATH = os.environ.get('LLR_REF', REF_DEFAULT)
# ╚══════════════════════════════════════════════════════════════════╝
# 2026-09-15 [세대 관문 · REF] CAND 만 보던 관문이 참조표는 안 봤다.
#   j_llrref (33,903창) 와 j_llrref15 (22,685창) 는 값이 다른 표다.
#   REF_DEFAULT 는 15 가 아니므로, 러너가 LLR_REF 를 안 내보내면 조용히 옛 표를 쓴다.
#   러너가 내보내고 있지만, 관문이 그것을 증명해야 한다.
_wref = 'j_llrref' + _G
# 2026-09-29: 위 CAND 관문과 같은 이유로 VERSION 조건을 없앤다. 여기가 더 위험하다 —
#   j_llrref(33,903창) 와 j_llrref15(22,685창) 는 「파일 이름이 같다」.
if os.path.basename(os.path.dirname(REF_PATH)) != _wref:
    sys.exit('중단: 세대 불일치 — VERSION=%s 인데 LLR_REF=%s\n'
             '  러너가 LLR_REF 를 내보내는지 보십시오 (3_validate/steps/run_one_sample.sh).'
             % (C.VERSION, REF_PATH))

_self_proof()   # 경로·설정이 다 정해진 뒤에 찍는다
_a = sys.argv
LIMIT = int(_a[_a.index('--limit') + 1]) if '--limit' in _a else 0


def build_windows(keep_blocks):
    """창을 03_jsdcount 와 같은 순서로 만든다. wid 가 참조표와 맞아야 한다."""
    P = pd.read_parquet(CAND + '/blocks_cpg.parquet').sort_values(['chr', 'blk', 'pos'])
    gp = P.groupby(['chr', 'blk'], observed=True)['pos'].apply(list)
    blkwin = {}
    nwin = 0
    for (ch, bk), pos in gp.items():
        for s in range(len(pos) - K + 1):
            w = tuple(pos[s:s + K])
            if (str(ch), int(bk)) in keep_blocks:
                blkwin.setdefault((ch, bk), []).append((nwin, w))
            nwin += 1
    return blkwin, nwin


def load_reference(nwin, blkwin):
    """참조표를 읽어 (nwin, 2^K) 배열 둘과 have 마스크로 편다."""
    if not os.path.exists(REF_PATH):
        sys.exit('참조표가 없습니다: %s\n  _dmr/05_llrref.py 를 먼저 돌리세요.' % REF_PATH)
    R = pd.read_parquet(REF_PATH)
    need = ['wid', 'pos1', 'pos2', 'pos3'] + \
           ['n%d' % i for i in range(M_PAT)] + ['g%d' % i for i in range(M_PAT)]
    lack = [c for c in need if c not in R.columns]
    if lack:
        sys.exit('참조표에 없는 열: %s' % lack)

    wid = R['wid'].to_numpy()
    if wid.max() >= nwin:
        sys.exit('참조표의 wid 최댓값 %d 이 이 판의 창 수 %d 를 넘습니다.\n'
                 '  blocks_cpg.parquet 가 참조표를 만들 때와 다릅니다.' % (wid.max(), nwin))

    PN = np.zeros((nwin, M_PAT))
    PG = np.zeros((nwin, M_PAT))
    have = np.zeros(nwin, dtype=bool)
    PN[wid] = R[['n%d' % i for i in range(M_PAT)]].to_numpy(dtype=np.float64)
    PG[wid] = R[['g%d' % i for i in range(M_PAT)]].to_numpy(dtype=np.float64)
    have[wid] = True

    # ── 정합성 검사. 재구성한 창 좌표와 참조표의 pos 가 같아야 한다.
    ref_pos = {int(w): (int(a), int(b), int(c)) for w, a, b, c
               in zip(R['wid'], R['pos1'], R['pos2'], R['pos3'])} if K == 3 else {}
    if ref_pos:
        chk = bad = 0
        for lst in blkwin.values():
            for w, coords in lst:
                r = ref_pos.get(int(w))
                if r is None:
                    continue
                chk += 1
                if tuple(int(x) for x in coords) != r:
                    bad += 1
                    if bad <= 3:
                        print('  !! wid %d 좌표 불일치: 재구성 %s / 참조 %s'
                              % (w, tuple(int(x) for x in coords), r))
        if bad:
            sys.exit('창 번호가 참조표와 어긋납니다 (%d/%d). blocks_cpg 가 바뀐 것입니다.'
                     % (bad, chk))
        print('  창 좌표 정합성 %s개 확인: 모두 일치' % format(chk, ','))

    # ── 평활. 참조표는 날 확률이므로 여기서 넣는다.
    if EPS > 0:
        PN = (1.0 - EPS) * PN + EPS / M_PAT
        PG = (1.0 - EPS) * PG + EPS / M_PAT
        PN[~have] = 0.0
        PG[~have] = 0.0
    print('  참조표 %s  창 %s개 · 평활 eps=%g'
          % (os.path.basename(REF_PATH), format(int(have.sum()), ','), EPS))
    return PN, PG, have


def scan(bam, blkwin, nwin):
    """무늬 세기. 03_jsdcount · 2_train/steps/06 의 scan 과 같은 규칙.

    덤으로 **온타겟 조각 수와 고유 조각 수**를 센다 (2026-09-03 추가).
      이 파이프라인은 PCR 중복을 어디에서도 제거하지 않는다
      (is_duplicate · markdup · rmdup 어디에도 없다 — 실측 확인).
      중복은 새 정보가 아니므로 LLR 의 리드 합을 부풀린다.
      조각의 (chr, 시작, 끝) 이 같으면 같은 분자로 보고 고유 수를 센다.
      중복 배수 = 조각수 / 고유수. 1.0 에 가까우면 이 걱정은 닫힌다.
    """
    Cm = np.zeros((nwin, M_PAT), dtype=np.int32)
    seen = set()
    nfrag = [0]

    def flush(frag, span):
        if span is not None:
            nfrag[0] += 1
            seen.add(span)
        for k, d in frag.items():
            for wid, w in blkwin.get(k, ()):
                v = 0
                for p in w:
                    b = d.get(p)
                    if b is None:
                        v = -1
                        break
                    v = (v << 1) | b
                if v >= 0:
                    Cm[wid, v] += 1

    f = pysam.AlignmentFile(bam, 'rb', check_sq=False)
    last = None
    frag = {}
    span = None
    for r in f:
        if r.is_unmapped or r.is_secondary or r.is_supplementary:
            continue
        if r.query_name != last:
            if frag or span is not None:
                flush(frag, span)
            frag = {}
            span = None
            last = r.query_name
        ch = r.reference_name
        e = r.reference_end
        if e is None:
            continue
        lo = (r.reference_start + 1) // BLOCK
        hi = e // BLOCK + 1
        if not any((ch, b) in blkwin for b in range(lo, hi)):
            continue
        # 온타겟으로 확인된 뒤에만 조각 범위를 넓힌다 (짝 두 개를 하나로 본다)
        span = ((ch, r.reference_start, e) if span is None
                else (span[0], min(span[1], r.reference_start), max(span[2], e)))
        try:
            xm = r.get_tag('XM')
        except KeyError:
            continue
        for qp, rp in r.get_aligned_pairs(matches_only=True):
            c = xm[qp]
            if c != 'z' and c != 'Z':
                continue
            p = rp + 1
            k = (ch, p // BLOCK)
            if k in blkwin:
                frag.setdefault(k, {}).setdefault(p, 1 if c == 'Z' else 0)
    if frag or span is not None:
        flush(frag, span)
    f.close()
    return Cm, nfrag[0], len(seen)


def llr_row(Cm, PN, PG, have):
    """θ 마다 LLR 합. 창·무늬를 전부 더한다.

    ratio = P_GBM / P_Normal 을 미리 만들어 두면
        ln(P_mix/P_N) = ln(1 + θ·(ratio − 1))
    로 θ 마다 곱셈 한 번이면 된다. θ 별 참조를 저장할 이유가 없는 까닭이다.

    P_Normal 이 0 인 무늬는 그 항을 버린다(개수를 세어 보고한다).
    LLR 은 리드 단위 합이므로 창당 최소 리드 문턱을 두지 않는다 — 이게 JSD 와 다른 점이다.
    """
    use = have & (Cm.sum(1) > 0)
    out = {}
    if not use.any():
        for t in THETAS:
            out['LLR_%g' % t] = 0.0
            out['LLR_%g_pr' % t] = 0.0
        out.update({'LLR_max': 0.0, 'LLR_max_pr': 0.0, 'best_theta': 0.0,
                    '_nreads': 0, '_drop': 0})
        return out

    n = Cm[use].astype(np.float64)
    pn = PN[use]
    pg = PG[use]
    ok = pn > 0
    drop = int(n[~ok].sum())          # P_Normal=0 인데 리드가 있는 칸
    n = np.where(ok, n, 0.0)
    ratio = np.divide(pg, pn, out=np.zeros_like(pg), where=ok)
    d = np.where(ok, ratio - 1.0, 0.0)

    # ── 리드 수. LLR 은 리드 하나하나의 합이라 리드가 2배면 LLR 도 2배다.
    #   검체마다 패널에 걸리는 리드 수가 다르면 LLR 스케일이 검체마다 달라지고,
    #   그러면 종양 분획이 아니라 **패널 커버리지**를 먼저 재게 된다.
    #   그래서 리드당 값(_pr)을 함께 낸다. 어느 쪽이 쓸모 있는지는 2_train/steps/08 가 고른다.
    #   (JSD 는 창마다 합=1 로 정규화해 이 문제에 면역이지만 LLR 은 아니다.)
    nreads = float(n.sum())
    den = nreads if nreads > 0 else 1.0

    vals = []
    for t in THETAS:
        v = float((n * np.log1p(t * d)).sum())
        out['LLR_%g' % t] = v
        out['LLR_%g_pr' % t] = v / den
        vals.append(v)

    # θ=0 을 후보에 넣는다. LLR(0) ≡ 0 이므로 LLR_max 는 항상 0 이상이 되고,
    # 이것이 일반화우도비 검정(GLRT) 통계량 그 자체다.
    #   다만 그 때문에 정상도 LLR_max ≥ 0 이라 양수 쪽으로 몰린다.
    #   개별 LLR_<θ> 열이 부호를 그대로 갖고 있으므로 판별은 그쪽이 맡는다.
    cand = [0.0] + list(THETAS)
    cvals = [0.0] + vals
    i = int(np.argmax(cvals))
    out['LLR_max'] = float(cvals[i])
    out['LLR_max_pr'] = float(cvals[i]) / den
    out['best_theta'] = float(cand[i])
    out['_nreads'] = int(nreads)
    out['_drop'] = drop
    return out


def main():
    """theta 별 LLR 과 LLR_max·best_theta 를 만든다. 163줄.

    받는 것
      REF_PATH   j_llrref<GEN>/llrref_K3_minr<MINR>.parquet  (정상 기준 무늬 분포)
      CAND       j_candidates<GEN>  (창 목록 — 03_jsdcount 와 같은 순서여야 한다)
      BISMARK_OUT/mut_<n>_reads/  섞어 놓은 BAM

    내는 것
      ML_OUT/<비율>_llr_<판>.csv    열 = <chr>_<blk>_llr_theta<θ> · _llrmax · _besttheta

    구역
      1) 참조표 읽기   세대 관문 두 개가 여기서 돈다 (CAND · REF_PATH)
      2) 창 만들기     06_jsdfeat 와 같은 규칙. wid 가 참조표와 맞아야 한다
      3) 무늬 세기     BAM 을 훑어 K=3 무늬 빈도를 센다
      4) LLR 계산      theta 마다 로그우도비. theta=0 은 후보로만 참여
      5) 저장

    ※ LLR(0) ≡ 0 이므로 LLR_max 는 0 에서 바닥을 친다. 중앙값 0.00 은
      「없음」 이 아니라 「바닥」 이다.
    """
    C.guard_not_original()
    C.banner('2_train/steps/07  Multi-theta LLR   [VERSION=%s . PANEL=%s]' % (C.VERSION, C.PANEL))
    out_dir = C.ML_OUT
    save_dir = os.environ.get('LLR_OUT') or out_dir
    if save_dir != out_dir:
        os.makedirs(save_dir, exist_ok=True)
        print('  시험 모드: 결과를 %s 에 씁니다' % save_dir)

    means = sorted(glob.glob(out_dir + '/*_mean_%s.csv' % C.VERSION),
                   key=lambda p: int(re.search(r'[/\\](\d+)_mean_', p).group(1)))
    if not means:
        sys.exit('mean CSV 가 없습니다. step04_1 을 먼저 돌리세요: ' + out_dir)

    ref = pd.read_csv(means[0], index_col=0, nrows=1)
    keep = set()
    for c in ref.columns:
        ch, _, bk = str(c).rpartition('_')
        if bk.isdigit() and ch:
            keep.add((ch, int(bk)))
    if not keep:
        sys.exit('열 이름을 chr_blk 로 못 읽었습니다: %s' % list(ref.columns)[:3])

    print('  칸 %d개 · 창 K=%d · θ %d개' % (len(keep), K, len(THETAS)))
    print('  θ = %s' % ' · '.join('%g' % t for t in THETAS))

    blkwin, nwin = build_windows(keep)
    PN, PG, have = load_reference(nwin, blkwin)
    inpanel = sum(len(v) for v in blkwin.values())
    usable = int(have[[w for v in blkwin.values() for w, _ in v]].sum()) if inpanel else 0
    print('  패널 안 창 %s개 · 그중 참조 있는 것 %s개'
          % (format(inpanel, ','), format(usable, ',')))
    if usable == 0:
        sys.exit('패널 창 중 참조표에 있는 것이 없습니다. 패널·참조 조합을 확인하세요.')

    # ── 행 이름 ↔ BAM 번호 대응 검사.
    #   GBM{i} 의 i 는 step04_1 이 bedGraph 를 numeric 정렬해 붙인 번호다.
    #   여기서는 그 i 를 sampled_reads_{i}.bam 번호로 그대로 쓴다.
    #   복제본이 하나라도 빠지면 뒤 번호가 당겨져 **조용히 어긋난다** — 오류 없이 값만 틀린다.
    #   step04_1 이 남긴 매핑 파일로 대조한다 (:86 sample_file_mapping).
    mp = '%s/sample_file_mapping_%s.csv' % (out_dir, C.VERSION)
    if os.path.exists(mp):
        MP = pd.read_csv(mp)
        bad = []
        for sid, fn in zip(MP['sample_id'].astype(str), MP['file'].astype(str)):
            g = re.match(r'^(Normal|GBM)(\d+)$', sid)
            h = re.search(r'sampled_reads_(\d+)\.', fn)
            if g and h and int(g.group(2)) != int(h.group(1)):
                bad.append((sid, fn))
        if bad:
            for s, f in bad[:5]:
                print('  !! %s 가 %s 에 붙어 있습니다' % (s, f))
            sys.exit('행 이름과 BAM 번호가 어긋납니다 (%d건). '
                     'step04_1 의 sample_file_mapping 을 보고 맞춘 뒤 다시 도세요.' % len(bad))
        print('  행↔BAM 번호 대응 %s건 확인: 모두 일치' % format(len(MP), ','))
    else:
        print('  ⚠ 매핑 파일이 없어 대응을 검사하지 못했습니다: %s' % os.path.basename(mp))

    # 원값 · 리드당(_pr) · 고유조각당(_pu) 셋을 다 낸다.
    #   원값   실제 GLRT 통계량. 학습은 행마다 깊이가 같아 이쪽이 맞다
    #   _pr    검체마다 커버리지가 다를 때 (검증에서 필요)
    #   _pu    PCR 중복까지 보정. 중복 배수가 1 이면 _pr 과 같아진다
    #   어느 것이 뽑혔는지 반드시 보고할 것: 원값이 뽑히면 커버리지를 잰 것일 수 있다
    cols = ([c for t in THETAS
             for c in ('LLR_%g' % t, 'LLR_%g_pr' % t, 'LLR_%g_pu' % t)]
            + ['LLR_max', 'LLR_max_pr', 'LLR_max_pu', 'best_theta'])
    mix = C.MIX_OUT
    folders = os.listdir(mix)
    zero = C.folder_name(0.0)
    cache = {}
    diag = []

    def one(folder, i):
        bam = '%s/%s/sampled_reads_%d.bam' % (mix, folder, i)
        if not os.path.exists(bam):
            return None
        Cm, nfr, nuq = scan(bam, blkwin, nwin)
        v = llr_row(Cm, PN, PG, have)
        v['_nfrag'] = nfr
        v['_nuniq'] = nuq
        # 고유 조각으로 나눈 값. 중복이 있으면 이것만이 해석 가능하다.
        den = float(nuq) if nuq > 0 else 1.0
        for t in THETAS:
            v['LLR_%g_pu' % t] = v['LLR_%g' % t] / den
        v['LLR_max_pu'] = v['LLR_max'] / den
        return v

    for m in means:
        n = int(re.search(r'[/\\](\d+)_mean_', m).group(1))
        M = pd.read_csv(m, index_col=0)
        hit = [f for f in folders if f == 'mut_%d_reads' % n]
        rows, names = [], []
        for name in M.index:
            g = re.match(r'^(Normal|GBM)(\d+)$', str(name))
            if not g:
                sys.exit('행 이름을 못 읽었습니다: %s' % name)
            kind, i = g.group(1), int(g.group(2))
            if kind == 'Normal':
                if i not in cache:
                    cache[i] = one(zero, i)
                v = cache[i]
            else:
                v = one(hit[0], i) if hit else None
            rows.append(v)
            names.append(name)
            if LIMIT and len(rows) >= LIMIT:
                break

        F = pd.DataFrame(index=names, columns=cols, dtype=float)
        nr, dp, nf, nu = [], [], [], []
        for r_i, v in enumerate(rows):
            if not v:
                nr.append(0)
                dp.append(0)
                nf.append(0)
                nu.append(0)
                continue
            for c in cols:
                F.iat[r_i, F.columns.get_loc(c)] = v[c]
            nr.append(v['_nreads'])
            dp.append(v['_drop'])
            nf.append(v.get('_nfrag', 0))
            nu.append(v.get('_nuniq', 0))

        if 'Type' not in M.columns:
            sys.exit('mean CSV 에 Type 열이 없습니다: %s' % m)
        F['Type'] = M['Type'].reindex(F.index)

        fp = '%s/%d_llr_%s.csv' % (save_dir, n, C.VERSION)
        F.to_csv(fp)
        miss = 100.0 * float(F[cols].isna().mean().mean())
        _sf, _su = float(sum(nf)), float(sum(nu))
        dup = (_sf / _su) if _su > 0 else float('nan')
        print('    저장 %s  (%d행 · 결측 %.1f%% · 창리드 중앙값 %s · 중복배수 %.3f · P_N=0 버림 %s)'
              % (os.path.basename(fp), F.shape[0], miss,
                 format(int(np.median(nr)) if nr else 0, ','), dup,
                 format(int(sum(dp)), ',')))
        diag.append({'ratio_reads': n, 'rows': F.shape[0], 'na_pct': miss,
                     'reads_median': int(np.median(nr)) if nr else 0,
                     'frag_median': int(np.median(nf)) if nf else 0,
                     'uniq_median': int(np.median(nu)) if nu else 0,
                     'dup_ratio': round(dup, 4),
                     'dropped_reads': int(sum(dp)),
                     'best_theta_at_top_pct':
                         100.0 * float((F['best_theta'] >= max(THETAS)).mean()),
                     'best_theta_zero_pct':
                         100.0 * float((F['best_theta'] <= 0).mean())})

    D = pd.DataFrame(diag)
    dp_path = '%s/llr_diag_%s.csv' % (save_dir, C.VERSION)
    D.to_csv(dp_path, index=False)
    print('')
    print('  진단 저장 %s' % os.path.basename(dp_path))
    top = D['best_theta_at_top_pct'].max() if len(D) else 0.0
    if top >= 20.0:
        print('  ⚠ best_theta 가 격자 최댓값(%g)에 %.0f%% 몰렸습니다 — LLR_THETAS 로 넓히세요.'
              % (max(THETAS), top))
    if len(D) and D['dropped_reads'].sum() > 0:
        print('  ⚠ P_Normal=0 인 무늬의 리드 %s개를 버렸습니다. LLR_EPS 로 평활을 넣으면 살릴 수 있습니다.'
              % format(int(D['dropped_reads'].sum()), ','))

    print('')
    print('2_train/steps/07 완료.')


if __name__ == '__main__':
    sys.exit(main())
