#!/usr/bin/env python
# 2026-09-26: 09-23 사고에서 가드 없는 이 부류가 importlib 로 «확인» 되다 실제로 돌았다.
#   82개 파일이 덮였다(내용은 결정적이라 같았지만 그건 운이었다).
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
    # 검증이 이 파일을 검체 폴더로 «복사» 해 돌린다. 그때는 0_setup 이 없다.
    #   조용히 넘어가면 어느 뿌리를 썼는지 로그에 안 남는다. 알리고 넘어간다.
    _M = _os.environ.get('METH_ROOT', '/ssd_data/Methylation').rstrip('/')
    _D = _os.environ.get('DATASET_ROOT', '/ssd_data/dataset').rstrip('/')
    _sys.stderr.write('  [설정 없음] 환경변수/기본값으로 진행 — METH_ROOT=%s\n' % _M)


# coding: utf-8
"""2_train/steps/06 - JSD 피처. 회의 표의 JSD-2 에 쓰는 것.

무엇을 만드나
  칸마다 "이 판의 리드 무늬가 정상 기준에서 얼마나 떨어졌나" 를 숫자 하나로.
  mean · entropy 와 같은 모양(판 200개 x 칸 N개)이라 2_train/steps/08 가 그대로 읽는다.

    JSD(Q, P_N)
      Q    이 판의 창별 무늬 분포      믹싱 BAM 을 읽어 새로 센다
      P_N  정상 20명의 기준 무늬 분포   03_jsdcount 가 세어둔 것을 재사용

왜 BAM 을 다시 읽나
  cov 파일은 자리마다 "메틸 몇 / 비메틸 몇" 만 담는다. 리드 하나가 어떤
  무늬였는지(000 인지 101 인지)는 사라진다. 무늬는 BAM 에만 있다.

무늬 규칙은 03_jsdcount.py 와 완전히 동일하게 맞춘다.
  XM 태그 'Z'=1 · 'z'=0 · 자리는 rp+1 · 첫 자리가 최상위 비트
  query_name 이 같은 연속 리드를 한 조각으로 (PE 두 짝 병합)
  창의 K 자리를 다 덮은 조각만 센다
어긋나면 같은 판인데도 JSD 가 0 이 아니게 나와 값 전체가 무의미해진다.

출력 형식은 같은 판의 <N>_mean_<VERSION>.csv 를 그대로 따른다.
행 이름(Normal1..100 · GBM1..100)과 열(칸)을 복사해 쓰므로 어긋날 수 없다.

명령줄:  --limit N (판마다 앞 N개만, 시험용)
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
#   배치가 다르면 «여기만» 고칩니다. 0_setup/config.conf 에는 뿌리 둘
#   (METH_ROOT=_M · DATASET_ROOT=_D)과 실험 설정만 있습니다.
#   환경변수나 명령줄 인자로 그때그때 덮어쓸 수도 있습니다(아래 참조).
# ==== [해석 블록] 경로를 정한다. 검사는 아래 블록에서 한 번에 한다. ====
#   2026-09-14: 앞서 JCOUNT 가 세대 검사 **뒤에** 해석되어 검사를 비켜 갔다.
#   CAND 때와 같은 결함이다. 해석과 검사를 갈라 사이에 아무것도 못 끼게 한다.
CAND = _os.environ.get('METH_ROOT', _M) + '/results/dmr/j_candidates'
if '--cand' in sys.argv:
    CAND = sys.argv[sys.argv.index('--cand') + 1]
CAND = os.environ.get('JSD_CAND', CAND)

JCOUNT = _os.environ.get('METH_ROOT', _M) + '/results/dmr/j_jsdcount'
if '--jcount' in sys.argv:
    JCOUNT = sys.argv[sys.argv.index('--jcount') + 1]
JCOUNT = os.environ.get('JSD_JCOUNT', JCOUNT)
# ╚══════════════════════════════════════════════════════════════════╝

# ==== [검사 블록] 여기 위로 해석을 더 넣지 않는다. ====
#   인자를 빠뜨리면 **조용히 옛 풀**로 돈다. 실제로 복구 스크립트가 --cand 를
#   빠뜨려 정상 20명 기준으로 돌기 시작했다. 호출하는 쪽을 전부 맞추는 것은
#   계속 실패한다. 지나가는 잎에서 한 번 막는다.
# 2026-09-26: 관문이 'j15' 로 못 박혀 있어 **다음 세대에 스스로 꺼졌다.**
#   GEN=16 이면 VERSION 이 j16… 이 되고 이 if 가 아예 안 돈다. 보호 장치가 자기를 끈다.
#   또 «포함»(in)이라 j_candidates150 같은 것도 통과했다. «같음»으로 바꾼다.
# 2026-09-29: 조건을 없앤다. 'j'+GEN 으로 시작하는 VERSION 에서만 검사하면
#   이름을 다르게 지은 판(panelA_5k_uni 등)에서 관문이 또 스스로 꺼진다 —
#   09-26 에 'j15' 못 박기를 고치면서 남겨 둔 같은 모양의 구멍이다.
#   검사의 뜻은 «풀이 설정된 GEN 과 맞는가» 이고, 그건 VERSION 이름과 무관하다.
#   GEN 이 비면 'j_candidates' + '' = 'j_candidates' 라 옛 세대도 그대로 통과한다.
_G = getattr(C, 'GEN')
for _name, _val, _want in (('CAND', CAND, 'j_candidates' + _G),
                           ('JCOUNT', JCOUNT, 'j_jsdcount' + _G)):
    if os.path.basename(_val.rstrip('/')) != _want:
        sys.exit('중단: 세대 불일치 — VERSION=%s · GEN=%s 인데 %s=%s (기대 %s)'
                 % (C.VERSION, _G, _name, _val, _want)
                 + '\n  검증 러너가 JSD_CAND/JSD_JCOUNT 를 내보내는지 보십시오'
                   ' (3_validate/steps/run_one_sample.sh).')
# 2026-09-26: getattr 의 셋째 인자를 지웠다. 기본값이 곧 조용한 오답이다 —
#   MIN_READS_PER_WINDOW 의 기본값 4 는 09-14 에 6 으로 «철회된» 값이었다.
#   복사본이 config 를 못 읽으면 06 은 4 로 세고 07 은 minr6 참조표를 읽어
#   두 피처가 다른 문턱으로 계산된다. 파일이 따로라 아무 에러도 안 난다.
K = getattr(C, 'PATTERN_K')
# 08-21: 환경변수 JSD_MINR 로 덮어쓸 수 있게 했다. 안 주면 config 그대로다.
#        결측 70% 를 낮출 수 있는지 시험하려는 것 — 되돌릴 필요가 없는 방식.
MINR = int(os.environ.get('JSD_MINR', 0)) or getattr(C, 'MIN_READS_PER_WINDOW')
# 08-21 소영 확정: 결측이 이 값 이상인 칸은 버린다.
#   실측 — 200칸 중 136칸이 80~100% 결측이고 61칸은 10% 아래다. 가운데가 비어 있다.
#   버린 칸은 원래도 값이 없던 칸이라, 없는 것을 없다고 인정하는 것이다.
#   JSD_MAXNA=1 을 주면 안 거른다.
MAXNA = float(os.environ.get('JSD_MAXNA', 0.30))
# 08-24: JSD_COLS 에 학습이 고른 칸 목록 파일을 주면 그대로 쓴다.
#   검증에서 칸을 다시 고르면 학습과 열이 달라져 모델을 못 쓴다.
#   규칙 — 학습이 고른 것을 그대로 들고 간다. 검증 자료를 보고 정하지 않는다.
COLS_FILE = os.environ.get('JSD_COLS', '')
BLOCK = C.BLOCK_SIZE
_a = sys.argv
LIMIT = int(_a[_a.index('--limit') + 1]) if '--limit' in _a else 0

# 2026-09-14 [자기 증명] 실행 중인 코드가 자기 소스를 해시해 찍는다.
#   이 줄은 **이 코드를 실제로 돌려야만** 로그에 생긴다. 사후에 못 만든다.
#   해석된 입력 경로도 같이 찍는다 — 참조표 두 개가 파일 이름이 같아 로그로
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


_self_proof()



def _save_cols(save_dir, keep_cols, why):
    """학습이 쓴 칸 목록을 파일로 남긴다. 검증이 이것을 읽어 같은 칸을 쓴다.

    2026-09-29: 이 파일을 만드는 코드가 트리에 «없었다». 없으면 검증이
      JSD_COLS="" 로 내려가 검증 검체의 결측을 보고 칸을 다시 골랐고,
      BEST.joblib 의 열과 어긋난 자리를 채점이 학습 평균으로 메워 AUC 0.500 이 났다.
      이제 검증은 이 파일이 없으면 «멈춘다» — 그래서 거르든 안 거르든 남겨야 한다.
    """
    # 2026-09-29: 빈 목록을 쓰지 않는다. 0바이트 파일은 검증 관문을 통과해
    #   피처가 빈 채 완주하게 만든다 — 없느니만 못하다. 멈춰서 알린다.
    if not keep_cols:
        sys.exit('중단: 남은 칸이 0개입니다 (MAXNA=%s). 커버리지나 JSD_MAXNA 를 '
                 '보십시오 — 빈 칸 목록을 쓰면 검증이 빈 피처로 완주합니다.' % MAXNA)
    cf = '%s/jsd_cols_%s.txt' % (save_dir, C.VERSION)
    try:
        os.makedirs(save_dir, exist_ok=True)
        with open(cf, 'w', encoding='utf-8') as fh:
            fh.write(chr(10).join(keep_cols) + chr(10))
        print('    칸 목록 저장 (%s, %d칸): %s' % (why, len(keep_cols), cf))
    except OSError as e:
        print('    !! 칸 목록을 못 썼습니다 (%s) — 검증이 멈춥니다' % e)

def build_windows(keep_blocks):
    """창을 03_jsdcount 와 같은 순서로 만든다. wid 가 count 파일과 맞아야 한다."""
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


def build_pn(nwin):
    """정상 기준 분포. 검체마다 합=1 로 만든 뒤 평균 (한 검체 = 한 표)."""
    fs = [f for f in sorted(glob.glob(JCOUNT + '/count_*_K%d.parquet' % K))
          if 'SNU' not in os.path.basename(f)]
    if not fs:
        sys.exit('정상 count 파일을 못 찾았습니다: ' + JCOUNT)
    S = np.zeros((nwin, 1 << K))
    N = np.zeros(nwin)
    for f in fs:
        d = pd.read_parquet(f)
        wid = d['wid'].to_numpy()
        cnt = d[['p%d' % i for i in range(1 << K)]].to_numpy(dtype=np.float64)
        tot = cnt.sum(1)
        ok = tot >= MINR
        S[wid[ok]] += cnt[ok] / tot[ok, None]
        N[wid[ok]] += 1
    have = N > 0
    P = np.zeros_like(S)
    P[have] = S[have] / N[have, None]
    print('  P_N : 정상 %d명 · 쓸 수 있는 창 %s개' % (len(fs), format(int(have.sum()), ',')))
    return P, have


def scan(bam, blkwin, nwin):
    """무늬 세기. 03_jsdcount 의 scan 과 같은 규칙."""
    Cm = np.zeros((nwin, 1 << K), dtype=np.int32)

    def flush(frag):
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
    for r in f:
        if r.is_unmapped or r.is_secondary or r.is_supplementary:
            continue
        if r.query_name != last:
            if frag:
                flush(frag)
            frag = {}
            last = r.query_name
        ch = r.reference_name
        e = r.reference_end
        if e is None:
            continue
        lo = (r.reference_start + 1) // BLOCK
        hi = e // BLOCK + 1
        if not any((ch, b) in blkwin for b in range(lo, hi)):
            continue
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
    if frag:
        flush(frag)
    f.close()
    return Cm


def jsd_rows(Q, P):
    """행마다 JSD(Q, P). 밑 2 이므로 0~1."""
    M = 0.5 * (Q + P)
    qs = np.where(Q > 0, Q, 1.0)
    ps = np.where(P > 0, P, 1.0)
    ms = np.where(M > 0, M, 1.0)
    x = np.where(Q > 0, Q * np.log2(qs / ms), 0.0)
    y = np.where(P > 0, P * np.log2(ps / ms), 0.0)
    return 0.5 * x.sum(1) + 0.5 * y.sum(1)


def sample_jsd(bam, blkwin, nwin, P, have, wid2blk, blocks):
    Cm = scan(bam, blkwin, nwin)
    tot = Cm.sum(1)
    use = (tot >= MINR) & have
    if not use.any():
        return {}, {}
    Q = Cm[use] / tot[use, None]
    j = jsd_rows(Q, P[use])
    # 2026-09-14 PDR (Landau 2014) — 같은 Cm · 같은 use · 같은 창(K=3).
    #   학습쪽 2_train/steps/06 와 정의를 맞춘다. BAM 을 다시 읽지 않는다.
    d = 1.0 - (Cm[use, 0] + Cm[use, (1 << K) - 1]) / tot[use].astype(float)
    acc, acd = {}, {}
    for wid, val, dv in zip(np.nonzero(use)[0], j, d):
        b = wid2blk.get(int(wid))
        if b is not None:
            acc.setdefault(b, []).append(val)
            acd.setdefault(b, []).append(dv)
    return ({b: float(np.mean(v)) for b, v in acc.items()},
            {b: float(np.mean(v)) for b, v in acd.items()})


def main():
    C.guard_not_original()
    C.banner('2_train/steps/06  JSD feature   [VERSION=%s . PANEL=%s]' % (C.VERSION, C.PANEL))
    out_dir = C.ML_OUT
    # 08-21: JSD_OUT 을 주면 결과만 그 폴더에 쓴다. mean CSV 는 원래 자리에서 읽는다.
    #        시험용 — 진짜 파일을 덮어쓰지 않으려는 것.
    save_dir = os.environ.get('JSD_OUT') or out_dir
    if save_dir != out_dir:
        os.makedirs(save_dir, exist_ok=True)
        print('  시험 모드 — 결과를 %s 에 씁니다' % save_dir)

    means = sorted(glob.glob(out_dir + '/*_mean_%s.csv' % C.VERSION),
                   key=lambda p: int(re.search(r'[/\\](\d+)_mean_', p).group(1)))
    if not means:
        sys.exit('mean CSV 가 없습니다. step04_1 을 먼저 돌리세요: ' + out_dir)

    ref = pd.read_csv(means[0], index_col=0, nrows=1)
    cols = list(ref.columns)
    colblk = []
    keep = set()
    for c in cols:
        ch, _, bk = str(c).rpartition('_')
        if bk.isdigit() and ch:
            colblk.append((ch, int(bk)))
            keep.add((ch, int(bk)))
        else:
            colblk.append(None)
    if not keep:
        sys.exit('열 이름을 chr_blk 로 못 읽었습니다: %s' % cols[:3])
    print('  칸 %d개 · 창 K=%d · 창당 최소 리드 %d' % (len(keep), K, MINR))

    blkwin, nwin = build_windows(keep)
    P, have = build_pn(nwin)
    wid2blk = {}
    for k, lst in blkwin.items():
        for wid, _ in lst:
            wid2blk[wid] = (str(k[0]), int(k[1]))

    mix = C.MIX_OUT
    folders = os.listdir(mix)
    zero = C.folder_name(0.0)
    cache = {}
    keep_cols = None      # 첫 비율에서 정하고 나머지 비율에 그대로 쓴다
    if COLS_FILE:
        keep_cols = [x.strip() for x in open(COLS_FILE, encoding='utf-8') if x.strip()]
        print('  칸 목록을 파일에서 읽었습니다 — %d칸 (%s)' % (len(keep_cols), COLS_FILE))

    def one(folder, i):
        bam = '%s/%s/sampled_reads_%d.bam' % (mix, folder, i)
        if not os.path.exists(bam):
            return {}, {}
        return sample_jsd(bam, blkwin, nwin, P, have, wid2blk, keep)

    for m in means:
        n = int(re.search(r'[/\\](\d+)_mean_', m).group(1))
        M = pd.read_csv(m, index_col=0)
        idx = M.index
        hit = [f for f in folders if f == 'mut_%d_reads' % n]
        rows = []
        names = []
        for name in idx:
            g = re.match(r'^(Normal|GBM)(\d+)$', str(name))
            if not g:
                sys.exit('행 이름을 못 읽었습니다: %s' % name)
            kind = g.group(1)
            i = int(g.group(2))
            if kind == 'Normal':
                if i not in cache:
                    cache[i] = one(zero, i)
                v = cache[i]
            else:
                v = one(hit[0], i) if hit else {}
            rows.append(v)
            names.append(name)
            if LIMIT and len(rows) >= LIMIT:
                break
        # 2026-09-14: jsd 와 pdr 두 벌을 같은 규칙으로 낸다.
        #   keep_cols 는 jsd 에서 정해지고 pdr 이 그대로 쓴다 — 두 벌의 칸이 같아야 한다.
        for _slot, _tag in ((0, 'jsd'), (1, 'pdr')):
            F = pd.DataFrame(index=names, columns=cols, dtype=float)
            for r_i, _vv in enumerate(rows):
                v = _vv[_slot] if _vv else None
                if not v:
                    continue
                for c_i, b in enumerate(colblk):
                    if b is not None and b in v:
                        F.iat[r_i, c_i] = v[b]
            # 08-21: 블록이 아닌 열(Type 등)은 mean CSV 에서 값째로 옮긴다.
            #        이름만 복사하면 Type 이 비어 2_train/steps/08 의 astype(int) 가 터진다.
            for _c, _b in zip(cols, colblk):
                if _b is None and _c in M.columns:
                    F[_c] = M[_c].reindex(F.index)
            if MAXNA < 1.0:
                _blk = [x for x, y in zip(cols, colblk) if y is not None]
                _oth = [x for x, y in zip(cols, colblk) if y is None]
                if keep_cols is None:
                    _na = F[_blk].isna().mean()
                    keep_cols = [x for x in _blk if _na[x] < MAXNA]
                    print('    칸 거르기 — 결측 %.0f%% 미만 %d / %d 칸만 남깁니다'
                          % (100 * MAXNA, len(keep_cols), len(_blk)))
                    if not COLS_FILE:
                        _save_cols(save_dir, keep_cols, '결측 %.0f%% 미만' % (100 * MAXNA))
                # 08-24: 학습이 고른 칸이 검증 검체에 없을 수 있다 (커버리지가 얕다).
                #   F[...] 는 KeyError 를 내지만 reindex 는 없는 칸을 빈 값으로 만든다.
                #   빈 값은 뒤에서 학습셋 평균(train_mean)으로 채운다 — 검증 자료 평균이 아니다.
                _want = keep_cols + _oth
                _lack = [x for x in _want if x not in F.columns]
                if _lack:
                    print('    학습 칸 중 이 검체에 없는 것 %d개 — 빈 값으로 둡니다' % len(_lack))
                F = F.reindex(columns=_want)
            elif keep_cols is not None:
                # 2026-09-29: JSD_MAXNA=1 일 때도 «학습이 고른 칸» 을 적용한다.
                #   예전에는 적용(_want/reindex)이 `if MAXNA < 1.0` 안에만 있어서,
                #   JSD_COLS 를 읽어 놓고 쓰지 않았다 — 검증 열이 학습과 어긋나
                #   채점이 학습 평균으로 메우고 AUC 0.500 이 났다.
                _oth = [x for x, y in zip(cols, colblk) if y is None]
                _want = keep_cols + _oth
                _lack = [x for x in _want if x not in F.columns]
                if _lack:
                    print('    학습 칸 중 이 검체에 없는 것 %d개 — 빈 값으로 둡니다' % len(_lack))
                F = F.reindex(columns=_want)
            elif not COLS_FILE:
                # 2026-09-29: JSD_MAXNA=1 (거르지 않음)에서도 남긴다. 위 쓰기는
                #   `if MAXNA < 1.0` 안에 있어서, 안 거르면 파일이 안 생겼다 —
                #   그런데 검증은 이제 그 파일이 없으면 멈춘다. 문서에 있는
                #   JSD_MAXNA=1 을 쓰면 학습은 되고 검증은 못 도는 상태가 됐다.
                keep_cols = [x for x, y in zip(cols, colblk) if y is not None]
                _save_cols(save_dir, keep_cols, '거르지 않음')
            fp = '%s/%d_%s_%s.csv' % (save_dir, n, _tag, C.VERSION)
            F.to_csv(fp)
            print('    저장 %s  (%d행 x %d칸 · 결측 %.1f%%)'
                  % (os.path.basename(fp), F.shape[0], F.shape[1],
                     100.0 * float(F.isna().mean().mean())))

    print('')
    print('2_train/steps/06 완료.')


if __name__ == '__main__':
    sys.exit(main())
