#!/usr/bin/env python
# 스크립트 전용 - import 하면 본체가 돌아 자료를 덮어씀
if __name__ != '__main__':
    raise ImportError(__file__ + ' 은(는) 스크립트다 — import 하지 않는다')

import os as _os, sys as _sys
# 0_setup 찾기: 위 폴더로 올라가며 → METH_CONF_DIR → 없으면 기본값 (검증 복사본용)
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
    # 설정 못 찾으면 환경변수/기본값 사용, 어느 뿌리 썼는지 로그에 남김
    _M = _os.environ.get('METH_ROOT', '/ssd_data/Methylation').rstrip('/')
    _D = _os.environ.get('DATASET_ROOT', '/ssd_data/dataset').rstrip('/')
    _sys.stderr.write('  [설정 없음] 환경변수/기본값으로 진행 — METH_ROOT=%s\n' % _M)


"""2단계 06: JSD 피처
  칸마다 이 판의 리드 무늬 분포가 정상 기준 분포에서 얼마나 떨어졌는지 (JSD 하나)
  Q   = 이 판의 창별 무늬 분포 (섞은 BAM 을 다시 읽어 셈, cov 에는 무늬가 없어서)
  P_N = 학습 정상 기준 분포 (03_jsdcount 결과 재사용)
  무늬 규칙은 03_jsdcount.py 와 같아야 함 (XM 'Z'=1, 'z'=0, PE 두 짝은 한 조각, 창 K 자리 다 덮은 조각만)
  출력은 같은 판의 <N>_mean_<VERSION>.csv 와 같은 행 / 열
사용: python 06_jsdfeat.py [--cand <후보 폴더>] [--jcount <jsdcount 폴더>] [--limit N]   (--limit: 판마다 앞 N개만, 시험용)
"""
import os
import sys
import glob
import re
import numpy as np
import pandas as pd
import pysam
import config as C

# 이 단계 경로. 폴더 구조 다르면 여기만 수정
# 명령줄 인자, 환경변수(JSD_CAND, JSD_JCOUNT) 순으로 덮어씀
CAND = _os.environ.get('METH_ROOT', _M) + '/results/dmr/j_candidates'
if '--cand' in sys.argv:
    CAND = sys.argv[sys.argv.index('--cand') + 1]
CAND = os.environ.get('JSD_CAND', CAND)

JCOUNT = _os.environ.get('METH_ROOT', _M) + '/results/dmr/j_jsdcount'
if '--jcount' in sys.argv:
    JCOUNT = sys.argv[sys.argv.index('--jcount') + 1]
JCOUNT = os.environ.get('JSD_JCOUNT', JCOUNT)

# 세대 확인: 경로 끝 폴더 이름이 정확히 j_candidates<GEN>, j_jsdcount<GEN> 이어야 함
# (경로 해석은 이 위에서 끝냄. 경로 안 넘기면 세대 없는 기본 경로라 여기서 중단)
_G = getattr(C, 'GEN')
for _name, _val, _want in (('CAND', CAND, 'j_candidates' + _G),
                           ('JCOUNT', JCOUNT, 'j_jsdcount' + _G)):
    if os.path.basename(_val.rstrip('/')) != _want:
        sys.exit('중단: 세대 불일치 — VERSION=%s · GEN=%s 인데 %s=%s (기대 %s)'
                 % (C.VERSION, _G, _name, _val, _want)
                 + '\n  검증 러너가 JSD_CAND/JSD_JCOUNT 를 내보내는지 보십시오'
                   ' (3_validate/steps/run_one_sample.sh).')
# 기본값 없이 config 값 사용 - 07 의 LLR 참조표와 문턱이 같아야 함
K = getattr(C, 'PATTERN_K')
# JSD_MINR 로 덮어쓰기 가능 (시험용), 안 주면 config 값
MINR = int(os.environ.get('JSD_MINR', 0)) or getattr(C, 'MIN_READS_PER_WINDOW')
# 결측 비율이 이 값 이상인 블록은 버림 (JSD_MAXNA=1 이면 안 거름)
MAXNA = float(os.environ.get('JSD_MAXNA', 0.30))
# JSD_COLS: 학습이 고른 블록 목록 파일 (검증은 이걸 그대로 씀, 다시 고르면 모델 열과 어긋남)
COLS_FILE = os.environ.get('JSD_COLS', '')
BLOCK = C.BLOCK_SIZE
_a = sys.argv
LIMIT = int(_a[_a.index('--limit') + 1]) if '--limit' in _a else 0

# 실행 중인 코드 md5, 입력 경로, 설정값을 로그에 찍음
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
    """학습이 쓴 블록 목록 저장. 검증이 읽어 같은 블록 사용 (없으면 검증 중단)"""
    # 빈 목록이면 중단 (빈 파일 쓰면 검증이 빈 피처로 끝까지 돎)
    if not keep_cols:
        sys.exit('중단: 남은 칸이 0개입니다 (MAXNA=%s). 커버리지나 JSD_MAXNA 를 '
                 '보십시오: 빈 칸 목록을 쓰면 검증이 빈 피처로 완주합니다.' % MAXNA)
    cf = '%s/jsd_cols_%s.txt' % (save_dir, C.VERSION)
    try:
        os.makedirs(save_dir, exist_ok=True)
        with open(cf, 'w', encoding='utf-8') as fh:
            fh.write(chr(10).join(keep_cols) + chr(10))
        print('    칸 목록 저장 (%s, %d칸): %s' % (why, len(keep_cols), cf))
    except OSError as e:
        print('    !! 칸 목록을 못 썼습니다 (%s): 검증이 멈춥니다' % e)

def build_windows(keep_blocks):
    """창을 03_jsdcount 와 같은 순서로 만듦 (wid 가 count 파일과 맞아야 함)"""
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
    """정상 기준 분포 P_N. 사람마다 합=1 로 만든 뒤 평균"""
    fs = sorted(glob.glob(JCOUNT + '/count_*_K%d.parquet' % K))
    # 정상, 암 구분은 samples.csv 의 type 열로 (파일 이름으로 추측 안 함)
    _sc = CAND + '/samples.csv'
    if not os.path.exists(_sc):
        sys.exit('samples.csv 가 없습니다 — 1단계 01_candidates 를 먼저 돌리세요: ' + _sc)
    _meta = pd.read_csv(_sc)
    for _c in ('sample', 'type'):
        if _c not in _meta.columns:
            sys.exit('samples.csv 에 %s 열이 없습니다: %s' % (_c, _sc))
    _norm = set(_meta.loc[_meta['type'] == 'Normal', 'sample'].astype(str))
    if not _norm:
        sys.exit('samples.csv 에 type=Normal 인 검체가 없습니다: ' + _sc)
    # held-out 이면 참조군 정상만 사용 (LLR 참조표와 같게, 섞기에 쓴 정상이 기준에 안 들어가게)
    if int(getattr(_cfg, 'LLR_HELDOUT', 0)):
        _ref, _ = _cfg.split_ref_mix(sorted(_norm))
        print('  P_N : held-out — 참조군 정상 %d/%d 명만 씀' % (len(_ref), len(_norm)))
        _norm = set(_ref)
    _suf = '_K%d.parquet' % K
    _nm = lambda f: os.path.basename(f)[len('count_'):-len(_suf)]
    _drop = [f for f in fs if _nm(f) not in _norm]
    fs = [f for f in fs if _nm(f) in _norm]
    if not fs:
        sys.exit('정상 count 파일을 못 찾았습니다 — samples.csv 의 정상 %d명과 이름이 겹치는 것이 없습니다: %s'
                 % (len(_norm), JCOUNT))
    if _drop:
        print('  P_N : 정상 아닌 %d개 제외' % len(_drop))
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
    print('  P_N : 정상 %d명 · 쓸 수 있는 창 %s개 · 빈 창 %s개(정상 누구도 %d 읽기 이상 없음)'
          % (len(fs), format(int(have.sum()), ','), format(int(nwin - have.sum()), ','), MINR))
    return P, have


def scan(bam, blkwin, nwin):
    """창별 무늬 세기 (03_jsdcount 의 scan 과 같은 규칙)"""
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

    # PE 두 짝을 한 조각으로 묶음 (03_jsdcount 와 같은 meth_config.fragments)
    f = pysam.AlignmentFile(bam, 'rb', check_sq=False)
    st = {}
    for grp in _cfg.fragments(f, st):
        frag = {}
        for r in grp:
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
    _cfg.dup_check(st.get('dup', 0), C.ALLOW_REPLACEMENT, '06_jsdfeat')
    return Cm


def jsd_rows(Q, P):
    """행마다 JSD(Q, P). 밑 2 라 0~1"""
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
    # PDR (Landau 2014): 같은 창 무늬 수에서 전부 메틸/전부 비메틸 아닌 비율
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
    # JSD_OUT 주면 결과만 그 폴더에 씀 (시험용, mean CSV 는 원래 자리에서 읽음)
    save_dir = os.environ.get('JSD_OUT') or out_dir
    if save_dir != out_dir:
        os.makedirs(save_dir, exist_ok=True)
        print('  시험 모드: 결과를 %s 에 씁니다' % save_dir)

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
    keep_cols = None      # 첫 비율에서 정하고 나머지 비율에 그대로 씀
    if COLS_FILE:
        keep_cols = [x.strip() for x in open(COLS_FILE, encoding='utf-8') if x.strip()]
        print('  칸 목록을 파일에서 읽었습니다. %d칸 (%s)' % (len(keep_cols), COLS_FILE))

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
        # jsd, pdr 두 벌 저장. 블록 목록은 jsd 에서 정하고 pdr 도 같은 것 사용
        for _slot, _tag in ((0, 'jsd'), (1, 'pdr')):
            F = pd.DataFrame(index=names, columns=cols, dtype=float)
            for r_i, _vv in enumerate(rows):
                v = _vv[_slot] if _vv else None
                if not v:
                    continue
                for c_i, b in enumerate(colblk):
                    if b is not None and b in v:
                        F.iat[r_i, c_i] = v[b]
            # 블록 아닌 열(Type 등)은 mean CSV 에서 값까지 복사 (비면 08 에서 오류)
            for _c, _b in zip(cols, colblk):
                if _b is None and _c in M.columns:
                    F[_c] = M[_c].reindex(F.index)
            if MAXNA < 1.0:
                _blk = [x for x, y in zip(cols, colblk) if y is not None]
                _oth = [x for x, y in zip(cols, colblk) if y is None]
                if keep_cols is None:
                    # 결측률 재는 행: all 은 전 행, normal 은 정상 행만
                    _normal = getattr(C, 'COLSET_FROM', 'all') == 'normal'
                    _isn = F.index.astype(str).str.startswith('Normal')
                    if _normal and not _isn.any():
                        sys.exit('COLSET_FROM=normal 인데 정상 행이 0개입니다 — 행 이름이 Normal 로 시작해야 합니다')
                    _na = F.loc[_isn if _normal else [True] * len(F), _blk].isna().mean()
                    keep_cols = [x for x in _blk if _na[x] < MAXNA]
                    if _normal:
                        # 정상은 덮이고 암 행에서 많이 비는 블록 수 (05 와 같은 기록)
                        _tna = F.loc[~_isn, keep_cols].isna().mean() if (~_isn).any() else None
                        print('    칸 거르기(normal): 암 행 결측 %.0f%% 이상인 남은 칸 %d개'
                              % (100 * MAXNA, 0 if _tna is None else int((_tna >= MAXNA).sum())))
                    print('    칸 거르기: 결측 %.0f%% 미만 %d / %d 칸만 남깁니다'
                          % (100 * MAXNA, len(keep_cols), len(_blk)))
                    if not COLS_FILE:
                        _save_cols(save_dir, keep_cols, '결측 %.0f%% 미만' % (100 * MAXNA))
                # 학습 블록이 검증 샘플에 없을 수 있어 reindex 로 빈 값 처리 (나중에 학습 평균으로 채움)
                _want = keep_cols + _oth
                _lack = [x for x in _want if x not in F.columns]
                if _lack:
                    print('    학습 칸 중 이 검체에 없는 것 %d개: 빈 값으로 둡니다' % len(_lack))
                F = F.reindex(columns=_want)
            elif keep_cols is not None:
                # JSD_MAXNA=1 이어도 JSD_COLS 의 학습 블록 목록 적용
                _oth = [x for x, y in zip(cols, colblk) if y is None]
                _want = keep_cols + _oth
                _lack = [x for x in _want if x not in F.columns]
                if _lack:
                    print('    학습 칸 중 이 검체에 없는 것 %d개: 빈 값으로 둡니다' % len(_lack))
                F = F.reindex(columns=_want)
            elif not COLS_FILE:
                # 안 걸러도 블록 목록은 저장 (없으면 검증 중단)
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
