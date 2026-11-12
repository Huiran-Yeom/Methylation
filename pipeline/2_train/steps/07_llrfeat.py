#!/usr/bin/env python
# 스크립트 전용 - import 하면 본체가 돌아 자료를 덮어씀
if __name__ != '__main__':
    raise ImportError(__file__ + ' 은(는) 스크립트다. import 하지 않는다')

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
    _sys.stderr.write('  [설정 없음] 환경변수/기본값으로 진행. METH_ROOT=%s\n' % _M)


"""2단계 07: 여러 θ 의 LLR 피처
  P_mix(θ)(x) = (1 - θ) P_Normal(x) + θ P_GBM(x)
  LLR(θ)      = 리드 합 ln[P_mix(θ)(x) / P_Normal(x)]   (x = 창 K 자리 메틸 무늬)
  열: LLR_<θ>, LLR_max (θ=0 포함 최댓값이라 0 이상), best_theta, Type (mean CSV 에서 복사)
  참조표는 1단계 05_llrref.py 결과를 그대로 읽음
    <METH_ROOT>/results/dmr/j_llrref/llrref_K3_minr6.parquet
  평활 없음이 기본 (LLR_EPS=0) - P_Normal(x)=0 인 무늬의 리드는 빼고, 뺀 수는 진단 파일에 남김
  무늬 규칙은 03_jsdcount.py, 06_jsdfeat.py 와 같음
  창 번호(wid)가 참조표 좌표와 다르면 바로 중단
환경변수: LLR_REF (참조 경로), LLR_THETAS (쉼표 목록), LLR_EPS (평활 상수), LLR_OUT (출력 폴더, 시험용)
사용: python 07_llrfeat.py [--cand <후보 폴더>] [--limit N]   (--limit: 비율마다 앞 N개만, 시험용)
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
CAND = _os.environ.get('METH_ROOT', _M) + '/results/dmr/j_candidates'
if '--cand' in sys.argv:            # 판마다 후보 풀이 다름
    CAND = sys.argv[sys.argv.index('--cand') + 1]
# 검증은 인자 없이 부르니 환경변수 LLR_CAND 로도 받음
CAND = os.environ.get('LLR_CAND', CAND)
# 세대 확인: 후보 폴더 이름이 정확히 j_candidates<GEN> 이어야 함
_G = getattr(C, 'GEN')
_want = 'j_candidates' + _G
# VERSION 이름과 상관없이 항상 확인
if os.path.basename(CAND.rstrip('/')) != _want:
    sys.exit('중단: 세대 불일치. VERSION=%s · GEN=%s 인데 CAND=%s (기대 %s)'
             % (C.VERSION, _G, CAND, _want))
REF_DEFAULT = _os.environ.get('METH_ROOT', _M) + '/results/dmr/j_llrref/llrref_K3_minr6.parquet'

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



# 기본값 없이 config 값 사용 - 06 과 같은 값이어야 함
K = getattr(C, 'PATTERN_K')
BLOCK = C.BLOCK_SIZE
M_PAT = 1 << K

# θ 목록 (비율 5% 까지 있어서 0.05 까지)
# θ=0 은 LLR 이 늘 0 이라 열로는 안 만들고 LLR_max, best_theta 후보로만 씀
THETAS_DEFAULT = [0.0005, 0.001, 0.002, 0.005, 0.01, 0.02, 0.05]

_env = os.environ.get('LLR_THETAS', '').strip()
THETAS = ([float(x) for x in _env.replace(' ', '').split(',') if x]
          if _env else list(THETAS_DEFAULT))
EPS = float(os.environ.get('LLR_EPS', 0.0))
REF_PATH = os.environ.get('LLR_REF', REF_DEFAULT)
# 참조표도 세대 확인 - 세대마다 파일 이름이 같아서 폴더 이름(j_llrref<GEN>)으로 봄
_wref = 'j_llrref' + _G
# VERSION 이름과 상관없이 항상 확인
if os.path.basename(os.path.dirname(REF_PATH)) != _wref:
    sys.exit('중단: 세대 불일치. VERSION=%s 인데 LLR_REF=%s\n'
             '  러너가 LLR_REF 를 내보내는지 보십시오 (3_validate/steps/run_one_sample.sh).'
             % (C.VERSION, REF_PATH))

_self_proof()   # 경로, 설정 다 정한 뒤에 찍음
_a = sys.argv
LIMIT = int(_a[_a.index('--limit') + 1]) if '--limit' in _a else 0


def build_windows(keep_blocks):
    """창을 03_jsdcount 와 같은 순서로 만듦 (wid 가 참조표와 맞아야 함)"""
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
    """참조표 읽어 P_Normal, P_GBM (nwin, 2^K) 배열과 have 마스크로 만듦"""
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

    # 창 좌표 확인: 다시 만든 창 좌표와 참조표 pos 가 같아야 함
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

    # 평활 (EPS > 0 일 때만): P = (1 - eps) * P + eps / 2^K
    if EPS > 0:
        PN = (1.0 - EPS) * PN + EPS / M_PAT
        PG = (1.0 - EPS) * PG + EPS / M_PAT
        PN[~have] = 0.0
        PG[~have] = 0.0
    print('  참조표 %s  창 %s개 · 평활 eps=%g'
          % (os.path.basename(REF_PATH), format(int(have.sum()), ','), EPS))
    return PN, PG, have


def scan(bam, blkwin, nwin):
    """창별 무늬 세기 (03_jsdcount, 06 과 같은 규칙)

    온타겟 조각 수와 고유 조각 수도 셈 - PCR 중복은 따로 안 지워서 LLR 합을 부풀릴 수 있음
    (chr, 시작, 끝) 같으면 같은 분자로 봄. 중복 배수 = 조각 수 / 고유 수
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

    # PE 두 짝을 한 조각으로 묶음 (03_jsdcount, 06 과 같은 meth_config.fragments)
    f = pysam.AlignmentFile(bam, 'rb', check_sq=False)
    st = {}
    for grp in _cfg.fragments(f, st):
        frag = {}
        span = None
        for r in grp:
            ch = r.reference_name
            e = r.reference_end
            if e is None:
                continue
            lo = (r.reference_start + 1) // BLOCK
            hi = e // BLOCK + 1
            if not any((ch, b) in blkwin for b in range(lo, hi)):
                continue
            # 온타겟일 때만 조각 범위 넓힘 (두 짝을 한 조각으로)
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
    _cfg.dup_check(st.get('dup', 0), C.ALLOW_REPLACEMENT, '07_llrfeat')
    return Cm, nfrag[0], len(seen)


def llr_row(Cm, PN, PG, have):
    """θ 마다 LLR 합 (모든 창, 무늬)

    ln(P_mix/P_N) = ln(1 + θ(ratio - 1)), ratio = P_GBM / P_Normal
    P_Normal 이 0 인 무늬는 빼고 개수만 셈. 리드 단위 합이라 창당 최소 리드 문턱 없음
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
    drop = int(n[~ok].sum())          # P_Normal=0 인데 리드 있는 무늬
    n = np.where(ok, n, 0.0)
    ratio = np.divide(pg, pn, out=np.zeros_like(pg), where=ok)
    d = np.where(ok, ratio - 1.0, 0.0)

    # LLR 은 리드 수에 비례해서 리드당 값(_pr)도 같이 냄 (샘플마다 커버리지 다를 때용)
    nreads = float(n.sum())
    den = nreads if nreads > 0 else 1.0

    vals = []
    for t in THETAS:
        v = float((n * np.log1p(t * d)).sum())
        out['LLR_%g' % t] = v
        out['LLR_%g_pr' % t] = v / den
        vals.append(v)

    # θ=0 도 후보에 넣음 → LLR_max 는 늘 0 이상 (GLRT 통계량)
    # 정상도 0 이상으로 몰려서 부호는 개별 LLR_<θ> 열에서 봄
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
    """θ 별 LLR 과 LLR_max, best_theta 행렬 만듦

    입력
      REF_PATH   j_llrref<GEN>/llrref_K3_minr<MINR>.parquet  (정상, 암 기준 무늬 분포)
      CAND       j_candidates<GEN>  (창 목록)
      MIX_OUT/mut_<n>_reads/  섞은 BAM

    출력
      ML_OUT/<n>_llr_<판>.csv   열 = LLR_<θ>, _pr, _pu, LLR_max, best_theta, Type
      ML_OUT/llr_diag_<판>.csv  비율별 진단

    LLR_max 는 0 이 바닥이라 중앙값 0 은 신호 없음이 아니라 바닥값
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
    # 검증은 블록 집합을 학습 mean 표에서 가져옴 (LLR_KEEP_FROM, 검증 mean 표는 블록을 안 거름)
    _kf = os.environ.get('LLR_KEEP_FROM')
    if _kf:
        _tm = sorted(glob.glob(_kf + '/*_mean_%s.csv' % C.VERSION),
                     key=lambda p: int(re.search(r'[/\\](\d+)_mean_', p).group(1)))
        if not _tm:
            sys.exit('LLR_KEEP_FROM 에 학습 mean CSV 가 없습니다: ' + _kf)
        print('  블록 집합: 학습 mean 표 %s' % os.path.basename(_tm[0]))
    keep_src = _tm[0] if _kf else means[0]     # 블록 집합만 학습에서, 행렬은 이 샘플 것

    ref = pd.read_csv(keep_src, index_col=0, nrows=1)
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

    # 행 이름 번호와 BAM 번호가 맞는지 확인 (sample_file_mapping 으로)
    # 반복 하나라도 빠지면 05 의 번호가 당겨져 조용히 어긋남
    mp = '%s/sample_file_mapping_%s.csv' % (out_dir, C.VERSION)
    if os.path.exists(mp):
        MP = pd.read_csv(mp)
        bad = []
        for sid, fn in zip(MP['sample_id'].astype(str), MP['file'].astype(str)):
            # Normal/GBM 은 05 가 붙인 행 이름 접두사 (암종 이름 아님)
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
        print('  !! 매핑 파일이 없어 대응을 검사하지 못했습니다: %s' % os.path.basename(mp))

    # 원값, 리드당(_pr), 고유 조각당(_pu) 셋 다 냄
    #   원값  GLRT 통계량 (학습은 깊이가 같아 이걸로 충분)
    #   _pr   샘플마다 커버리지 다를 때
    #   _pu   PCR 중복까지 보정
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
        # 고유 조각 수로 나눈 값 (중복 있을 때 사용)
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
        print('  !! best_theta 가 격자 최댓값(%g)에 %.0f%% 몰렸습니다. LLR_THETAS 로 넓히세요.'
              % (max(THETAS), top))
    if len(D) and D['dropped_reads'].sum() > 0:
        print('  !! P_Normal=0 인 무늬의 리드 %s개를 버렸습니다. LLR_EPS 로 평활을 넣으면 살릴 수 있습니다.'
              % format(int(D['dropped_reads'].sum()), ','))

    print('')
    print('2_train/steps/07 완료.')


if __name__ == '__main__':
    sys.exit(main())
