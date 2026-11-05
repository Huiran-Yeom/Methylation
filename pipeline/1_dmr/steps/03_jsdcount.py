# -*- coding: utf-8 -*-
"""창별 리드 무늬 세기 (라벨 안 씀)
02 의 표준 CpG 자리로 창을 만들고 BAM 전체를 훑어 창을 다 덮은 리드의 무늬를 셈
창은 모든 검체 공통, 최소 리드 수 문턱은 여기서 안 걺, PE 두 짝은 한 조각으로 합침

사용: python 03_jsdcount.py [--nd <정상 cov 폴더>] [--cand ..] [--out ..] [--sample <검체>] [--limit N]
"""

# 스크립트 전용. import 하면 본체가 돌아 자료를 덮어씀
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


import os, re, sys, glob
import numpy as np, pandas as pd, pysam

# 산출 폴더 이름에 GEN(패널 세대) 붙임. 못 읽으면 중단
_GEN = _os.environ.get('GEN', '')
if not _GEN:
    try:
        _GEN = _cfg.GEN
    except Exception:
        raise SystemExit('GEN(패널 세대)을 못 읽었습니다 — config.conf 의 GEN 을 채우거나 '
                         '환경변수로 주십시오. 비워 두면 옛 세대 폴더를 가리킵니다.')
# 이 단계 경로. 폴더 구조 다르면 여기만 수정 (환경변수, 인자로 덮어쓰기 가능)
#   NORMAL_SET      정상 cov (검체 목록용)
#   CANCER_COV_DIR  암 cov (검체 목록용 · 옛 이름 GBM_COV_DIR)
#   NORMAL_BAM_DIR  정상 BAM (실제로 읽는 곳)
#   CANCER_BAM_DIR  암 BAM (실제로 읽는 곳 · 옛 이름 GBM_BAM_DIR)
CAND = _M + '/results/dmr/j_candidates' + _GEN
OUT  = _M + '/results/dmr/j_jsdcount' + _GEN
NB   = _os.environ.get('NORMAL_BAM_DIR') or _D + '/for_in_silico_test/normal_cfDNA_public/aligned_bam'
GB   = _cfg.cancer_env('CANCER_BAM_DIR', 'GBM_BAM_DIR') or _D + '/training/GBM_cell-line/aligned_bam'
ND   = _os.environ.get('NORMAL_SET')     or _M + '/sample_data/_input/normal_pub15'
GD   = _cfg.cancer_env('CANCER_COV_DIR', 'GBM_COV_DIR') or _D + '/training/GBM_cell-line/cov'
BLOCK, KS = 100, (3, 4, 5)
try:   # config.conf 의 BLOCK_SIZE 사용 (01, 02, 11 과 같은 격자)
    BLOCK = int(_cfg.BLOCK_SIZE)
except Exception:
    pass

a = sys.argv
if '--nd'   in a: ND   = a[a.index('--nd')+1]        # 학습 정상 폴더를 CLI 로
if '--cand' in a: CAND = a[a.index('--cand')+1]
if '--out'  in a: OUT  = a[a.index('--out')+1]
ONLY  = a[a.index('--sample')+1] if '--sample' in a else None
LIMIT = int(a[a.index('--limit')+1]) if '--limit' in a else 0

# 검체 이름, 단위 규칙은 sample_id.py (SID_RE, GROUP_RE)
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import sample_id as _sid
sid = _sid.sid

P = pd.read_parquet(CAND + '/blocks_cpg.parquet').sort_values(['chr','blk','pos'])
gp = P.groupby(['chr','blk'], observed=True)['pos'].apply(list)
wins   = {K: [] for K in KS}
blkwin = {}
for (ch, bk), pos in gp.items():
    for K in KS:
        for s in range(len(pos) - K + 1):
            w = tuple(pos[s:s+K])
            blkwin.setdefault((ch,bk), {}).setdefault(K, []).append((len(wins[K]), w))
            wins[K].append((ch, bk, s, w))
print('칸 %d개 · 창 %s' % (len(gp), ' · '.join('K%d %d' % (K, len(wins[K])) for K in KS)))
C = {K: np.zeros((len(wins[K]), 1 << K), dtype=np.int32) for K in KS}

def flush(frag):
    for k, d in frag.items():
        for K, lst in blkwin[k].items():
            arr = C[K]
            for wid, w in lst:
                v = 0
                for p in w:
                    b = d.get(p)
                    if b is None: v = -1; break
                    v = (v << 1) | b
                if v >= 0: arr[wid, v] += 1

def scan(bam):
    for K in KS: C[K][:] = 0
    f = pysam.AlignmentFile(bam, 'rb', check_sq=False)
    nr, nh = 0, 0

    def _limited(it):            # LIMIT 은 1차 정렬 리드 수 기준
        nonlocal nr
        for r in it:
            if r.is_unmapped or r.is_secondary or r.is_supplementary: continue
            nr += 1
            if LIMIT and nr > LIMIT: return
            yield r

    # 두 짝 묶기는 meth_config.fragments 사용 (좌표 정렬 BAM 에서도 짝이 묶임)
    st = {}
    for grp in _cfg.fragments(_limited(f), st):
        frag = {}
        for r in grp:
            ch, e = r.reference_name, r.reference_end
            if e is None: continue
            if not any((ch, b) in blkwin for b in range((r.reference_start+1)//BLOCK, e//BLOCK + 1)):
                continue
            nh += 1
            try: xm = r.get_tag('XM')
            except KeyError: continue
            for qp, rp in r.get_aligned_pairs(matches_only=True):
                c = xm[qp]
                if c != 'z' and c != 'Z': continue
                p = rp + 1
                k = (ch, p // BLOCK)
                if k in blkwin: frag.setdefault(k, {}).setdefault(p, 1 if c == 'Z' else 0)
        if frag: flush(frag)
    f.close()
    _cfg.dup_check(st.get('dup', 0), True, '03_jsdcount')
    return nr, nh

_CANCER = getattr(_cfg, 'CANCER_LABEL', 'GBM')   # samples.csv 의 type 값

# cov, BAM 짝은 아무것도 쓰기 전에 맞춰 봄 (중간에 멈추면 덮어쓰기 관문에 막힘)
jobs = [(sid(f), NB, 'Normal') for f in sorted(glob.glob(ND+'/*.cov.gz'))] + \
       [(sid(f), GB, _CANCER)    for f in sorted(glob.glob(GD+'/*.cov.gz'))]
if ONLY: jobs = [j for j in jobs if j[0] == ONLY]
if not jobs: sys.exit('대상 검체 없음')

# 검체마다 짝 BAM 찾기
_pair, _miss = [], []
for s, bd, ty in jobs:
    hit = [f for f in glob.glob(bd+'/*.bam') if sid(f) == s]
    if len(hit) == 1: _pair.append((s, hit[0], ty))
    else:             _miss.append((s, ty, len(hit), bd))

_have = {}
for _s, _f, _t in _pair: _have[_t] = _have.get(_t, 0) + 1
for _t in sorted({t for _, _, t in jobs}):
    _m = sum(1 for x in _miss if x[1] == _t)
    print('  %-8s 짝 %d개 · 못 찾음 %d개' % (_t, _have.get(_t, 0), _m))

if _miss:
    for _s, _t, _n, _bd in _miss[:20]:
        print('    %-16s %s BAM %d개  (%s)' % (_s, _t, _n, _bd))
    if len(_miss) > 20: print('    … 외 %d개' % (len(_miss) - 20))
    # 짝 없는 검체 있으면 중단 (04 가 samples.csv 기준이라 빠지면 다른 패널이 나옴)
    # 알고 넘기려면 ALLOW_SKIP=1
    if os.environ.get('ALLOW_SKIP') != '1':
        sys.exit('중단: 짝 BAM 을 못 찾은 검체가 %d개 있습니다.' % len(_miss) + chr(10) +
                 '  cov 이름과 BAM 이름이 첫 _ 앞까지 같아야 짝이 맞습니다.' + chr(10) +
                 '  정상 BAM: ' + NB + chr(10) +
                 '  암   BAM: ' + GB + chr(10) +
                 '  폴더가 다르면 NORMAL_BAM_DIR · CANCER_BAM_DIR 로 주십시오.' + chr(10) +
                 '  빠진 채로 진행하려면 ALLOW_SKIP=1 (04 의 정상 기준이 줄어듭니다).')
    # ALLOW_SKIP 이어도 정상 또는 암이 통째로 비면 중단
    _gone = [t for t in sorted({j[2] for j in jobs}) if not any(x[2] == t for x in _pair)]
    if _gone:
        sys.exit(('중단: %s 을 한 개도 짝짓지 못했습니다 (ALLOW_SKIP 으로도 넘길 수 없습니다).'
                  % ' · '.join(_gone)) + chr(10) +
                 '  정상 BAM: ' + NB + chr(10) + '  암   BAM: ' + GB + chr(10) +
                 '  폴더가 다르면 NORMAL_BAM_DIR · CANCER_BAM_DIR 로 주십시오.')
    print('  ALLOW_SKIP=1 — 빠진 채로 진행합니다. 04 의 정상 기준이 줄어듭니다.')

if not _pair: sys.exit('중단: 짝지은 검체가 하나도 없습니다.')

os.makedirs(OUT, exist_ok=True)
for K in KS:
    fp = OUT + '/windows_K%d.parquet' % K
    if not os.path.exists(fp):
        d = pd.DataFrame({'chr':[x[0] for x in wins[K]], 'blk':[x[1] for x in wins[K]],
                          'wi':[x[2] for x in wins[K]]})
        for j in range(K): d['pos%d' % (j+1)] = [x[3][j] for x in wins[K]]
        d.to_parquet(fp, index=False); print('창 표 저장 K%d' % K)


for i, (s, bam, ty) in enumerate(_pair, 1):
    nr, nh = scan(bam)
    msg = []
    for K in KS:
        tot = C[K].sum(1)
        nz  = int((tot > 0).sum())
        med = int(np.median(tot[tot > 0])) if nz else 0
        msg.append('K%d %d창(%.0f%%) 중앙%d' % (K, nz, 100.0*nz/max(1,len(tot)), med))
        cd = pd.DataFrame(C[K], columns=['p%d' % j for j in range(1 << K)])
        cd.insert(0, 'wid', np.arange(len(cd), dtype=np.int32))
        cd = cd[tot > 0]
        cd.to_parquet(OUT + '/count_%s_K%d.parquet' % (s, K), index=False)
    print('  [%2d/%d] %-16s 리드 %10d · 후보걸침 %9d · %s'
          % (i, len(_pair), s, nr, nh, ' | '.join(msg)), flush=True)
print('\n저장: ' + OUT)
