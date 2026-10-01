# -*- coding: utf-8 -*-
"""S0a 후보 칸의 표준 CpG 자리 (2026-08-20) — 라벨을 쓰지 않는다.

JSD 창은 모든 검체에서 같아야 한다. 검체마다 관측된 자리로 창을 만들면
서로 다른 자리의 분포를 맞대게 되어 JSD 가 성립하지 않는다.
그래서 칸마다 '충분한 검체가 함께 보는 자리'를 먼저 확정한다.
기준은 후보 칸과 같다. 정상 60% AND 암 60%(종당 파일 MINREP개 이상).
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


import os, re, sys, glob, math
import numpy as np, pandas as pd

# ╔══ 이 단계가 쓰는 경로 ═══════════════════════════════════════════╗
#   배치가 다르면 「여기만」 고친다. 0_setup/config.conf 에는 뿌리 둘
#   (METH_ROOT=_M · DATASET_ROOT=_D)과 실험 설정만 둔다.
#   환경변수나 명령줄 인자로 그때그때 덮어쓸 수 있다(아래 참조).
#   이 단계가 읽는 자료 폴더는 환경변수로 옮길 수 있다.
#     NORMAL_SET      정상 cov (*.cov.gz)
#     GBM_COV_DIR     암 cov
ND   = _os.environ.get('NORMAL_SET')  or _M + '/sample_data/_input/normal_pub15'
GD   = _os.environ.get('GBM_COV_DIR') or _D + '/training/GBM_cell-line/cov'
# 산출 폴더는 GEN(패널 세대)을 따른다. 비우면 세대 없는 옛 폴더를 가리킨다.
_GEN = _os.environ.get('GEN', '')
if not _GEN:
    try:
        _GEN = _cfg.GEN
    except Exception:
        raise SystemExit('GEN(패널 세대)을 못 읽었습니다 — config.conf 의 GEN 을 채우거나 '
                         '환경변수로 주십시오. 비워 두면 옛 세대 폴더를 가리킵니다.')
CAND = _M + '/results/dmr/j_candidates' + _GEN
BLOCK, MINCOV, FRAC, MINREP = 100, 9, 0.60, 2
try:   # 2026-09-29: BLOCK_SIZE 가 config.conf 에 있는데 1단계는 100 을 박고 있었다.
    BLOCK = int(_cfg.BLOCK_SIZE)
except Exception:
    pass

a = sys.argv
if '--nd'   in a: ND   = a[a.index('--nd')+1]        # 2026-09-10 학습 정상 폴더를 CLI 로
if '--cand' in a: CAND = a[a.index('--cand')+1]
if '--depth' in a: MINCOV = int(a[a.index('--depth')+1]) - 1
OUT = a[a.index('--out')+1] if '--out' in a else CAND
# ╚══════════════════════════════════════════════════════════════════╝

# 파일 이름에서 검체 이름을 뽑는다. 기본은 첫 _ 앞까지 — cov 와 BAM 의 뒤쪽
#   꼬리표가 달라도 짝이 맞게 하기 위해서다.
#   자료마다 이름 규칙이 다르다. 안 맞으면 «자료를 고치지 말고» SID_RE 를 준다.
#   첫 괄호가 검체 이름이 되고, 안 걸리면 기본 규칙으로 떨어진다. 예:
#     SID_RE='(.+?)__'   H1876_0d__SRR…  -> H1876_0d   (밑줄 둘로 가르는 자료)
#                        H1876__SRR…     -> H1876
_SRE = _os.environ.get('SID_RE')
_srx = re.compile(_SRE) if _SRE else None
def sid(p):
    b = os.path.basename(p)
    for suf in ('.bismark.cov.gz', '.cov.gz'):
        if b.endswith(suf):
            b = b[:-len(suf)]; break
    if _srx is not None:
        m = _srx.match(b)
        if m and m.groups(): return m.group(1)
    return b.split('_')[0]
# 반복본을 독립단위(세포주 한 종)로 묶는 규칙. 커버 조건을 파일이 아니라
#   독립단위로 세기 때문에 이 규칙이 need_g 와 MINREP 에 바로 들어간다.
#   기본값은 이 연구의 GBM 자료 이름 규칙(SNU<n>-<반복>)이다. 다른 자료는
#   GROUP_RE 로 준다 — 첫 괄호가 독립단위 이름이 된다. 예:
#     GROUP_RE='([A-Za-z0-9]+)'      H1876_0d -> H1876 · H1876 -> H1876
#   안 주면 이름이 그대로 독립단위가 되어, 반복본이 서로 다른 종으로 세어진다.
_GRE = _os.environ.get('GROUP_RE')
_grx = re.compile(_GRE) if _GRE else None


def grp(s, t='GBM'):
    # GROUP_RE 는 «암 쪽에만» 쓴다. 반복본 묶기는 세포주 개념이고, 정상은
    #   한 사람당 한 파일이 전제다(아래 assert u_norm == n_norm 이 그것을 지킨다).
    #   정상까지 묶으면 NC-P-1·NC-P-3·… 아홉이 'NC' 하나가 되어 커버 기준이
    #   9명에서 5명으로 느슨해진다.
    if _grx is not None and t == 'GBM':
        m = _grx.match(s)
        return m.group(1) if (m and m.groups()) else s
    return s.split('-')[0] if s.upper().startswith('SNU') and '-' in s else s

files = [(f,'Normal') for f in sorted(glob.glob(ND+'/*.cov.gz'))] + \
        [(f,'GBM')    for f in sorted(glob.glob(GD+'/*.cov.gz'))]
if not files: sys.exit('cov 파일 없음')
# 같은 sid 가 둘 이상이면 멈춘다. sid 는 검체 이름이자 표의 열 이름이라,
#   겹치면 pivot_table(aggfunc='first') 이 뒤엣것을 «조용히 버린다».
_byid = {}
for _f, _t in files: _byid.setdefault(sid(_f), []).append(_os.path.basename(_f))
_dup = {k: v for k, v in _byid.items() if len(v) > 1}
if _dup:
    _msg = ['중단: 검체 이름이 겹칩니다 (첫 _ 앞까지가 이름입니다).']
    for _k, _v in sorted(_dup.items())[:10]:
        _msg.append('  %s : %s' % (_k, ' · '.join(_v)))
    if len(_dup) > 10: _msg.append('  … 외 %d개' % (len(_dup) - 10))
    _msg.append('  겹치면 표의 열이 덮여 한 쪽이 조용히 사라집니다.')
    _msg.append('  파일 이름의 첫 _ 앞을 서로 다르게 바꾸십시오.')
    sys.exit(chr(10).join(_msg))
_nn = sum(1 for f, _t in files if _t == 'Normal')
if _nn == 0:          sys.exit('정상 cov 가 없다 — NORMAL_SET 또는 --nd 를 확인하라: ' + ND)
if _nn == len(files): sys.exit('암 cov 가 없다 — GBM_COV_DIR 을 확인하라: ' + GD)

u_norm = len({grp(sid(f), 'Normal') for f,t in files if t=='Normal'})
u_gbm  = len({grp(sid(f)) for f,t in files if t=='GBM'})
need_n = int(math.ceil(FRAC*u_norm))
need_g = int(math.ceil(FRAC*u_gbm))
_nrep = {}
for _f,_t in files:
    if _t == 'GBM':
        _g = grp(sid(_f)); _nrep[_g] = _nrep.get(_g,0)+1
_req = {g: min(MINREP,n) for g,n in _nrep.items()}
print('정상 %d명 · 암 %d종 · 자리 기준 정상 %d AND 암 %d'
      % (u_norm, u_gbm, need_n, need_g))

B  = pd.read_csv(CAND + '/blocks_frac%03d.csv' % round(FRAC*100), dtype={'chr':str})
BM = pd.MultiIndex.from_arrays([B['chr'].values, B['blk'].values])
print('후보 칸 %d개' % len(BM))

rows = []
for i,(p,t) in enumerate(files,1):
    s = sid(p)
    d = pd.read_csv(p, sep='\t', header=None, usecols=[0,1,4,5],
                    names=['chr','start','tn','cn'], compression='gzip',
                    dtype={'chr':str,'start':np.int64,'tn':np.int32,'cn':np.int32})
    d = d[(d.tn + d.cn) > MINCOV]
    d['blk'] = d.start // BLOCK
    d = d[pd.MultiIndex.from_arrays([d['chr'].values, d['blk'].values]).isin(BM)]
    print('  [%2d/%d] %-16s %8d' % (i,len(files),s,len(d)))
    if d.empty: continue
    d = d[['chr','blk','start']].copy()
    d['grp'] = grp(s); d['ty'] = t
    rows.append(d)

A = pd.concat(rows, ignore_index=True)
del rows

N  = A[A.ty=='Normal'].groupby(['chr','blk','start'], observed=True)['grp'].nunique()
G  = A[A.ty=='GBM'   ].groupby(['chr','blk','start','grp'], observed=True).size()
Gv = np.array([_req[g] for g in G.index.get_level_values('grp')])
Go = pd.Series((G.values >= Gv).astype(np.int32),
               index=G.index).groupby(level=[0,1,2]).sum()
del A, G, Gv

ix = N.index.union(Go.index)
N  = N.reindex(ix, fill_value=0)
Go = Go.reindex(ix, fill_value=0)
keep = ix[(N.values >= need_n) & (Go.values >= need_g)]
print('\n표준 자리 %d개 / 관측 자리 %d개' % (len(keep), len(ix)))

P = pd.DataFrame({'chr': keep.get_level_values(0),
                  'blk': keep.get_level_values(1),
                  'pos': keep.get_level_values(2)}).sort_values(['chr','blk','pos'])
c = P.groupby(['chr','blk'], observed=True).size()
print('자리가 남은 칸 %d개 / 후보 %d개' % (len(c), len(BM)))
print('칸당 자리 수  최소 %d · 중앙 %d · 평균 %.1f · 최대 %d'
      % (c.min(), int(c.median()), c.mean(), c.max()))
print()
for K in (3,4,5):
    ok = c[c >= K]
    print('  K=%d   칸 %7d (%5.1f%%)   창 %9d'
          % (K, len(ok), 100.0*len(ok)/len(BM), int((ok-K+1).sum())))

os.makedirs(OUT, exist_ok=True)
P.to_parquet(OUT+'/blocks_cpg.parquet', index=False)
txt = ['표준 자리 %d / 칸 %d' % (len(P), len(c)),
       '기준 정상 %d명 AND 암 %d종 (종당 파일 %d개 이상)' % (need_n, need_g, MINREP),
       '칸당 자리 최소%d 중앙%d 평균%.1f 최대%d' % (c.min(), int(c.median()), c.mean(), c.max())]
for K in (3,4,5):
    ok = c[c >= K]
    txt.append('K=%d 칸 %d 창 %d' % (K, len(ok), int((ok-K+1).sum())))
open(OUT+'/blocks_cpg_summary.txt','w',encoding='utf-8').write('\n'.join(txt)+'\n')
print('\n저장: '+OUT+'/blocks_cpg.parquet')
