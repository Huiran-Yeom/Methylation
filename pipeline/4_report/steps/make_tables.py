
# 스크립트 전용 (import 하면 본체가 돎). python <파일> 로만 실행
if __name__ != '__main__':
    raise ImportError(__file__ + ' 은(는) 스크립트다 — import 하지 않는다')

import os as _os, sys as _sys
# 0_setup 찾기: 위 폴더로 올라가며 찾고, 없으면 METH_CONF_DIR
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


# 외부검증 표: 비율마다 그 비율의 양성, 음성만
#   검체 단위: 행 그대로
#   사람 단위: 양성만 사람으로 묶어 평균, 음성은 검체 단위 (aggregate_scores.py 와 같은 규칙)
import pandas as pd, numpy as np, glob, os, re, sys
# COHORT 기본값 없음 - 없는 코호트 읽으면 빈 표로 정상 종료돼서
COHORT = os.environ.get('COHORT') or ''
if not COHORT:
    raise SystemExit('COHORT 를 주십시오: 3단계에서 쓴 것과 같아야 합니다.'
                     '  예: COHORT=<your_cohort> python ' + os.path.basename(__file__))
V = _M + '/%s_val' % COHORT
NTOT = len(glob.glob(V + '/*/results'))   # 분모는 폴더 수에서 셈
# 판 이름 기본값 없음, 인자 필수
if len(sys.argv) < 2:
    sys.exit('사용: python make_tables.py <판이름>   예) python make_tables.py %s'
             '   ·  판 명단은 config.conf 의 PANELS · NULL_PANELS 에 있습니다: %s'
             % ('j' + _cfg.GEN + _cfg.PANELS[0],
                ', '.join(_cfg.PANELS + _cfg.NULL_PANELS)))
PAN = sys.argv[1]
COMBOS = ['mean|entropy', 'mean|entropy|jsd', 'llr', 'mean|entropy|llr',
          'mean|entropy|jsd|llr', 'pdr', 'mean|entropy|pdr',
          'mean|entropy|jsd|pdr', 'mean|entropy|jsd|llr|pdr']
SHORT = {'mean|entropy': 'ME', 'mean|entropy|jsd': 'ME+JSD', 'llr': 'LLR',
         'mean|entropy|llr': 'ME+LLR', 'mean|entropy|jsd|llr': 'ME+JSD+LLR', 'pdr': 'PDR',
         'mean|entropy|pdr': 'ME+PDR', 'mean|entropy|jsd|pdr': 'ME+JSD+PDR',
         'mean|entropy|jsd|llr|pdr': 'ME+JSD+LLR+PDR'}
# 읽을 깊이 (다른 자료면 DEPTHS='5k 50k' 처럼)
DEPTHS = tuple((os.environ.get('DEPTHS') or '5k 15k 50k 100k').split())

# 사람 뽑기 (PERSON_RE), all_combos 파일 찾기는 report_common.py
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from report_common import person, person_warn, find_allcombo   # noqa: E402


_GEN_SEEN = {}      # 판_깊이 -> (이름세대, 조합수)


def _width(f):
    """보고 9조합(COMBOS) 중 그 파일에 있는 수 (보고 조합 칸이 빠진 판 확인용)"""
    try:
        have = set(pd.read_csv(f, usecols=['feature_set']).feature_set.unique())
    except Exception:
        return -1
    return len([c for c in COMBOS if c in have])


def _allcombo(prefix):
    """all_combos 채점 파일 찾기 (옛 이름 _모든조합.csv 도)"""
    return find_allcombo(prefix, _width, _GEN_SEEN)


def gen_warn():
    """칸마다 읽은 산출물의 이름세대나 조합 수가 다르면 표 밑에 경고"""
    if len(_GEN_SEEN) < 2:
        return
    gens = {g for g, _ in _GEN_SEEN.values()}
    wids = {w for _, w in _GEN_SEEN.values()}
    if len(gens) < 2 and len(wids) < 2:
        return
    if len(wids) > 1:
        print('   !! 조합 집합 비대칭 — 판마다 보고 %d조합 중 가진 수가 다릅니다.'
              ' 빈 칸이 생기는 열은 판마다 다른 조합을 나란히 놓은 것입니다.'
              ' 그 열로 판 순위를 매기지 마십시오.' % len(COMBOS))
    else:
        print('   !! 이름세대 혼합 — 보고 조합은 %d개로 같으나 읽은 산출물의 세대가'
              ' 다릅니다. 값은 비교할 수 있으나 재현 경로가 다릅니다.' % wids.pop())
    for _n in sorted(_GEN_SEEN):
        _g, _w = _GEN_SEEN[_n]
        print('        %-28s %-11s 보고조합 %d/%d' % (_n, _g, _w, len(COMBOS)))

def auc(y, s):
    y, s = np.asarray(y, float), np.asarray(s, float)
    npos, nneg = (y == 1).sum(), (y == 0).sum()
    if not npos or not nneg:
        return np.nan
    r = pd.Series(s).rank().values
    return (r[y == 1].sum() - npos * (npos + 1) / 2) / (npos * nneg)


def tie_width(y, s):
    """동점 때문에 안 정해지는 AUC 폭 → (최소, 최대, 고유값 수)

    RandomForest 점수는 k/100 으로 양자화돼서 동점이 많음
    폭이 판 사이 차이보다 크면 순위 주장 불가라 같이 출력
    """
    y, s = np.asarray(y, float), np.asarray(s, float)
    P, N = (y == 1).sum(), (y == 0).sum()
    if not P or not N:
        return np.nan, np.nan, 0
    o = np.argsort(s, kind='mergesort')
    ys, ss = y[o], s[o]
    def w(t):                      # t=0 동점을 양성에 불리하게, 1 유리하게
        U = 0.0; seen = 0; i = 0
        while i < len(ss):
            j = i
            while j < len(ss) and ss[j] == ss[i]:
                j += 1
            blk = ys[i:j]; p = (blk == 1).sum(); n = (blk == 0).sum()
            U += p * seen + p * n * t
            seen += n; i = j
        return U / (P * N)
    return w(0.0), w(1.0), len(set(ss.tolist()))



def baseline_lines(D, depth, ratios):
    """REPORT_BASELINE=on 일 때 패널 표 밑에 붙는 대조 줄 (판정에 안 씀)

    값: 행(검체, 비율, 반복)마다 같은 뽑기의 mean 특징 표(<비율>_mean_<판>.csv) 블록 평균을 다시 평균
        (가중 없음, 빈 블록 제외)
    전역 평균 메틸화 아님 - 패널 블록만의 평균이라 패널 점수가 블록 평균 하나보다 나은지 보는 용도
    행: 교집합 거친 D 의 1차 조합(COMBOS 마지막) 행
    방향: 같은 판 학습 mean 표(비율 0 제외)에서 정함 (검증 라벨로 정하면 기준선이 유리해짐)
    """
    v = '%s_%s_uni' % (PAN, depth)
    base = D[D.feature_set == COMBOS[-1]]
    if base.empty:
        print('%-16s %-5s   (기준선: 1차 조합 행이 없어 내지 않는다)' % ('', ''))
        return
    parts, lack = [], 0
    for (s, p), g in base.groupby(['sample', 'mut_reads']):
        f = '%s/%s/results/%s/step04_ML_classifier/%d_mean_%s.csv' % (V, s, v, int(p), v)
        if not os.path.exists(f):
            lack += len(g)
            continue
        m = pd.read_csv(f, index_col=0)
        bm = m.drop(columns=['Type'], errors='ignore').mean(axis=1)
        parts.append(g.assign(_b=g['rep'].map(bm)))
    if not parts:
        print('%-16s %-5s   (기준선: mean 특징 표를 못 찾았다)' % ('', ''))
        return
    trd = '%s/results/panels/%s/%s/step04_ML_classifier' % (_M, PAN, v)
    tr = [f for f in glob.glob('%s/*_mean_%s.csv' % (trd, v))
          if not os.path.basename(f).startswith('0_')]
    if not tr:
        print('%-16s %-5s   (기준선: 학습 mean 표가 없어 방향을 못 정한다)' % ('', ''))
        return
    T = pd.concat([pd.read_csv(f, index_col=0) for f in tr])
    tb = T.drop(columns=['Type']).mean(axis=1)
    up = tb[T['Type'] == 1].mean() >= tb[T['Type'] == 0].mean()
    B = pd.concat(parts, ignore_index=True)
    nan = int(B['_b'].isna().sum()) + lack
    B = B.dropna(subset=['_b'])
    B['_b'] = B['_b'] if up else -B['_b']
    rs, rp = [], []
    for r in ratios:
        x = B[B.fraction_pct == r]
        rs.append(auc(x.y_true.values, x._b.values))
        gp = x[x.y_true == 1].groupby('person')['_b'].mean()
        gn = x[x.y_true == 0].groupby('sample')['_b'].mean()
        rp.append(auc(np.r_[np.ones(len(gp)), np.zeros(len(gn))], np.r_[gp.values, gn.values]))
    print('%-16s %-5s %s   검체' % ('기준선 블록평균', '-', ''.join('%9.3f' % v for v in rs)))
    print('%-16s %-5s %s   사람' % ('', '', ''.join('%9.3f' % v for v in rp)))
    print('%-16s %-5s   ↑ 대조 · 판정 아님: 같은 뽑기의 패널 블록 평균 메틸화 하나(전역 아님) · 방향 %s (학습 표로 정함)%s'
          % ('', '', '높을수록 암' if up else '낮을수록 암', (' · 값이 없어 뺀 행 %d' % nan) if nan else ''))
    # 방향 정한 학습 표 경로 출력 (다른 트리 모델로 채점했을 때 확인용)
    print('%-16s %-5s   방향을 정한 학습 표: %s (%d개)' % ('', '', trd, len(tr)))


for depth in DEPTHS:
    # all_combos 파일만 읽음 (단일모델 파일과 안 섞이게)
    fs = _allcombo('%s/*/results/%s_%s_uni/step04_ML_classifier/y_prob_all_%s_%s_uni'
                   % (V, PAN, depth, PAN, depth))
    if not fs:
        print('\n== %s 리드: all_combos 채점 파일 없음' % depth); continue
    D = pd.concat([pd.read_csv(f) for f in fs], ignore_index=True)
    D, _nu = _cfg.unify_score_unit(D)          # SCORE_UNIT=consistent 일 때만 바뀜
    if _nu:
        print('  [점수 단위] 확률과 decision_function 이 섞인 조합 %d개를 decision_function 으로 맞췄다' % _nu)
    D['person'] = D['sample'].map(person)
    person_warn()
    # 비율마다 세 판 모두 값 낸 검체만 사용
    #   깊이, 비율 높으면 리드 못 채워 빠지는 검체 있음 (무작위 아님, 판마다 집단 달라짐)
    _ISEC = {}      # 비율 -> 공통 검체 집합
    _HAS = {}       # 비율 -> 그 비율이 하나라도 있는 판 수
    _NPAN = 0
    _UNI = {}
    for _p in ['j' + _cfg.GEN + _x for _x in _cfg.PANELS]:
        _byr = {}
        for _f in _allcombo('%s/*/results/%s_%s_uni/step04_ML_classifier/y_prob_all_%s_%s_uni'
                            % (V, _p, depth, _p, depth)):
            try:
                _d = pd.read_csv(_f, usecols=['sample', 'feature_set', 'fraction_pct'])
                _d = _d[_d.feature_set == COMBOS[-1]]
                if not len(_d):
                    continue
                _sm = _d['sample'].iloc[0]
                for _r in _d['fraction_pct'].unique():
                    _byr.setdefault(round(float(_r), 4), set()).add(_sm)
            except Exception:
                pass
        for _r, _s in _byr.items():
            _UNI[_r] = _UNI.get(_r, set()) | _s
            _ISEC[_r] = _s if _r not in _ISEC else (_ISEC[_r] & _s)
        if _byr:                # 채점 파일 있는 판만 셈 (그 깊이에 안 돌린 판은 제외)
            _NPAN += 1
        for _r in _byr:
            _HAS[_r] = _HAS.get(_r, 0) + 1
    # 어떤 판에 그 비율이 통째로 없으면 공통 검체 없음으로 처리
    for _r in list(_ISEC):
        if _HAS.get(_r, 0) < _NPAN:
            _ISEC[_r] = set()
    _cutinfo = []
    if _ISEC:
        _keep = []
        for _r, _s in _ISEC.items():
            if len(_s) < len(_UNI.get(_r, _s)):
                _cutinfo.append('%g%% %d/%d' % (_r, len(_s), len(_UNI[_r])))
            _keep.append(D[(D.fraction_pct.round(4) == _r) & (D['sample'].isin(_s))])
        if _keep:
            D = pd.concat(_keep, ignore_index=True)
    ratios = sorted(D['fraction_pct'].unique())
    n1 = len(D[(D.feature_set == COMBOS[0]) & (D.y_true == 1) & (D.fraction_pct == ratios[0])])
    print('\n== %s 리드 · %s · %s · 채점 %d/%d · 비율당 양성/음성 각 %d'
          % (depth, PAN, COHORT, len(fs), NTOT, n1))
    if _cutinfo:
        print('   [교집합] 비율마다 세 판 공통만 씀 (공통/합집합) : %s'
              % ' · '.join(_cutinfo))
    print('%-16s %-5s %s' % ('피처', '모델', ''.join('%9s' % ('%g%%' % r) for r in ratios)))
    for c in COMBOS:
        sub = D[D.feature_set == c]
        if sub.empty:
            continue
        mdl = sub['model'].mode()
        rs, rp, tw = [], [], []
        for r in ratios:
            x = sub[sub.fraction_pct == r]
            rs.append(auc(x.y_true.values, x.y_score.values))
            _lo, _hi, _u = tie_width(x.y_true.values, x.y_score.values)
            tw.append((_lo, _hi, _u))
            gp = x[x.y_true == 1].groupby('person')['y_score'].mean()     # 양성만 사람으로
            gn = x[x.y_true == 0].groupby('sample')['y_score'].mean()     # 음성은 검체로
            rp.append(auc(np.r_[np.ones(len(gp)), np.zeros(len(gn))], np.r_[gp.values, gn.values]))
        print('%-16s %-5s %s   검체' % (SHORT.get(c, c), (mdl.iloc[0][:3] if len(mdl) else '?'),
                                        ''.join('%9.3f' % v for v in rs)))
        print('%-16s %-5s %s   사람' % ('', '', ''.join('%9.3f' % v for v in rp)))
        # 상수 예측 경고 - 고유값 1개면 AUC 0.5 지만 성능 아니라 구별 실패
        #   (동점 폭 0 이라 아래 동점폭 경고에는 안 걸림)
        _const = [i for i, (_, _, u) in enumerate(tw) if u == 1]
        if _const:
            print('%-16s %-5s   !! 상수 예측: 비율 %s 에서 점수 고유값이 1개다.'
                  '  이 줄의 AUC 는 성능이 아니라 구별 실패다.'
                  % ('', '', ','.join('%g%%' % ratios[i] for i in _const)))

        # 동점 폭 0.02 넘는 칸 있으면 그 줄 밑에 폭 출력 (판 사이 차이 최소 단위 정도)
        if any((h - l) > 0.02 for l, h, _ in tw if l == l):
            # 폭만 출력 (구간 두 수는 9칸에 안 들어감)
            print('%-16s %-5s %s   ↑동점폭' % ('', '',
                  ''.join(('%9.3f' % (h - l)) if (h - l) > 0.02 else '%9s' % '·'
                          for l, h, _ in tw)))
            _mu = min(u for _, _, u in tw if u)
            print('%-16s %-5s   고유값 최소 %d개: 이 줄은 동점이 많아 순위를 주장할 수 없다'
                  % ('', '', _mu))

    if getattr(_cfg, 'REPORT_BASELINE', 'off') == 'on':
        baseline_lines(D, depth, ratios)

gen_warn()
