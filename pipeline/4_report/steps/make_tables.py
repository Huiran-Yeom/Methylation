
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


# 외부검증 표: 비율마다 그 비율의 양성/음성만 (1,780 대 1,780).
#   검체 단위 : 행 그대로
#   사람 단위 : 양성만 person 으로 묶는다(평균). 음성은 검체 단위로 묶는다.
#               aggregate_scores.py 의 [사람 단위] 집계와 같은 규칙.
#               (음성을 공여자로 묶으면 분포가 뭉개지므로 묶지 않는다)
import pandas as pd, numpy as np, glob, os, re, sys
# 기본값을 우리 내부 코호트 이름으로 두지 않는다. 두면 남이 받아
#   돌릴 때 「없는 코호트」 를 읽어 표가 다 빈 채 정상 종료한다.
COHORT = os.environ.get('COHORT') or ''
if not COHORT:
    raise SystemExit('COHORT 를 주십시오: 3단계에서 쓴 것과 같아야 합니다.'
                     '  예: COHORT=<your_cohort> python ' + os.path.basename(__file__))
V = _M + '/%s_val' % COHORT
NTOT = len(glob.glob(V + '/*/results'))   # 분모를 178 로 박지 않는다
# 판 이름에 기본값을 두지 않는다. 두면 인자 없이 돌릴 때 조용히 특정
#   세대·특정 판을 본다.
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
# 읽을 깊이. 자료를 바꾸면 DEPTHS='5k 50k' 처럼 준다.
DEPTHS = tuple((os.environ.get('DEPTHS') or '5k 15k 50k 100k').split())

# 검체 이름에서 「사람」을 뽑는다. 사람 단위로 묶어 한 사람이 표를 여럿 쥐지 않게 한다.
#   기본값은 이 논문의 이름 규칙(pat_1_... · ind_3_...)이다.
#   다른 자료를 쓰면 PERSON_RE 로 바꾼다. 괄호 하나가 사람 이름이 된다.
#   예: PERSON_RE='^([A-Z0-9]+)-'
_PRE = re.compile(os.environ.get('PERSON_RE') or r'^((pat|ind)_\d+)_')
_pmiss = [0, 0]


def person(s):
    m = _PRE.match(str(s))
    if m and m.groups():
        _pmiss[1] += 1
        return m.group(1)
    _pmiss[0] += 1
    return str(s)


def person_warn():
    """하나도 못 맞추면 사람 묶기가 검체 묶기와 같아진다. 조용히 넘기지 않는다."""
    if _pmiss[1] == 0 and _pmiss[0]:
        print('  ! PERSON_RE 가 검체 %d개 중 하나도 못 맞췄습니다 — 사람 단위 묶기가 검체 단위와'
              ' 같아집니다. 환경변수 PERSON_RE 를 자료의 이름 규칙에 맞추세요.' % _pmiss[0])


_GEN_SEEN = {}      # 판_깊이 -> (이름세대, 조합수)


def _width(f):
    """이 표가 쓰는 보고 9조합 중 그 파일에 몇 개가 있는지.

    담긴 조합 총수(31 이냐 9 냐)는 이 표를 치우치게 하지 않는다 — 이 표는
    COMBOS 를 한 줄씩 내고 조합 최대값을 쓰지 않는다. 해로운 것은 보고
    조합인데 칸이 아예 없는 경우다. 그러면 그 열은 판마다 다른 조합 집합을
    나란히 놓은 것이 된다.
    """
    try:
        have = set(pd.read_csv(f, usecols=['feature_set']).feature_set.unique())
    except Exception:
        return -1
    return len([c for c in COMBOS if c in have])


def _allcombo(prefix):
    """all_combos 채점 파일을 찾는다.

    옛 산출물은 꼬리가 _모든조합.csv 다. 새 이름만 찾으면 옛 트리에서는
    한 칸도 못 읽는다. 새 이름이 하나도 없을 때만 옛 이름을 본다. 끝을 맞춰 글롭하므로
    _모든조합_포화전.csv (중간 산출물) 는 안 걸린다.
    """
    _k = prefix.split('/')[-1][len('y_prob_all_'):]
    fs = glob.glob(prefix + '_all_combos.csv')
    if fs:
        _GEN_SEEN[_k] = ('all_combos', _width(fs[0]))
        return fs
    fs = glob.glob(prefix + '_모든조합.csv')
    if fs:
        # [기록] 이름만이 아니라 보고 조합이 몇 개 있는지까지 남긴다. 이름이
        #   다른 것은 그 자체로 해롭지 않다. 해로운 것은 보고 조합인데 칸이
        #   없는 것이다. 이름만 찍으면 읽는 사람이 「이름만 다르겠지」 로
        #   넘긴다.
        _GEN_SEEN[_k] = ('모든조합', _width(fs[0]))
        print('   [옛 이름으로 읽음] %s_모든조합.csv (%d개)' % (prefix.split('/')[-1], len(fs)))
    return fs


def gen_warn():
    """칸마다 읽은 산출물의 이름세대나 조합 수가 다르면 표 밑에 찍는다.
       막을 수는 없어도 모르고 넘어가는 것은 막는다."""
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
    """동점 때문에 정해지지 않는 AUC 폭. (최소, 최대, 고유값 수)

    RandomForest 를 고른 칸에서 predict_proba 가 k/100 으로 양자화되고,
      학습이 양성 90% 라 표가 몰려 3,300행에 고유값 10개 까지 떨어진다.
      같은 점수끼리는 누가 위인지 모르므로 AUC 가 한 점이 아니라 구간이 된다.
      보고값(auc)은 동점을 반반으로 센 구간 안의 한 점일 뿐이다.
      폭이 판 사이 격차보다 크면 순위를 주장할 수 없다. 그래서 같이 낸다.
      RandomForest 에는 decision_function 이 없어 df 전환으로 못 고친다.
        값을 바꾸는 대신 폭을 함께 보고한다(사전등록 불변).
    """
    y, s = np.asarray(y, float), np.asarray(s, float)
    P, N = (y == 1).sum(), (y == 0).sum()
    if not P or not N:
        return np.nan, np.nan, 0
    o = np.argsort(s, kind='mergesort')
    ys, ss = y[o], s[o]
    def w(t):                      # t=0 동점을 양성에 불리 · 1 유리
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

RAT_LABEL = ['0.1%', '0.5%', '1%', '2%', '2.5%', '3%', '5%']


def baseline_lines(D, depth, ratios):
    """REPORT_BASELINE=on 일 때 패널 표 밑에 붙는 대조 줄. 판정에 쓰지 않는다.

    값: 행(검체 · 비율 · 반복)마다, 패널 점수와 같은 깊이의 같은 뽑기로 만든 mean 특징 표
        (<비율>_mean_<판>.csv)의 블록 평균 메틸화를 다시 평균한 숫자 하나. 리드나 CpG 로 가중하지 않고,
        빈 블록은 행마다 빠진다.
    주의: 이것은 «전역» 평균 메틸화가 아니다. DMR 로 고른 패널 블록만의 평균이라 암 신호가 이미 들어 있다.
        사전등록 74-1 이 요구한 «같은 뽑기의 전역 평균» 은 뽑은 리드의 cov 가 남지 않아 낼 수 없어 이것으로
        대신했다. 그래서 플랫폼 · 전역 메틸화 교란을 재는 대조가 아니라, «패널 점수가 블록 평균 하나보다
        나은가» 를 보는 대조다.
    행: 1차 조합(COMBOS 의 마지막)의 행과 같다 — 교집합 거르기를 이미 거친 D 를 쓴다.
    방향: 같은 판 학습 mean 표(비율 0 을 뺀 전부)에서 한 번 정한다 — 암 행 평균 ≥ 정상 행 평균이면 높을수록 암.
        패널 점수의 방향도 학습 모델이 정하므로 같은 조건이다. 검증 라벨로 고르면 기준선이 유리해진다(사전등록 87-4).
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
    tr = [f for f in glob.glob('%s/results/panels/%s/%s/step04_ML_classifier/*_mean_%s.csv' % (_M, PAN, v, v))
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


for depth in DEPTHS:
    # 태그 없는 단일모델 파일과 섞이지 않게 all_combos만 읽는다
    fs = _allcombo('%s/*/results/%s_%s_uni/step04_ML_classifier/y_prob_all_%s_%s_uni'
                   % (V, PAN, depth, PAN, depth))
    if not fs:
        print('\n== %s 리드: all_combos 채점 파일 없음' % depth); continue
    D = pd.concat([pd.read_csv(f) for f in fs], ignore_index=True)
    D, _nu = _cfg.unify_score_unit(D)          # SCORE_UNIT=consistent 일 때만 바뀐다
    if _nu:
        print('  [점수 단위] 확률과 decision_function 이 섞인 조합 %d개를 decision_function 으로 맞췄다' % _nu)
    D['person'] = D['sample'].map(person)
    person_warn()
    # [교집합 · 비율별] 사전등록 2-2(94·256줄): 「세 패널이 모두 값을 낸 검증 검체만」.
    #   실측: 빠짐은 검체 전체가 아니라 (검체,비율) 칸 단위로 일어난다.
    #     5k   모든 비율 178/178/178            깨끗
    #     50k  2% 부터 갈림 · 5% 164/165/165
    #     100k 1% 부터 갈림 · 5% 110/111/120 → 교집합 99 (178의 56%)
    #   높은 비율에서 종양 리드를 채우려면 원본에 그만큼 리드가 있어야 하는데 못 채우는 검체가 있다.
    #   빠짐이 커버리지와 상관하므로 무작위가 아니다 → 전수로 내면 판마다 다른 환자 집단을 비교하게 된다.
    #   그래서 비율마다 세 판 공통 검체로 거른다.
    _ISEC = {}      # 비율 -> 공통 검체 집합
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
        print('   [교집합] 비율마다 세 판 공통만 씀 (공통/합집합) : %s — 사전등록 2-2'
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
        # [관문] 모델이 상수를 내뱉으면 동점 폭이 0 이라 아래 경고를 비켜 간다.
        #   고유값 1개면 양성·음성에 같은 점수를 줬다는 뜻이고 AUC 는 정확히 0.5 로
        #   「성능이 낮다」 처럼 보인다. 실제로 폐 100k 의 llr 단독에서 XGBoost 가
        #   840행 전부 0.9995 를 냈다. 낮은 성능이 아니라 못 쓰는 칸이다.
        _const = [i for i, (_, _, u) in enumerate(tw) if u == 1]
        if _const:
            print('%-16s %-5s   !! 상수 예측: 비율 %s 에서 점수 고유값이 1개다.'
                  '  이 줄의 AUC 는 성능이 아니라 구별 실패다.'
                  % ('', '', ','.join(RAT_LABEL[i] for i in _const)))

        # 동점 폭이 0.02 를 넘는 칸이 있으면 그 줄 밑에 구간을 적는다.
        #   0.02 는 「세 판 사이 격차의 최소 단위」 정도. 이보다 크면 순위 주장이 흔들린다.
        if any((h - l) > 0.02 for l, h, _ in tw if l == l):
            # 폭만 낸다. 구간 두 수를 9칸에 넣으면 붙어서 안 읽힌다
            print('%-16s %-5s %s   ↑동점폭' % ('', '',
                  ''.join(('%9.3f' % (h - l)) if (h - l) > 0.02 else '%9s' % '·'
                          for l, h, _ in tw)))
            _mu = min(u for _, _, u in tw if u)
            print('%-16s %-5s   고유값 최소 %d개: 이 줄은 동점이 많아 순위를 주장할 수 없다'
                  % ('', '', _mu))

    if getattr(_cfg, 'REPORT_BASELINE', 'off') == 'on':
        baseline_lines(D, depth, ratios)

gen_warn()
