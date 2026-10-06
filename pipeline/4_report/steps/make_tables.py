
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
#               (앞판은 음성을 공여자로 묶어 분포를 뭉갰다. 지금은 안 묶는다)
import pandas as pd, numpy as np, glob, os, re, sys
import os
# 기본값을 우리 내부 코호트 이름으로 두지 않는다. 안 주면 남이 받아
#   돌릴 때 「없는 코호트」 를 읽어 표가 다 비는데도 정상 종료했다.
COHORT = os.environ.get('COHORT') or ''
if not COHORT:
    raise SystemExit('COHORT 를 주십시오: 3단계에서 쓴 것과 같아야 합니다.'
                     '  예: COHORT=<your_cohort> python ' + os.path.basename(__file__))
V = _M + '/%s_val' % COHORT
NTOT = len(glob.glob(V + '/*/results'))   # 분모를 178 로 박지 않는다
# 기본값 'j15bl200' 을 지웠다. 인자 없이 돌리면 조용히 특정
#   세대·특정 판을 본다. 07_panel 의 m_candidates 기본값과 같은 모양이다.
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

# 검체 이름에서 «사람»을 뽑는다. 사람 단위로 묶어 한 사람이 표를 여럿 쥐지 않게 한다.
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


def _allcombo(prefix):
    """all_combos 채점 파일을 찾는다.

    2026-09 까지 이 꼬리가 _모든조합.csv 였다. 한글->영문 이름바꿈을 코드만 하고
    저장된 산출물 쪽 짝을 안 지어서, 옛 트리에 대고 돌리면 한 칸도 못 읽는다.
    새 이름이 하나도 없을 때만 옛 이름을 본다. 끝을 맞춰 글롭하므로
    _모든조합_포화전.csv (중간 산출물) 는 안 걸린다.
    """
    fs = glob.glob(prefix + '_all_combos.csv')
    return fs or glob.glob(prefix + '_모든조합.csv')

def auc(y, s):
    y, s = np.asarray(y, float), np.asarray(s, float)
    npos, nneg = (y == 1).sum(), (y == 0).sum()
    if not npos or not nneg:
        return np.nan
    r = pd.Series(s).rank().values
    return (r[y == 1].sum() - npos * (npos + 1) / 2) / (npos * nneg)


def tie_width(y, s):
    """동점 때문에 정해지지 않는 AUC 폭. (최소, 최대, 고유값 수)

    2026-09-22: RandomForest 를 고른 칸에서 predict_proba 가 k/100 으로 양자화되고,
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

for depth in DEPTHS:
    # 태그 없는 단일모델 파일과 섞이지 않게 all_combos만 읽는다
    fs = _allcombo('%s/*/results/%s_%s_uni/step04_ML_classifier/y_prob_all_%s_%s_uni'
                   % (V, PAN, depth, PAN, depth))
    if not fs:
        print('\n== %s 리드: all_combos 채점 파일 없음' % depth); continue
    D = pd.concat([pd.read_csv(f) for f in fs], ignore_index=True)
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
