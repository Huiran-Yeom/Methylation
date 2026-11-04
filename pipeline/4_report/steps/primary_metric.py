# -*- coding: utf-8 -*-
"""1차 지표 계산 (규칙은 결과 보기 전에 정함)

    지표    외부검증 검체 단위 AUC (검체, 복제본 전체)
    모델    깊이 5k 판의 pooled 모델, 채점은 0.5% 양성
    조합    mean|entropy|jsd|llr|pdr 고정 (+ llr 뺀 조합 병기). 최고 조합 고르지 않음
    교집합  세 판 모두 값 낸 행만
    판정    짝지은 부트스트랩 2000회, Bonferroni 0.05/3 → 98.33% 구간
            복원추출로 여러 번 뽑힌 사람은 뽑힌 횟수만큼 가중

    원점수는 make_tables.py 와 같은 검증 산출 트리에서 읽음
    혼합판은 음성 rep 번호가 공여자 아니라 반복본 번호라서 음성 고정, 양성 사람만 재표본
    (구간에 음성 쪽 불확실성 없음, 특이도도 합본만)

사용법:  primary_metric.py [--depth 5k] [--ratio 0.005] [--combo "..."] [--nboot 2000]
         COHORT 환경변수 필수 (기본값 없음)
"""
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


import os, sys, re, glob, io, itertools
import numpy as np, pandas as pd
from sklearn.metrics import roc_auc_score
from collections import Counter

a = sys.argv
DEPTH = a[a.index('--depth') + 1] if '--depth' in a else '5k'
RATIO = float(a[a.index('--ratio') + 1]) if '--ratio' in a else 0.005
COMBO = a[a.index('--combo') + 1] if '--combo' in a else 'mean|entropy|jsd|llr|pdr'
NBOOT = int(a[a.index('--nboot') + 1]) if '--nboot' in a else 2000
SEED = 42
try:   # config.conf 의 SEED 사용, 못 읽으면 42
    SEED = int(_cfg.SEED)
except Exception:
    pass
# COHORT 기본값 없음 - 없는 코호트 읽으면 빈 표로 정상 종료돼서
COHORT = os.environ.get('COHORT') or ''
if not COHORT:
    raise SystemExit('COHORT 를 주십시오: 3단계에서 쓴 것과 같아야 합니다.'
                     '  예: COHORT=<your_cohort> python ' + os.path.basename(__file__))
LIVE = _M + '/%s_val' % COHORT
# 판 명단은 config.conf 에서만 읽음
PANELS = _cfg.PANELS

# 플랫폼 구분용 배경 번호. 기본 1 2 3 = 이 GBM 자료의 EGAD 배경 (나머지는 GSE)
# 다른 자료는 EGAD_REPS 로 지정, 빈 값이면 플랫폼 구분 안 함 (특이도는 합본만)
# 혼합판에서는 안 씀
_er = os.environ.get('EGAD_REPS')
EGAD_REPS = (1, 2, 3) if _er is None else tuple(int(x) for x in _er.replace(',', ' ').split())


# 검체 이름에서 사람 뽑기 (한 사람이 표 여럿 쥐지 않게)
#   기본은 pat_1_..., ind_3_... 규칙. 다른 자료면 PERSON_RE 로, 괄호 하나 = 사람 이름
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
    """PERSON_RE 가 하나도 못 맞추면 경고 (사람 단위 = 검체 단위가 돼서)"""
    if _pmiss[1] == 0 and _pmiss[0]:
        print('  ! PERSON_RE 가 검체 %d개 중 하나도 못 맞췄습니다 — 사람 단위 묶기가 검체 단위와'
              ' 같아집니다. 환경변수 PERSON_RE 를 자료의 이름 규칙에 맞추세요.' % _pmiss[0])


def donor(rep):
    """음성 행 rep 'Normal7' → 7. 혼합 이전 판에서만 공여자 번호 (혼합판은 반복본 번호)"""
    m = re.search(r'(\d+)$', str(rep))
    return int(m.group(1)) if m else 0


def val_dir(panel):
    return '%s/*/results/%s' % (LIVE, _cfg.version_name(panel, DEPTH))


def is_mixed():
    """혼합판이면 참 (background_map CSV 에 take 열 있으면 혼합판)
       못 찾으면 참으로 봄 (모르고 층화하는 쪽이 더 나빠서)"""
    for p in PANELS:
        fs = glob.glob(val_dir(p) + '/step03_Preprocessing_for_ML/*/background_map*.csv')
        if fs:
            try:
                return 'take' in io.open(fs[0], encoding='utf-8').readline()
            except Exception:
                return True
    return True



_GEN_SEEN = {}      # 판_깊이 -> (이름세대, 1차 조합 유무 1/0/-1)


def _width(f):
    """1차 조합이 그 파일에 있으면 1, 없으면 0 (못 읽으면 -1)"""
    try:
        return int(COMBO in set(pd.read_csv(f, usecols=['feature_set']).feature_set.unique()))
    except Exception:
        return -1


def _allcombo(prefix):
    """all_combos 채점 파일 찾기

    없으면 옛 이름 _모든조합.csv 로 찾음 (_모든조합_포화전.csv 중간 산출은 안 걸림)
    """
    _k = prefix.split('/')[-1][len('y_prob_all_'):]
    fs = glob.glob(prefix + '_all_combos.csv')
    if fs:
        _GEN_SEEN[_k] = ('all_combos', _width(fs[0]))
        return fs
    fs = glob.glob(prefix + '_모든조합.csv')
    if fs:
        # 이름과 함께 1차 조합 유무도 기록 (gen_warn 에서 출력)
        _GEN_SEEN[_k] = ('모든조합', _width(fs[0]))
        print('   [옛 이름으로 읽음] %s_모든조합.csv (%d개)' % (prefix.split('/')[-1], len(fs)))
    return fs


def gen_warn():
    """칸마다 읽은 산출물의 이름세대나 1차 조합 유무가 다르면 표 밑에 경고"""
    if len(_GEN_SEEN) < 2:
        return
    gens = {g for g, _ in _GEN_SEEN.values()}
    wids = {w for _, w in _GEN_SEEN.values()}
    if len(gens) < 2 and len(wids) < 2:
        return
    if 0 in wids or -1 in wids:
        print('   !! 1차 조합이 없는 판이 있습니다 (%s). 그 판은 1차 지표를'
              ' 낼 수 없습니다 — 빠진 채로 순위를 매기지 마십시오.' % COMBO)
    else:
        print('   !! 이름세대 혼합 — 1차 조합(%s)은 모든 판에 있으나 읽은'
              ' 산출물의 세대가 다릅니다. 값은 비교할 수 있으나 재현 경로가'
              ' 다릅니다.' % COMBO)
    for _n in sorted(_GEN_SEEN):
        _g, _w = _GEN_SEEN[_n]
        print('        %-28s %-11s 1차조합 %s' % (_n, _g, '있음' if _w == 1 else '없음'))

def load(panel, combo):
    """make_tables.py 와 같은 원본 읽기 (검체마다 파일 하나)"""
    v = _cfg.version_name(panel, DEPTH)
    pat = '%s/step04_ML_classifier/y_prob_all_%s' % (val_dir(panel), v)
    fs = sorted(_allcombo(pat))
    if not fs:
        return None, '채점 파일 없음: %s' % pat
    d = pd.concat([pd.read_csv(f) for f in fs], ignore_index=True)
    d, _nu = _cfg.unify_score_unit(d)          # SCORE_UNIT=consistent 일 때만 바뀜
    d = d[d.feature_set == combo]
    if not len(d):
        return None, '조합 %s 가 %d개 파일에 없음' % (combo, len(fs))
    d = d.copy()
    d['person'] = d['sample'].map(person)
    person_warn()
    d['donor'] = d['rep'].map(donor) if 'rep' in d.columns else 0
    return d, None


def a1w(pos, neg, wp, wn):
    """a1 에 뽑힌 횟수로 가중 (복원추출 부트스트랩용)"""
    mp, mn = wp > 0, wn > 0
    if not mp.any() or not mn.any():
        return float('nan')
    y = np.r_[np.zeros(mn.sum()), np.ones(mp.sum())]
    sc = np.r_[neg['y_score'].values[mn], pos['y_score'].values[mp]]
    return float(roc_auc_score(y, sc, sample_weight=np.r_[wn[mn], wp[mp]]))


def a1(pos, neg):
    """a1: 검체, 복제본 전체 AUC"""
    if not len(pos) or not len(neg):
        return float('nan')
    y = np.r_[np.zeros(len(neg)), np.ones(len(pos))]
    s = np.r_[neg['y_score'].values, pos['y_score'].values]
    return float(roc_auc_score(y, s))


GATE = 0.70   # 판정 문턱: 셋 다 못 넘으면 순위 안 매김


def run(combo, mixed):
    """조합 하나의 판별 AUC, 짝지은 부트스트랩 구간, 특이도 출력 (표준출력만)

      combo   '|' 로 이은 피처 조합 (예: mean|entropy|jsd|llr|pdr)
      mixed   혼합판 여부, 참이면 음성 공여자 층화 안 함
    """
    print('')
    print('=' * 74)
    print('조합 %s   (깊이 %s · 비율 %g%%)' % (combo, DEPTH, RATIO * 100))
    print('=' * 74)
    D = {}
    for p in PANELS:
        d, err = load(p, combo)
        if err:
            print('  !! %s: %s' % (p, err))
            return 1
        D[p] = d
    # (검체, 비율, 반복, 라벨) 단위로 세 판 교집합 (검체 단위로만 맞추면 판마다 행 수 달라짐)
    _keys = ['sample', 'ratio', 'rep', 'y_true']
    _rows = set.intersection(*[set(d[_keys].itertuples(index=False, name=None)) for d in D.values()])
    common = {k[0] for k in _rows}
    _drop = {p: len(d) - sum(1 for k in d[_keys].itertuples(index=False, name=None) if k in _rows)
             for p, d in D.items()}
    print('  교집합 검체 %d개 · 사람 %d명 · 세 판 공통 행만 씀(판마다 뺀 행 %s)'
          % (len(common), len({person(s) for s in common}),
             ' · '.join('%s %d' % (p, n) for p, n in _drop.items())))

    out = {}
    for p in PANELS:
        d = D[p]
        d = d[[k in _rows for k in d[_keys].itertuples(index=False, name=None)]]
        pos = d[(d.y_true == 1) & (np.isclose(d.ratio, RATIO))]
        neg = d[d.y_true == 0]
        out[p] = (pos, neg, a1(pos, neg))
        print('  %-8s AUC %.4f  (양성 %d행 · 음성 %d행)'
              % (p, out[p][2], len(pos), len(neg)))

    # 판정 문턱(0.70) 확인
    aucs = [out[p][2] for p in PANELS]
    if all(v == v for v in aucs) and max(aucs) < GATE:
        print('')
        print('  순위 매기지 않음: 셋 중 어느 것도 %.2f 을 못 넘었다 (최고 %.4f).'
              % (GATE, max(aucs)))
        print('    이 기준은 결과를 보기 전에 정했다.')
        print('    아래 구간은 순위를 내려는 것이 아니라 구간을 보고하려고 낸다.')

    # 음성 공여자 재표본 가능 여부
    if mixed:
        print('')
        print('  [혼합판] 음성 한 검체에 배경 전원의 리드가 섞여 들어간다(배경 구성은 코호트 설정의 NORMALS*).')
        print('           음성 행의 rep(Normal1..10)은 반복본 번호이지 공여자가 아니다.')
        print('           -> 음성 공여자 층화를 하지 않는다. 양성만 사람 단위로 재표본한다.')
        dn, eg, gs = [], [], []
    else:
        dn = sorted({int(x) for x in pd.concat([out[p][1]['donor'] for p in PANELS]) if x})
        eg = [x for x in dn if x in EGAD_REPS]
        gs = [x for x in dn if x not in EGAD_REPS]
        print('  음성 공여자 %d명: EGAD %s · GSE %s' % (len(dn), eg, gs))
        if not eg or not gs:
            print('  (알림) 한쪽 플랫폼이 없습니다. 층화 없이 공여자 전체를 재표본합니다.')

    rs = np.random.RandomState(SEED)
    people = sorted({person(s) for s in common})
    print('')
    print('  쌍별 차이 (짝지은 부트스트랩 %d회 · 재표본 단위 %s · Bonferroni 98.33%%)'
          % (NBOOT, '양성 사람만 (음성 고정)' if mixed else '양성 사람 · 음성 공여자 층화'))
    for x, y in itertools.combinations(PANELS, 2):
        diffs = []
        for _ in range(NBOOT):
            keep = Counter(rs.choice(people, len(people), replace=True))
            kd = Counter()
            if not mixed:
                if eg and gs:
                    kd = Counter(rs.choice(eg, len(eg), replace=True)) + \
                         Counter(rs.choice(gs, len(gs), replace=True))
                elif dn:
                    kd = Counter(rs.choice(dn, len(dn), replace=True))
            dd = []
            for p in (x, y):
                pos, neg, _ = out[p]
                wp = pos['person'].map(keep).fillna(0).values
                wn = np.ones(len(neg)) if mixed else neg['person'].map(keep).fillna(0).values
                if kd:
                    wn = wn * neg['donor'].map(kd).fillna(0).values
                dd.append(a1w(pos, neg, wp, wn))
            if all(v == v for v in dd):
                diffs.append(dd[0] - dd[1])
        if not diffs:
            print('    %-8s - %-8s  계산 불가' % (x, y)); continue
        q = np.percentile(diffs, [0.835, 99.165])
        star = '유의' if (q[0] > 0 or q[1] < 0) else '0 포함'
        print('    %-8s - %-8s  %+.4f  [%+.4f, %+.4f]  %s'
              % (x, y, out[x][2] - out[y][2], q[0], q[1], star))
    if mixed:
        print('    ※ 각주(표에 같이 박는다): 음성을 고정했으므로 이 구간에는')
        print('       음성 축 불확실성이 0 으로 들어간다. 새 음성에서의')
        print('       위양성 변동은 이 구간에 없다.')

    # 특이도
    print('')
    if mixed:
        print('  특이도 (합본만: 혼합판이라 플랫폼별로 가를 수 없다)')
    elif not EGAD_REPS:
        print('  특이도 (합본만: EGAD_REPS 가 비어 플랫폼으로 가르지 않는다)')
    else:
        print('  특이도 (플랫폼별 · 낮은 쪽을 주장값으로 쓴다 — 조항 (2))')
    for p in PANELS:
        neg = out[p][1]
        if 'y_pred' not in neg.columns:
            print('    %-8s y_pred 칸이 없어 못 잽니다' % p); continue
        def sp(sub):
            return float('nan') if not len(sub) else 1.0 - float(sub['y_pred'].mean())
        nvalid = neg[neg['y_pred'] >= 0]
        if mixed or not EGAD_REPS:
            print('    %-8s 합본 %.4f  (위양성 %d개 / 음성 %d행%s)'
                  % (p, sp(nvalid), int((nvalid['y_pred'] == 1).sum()), len(nvalid),
                     '' if len(nvalid) == len(neg)
                     else ' · y_pred 없는 %d행 제외' % (len(neg) - len(nvalid))))
            continue
        se = sp(nvalid[nvalid['donor'].isin(EGAD_REPS)])
        sg = sp(nvalid[(~nvalid['donor'].isin(EGAD_REPS)) & (nvalid['donor'] > 0)])
        sa = sp(nvalid)
        lo = np.nanmin([se, sg])
        fp = nvalid[nvalid['y_pred'] == 1]
        sh = float('nan')
        if len(fp):
            sh = max(float(fp['donor'].isin(EGAD_REPS).mean()),
                     float((~fp['donor'].isin(EGAD_REPS)).mean()))
        print('    %-8s EGAD %.4f · GSE %.4f · 합본 %.4f  -> 주장 %.4f  '
              '(위양성 %d개 · 한쪽 쏠림 %.0f%%%s)'
              % (p, se, sg, sa, lo, len(fp), sh * 100 if sh == sh else float('nan'),
                 '  70% 초과: 합본으로 주장하지 않는다' if (sh == sh and sh > 0.70) else ''))
    return 0


def main():
    mixed = is_mixed()
    print('1차 지표: 코호트 %s · 깊이 %s · 비율 %g%% · %s'
          % (COHORT, DEPTH, RATIO * 100, '혼합판' if mixed else '혼합 이전 판'))
    rc = run(COMBO, mixed)
    if rc:
        return rc
    alt = '|'.join([c for c in COMBO.split('|') if c != 'llr'])
    if alt != COMBO:
        print('')
        print('[플랫폼 조항 (1)] llr 을 뺀 한 벌을 병기한다.')
        print('  둘이 크게 다르면 llr 의 기여가 크다는 뜻이다(플랫폼 · 대조 종류 · 참조 누수 · 신호를 가르지 않는다).')
        print('  대표는 위 조합이다.')
        run(alt, mixed)
    gen_warn()
    return 0


if __name__ == '__main__':
    sys.exit(main())
