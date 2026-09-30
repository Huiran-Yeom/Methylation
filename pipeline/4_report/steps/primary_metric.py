# -*- coding: utf-8 -*-
"""사전등록 2절·3절의 1차 지표를 계산한다. 2026-09-14 · 2026-09-22 대폭 수정

  규칙은 사전등록 문서(이 저장소에 없음) 에 결과 보기 전에 적혀 있다.
    지표   외부 검증 검체 단위 AUC = primary_auc.py 의 a1 (검체·복제본 전체)
    검체   EGAD 환자 178검체 / 101명(양성) · 음성 배경 정상 10명 혼합
    모델   깊이 5k · 비율 0.5% 로 학습된 판의 pooled 모델
    조합   mean|entropy|jsd|llr|pdr  (5종 고정)
           + 병기 mean|entropy|jsd|pdr  (llr 뺀 것 — 플랫폼 조항 (1))
    교집합 세 패널이 모두 값을 낸 검증 검체만
    판정   짝지은 부트스트랩 2,000회 · 쌍이 셋이라 Bonferroni a=0.05/3 -> 98.33% 구간

   최고 조합을 고르지 않는다.  전수 표는 부차이고 여기서는 박은 조합만 본다.

  [2026-09-22 세 곳을 고쳤다. 앞의 둘은 조용히 틀린 값을 내고 있었다]

  (1) 읽는 곳이 죽은 경로였다. `~/tmp/results/` 는 09-16 혼합 이전 세대의 자리다.
      지금 원점수는 make_tables.py 와 같은 곳: 검증 산출 트리 안에 있다.
      같은 원본을 읽어야 표의 숫자와 여기 숫자가 어긋나지 않는다.

  (2) 음성 「공여자 층화」 가 유령 축이었다. 09-16 혼합 뒤로 음성 한 검체에
      정상 10명이 리드 수준에서 다 들어간다. 그런데 음성 행의 `rep` 은 여전히
      `Normal1`..`Normal10` 이고 이것은 공여자가 아니라 반복본 번호다.
      옛 코드는 그 숫자를 공여자로 읽어 EGAD=[1,2,3] · GSE=[4..10] 로 「층화」 했다.
      둘 다 안 비므로 경고도 안 떴다. 없는 축으로 재표본하고 층화했다고 적었다.
      09-22 에 넣은 가드도 `COHORT` 가 정의된 적이 없어 NameError -> except 로
      삼켜져 한 번도 발동하지 않았다. 문법이 맞는 것과 도는 것은 다르다.

      고침: 혼합판이면 음성 축을 아예 만들지 않는다. 양성만 사람 단위로
      재표본하고 음성은 고정한다(사전등록 3절 09-16 21:10 사용자 결정).
      그 대가를 출력에 매번 박는다: 산출 구간에 음성 축 불확실성이 0 이므로
      새 정상인에서의 위양성 변동은 이 구간에 없다.

  (3) 특이도의 EGAD/GSE 분리도 같은 유령 축이었다. 혼합판에서는 합본만 낸다.

사용:  primary_metric.py [--depth 5k] [--ratio 0.005] [--combo "..."] [--nboot 2000]
       COHORT 는 환경변수이고 「반드시」 주어야 한다 (기본값을 두지 않는다).
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


import os, sys, re, glob, io, itertools
import numpy as np, pandas as pd
from sklearn.metrics import roc_auc_score

a = sys.argv
DEPTH = a[a.index('--depth') + 1] if '--depth' in a else '5k'
RATIO = float(a[a.index('--ratio') + 1]) if '--ratio' in a else 0.005
COMBO = a[a.index('--combo') + 1] if '--combo' in a else 'mean|entropy|jsd|llr|pdr'
NBOOT = int(a[a.index('--nboot') + 1]) if '--nboot' in a else 2000
SEED = 42
try:   # 2026-09-29: SEED 가 config.conf 에 있는데 여기서만 박혀 있었다.
    SEED = int(_cfg.SEED)
except Exception:
    pass
# 기본값을 우리 내부 코호트 이름으로 두지 않는다. 안 주면 남이 받아
#   돌릴 때 「없는 코호트」 를 읽어 표가 다 비는데도 정상 종료했다.
COHORT = os.environ.get('COHORT') or ''
if not COHORT:
    raise SystemExit('COHORT 를 주십시오: 3단계에서 쓴 것과 같아야 합니다.'
                     '  예: COHORT=<your_cohort> python ' + os.path.basename(__file__))
LIVE = _M + '/%s_val' % COHORT
# 판 명단을 config.conf 한 곳에서 읽는다. 여기 복사본이 있으면
#   표·감시·귀무사슬이 서로 어긋날 수 있다.
PANELS = _cfg.PANELS

# 플랫폼: 배경 1~3 = EGAD(페어드엔드 48bp) · 4~10 = GSE301266(단일말단 80bp)
#  혼합판에서는 이 축이 존재하지 않는다. MIXED 가 참이면 쓰지 않는다. 
EGAD_REPS = (1, 2, 3)


def person(s):
    m = re.match(r'^((pat|ind)_\d+)_', str(s))
    return m.group(1) if m else str(s)


def donor(rep):
    """혼합 이전 판에서만 뜻이 있다. 음성 행의 rep 'Normal7' -> 7.
       혼합판에서는 이 숫자가 반복본 번호이지 공여자가 아니다. is_mixed() 로 가른다."""
    m = re.search(r'(\d+)$', str(rep))
    return int(m.group(1)) if m else 0


def val_dir(panel):
    return '%s/*/results/%s' % (LIVE, _cfg.version_name(panel, DEPTH))


def is_mixed():
    """혼합판이면 참. background_map CSV 의 `take` 열이 NEW_MIX 만 내는 서명이다.
       못 찾으면 참으로 본다: 모르는데 층화하는 쪽이 훨씬 나쁘다."""
    for p in PANELS:
        fs = glob.glob(val_dir(p) + '/step03_Preprocessing_for_ML/*/background_map*.csv')
        if fs:
            try:
                return 'take' in io.open(fs[0], encoding='utf-8').readline()
            except Exception:
                return True
    return True


def load(panel, combo):
    """make_tables.py 와 같은 원본을 읽는다. 검체마다 파일 하나."""
    v = _cfg.version_name(panel, DEPTH)
    pat = '%s/step04_ML_classifier/y_prob_all_%s_all_combos.csv' % (val_dir(panel), v)
    fs = sorted(glob.glob(pat))
    if not fs:
        return None, '채점 파일 없음: %s' % pat
    d = pd.concat([pd.read_csv(f) for f in fs], ignore_index=True)
    d = d[d.feature_set == combo]
    if not len(d):
        return None, '조합 %s 가 %d개 파일에 없음' % (combo, len(fs))
    d = d.copy()
    d['person'] = d['sample'].map(person)
    d['donor'] = d['rep'].map(donor) if 'rep' in d.columns else 0
    return d, None


def a1(pos, neg):
    """사전등록의 a1: 검체·복제본 전체."""
    if not len(pos) or not len(neg):
        return float('nan')
    y = np.r_[np.zeros(len(neg)), np.ones(len(pos))]
    s = np.r_[neg['y_score'].values, pos['y_score'].values]
    return float(roc_auc_score(y, s))


GATE = 0.70   # 사전등록 되돌림 ㉠: 셋 다 이것을 못 넘으면 순위를 매기지 않는다


def run(combo, mixed):
    """한 피처 조합에 대해 판별 표와 짝지은 부트스트랩 구간을 찍는다. 120줄.

    받는 것
      combo   '|' 로 이은 피처 조합 (예: mean|entropy|jsd|llr|pdr)
      mixed   혼합판인가. 참이면 공여자 층화를 하지 않는다
      전역    DEPTH · RATIO · NBOOT · SEED · PANELS (config.conf)

    읽는 것
      <뿌리>/<코호트>_val/*/results/<판>/step04_ML_classifier/
        y_prob_all_<판>_all_combos.csv      <- make_tables.py 와 「같은 원본」

    내는 것 (파일 없음. 표준출력만)
      판별 AUC 표 · 동점폭 · 쌍별 차이 구간 · 특이도

    구역
      1) 자료 읽기    판마다 load() · 없으면 이유를 적고 넘어간다
      2) 검체 맞추기  비율마다 세 판 공통 검체로 거른다 (사전등록 2-2)
      3) 지표         검체 단위 · 사람 단위 AUC · 동점폭
      4) 부트스트랩   양성 사람만 재표본 (음성 고정) · Bonferroni
      5) 특이도       합본만. 혼합판이라 플랫폼별로 못 가른다

    ※ 음성을 고정하므로 구간에 음성 축 불확실성이 0 으로 들어간다.
      새 정상인에서의 위양성 변동은 이 구간에 없다. 각주로 표에 박는다.
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
    common = set.intersection(*[set(d['sample']) for d in D.values()])
    print('  교집합 검체 %d개 · 사람 %d명 (사전등록 2-2)'
          % (len(common), len({person(s) for s in common})))

    out = {}
    for p in PANELS:
        d = D[p]
        d = d[d['sample'].isin(common)]
        pos = d[(d.y_true == 1) & (np.isclose(d.ratio, RATIO))]
        neg = d[d.y_true == 0]
        out[p] = (pos, neg, a1(pos, neg))
        print('  %-8s AUC %.4f  (양성 %d행 · 음성 %d행)'
              % (p, out[p][2], len(pos), len(neg)))

    # ── 사전등록 되돌림 ㉠ 을 여기서 판정한다 ───────────────────────────
    aucs = [out[p][2] for p in PANELS]
    if all(v == v for v in aucs) and max(aucs) < GATE:
        print('')
        print('  되돌림 ㉠ 발동: 셋 중 어느 것도 %.2f 을 못 넘었다 (최고 %.4f).'
              % (GATE, max(aucs)))
        print('    사전등록은 이 조건에서 순위를 매기지 않는다 고 결과 전에 정했다.')
        print('    아래 구간은 순위를 내려는 것이 아니라 구간을 보고하려고 낸다.')

    # ── 음성 재추출 축이 있는가 ──────────────────────────────────────────
    if mixed:
        print('')
        print('  [혼합판] 음성 한 검체에 정상 10명이 리드 수준에서 다 들어간다.')
        print('           음성 행의 rep(Normal1..10)은 반복본 번호이지 공여자가 아니다.')
        print('           -> 음성 공여자 층화를 하지 않는다. 양성만 사람 단위로 재표본한다.')
        print('           (사전등록 3절 09-16 21:10 사용자 결정 · Z절 09-22 고정)')
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
            keep = set(rs.choice(people, len(people), replace=True))
            kd = set()
            if not mixed:
                if eg and gs:
                    kd = set(rs.choice(eg, len(eg), replace=True)) | \
                         set(rs.choice(gs, len(gs), replace=True))
                elif dn:
                    kd = set(rs.choice(dn, len(dn), replace=True))
            dd = []
            for p in (x, y):
                pos, neg, _ = out[p]
                nn = neg if mixed else neg[neg['person'].isin(keep)]
                if kd:
                    nn = nn[nn['donor'].isin(kd)]
                dd.append(a1(pos[pos['person'].isin(keep)], nn))
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
        print('       음성 축 불확실성이 0 으로 들어간다. 새 정상인에서의')
        print('       위양성 변동은 이 구간에 없다.')

    # ── 특이도 ────────────────────────────────────────────────────────
    print('')
    if mixed:
        print('  특이도 (합본만: 혼합판이라 플랫폼별로 가를 수 없다)')
    else:
        print('  특이도 (플랫폼별 · 낮은 쪽을 주장값으로 쓴다 — 조항 (2))')
    for p in PANELS:
        neg = out[p][1]
        if 'y_pred' not in neg.columns:
            print('    %-8s y_pred 칸이 없어 못 잽니다' % p); continue
        def sp(sub):
            return float('nan') if not len(sub) else 1.0 - float(sub['y_pred'].mean())
        nvalid = neg[neg['y_pred'] >= 0]
        if mixed:
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
        print('  llr 은 플랫폼차 2.1702 대 암차 0.4816 (4.5배) 이고')
        print('  정상 EGAD -0.2989 · 정상 GSE +2.0087 로 부호가 갈린다.')
        print('  둘이 크게 다르면 그 차이가 llr 의 플랫폼 의존이다. 대표는 위 조합이다.')
        run(alt, mixed)
    return 0


if __name__ == '__main__':
    sys.exit(main())
