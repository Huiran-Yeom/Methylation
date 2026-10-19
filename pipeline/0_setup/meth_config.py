# -*- coding: utf-8 -*-
"""config.conf 를 읽는다. 고칠 곳은 옆의 config.conf 하나다.

보는 순서는 환경변수 → config.conf → 기본값이다.

    import meth_config as cfg
    cfg.METH_ROOT            '/ssd_data/Methylation'
    cfg.TARGET_DEPTH_READS   5000
    cfg.MUT_RATIOS           [0.0, 0.001, ...]
    cfg.snapshot('<output folder>')     해결된 값과 그 출처를 산출물 옆에 남긴다
"""
import datetime
import os
import re
import sys

# 코드가 있는 자리에서 기본값을 만든다. 저장소 안을 가리키면 내려받은 그대로 동작한다.
#   서버 고유 경로를 박아 두면 검증이 「학습 스크립트를 못 찾은 채」 돌아 그 피처만 빠진다.
_TREE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

_DEFAULTS = {
    'METH_ROOT': '/ssd_data/Methylation',
    'DATASET_ROOT': '/ssd_data/dataset',
    'WORK_ROOT': os.path.expanduser('~'),
    'TMP_ROOT': os.path.expanduser('~/tmp'),
    'TRAIN_SRC': os.path.join(_TREE, '2_train', 'steps'),
    'ANALYSIS_DIR': os.path.join(_TREE, '3_validate', 'steps'),
    'GEN': '15',
    # samples.csv 의 type 열에 적히는 「암」 쪽 이름. 쓰는 곳과 읽는 곳이 모두
    #   이 값을 본다. 바꾸면 1단계를 처음부터 다시 돌려야 한다.
    'CANCER_LABEL': 'GBM',
    # 빈 값이면 「어느 코호트든 허용」. 이름을 적으면 그 목록만 허용한다.
    'VAL_COHORTS': '',
    'PANELS': 'bl200,jsd200,jsdb200',
    'NULL_PANELS': 'rand1200,rand2200,rand3200,randb1200,randb2200,randb3200',
    'VERSION': 'j15bl200_5k_uni',
    'PANEL': 'bl200',
    'TARGET_DEPTH_READS': '5000',
    'DEPTH_MODE': 'total',
    'MUT_RATIOS': '0.0,0.001,0.005,0.01,0.02,0.025,0.03,0.05,0.1,1.0',
    'PATTERN_K': '3',
    'MIN_READS_PER_WINDOW': '6',
    'BLOCK_SIZE': '100',
    'MIN_SAMPLE_COVERAGE': '0.8',
    'COLSET_FROM': 'all',
    'BLAS_THREADS': '1',
    'SCORE_UNIT': 'as_is',
    'REPORT_BASELINE': 'off',
    'CONVERGENCE_GATE': 'off',
    'SEED': '42',
    'TEST_SIZE': '0.2',
    'N_JOBS': '1',
    # 비우면 모델 넷을 교차검증으로 겨뤄 1등.
    'FIXED_MODEL': '',
    # 0 이면 LLR 참조표와 학습 풀에 같은 사람이 들어간다.
    'LLR_HELDOUT': '0',
}
_CACHE = None
SOURCE = {}


def _conf():
    """config.conf 를 찾아 읽는다. 못 찾으면 멈춘다.

    조용히 기본값으로 도는 것을 막는다. 그러면 「어느 설정으로 돌았나」 를 알 방법이 없다.
    일부러 기본값만 쓰려면 METH_CONF_OPTIONAL=1 을 준다.
    """
    global _CACHE
    if _CACHE is not None:
        return _CACHE
    _CACHE = {}
    here = os.path.dirname(os.path.abspath(__file__))
    for _ in range(5):
        f = os.path.join(here, 'config.conf')
        if os.path.exists(f):
            for ln in open(f, encoding='utf-8'):
                ln = ln.split('#')[0].strip()
                if '=' in ln:
                    k, v = ln.split('=', 1)
                    _CACHE[k.strip()] = os.path.expandvars(v.strip().strip('"\''))
            return _CACHE
        here = os.path.dirname(here)
    if not os.environ.get('METH_CONF_OPTIONAL'):
        raise SystemExit(
            'config.conf 를 못 찾았다 (meth_config.py 위 5단계까지 봄).\n'
            '  기본값으로 조용히 도는 것을 막으려고 여기서 멈춘다.\n'
            '  일부러 기본값만 쓰려면 METH_CONF_OPTIONAL=1 을 준다.')
    return _CACHE


def get(key):
    """환경변수 -> conf -> 기본값. 출처를 출처[key] 에 적어 둔다.

    셸 헤더가 `set -a; . config.conf` 로 conf 를 환경에 export 하므로
      conf 값도 os.environ 에 올라온다. 그것을 「환경변수」 로 적으면 스냅샷의
      출처 열이 「항상」 거짓이 된다. 기록() 의 목적이 무너진다.
      그래서 두 값이 같으면 출처는 conf 로 적는다.
      한계: 사용자가 conf 와 「똑같은」 값을 환경변수로 줬다면 conf 로 적힌다.
      구별할 방법이 없고, 값이 같으므로 해석에 영향이 없다.
    """
    _e = os.environ.get(key)
    _c = _conf().get(key)
    if _e is not None and _e != '':
        v = _e
        origin = 'conf' if (_c is not None and _e == _c) else '환경변수'
    elif _c is not None and _c != '':
        v, origin = _c, 'conf'
    else:
        v, origin = _DEFAULTS.get(key, ''), '기본값'
    SOURCE[key] = origin
    return v


def get_path(key):
    return _expand(get(key))


# conf 안의 $METH_ROOT·$DATASET_ROOT 를 여기서 푼다.
#   셸은 `set -a; . config.conf` 로 풀지만, 파이썬이 conf 를 직접 읽을 때는
#   그 둘이 「환경변수가 아니라서」 expandvars 가 못 푼다. 문자 그대로 남는다.
#   뿌리 둘은 자기 자신을 참조할 수 없으므로 먼저 푼 값을 담아 두고 쓴다.
_ROOTS = {}


def _expand(v):
    for _k in ('METH_ROOT', 'DATASET_ROOT'):
        _val = _ROOTS.get(_k)
        if _val:
            v = v.replace('${%s}' % _k, _val).replace('$' + _k, _val)
    return os.path.expanduser(os.path.expandvars(v)).rstrip('/')


METH_ROOT = _ROOTS['METH_ROOT'] = get_path('METH_ROOT')
DATASET_ROOT = _ROOTS['DATASET_ROOT'] = get_path('DATASET_ROOT')
WORK_ROOT = get_path('WORK_ROOT')
TMP_ROOT = get_path('TMP_ROOT')


TRAIN_SRC = get_path('TRAIN_SRC')
ANALYSIS_DIR = get_path('ANALYSIS_DIR')

GEN = get('GEN')
CANCER_LABEL = get('CANCER_LABEL')
VAL_COHORTS = [x.strip() for x in get('VAL_COHORTS').split(',') if x.strip()]
PANELS = [x.strip() for x in get('PANELS').split(',') if x.strip()]
# STEP_MAP 도 다른 키처럼 모듈 속성으로 노출한다. 다루는 법을 하나로 맞춘다.
STEP_MAP = get('STEP_MAP')
NULL_PANELS = [x.strip() for x in get('NULL_PANELS').split(',') if x.strip()]

def version_name(panel, depth):
    """접두(세대) + 명단(자료) + 깊이 를 한 곳에서 조립한다."""
    return 'j%s%s_%s_uni' % (GEN, panel, depth)
VERSION = get('VERSION')
PANEL = get('PANEL')
DEPTH_MODE = get('DEPTH_MODE')

TARGET_DEPTH_READS = int(get('TARGET_DEPTH_READS'))
PATTERN_K = int(get('PATTERN_K'))
MIN_READS_PER_WINDOW = int(get('MIN_READS_PER_WINDOW'))
BLOCK_SIZE = int(get('BLOCK_SIZE'))
SEED = int(get('SEED'))
N_JOBS = int(get('N_JOBS'))
MIN_SAMPLE_COVERAGE = float(get('MIN_SAMPLE_COVERAGE'))
COLSET_FROM = get('COLSET_FROM').strip().lower()
BLAS_THREADS = get('BLAS_THREADS').strip().lower()
if BLAS_THREADS != 'none' and not (BLAS_THREADS.isdigit() and int(BLAS_THREADS) >= 1):
    raise SystemExit('BLAS_THREADS 는 none 또는 1 이상의 정수여야 합니다: %r' % BLAS_THREADS)
SCORE_UNIT = get('SCORE_UNIT').strip().lower()
if SCORE_UNIT not in ('as_is', 'consistent'):
    raise SystemExit('SCORE_UNIT 은 as_is 또는 consistent 이어야 합니다: %r' % SCORE_UNIT)
REPORT_BASELINE = get('REPORT_BASELINE').strip().lower()
if REPORT_BASELINE not in ('on', 'off'):
    raise SystemExit('REPORT_BASELINE 은 on 또는 off 이어야 합니다: %r' % REPORT_BASELINE)
CONVERGENCE_GATE = get('CONVERGENCE_GATE').strip().lower()
if CONVERGENCE_GATE not in ('on', 'off'):
    raise SystemExit('CONVERGENCE_GATE 는 on 또는 off 이어야 합니다: %r' % CONVERGENCE_GATE)
if COLSET_FROM not in ('all', 'normal'):
    raise SystemExit('COLSET_FROM 은 all 또는 normal 이어야 합니다: %r' % COLSET_FROM)
TEST_SIZE = float(get('TEST_SIZE'))
MUT_RATIOS = [float(x) for x in get('MUT_RATIOS').replace(' ', '').split(',') if x]

# 모델을 미리 하나로 고정한다. 이름은 08·09 의 모델 키와 같아야 한다.
#   결과를 보고 고르면 그 자체가 선택이라 귀무 바닥이 0.03~0.04 오른다(사전등록 41-2).
FIXED_MODEL = get('FIXED_MODEL').strip()
# LLR 참조표 공여자와 섞기 공여자를 가른다(사전등록 26-2). 1 이면 켠다.
LLR_HELDOUT = int(get('LLR_HELDOUT') or 0)


def split_ref_mix(units):
    """(참조군, 혼합군). 단위 이름을 바이트 순으로 줄 세워 번갈아 가른다 — 사전등록 26-2.

    바이트 순 = `LC_ALL=C sort`. 기본 로케일 정렬은 대소문자를 무시해 다른 명단을 낸다.
    1·3·5…번째가 참조군, 2·4·6…번째가 혼합군. 시드를 쓰지 않는다 — 고르는 행위가 선택이 된다.
    단위는 「사람·세포주」 다. 같은 세포주의 반복을 쪼개면 누수가 남는다.
    """
    u = sorted(set(units), key=lambda s: s.encode('utf-8'))
    return u[0::2], u[1::2]

# VERSION 이름 안의 깊이와 TARGET_DEPTH_READS 가 어긋나면 멈춘다.
#   둘은 conf 에 따로 적히므로 어긋날 수 있다. 적어 두는 것은 방어가 아니다.
_depth_match = re.search(r'_(\d+)k_', VERSION)
if _depth_match is None:
    # 못 찾은 것을 「통과」 로 두면 이름이 조금만 달라져도 가드가 조용히 사라진다.
    raise SystemExit(
        'VERSION(%s) 에서 깊이(_Nk_)를 못 찾았다.\n'
        '  j15bl200_5k_uni 처럼 이름 안에 깊이가 들어가야 한다.' % VERSION)
if int(_depth_match.group(1)) * 1000 != TARGET_DEPTH_READS:
    raise SystemExit(
        'VERSION(%s)의 깊이와 TARGET_DEPTH_READS(%d)가 어긋난다.\n'
        '  둘을 같이 고친다. 한쪽만 고치면 폴더 이름과 내용이 달라진다.'
        % (VERSION, TARGET_DEPTH_READS))


def snapshot(folder):
    """해결된 설정과 각 값의 출처를 산출 폴더에 「덧붙인다」.

    산출물 옆에 이것이 없으면 「어느 설정에서 나온 숫자인가」 를 나중에 못 댄다.
    환경변수로 한 번 돌린 실행은 특히 아무 데도 안 남는다.

    덮어쓰지 않고 이어쓰는 이유 —
      폴더 이름을 정하는 것은 ROOT 와 VERSION 뿐이다. MIN_READS_PER_WINDOW ·
      PATTERN_K · SEED · MUT_RATIOS · DEPTH_MODE 는 **경로를 안 바꾸면서
      산출물을 바꾼다.** 그래서 한 폴더에 서로 다른 설정의 산출물이 섞일 수 있다.
      MIN_READS_PER_WINDOW 를 4에서 6으로 바꿨을 때 실제로 섞였다.
      덮어쓰면 마지막 실행의 설정 하나만 남아 「이 폴더는 하나다」 라고
      단언한다. 없는 것보다 나쁘다. 이어쓰면 블록이 갈려 섞인 것이 보인다.
      폴더가 정말 한 설정이면 같은 블록이 반복될 뿐이라 잃는 것이 없다.
      「현재」 를 알고 싶으면 마지막 블록을 본다.
    """
    os.makedirs(folder, exist_ok=True)
    p = os.path.join(folder, 'config_snapshot.txt')
    who = os.path.basename(sys.argv[0]) if sys.argv and sys.argv[0] else '?'
    with open(p, 'a', encoding='utf-8') as f:
        f.write('\n# %s  %s\n'
                % (datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S'), who))
        for k in sorted(SOURCE):
            f.write('%-22s %-46s (%s)\n' % (k, globals().get(k, ''), SOURCE[k]))
        # 수치 연산이 실제로 본 스레드 값. 환경변수로 미리 준 값이 BLAS_THREADS 보다 이기므로 따로 적는다.
        for k in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS'):
            f.write('%-22s %-46s (%s)\n' % (k, os.environ.get(k, '없음(라이브러리 기본)'), 'env'))
    return p


if __name__ == '__main__':
    for k in ('METH_ROOT', 'DATASET_ROOT', 'WORK_ROOT', 'TMP_ROOT',
              'TRAIN_SRC', 'ANALYSIS_DIR', 'GEN', 'VAL_COHORTS', 'PANELS', 'NULL_PANELS', 'VERSION', 'PANEL',
              'TARGET_DEPTH_READS', 'DEPTH_MODE', 'MUT_RATIOS', 'PATTERN_K',
              'MIN_READS_PER_WINDOW', 'BLOCK_SIZE', 'MIN_SAMPLE_COVERAGE', 'COLSET_FROM', 'BLAS_THREADS', 'SCORE_UNIT', 'REPORT_BASELINE', 'CONVERGENCE_GATE',
              'SEED', 'TEST_SIZE', 'N_JOBS', 'FIXED_MODEL', 'LLR_HELDOUT'):
        print('%-22s %-44s (%s)' % (k, globals()[k], SOURCE.get(k, '-')))
    # 실제로 수치 연산이 본 스레드 값(진입점이 내보낸 값 · 없으면 라이브러리 기본)
    for k in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS'):
        print('%-22s %-44s (%s)' % (k, os.environ.get(k, '없음(라이브러리 기본)'), 'env'))


def count_convergence():
    """수렴 경고(ConvergenceWarning)를 세는 칸 [n] 을 돌려준다. 2단계 08 · 09 가 부른다.

    경고마다 세되 화면에는 처음 하나만 낸다. 다른 경고는 그대로 지나간다.
    'always' 는 판에 따라 같은 자리 경고를 한 번만 내는 것을 막는다(sklearn 1.6 은 적합마다 다시 낸다).
    이 프로세스에서 난 것만 센다 — joblib 일꾼(n_jobs > 1)에서 난 경고는 못 센다.
    """
    import warnings
    from sklearn.exceptions import ConvergenceWarning
    n, show = [0], warnings.showwarning

    def _show(message, category, *a, **k):
        if issubclass(category, ConvergenceWarning):
            n[0] += 1
            if n[0] > 1:
                return
        show(message, category, *a, **k)

    warnings.simplefilter('always', ConvergenceWarning)
    warnings.showwarning = _show
    return n


def convergence_check(n, failures, step):
    """수렴 경고 건수를 찍고, CONVERGENCE_GATE=on 이면 1건 이상일 때 failures 에 넣는다."""
    print('\n수렴 경고 %d건 (%s · CONVERGENCE_GATE=%s)' % (n[0], step, CONVERGENCE_GATE))
    if n[0] and CONVERGENCE_GATE == 'on':
        failures.append('수렴 경고 %d건 - CONVERGENCE_GATE=on' % n[0])


def unify_score_unit(D, keys=('feature_set',)):
    """SCORE_UNIT=consistent 이면 keys 묶음마다 점수 단위를 하나로 맞춘다.

    채점(score_one_sample.py)은 검체마다 점수 폭이 1e-6 미만으로 뭉치면 그 검체만
    decision_function 으로 바꾼다. 그러면 한 칸 안에 확률과 decision_function 이 섞여
    AUC 를 해석할 수 없다. consistent 에서는 채점이 y_df 열을 늘 같이 적고, 여기서
    decision_function 이 있는 모형(묶음 안 y_df 가 다 찬 묶음)은 묶음 전체를 y_df 로 바꾼다.
    전환이 없어도 바꾼다 — 폭은 넓지만 일부 행이 1.0 에 붙은 «부분 포화» 동률도 없앤다(사전등록 75).
    LR 은 시그모이드라 순서는 같고 동률만 풀린다. RF · XGB 처럼 decision_function 이 없으면 확률 그대로다.
    as_is 면 D 를 그대로 돌려준다(지금까지의 표를 낸 동작).
    돌려주는 값: (D, 단위를 바꾼 묶음 수)
    """
    if SCORE_UNIT != 'consistent' or not len(D):
        return D, 0
    if 'y_df' not in D.columns:
        raise SystemExit('SCORE_UNIT=consistent 인데 채점 파일에 y_df 열이 없습니다 — 같은 설정으로 다시 채점하십시오')
    K = [D[k] for k in keys]
    st = D['score_type'].astype(str)
    sw = st.str.contains('df')                                  # 포화로 전환된 행
    if D.loc[sw.groupby(K).transform('any'), 'y_df'].isna().any():
        raise SystemExit('전환된 행이 있는 묶음에 y_df 가 빈 행이 있습니다 — '
                         'consistent 로 다시 채점되지 않은 검체 파일이 섞였습니다')
    have = D['y_df'].notna().groupby(K)
    g = have.transform('all')                                     # decision_function 이 있는 모형의 묶음
    if (have.transform('any') & ~g).any():                        # 일부만 빔 — 옛 설정 파일이 섞였다
        raise SystemExit('한 묶음 안에서 y_df 가 일부 행만 비었습니다 — consistent 로 다시 채점되지 않은 '
                         '검체 파일이 섞였습니다')
    conv = g & st.str.startswith('predict_proba')               # 단위가 실제로 바뀌는 행 (SVM 은 이미 df)
    D = D.copy()
    D.loc[g, 'y_score'] = D.loc[g, 'y_df']
    # 확률 단위였던 행의 문턱(thr_spec95)은 이제 y_score 와 단위가 달라 비운다. y_pred 는 원래 단위로
    #   이미 판정된 값이라 그대로 둔다. 포화로 이미 decision_function 이던 행은 문턱도 z 공간이라 둔다.
    if 'thr_spec95' in D.columns:
        D.loc[conv & ~sw, 'thr_spec95'] = float('nan')
    D.loc[conv, 'score_type'] = 'decision_function(통일)'
    return D, int(D.loc[conv, list(keys)].drop_duplicates().shape[0])
