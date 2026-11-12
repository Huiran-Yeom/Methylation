# -*- coding: utf-8 -*-
"""config.conf 읽기 (고칠 곳은 옆의 config.conf 하나)
우선순위: 환경변수 → config.conf → 기본값

    import meth_config as cfg
    cfg.METH_ROOT            '/ssd_data/Methylation'
    cfg.TARGET_DEPTH_READS   5000
    cfg.MUT_RATIOS           [0.0, 0.001, ...]
    cfg.snapshot('<output folder>')     최종 값과 출처를 산출 폴더에 남김
"""
import datetime
import os
import re
import sys

# 기본 경로는 코드 위치 기준 (저장소 받은 그대로 동작)
_TREE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

_DEFAULTS = {
    'METH_ROOT': '/ssd_data/Methylation',
    'DATASET_ROOT': '/ssd_data/dataset',
    'WORK_ROOT': os.path.expanduser('~'),
    'TMP_ROOT': os.path.expanduser('~/tmp'),
    'TRAIN_SRC': os.path.join(_TREE, '2_train', 'steps'),
    'ANALYSIS_DIR': os.path.join(_TREE, '3_validate', 'steps'),
    'GEN': '15',
    # samples.csv type 열의 암 쪽 이름. 바꾸면 1단계부터 다시 돌림
    'CANCER_LABEL': 'GBM',
    # 비우면 아무 코호트나 허용, 적으면 그 목록만
    'VAL_COHORTS': '',
    'PANELS': 'bl200,jsd200,jsdb200',
    'NULL_PANELS': 'rand1200,rand2200,rand3200,randb1200,randb2200,randb3200',
    'VERSION': 'j15bl200_5k_uni',
    'PANEL': 'bl200',
    # 1단계 패널 블록 수. PANELS, NULL_PANELS, VERSION, PANEL 이름의 숫자도 같이 바꿈
    'PANEL_N': '200',
    # bl 안정성 빈도를 셀 때 중요도 상위 몇 개까지 볼지. N 이면 PANEL_N 과 같게
    'TOPCAP': '1000',
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
    # 비우면 모델 4개를 교차검증으로 겨뤄 1등 사용
    'FIXED_MODEL': '',
    # 0 이면 LLR 참조표와 학습 풀에 같은 사람이 들어감
    'LLR_HELDOUT': '0',
}
_CACHE = None
SOURCE = {}


def _conf():
    """config.conf 찾아 읽기. 못 찾으면 중단 (조용히 기본값으로 도는 것 방지)
    기본값만 쓰려면 METH_CONF_OPTIONAL=1
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
    """환경변수 → conf → 기본값. 출처는 SOURCE[key] 에 기록
    셸이 conf 를 export 해서 환경변수에도 같은 값이 올라옴 → 값이 conf 와 같으면 출처는 conf
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


# conf 안의 $METH_ROOT, $DATASET_ROOT 를 직접 풂
# 파이썬이 conf 를 읽을 땐 둘이 환경변수가 아니라 expandvars 로 안 풀림
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
# STEP_MAP 도 다른 키처럼 모듈 속성으로 노출
STEP_MAP = get('STEP_MAP')
NULL_PANELS = [x.strip() for x in get('NULL_PANELS').split(',') if x.strip()]
PANEL_N = get('PANEL_N').strip()
if not (PANEL_N.isdigit() and int(PANEL_N) >= 1):
    raise SystemExit('PANEL_N 은 1 이상의 정수여야 합니다: %r' % PANEL_N)
PANEL_N = int(PANEL_N)
TOPCAP = get('TOPCAP').strip()
if TOPCAP != 'N' and not (TOPCAP.isdigit() and int(TOPCAP) >= 1):
    raise SystemExit('TOPCAP 은 N 또는 1 이상의 정수여야 합니다: %r' % TOPCAP)
TOPCAP = PANEL_N if TOPCAP == 'N' else int(TOPCAP)


def cancer_env(new, old):
    """암 쪽 자료 위치. 새 이름(CANCER_*)과 옛 이름(GBM_*) 둘 다 받음. 둘 다 있고 다르면 중단"""
    a, b = os.environ.get(new, ''), os.environ.get(old, '')
    if a and b and a != b:
        raise SystemExit('%s 와 %s 가 다르다 (%s · %s). 하나만 주십시오' % (new, old, a, b))
    return a or b


def version_name(panel, depth):
    """판 이름 조립: 세대 + 패널 + 깊이"""
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

# 모델 고정. 이름은 2단계 08, 09 의 모델 키와 같아야 함
FIXED_MODEL = get('FIXED_MODEL').strip()
# 1 이면 LLR 참조표 공여자와 섞기 공여자를 가름
LLR_HELDOUT = int(get('LLR_HELDOUT') or 0)


def split_ref_mix(units):
    """(참조군, 혼합군) 반환. 단위 이름을 바이트 순(LC_ALL=C sort)으로 세워 번갈아 가름
    홀수 번째가 참조군, 짝수 번째가 혼합군. 시드 안 씀
    단위는 사람 또는 세포주 (같은 세포주 반복을 쪼개면 누수)
    """
    u = sorted(set(units), key=lambda s: s.encode('utf-8'))
    return u[0::2], u[1::2]

# VERSION 이름 안의 깊이와 TARGET_DEPTH_READS 가 어긋나면 중단
_depth_match = re.search(r'_(\d+)k_', VERSION)
if _depth_match is None:
    # 깊이를 못 찾아도 중단 (그냥 넘기면 가드가 꺼짐)
    raise SystemExit(
        'VERSION(%s) 에서 깊이(_Nk_)를 못 찾았다.\n'
        '  j15bl200_5k_uni 처럼 이름 안에 깊이가 들어가야 한다.' % VERSION)
if int(_depth_match.group(1)) * 1000 != TARGET_DEPTH_READS:
    raise SystemExit(
        'VERSION(%s)의 깊이와 TARGET_DEPTH_READS(%d)가 어긋난다.\n'
        '  둘을 같이 고친다. 한쪽만 고치면 폴더 이름과 내용이 달라진다.'
        % (VERSION, TARGET_DEPTH_READS))
# VERSION 의 세대, 판도 GEN, PANEL 과 맞아야 함 (어긋나면 다른 판 모델이 이 판 경로에 저장됨)
# VERSION, PANEL 을 환경변수로 준 실행(2, 3단계)만 검사. 1단계는 VERSION 안 씀
if '환경변수' in (SOURCE.get('VERSION'), SOURCE.get('PANEL')) and VERSION.startswith('j') and not VERSION.startswith('j%s%s_' % (GEN, PANEL)):
    raise SystemExit(
        'VERSION(%s) 이 GEN(%s) · PANEL(%s) 과 어긋난다. j<GEN><PANEL>_<깊이>k_<방식> 꼴이어야 한다.'
        % (VERSION, GEN, PANEL))


def snapshot(folder):
    """최종 설정과 출처를 산출 폴더의 config_snapshot.txt 에 이어씀
    덮어쓰지 않음 - 경로 안 바꾸는 설정(PATTERN_K, SEED 등)이 바뀌면 한 폴더에 섞여서
    지금 설정은 마지막 블록
    """
    os.makedirs(folder, exist_ok=True)
    p = os.path.join(folder, 'config_snapshot.txt')
    who = os.path.basename(sys.argv[0]) if sys.argv and sys.argv[0] else '?'
    with open(p, 'a', encoding='utf-8') as f:
        f.write('\n# %s  %s\n'
                % (datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S'), who))
        for k in sorted(SOURCE):
            f.write('%-22s %-46s (%s)\n' % (k, globals().get(k, ''), SOURCE[k]))
        # 실제 스레드 값도 기록 (환경변수가 BLAS_THREADS 보다 우선)
        for k in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS'):
            f.write('%-22s %-46s (%s)\n' % (k, os.environ.get(k, '없음(라이브러리 기본)'), 'env'))
    return p


if __name__ == '__main__':
    for k in ('METH_ROOT', 'DATASET_ROOT', 'WORK_ROOT', 'TMP_ROOT',
              'TRAIN_SRC', 'ANALYSIS_DIR', 'GEN', 'VAL_COHORTS', 'PANELS', 'NULL_PANELS', 'PANEL_N', 'TOPCAP', 'VERSION', 'PANEL',
              'TARGET_DEPTH_READS', 'DEPTH_MODE', 'MUT_RATIOS', 'PATTERN_K',
              'MIN_READS_PER_WINDOW', 'BLOCK_SIZE', 'MIN_SAMPLE_COVERAGE', 'COLSET_FROM', 'BLAS_THREADS', 'SCORE_UNIT', 'REPORT_BASELINE', 'CONVERGENCE_GATE',
              'SEED', 'TEST_SIZE', 'N_JOBS', 'FIXED_MODEL', 'LLR_HELDOUT'):
        print('%-22s %-44s (%s)' % (k, globals()[k], SOURCE.get(k, '-')))
    # 실제 스레드 값 (진입점이 내보낸 값, 없으면 라이브러리 기본)
    for k in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS'):
        print('%-22s %-44s (%s)' % (k, os.environ.get(k, '없음(라이브러리 기본)'), 'env'))


def check_same_refs(paths, what=''):
    """BAM 들의 참조 서열 목록(이름, 길이, 순서)이 같은지 확인, 다르면 중단
    리드를 다른 헤더로 옮기면 reference_id 가 그대로라 엉뚱한 염색체에 붙음
    """
    import pysam
    sig = None
    for p in paths:
        with pysam.AlignmentFile(p, 'rb', check_sq=False) as b:
            s = tuple(zip(b.references, b.lengths))
        if sig is None:
            sig, first = s, p
        elif s != sig:
            raise SystemExit('참조 서열 목록이 다르다%s: %s 와 %s. 같은 기준 유전체 · 같은 헤더로 정렬된 BAM 이어야 한다'
                             % (' (' + what + ')' if what else '', os.path.basename(first), os.path.basename(p)))


def fragments(reads, stats=None):
    """BAM 리드를 조각 단위로 묶어 내줌. 정렬 순서와 상관없이 PE 두 짝을 한 조각으로
    PE: read1, read2 다 올 때까지 이름으로 기다림 (끝까지 짝 없으면 마지막에 냄)
    짝 표시 없는 리드: 같은 이름이 바로 이어 오면 한 조각으로 합침
    정렬 안 됨, 2차, 보조 정렬은 건너뜀
    stats 주면 같은 이름, 같은 짝 번호 중복 수를 stats['dup'] 에 셈
    """
    pend, se = {}, []
    for r in reads:
        if r.is_unmapped or r.is_secondary or r.is_supplementary:
            continue
        if not r.is_paired or r.mate_is_unmapped:
            if se and se[0].query_name == r.query_name:
                se.append(r)
                continue
            if se:
                yield se
            se = [r]
            continue
        if se:
            yield se
            se = []
        k = 1 if r.is_read1 else 2
        e = pend.get(r.query_name)
        if e is None:
            e = pend[r.query_name] = ([], set())
        elif k in e[1] and stats is not None:
            stats['dup'] = stats.get('dup', 0) + 1
        e[0].append(r)
        e[1].add(k)
        if len(e[1]) == 2:
            del pend[r.query_name]
            yield e[0]
    if se:
        yield se
    for e in pend.values():
        yield e[0]


def dup_check(n, allow, what):
    """fragments 가 센 중복 리드 수 출력. allow 거짓이면 0 아닐 때 중단
    파일 이름은 안 찍음 (정상 BAM 이름에 사람 이름 있음)"""
    if not n:
        return
    m = '%s: 같은 이름 · 같은 짝 번호 리드 %d 개 (복원 추출 또는 이름 충돌)' % (what, n)
    if not allow:
        raise SystemExit(m + '. ALLOW_REPLACEMENT=False 인데 0 이 아니므로 멈춘다')
    print('  경고: ' + m)


def count_convergence():
    """수렴 경고(ConvergenceWarning) 세는 칸 [n] 반환. 2단계 08, 09 에서 사용
    다 세고 화면엔 첫 1건만 출력. 이 프로세스 것만 셈 (joblib 일꾼 경고는 못 셈)
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
    """수렴 경고 건수 출력. CONVERGENCE_GATE=on 이고 1건 이상이면 failures 에 추가"""
    print('\n수렴 경고 %d건 (%s · CONVERGENCE_GATE=%s)' % (n[0], step, CONVERGENCE_GATE))
    if n[0] and CONVERGENCE_GATE == 'on':
        failures.append('수렴 경고 %d건 - CONVERGENCE_GATE=on' % n[0])


def unify_score_unit(D, keys=('feature_set',)):
    """SCORE_UNIT=consistent 이면 keys 묶음마다 점수 단위를 하나로 맞춤
    채점은 점수가 뭉친 검체만 decision_function 으로 바꿔서 한 칸에 단위가 섞임
    y_df 가 다 찬 묶음은 전체를 y_df 로 바꿈 (RF, XGB 는 확률 그대로). as_is 면 그대로
    반환: (D, 단위 바꾼 묶음 수)
    """
    if SCORE_UNIT != 'consistent' or not len(D):
        return D, 0
    if 'y_df' not in D.columns:
        raise SystemExit('SCORE_UNIT=consistent 인데 채점 파일에 y_df 열이 없습니다. 같은 설정으로 다시 채점하십시오')
    K = [D[k] for k in keys]
    st = D['score_type'].astype(str)
    sw = st.str.contains('df')                                  # 포화로 전환된 행
    if D.loc[sw.groupby(K).transform('any'), 'y_df'].isna().any():
        raise SystemExit('전환된 행이 있는 묶음에 y_df 가 빈 행이 있습니다. '
                         'consistent 로 다시 채점되지 않은 검체 파일이 섞였습니다')
    have = D['y_df'].notna().groupby(K)
    g = have.transform('all')                                     # decision_function 이 있는 모형의 묶음
    if (have.transform('any') & ~g).any():                        # 일부만 빔 (다른 설정 채점 파일 섞임)
        raise SystemExit('한 묶음 안에서 y_df 가 일부 행만 비었습니다. consistent 로 다시 채점되지 않은 '
                         '검체 파일이 섞였습니다')
    conv = g & st.str.startswith('predict_proba')               # 단위가 실제로 바뀌는 행 (SVM 은 이미 df)
    D = D.copy()
    D.loc[g, 'y_score'] = D.loc[g, 'y_df']
    # 확률 단위였던 행의 문턱(thr_spec95)은 단위가 달라져서 비움. y_pred 는 그대로
    # 이미 decision_function 이던 행은 문턱도 같은 단위라 그대로
    if 'thr_spec95' in D.columns:
        D.loc[conv & ~sw, 'thr_spec95'] = float('nan')
    D.loc[conv, 'score_type'] = 'decision_function(통일)'
    return D, int(D.loc[conv, list(keys)].drop_duplicates().shape[0])
