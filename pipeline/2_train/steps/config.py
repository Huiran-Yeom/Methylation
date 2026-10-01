# -*- coding: utf-8 -*-
"""j 파이프라인 설정: _make_versions.py 가 이 틀에서 생성한다. 직접 고치지 말 것.

이름 규칙  j{패널}_{깊이}_{믹싱}
  패널  bl200 | bl300 | jsd200 | jsd300
        bl  = 기존 방식 (q<0.05 & |dbeta|>=30 -> 부트스트랩 순위)
        jsd = JSD 순위 상위 N
  깊이  5k | 50k | 100k     TARGET_DEPTH_READS = 샘플 하나의 총 리드 수
  믹싱  uni | cov           균일 무작위 | 커버리지 일치

m 판과 다른 곳 (2026-08-20)
  - 정상 풀이 EGAD 44파일(3명) -> 공개·랩 20명
  - 풀 BAM 을 패널 영역만 뽑아 만든다 (39G -> 28.8MB)
  - DMR 은 버전이 아니라 패널 이름으로 찾는다
"""
import os
import os as _os
import sys as _sys

# 설정은 0_setup/config.conf 에서만 읽는다. 다른 설정 기전을 두지 않는다.
# 검증은 이 파일을 작업폴더로 「복사」 해 돌린다 — 트리 밖이라
#   상위 폴더를 거슬러 올라가도 0_setup 이 없다. 그때는 METH_CONF_DIR 를 본다.
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
import meth_config as _cfg  # noqa: E402
_M, _D = _cfg.METH_ROOT, _cfg.DATASET_ROOT


ROOT = _M
GEN = _cfg.GEN        # 패널 세대 접미사. 06·07 의 세대 관문이 이걸 본다
VERSION = _cfg.VERSION
PANEL   = _cfg.PANEL

DMR_SOURCE_VERSION = None
SHARE_FROM = None          # 이 판이 2_train/steps/01 을 직접 돌린다
SEED = _cfg.SEED

BLOCK_SIZE = _cfg.BLOCK_SIZE
CPG_THRESHOLD = 4             # step02 전용 — j 는 step02 를 여기서 안 돌린다
CPG_MIN_COVERAGE = 4
DMR_CPG_SAMPLE_FRAC = 1.0
MIN_READS_PER_BLOCK = 5       # 2_train/steps/01 에서 실제로 쓰인다

DEPTH_MODE = _cfg.DEPTH_MODE
TARGET_DEPTH_READS = _cfg.TARGET_DEPTH_READS
ALLOW_REPLACEMENT = False

MUT_RATIOS = _cfg.MUT_RATIOS
NUM_FILES = 100
POOL_HEADROOM = 2             # cov 전용

SITE_ENTROPY_METHOD = 'histogram'   # 08-20 변경. value_counts 는 커버리지 대리지표
INPUT_FORMAT = 'bedgraph'
FILE_ORDER = 'numeric'
MIN_SAMPLE_COVERAGE = _cfg.MIN_SAMPLE_COVERAGE
TEST_SIZE = _cfg.TEST_SIZE
STRATIFY = True
SVM_PROBABILITY = False
SAVE_TRAIN_SCORES = True
N_JOBS = _cfg.N_JOBS

MERGE_STRANDS = True
PATTERN_K = _cfg.PATTERN_K                 # m200 학습·검증 모두 3 (로그 확인). m736 만 4였다
MIN_READS_PER_WINDOW = _cfg.MIN_READS_PER_WINDOW      # 2026-09-14 사용자 지시 4->6. DMR선정(04_jsd)·LLR참조(05_llrref)가
#                               둘 다 6이라 특징만 4로 느슨했다. 네 단계 문턱을 맞춘다.
#                               창당 4리드면 엔트로피가 가질 수 있는 값이 5개뿐이었다(최대 0.6667).
#                               6이면 11개. 믹싱 BAM 창당 리드 중앙 97이라 창은 거의 안 준다.
NORMALIZE_PATTERN = True

# ╔══ 이 단계가 쓰는 경로 ═══════════════════════════════════════════╗
#   배치가 다르면 「여기만」 고친다. config.conf 에는 뿌리 둘과 실험 설정만 둔다.
GENOME_FOLDER = ROOT + '/data/Bisulfite_Genome/'   # bismark 인덱스 (없으면 --genome_folder 없이 돈다)
# ╚══════════════════════════════════════════════════════════════════╝

# RAW_NORMAL_COV_DIR · RAW_GBM_COV_DIR 을 지웠다 — 트리의 어느 코드도
#   읽지 않는 죽은 값이었다. 정상·암 cov 자리는 1단계가 자기 경로 블록에서 정하고,
#   환경변수 NORMAL_SET·GBM_COV_DIR 로 옮긴다. 여기에 또 두면 둘이 어긋난다.
COV_SUFFIX = '.cov.gz'
TOP_N_SAMPLES = 0
DMR_PER_GROUP = 0
DMR_GROUP_RE = r'^([A-Za-z]+[0-9]*)[_-]'

# 판마다 「자기 풀」 을 본다. 예전에는 _cellline 하나로 박혀 있었다.
#   풀 BAM 은 패널 영역의 리드만 담는다. 귀무판은 무작위 블록이라 실판 합집합
#   영역 밖이고, cellline 풀에는 그 자리의 리드가 없다. 실측:
#     rand1200  통과 블록 8 / 200
#     rand3200  「풀 부족으로 중단」 (정상 2,678 보유 / 2,936 요청)
#   원본은 판마다 config.py 가 따로 있어 각자 자기 풀을 가리켰다:
#     j15rand*  -> j_pool<GEN>_rand    (무작위 패널 전용 풀)
#     j15randb* -> j_pool<GEN>_randb   (귀무B 전용 · bl200 CpG 분포 맞춤)
#     그 외     -> j_pool<GEN>_cellline (정상15 · 암30 · 세포주 단위)
#   1_dmr/steps/run_pool.sh 가 cellline|rand|randb 셋을 다 만든다.
#   randb 를 rand 보다 「먼저」 본다. 'randb1200' 은 'rand' 로도 시작한다.
if PANEL.startswith('randb'):  _POOLDIR = 'j_pool' + GEN + '_randb'
elif PANEL.startswith('rand'): _POOLDIR = 'j_pool' + GEN + '_rand'
else:                          _POOLDIR = 'j_pool' + GEN + '_cellline'
_POOL = os.environ.get('POOL_DIR') or (ROOT + '/results/dmr/' + _POOLDIR)
# 2026-09 까지 실판 풀 폴더 이름이 j_pool<GEN>_세포주 였다. 옛 트리에는 영문 이름이
#   없어 풀을 못 찾고 2단계 첫 걸음이 죽는다. 새 이름이 없으면 옛 이름을 본다.
if not os.path.isdir(_POOL) and _POOLDIR.endswith('_cellline'):
    _old = _POOL[:-len('_cellline')] + '_세포주'
    if os.path.isdir(_old):
        _POOL = _old
FULL_NORMAL_BAM = _POOL + '/pool_normal_all.bam'
FULL_GBM_BAM    = _POOL + '/pool_gbm_all.bam'
# NEGCTRL_BG_BAM · NEGCTRL_MIX_BAM 을 지웠다 — 읽는 코드가 없고,
#   가리키던 pool_normal_A/B.bam 을 만드는 코드도 없다(11_pool 은 *_all.bam 만 쓴다).

_R = ROOT.rstrip('/\\')
_P = VERSION.split('_')[0]          # 패널 폴더 (jbl200 · jjsd300 …)
SD = _R + '/sample_data/panels/' + _P + '/' + VERSION
RS = _R + '/results/panels/' + _P + '/' + VERSION
_S = SHARE_FROM or VERSION

DMR_DIR = _R + '/results/dmr/j_panel_dmr' + GEN + '_panels/' + PANEL
DMR_CSV = DMR_DIR + '/DMR_confirmed_' + PANEL + '.csv'
MIX_SD  = _R + '/sample_data/panels/' + _S.split('_')[0] + '/' + _S + '/step03_Preprocessing_for_ML'
MUT_SORTED = MIX_SD + '/mut_sorted_' + _S + '.bam'
WT_SORTED  = MIX_SD + '/wt_sorted_' + _S + '.bam'
MUT_NAMES  = MIX_SD + '/mut_names_' + _S + '.txt'
WT_NAMES   = MIX_SD + '/wt_names_' + _S + '.txt'

MIX_OUT      = RS + '/step03_Preprocessing_for_ML'
BISMARK_OUT  = SD + '/step04_ML_classifier'
ML_OUT       = RS + '/step04_ML_classifier'
READ_ENT_OUT = RS + '/step04_read_entropy'


def folder_name(ratio):
    return 'mut_%d_reads' % int(ratio * TARGET_DEPTH_READS)


MUT_READ_COUNTS = [int(r * TARGET_DEPTH_READS) for r in MUT_RATIOS]


def ensure_dirs():
    for d in (MIX_SD, MIX_OUT, BISMARK_OUT, ML_OUT):
        os.makedirs(d, exist_ok=True)
    # 산출 폴더를 만드는 자리에서 설정을 같이 남긴다.
    #   이것이 없으면 「어느 설정에서 나온 숫자인가」 를 나중에 못 댄다.
    #   import 시점이 아니라 여기서 부른다. import 부작용은 만들지 않는다.
    _cfg.snapshot(ML_OUT)


def guard_not_original():
    if not VERSION.startswith('j'):
        raise SystemExit('중단: VERSION=' + VERSION + ' — j 로 시작해야 합니다.')
    for out, parent in ((SD, _R + '/sample_data/panels/' + _P), (RS, _R + '/results/panels/' + _P)):
        a, b = os.path.abspath(out), os.path.abspath(parent)
        if a == b:
            raise SystemExit('중단: 출력이 원본 폴더 자체입니다. ' + a)
        if os.path.basename(a) != VERSION:
            raise SystemExit('중단: 출력 경로 끝이 VERSION 이 아닙니다. ' + a)


def banner(step):
    print('=' * 74)
    print('  %s   [VERSION=%s · PANEL=%s]' % (step, VERSION, PANEL))
    print('  깊이 %s · 복제본 %d · 비율 %d개'
          % (format(TARGET_DEPTH_READS, ','), NUM_FILES, len(MUT_RATIOS)))
    print('  DMR  %s' % DMR_CSV)
    print('  중간 %s' % SD)
    print('  결과 %s' % RS)
    print('=' * 74)


FORCE = '--force' in _sys.argv

