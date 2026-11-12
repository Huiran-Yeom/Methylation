# -*- coding: utf-8 -*-
"""2단계 공통 설정. _make_versions.py 가 이 틀로 만듦, 직접 수정 금지

이름 규칙  j{패널}_{깊이}_{섞기}
  패널  bl200 | bl300 | jsd200 | jsd300
        bl  = q<0.05 & |dbeta|>=30 -> 부트스트랩 순위
        jsd = JSD 순위 상위 N
  깊이  5k | 50k | 100k     TARGET_DEPTH_READS = 샘플 하나의 총 리드 수
  섞기  uni | cov           균일 무작위 | 커버리지 맞춤

풀 BAM 은 패널 영역 리드만 담음, DMR 은 패널 이름으로 찾음
"""
import os
import os as _os
import sys as _sys

# 설정은 0_setup/config.conf 에서만 읽음
# 검증 때는 작업폴더로 복사돼 돌아서 0_setup 이 안 보임 → METH_CONF_DIR 사용
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
GEN = _cfg.GEN        # 패널 세대 접미사 (06, 07 세대 확인에 사용)
VERSION = _cfg.VERSION
PANEL   = _cfg.PANEL

DMR_SOURCE_VERSION = None
SHARE_FROM = None          # 01 을 이 판에서 직접 돌림
SEED = _cfg.SEED
FIXED_MODEL = getattr(_cfg, 'FIXED_MODEL', '')   # 비면 모델 4개 비교
count_convergence, convergence_check = _cfg.count_convergence, _cfg.convergence_check   # 08, 09 수렴 경고 세기

BLOCK_SIZE = _cfg.BLOCK_SIZE
CPG_THRESHOLD = 4             # 여기서는 안 씀
CPG_MIN_COVERAGE = 4
DMR_CPG_SAMPLE_FRAC = 1.0
MIN_READS_PER_BLOCK = 5       # 01 에서 사용

DEPTH_MODE = _cfg.DEPTH_MODE
TARGET_DEPTH_READS = _cfg.TARGET_DEPTH_READS
ALLOW_REPLACEMENT = False

MUT_RATIOS = _cfg.MUT_RATIOS
NUM_FILES = 100
POOL_HEADROOM = 2             # cov 전용

SITE_ENTROPY_METHOD = 'histogram'   # value_counts 는 커버리지를 따라가서 안 씀
INPUT_FORMAT = 'bedgraph'
FILE_ORDER = 'numeric'
MIN_SAMPLE_COVERAGE = _cfg.MIN_SAMPLE_COVERAGE
COLSET_FROM = _cfg.COLSET_FROM
TEST_SIZE = _cfg.TEST_SIZE
STRATIFY = True
SVM_PROBABILITY = False
SAVE_TRAIN_SCORES = True
N_JOBS = _cfg.N_JOBS

MERGE_STRANDS = True
PATTERN_K = _cfg.PATTERN_K                 # 창 하나의 CpG 수
MIN_READS_PER_WINDOW = _cfg.MIN_READS_PER_WINDOW      # 창당 최소 리드 (1단계 04_jsd, 05_llrref 와 같은 값)
NORMALIZE_PATTERN = True

# 이 단계 경로. 폴더 구조 다르면 여기만 수정
GENOME_FOLDER = ROOT + '/data/Bisulfite_Genome/'   # bismark 인덱스 (없으면 --genome_folder 없이 실행)

# 정상, 암 cov 경로는 1단계가 정함 (NORMAL_SET, GBM_COV_DIR 환경변수)
COV_SUFFIX = '.cov.gz'
TOP_N_SAMPLES = 0
DMR_PER_GROUP = 0
DMR_GROUP_RE = r'^([A-Za-z]+[0-9]*)[_-]'

# 판마다 자기 풀 사용 - 풀은 자기 패널 영역 리드만 있어서 무작위 패널은 cellline 풀로 못 돎
#   j15rand*  -> j_pool<GEN>_rand
#   j15randb* -> j_pool<GEN>_randb  (bl200 CpG 분포 맞춤)
#   그 외     -> j_pool<GEN>_cellline
# randb 를 rand 보다 먼저 확인 ('randb...' 도 'rand' 로 시작)
if PANEL.startswith('randb'):  _POOLDIR = 'j_pool' + GEN + '_randb'
elif PANEL.startswith('rand'): _POOLDIR = 'j_pool' + GEN + '_rand'
else:                          _POOLDIR = 'j_pool' + GEN + '_cellline'
_POOL = os.environ.get('POOL_DIR') or (ROOT + '/results/dmr/' + _POOLDIR)
# _cellline 폴더 없으면 j_pool<GEN>_세포주 사용 (POOL_DIR 직접 주면 그대로)
if (not os.environ.get('POOL_DIR') and not os.path.isdir(_POOL)
        and _POOLDIR.endswith('_cellline')):
    _old = _POOL[:-len('_cellline')] + '_세포주'
    if os.path.isdir(_old):
        _POOL = _old
FULL_NORMAL_BAM = _POOL + '/pool_normal_all.bam'
FULL_GBM_BAM    = _POOL + '/pool_gbm_all.bam'

_R = ROOT.rstrip('/\\')
_P = VERSION.split('_')[0]          # 패널 폴더 (jbl200, jjsd300 등)
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
    # 산출 폴더에 설정 스냅샷 저장 (import 때가 아니라 여기서)
    _cfg.snapshot(ML_OUT)


def guard_not_original():
    if not VERSION.startswith('j'):
        raise SystemExit('중단: VERSION=' + VERSION + '. j 로 시작해야 합니다.')
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

