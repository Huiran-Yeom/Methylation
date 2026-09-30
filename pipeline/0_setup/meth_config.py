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

# 코드가 있는 자리에서 기본값을 만든다. 예전에는 우리 서버의
#   ~/Methylation_rev1/... 을 기본값으로 박아 두었는데, 받아서 돌리는 사람에게는
#   없는 경로다. 검증이 「학습 스크립트를 못 찾은 채」 돌아 그 피처만 빠진다.
#   저장소 안을 가리키면 내려받은 그대로 동작한다.
_TREE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

_DEFAULTS = {
    'METH_ROOT': '/ssd_data/Methylation',
    'DATASET_ROOT': '/ssd_data/dataset',
    'WORK_ROOT': os.path.expanduser('~'),
    'TMP_ROOT': os.path.expanduser('~/tmp'),
    'TRAIN_SRC': os.path.join(_TREE, '2_train', 'steps'),
    'ANALYSIS_DIR': os.path.join(_TREE, '3_validate', 'steps'),
    'GEN': '15',
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
    'SEED': '42',
    'TEST_SIZE': '0.2',
    'N_JOBS': '1',
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

    2026-09-28: 셸 헤더가 `set -a; . config.conf` 로 conf 를 환경에 export 하므로
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
    return get(key).rstrip('/')


METH_ROOT = get_path('METH_ROOT')
DATASET_ROOT = get_path('DATASET_ROOT')
WORK_ROOT = get_path('WORK_ROOT')
TMP_ROOT = get_path('TMP_ROOT')
# conf 안의 $METH_ROOT·$DATASET_ROOT 를 여기서 푼다.
#   셸은 `set -a; . config.conf` 로 풀지만, 파이썬이 conf 를 직접 읽을 때는
#   그 둘이 「환경변수가 아니라서」 expandvars 가 못 푼다 — 문자 그대로 남는다.
def _expand(v):
    for _k, _val in (('METH_ROOT', METH_ROOT), ('DATASET_ROOT', DATASET_ROOT)):
        v = v.replace('${%s}' % _k, _val).replace('$' + _k, _val)
    return os.path.expanduser(os.path.expandvars(v)).rstrip('/')


TRAIN_SRC = get_path('TRAIN_SRC')
ANALYSIS_DIR = get_path('ANALYSIS_DIR')

GEN = get('GEN')
VAL_COHORTS = [x.strip() for x in get('VAL_COHORTS').split(',') if x.strip()]
PANELS = [x.strip() for x in get('PANELS').split(',') if x.strip()]
# 자체점검이 「conf 에 있는데 모듈 속성으로 없다」 를 잡았다.
#   paths.py 는 값('STEP_MAP') 으로 읽어 동작은 했지만, 다른 키와 다루는 법이
#   달라 읽는 사람이 헷갈린다. 노출 방식을 하나로 맞춘다.
STEP_MAP = get('STEP_MAP')
NULL_PANELS = [x.strip() for x in get('NULL_PANELS').split(',') if x.strip()]

def version_name(panel, depth):
    """접두(세대) + 명단(자료) + 깊이 를 한 곳에서 조립한다."""
    return 'j%s%s_%s_uni' % (GEN, panel, depth)
          # 패널 산출물 폴더 접미사 (j_candidates<GEN> 등)
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
TEST_SIZE = float(get('TEST_SIZE'))
MUT_RATIOS = [float(x) for x in get('MUT_RATIOS').replace(' ', '').split(',') if x]

# VERSION 이름 안의 깊이와 TARGET_DEPTH_READS 가 어긋나면 멈춘다.
#   옛 구조는 config.py 를 「생성」 해 둘이 같이 박혔으므로 어긋날 수 없었다.
#   설정을 한 파일로 빼면서 어긋날 수 있게 됐다. 적어 두는 것은 방어가 아니다.
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
      실제로 섞인 적이 있다. 09-14 에 MIN_READS_PER_WINDOW 를 4에서 6으로
      바꿨고, 「DMR 과 LLR 참조표의 깊이가 달랐다」·「3-1 표의 LLR참조
      193·193·194 는 100k 값이다」가 그 흔적이다.
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
    return p


if __name__ == '__main__':
    for k in ('METH_ROOT', 'DATASET_ROOT', 'WORK_ROOT', 'TMP_ROOT',
              'TRAIN_SRC', 'ANALYSIS_DIR', 'GEN', 'VAL_COHORTS', 'PANELS', 'NULL_PANELS', 'VERSION', 'PANEL',
              'TARGET_DEPTH_READS', 'DEPTH_MODE', 'MUT_RATIOS', 'PATTERN_K',
              'MIN_READS_PER_WINDOW', 'BLOCK_SIZE', 'MIN_SAMPLE_COVERAGE',
              'SEED', 'TEST_SIZE', 'N_JOBS'):
        print('%-22s %-44s (%s)' % (k, globals()[k], SOURCE.get(k, '-')))
