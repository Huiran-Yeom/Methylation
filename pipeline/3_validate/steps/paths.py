"""버전 태그로 스크립트 폴더, DMR 파일 찾기 (찾는 규칙은 여기에만)

TRAIN_SRC 먼저 보고, 없으면 runs/<그룹>/<버전>, 보관함, scripts_<버전> 순서로 찾음
config.conf 에서 읽은 값(VERSION_CONF 등)도 여기서 내보냄
"""
import glob
import os
import re

import os as _os, sys as _sys
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
    ROOT = _os.path.dirname(_cfg.TRAIN_SRC.rstrip('/'))
    TRAIN_SRC = _cfg.TRAIN_SRC.rstrip('/')
    DATA = _cfg.METH_ROOT
    # 설정은 여기서만 읽고 make_val_config 에 내보냄
    VERSION_CONF = _cfg.VERSION
    # hasattr 이름은 실제 함수 이름(get)과 같아야 함 - 다르면 STEP_MAP 이 안 쓰임
    STEP_MAP_CONF = _cfg.get('STEP_MAP') if hasattr(_cfg, 'get') else ''
    WORK_ROOT_CONF = getattr(_cfg, 'WORK_ROOT', '') or ''
    # 코호트 정의 폴더
    COHORT_DIR_CONF = _os.environ.get('COHORT_DIR') or _os.path.expanduser('~/cohorts')
except ImportError:
    raise SystemExit(
        '설정을 못 찾았습니다. 0_setup/meth_config.py 가 있는 자리를 찾지 못했습니다.\n'
        '  이 파일은 정리본 안에서 도는 것이라 기본값으로 넘어가지 않습니다.\n'
        '  트리 밖에서 돌리려면 METH_CONF_DIR 로 0_setup 자리를 주십시오.')


def _first(pats, want_dir=True):
    """패턴 순서대로 보고 처음 맞는 것 하나 반환"""
    ok = os.path.isdir if want_dir else os.path.isfile
    for p in pats:
        for h in sorted(glob.glob(p)):
            if ok(h):
                return h
    return None


def script_dir_abs(v):
    """버전의 스크립트 폴더 절대경로. 없으면 None"""
    # TRAIN_SRC 에 config.py 가 있으면(평평 배치) 그 폴더 반환
    #   TRAIN_SRC 가 버전 폴더면 중단 - v 와 상관없이 같은 폴더를 돌려줌
    if os.path.isfile(os.path.join(TRAIN_SRC, 'config.py')):
        if re.match(r'^j.*_\d+k_', os.path.basename(TRAIN_SRC)):
            raise SystemExit(
                'TRAIN_SRC 가 버전 폴더다 (%s). '
                '그룹 폴더(runs/<그룹>) 나 평평한 2_train/steps/ 을 가리켜야 한다.' % TRAIN_SRC)
        return TRAIN_SRC
    h = _first([
        os.path.join(TRAIN_SRC, v),                               # 정리본/그룹 아래
        os.path.join(ROOT, 'runs', '*', v),                       # 그룹 폴더
        os.path.join(ROOT, 'runs', '*', '_archive_versions', v),  # 보관함
        os.path.join(ROOT, 'runs', '*', '*', v),                  # 한 겹 더 들어간 경우
        os.path.join(ROOT, 'scripts_%s' % v),                     # 옛 배치
    ])
    if h:
        return h
    c = os.path.join(ROOT, 'scripts')
    return c if os.path.isdir(c) else None


def script_dir(v):
    """script_dir_abs 의 ROOT 기준 상대경로"""
    h = script_dir_abs(v)
    return os.path.relpath(h, ROOT) if h else None


def dmr_csv(v):
    """DMR_confirmed_<버전>.csv 경로. 없으면 None"""
    n = 'DMR_confirmed_%s.csv' % v
    R = os.path.join(DATA, 'results')
    return _first([
        os.path.join(R, v, 'step02_DMR', n),                # results/<버전>
        os.path.join(R, '*', v, 'step02_DMR', n),           # results/<무엇>/<버전>
        os.path.join(R, '*', '*', v, 'step02_DMR', n),      # _archive/old_versions/<버전>
        os.path.join(R, '*', '*', '*', v, 'step02_DMR', n), # 한 겹 더
    ], want_dir=False)

