"""버전 태그 -> 스크립트 폴더 · DMR 파일 (찾는 규칙은 이 파일 한 곳에만).

2026-08-10: 스크립트가 scripts_<버전> 에서 runs/<그룹>/<버전> 으로 옮겨졌다.
2026-08-14: 폴더를 보관함으로 옮겼더니 검증이 끊겼다. 한 자리만 보던 것을
            여러 자리를 차례로 보도록 넓힌다. 앞으로 어디로 옮겨도 안 깨진다.
2026-09-25: 정리본에서는 학습 코드가 2_train/steps/ 한 폴더에 평평하게 있다.
            TRAIN_SRC 를 먼저 보고, 옛 배치는 그 뒤에 그대로 둔다.
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
    # 설정을 읽는 자리는 트리에서 여기 하나다. VERSION 도 같이 내보내
    #   부르는 쪽(make_val_config)이 또 다른 조회 경로를 만들지 않게 한다.
    VERSION_CONF = _cfg.VERSION
    # STEP_MAP 도 여기서 내보낸다 — 설정을 읽는 자리를 하나로 유지한다.
    # 관문이 hasattr(_cfg, '값') 이었다 — 식별자 영어화 때 호출은
    #   get() 으로 바꾸고 이 「이름 문자열」 은 안 바꿔서 영원히 False 였다.
    #   결과: STEP_MAP_CONF 가 항상 '' 이라 config.conf 의 STEP_MAP 이 죽어 있었다.
    STEP_MAP_CONF = _cfg.get('STEP_MAP') if hasattr(_cfg, 'get') else ''
    # WORK_ROOT 도 여기서 내보낸다. make_val_config 가 _cfg 를 직접
    #   쓰려다 「import 하지 않은 이름」 을 써서 NameError 가 났다 (WORK_ROOT 가
    #   빈 값일 때만 터지는 자리였다). 설정을 읽는 자리를 하나로 유지한다.
    WORK_ROOT_CONF = getattr(_cfg, 'WORK_ROOT', '') or ''
    # 코호트 정의 폴더도 여기서 내보낸다 (설정을 읽는 자리는 하나다).
    COHORT_DIR_CONF = _os.environ.get('COHORT_DIR') or _os.path.expanduser('~/cohorts')
except ImportError:
    raise SystemExit(
        '설정을 못 찾았습니다. 0_setup/meth_config.py 가 있는 자리를 찾지 못했습니다.\n'
        '  이 파일은 정리본 안에서 도는 것이라 기본값으로 넘어가지 않습니다.\n'
        '  트리 밖에서 돌리려면 METH_CONF_DIR 로 0_setup 자리를 주십시오.')


def _first(pats, want_dir=True):
    """패턴을 적은 순서대로 훑어, 처음 맞는 것 하나를 준다."""
    ok = os.path.isdir if want_dir else os.path.isfile
    for p in pats:
        for h in sorted(glob.glob(p)):
            if ok(h):
                return h
    return None


def script_dir_abs(v):
    """버전의 스크립트 폴더: 절대경로. 없으면 None."""
    # 정리본: 학습 코드가 한 폴더에 평평하게 있으면 그 폴더가 곧 답이다.
    #   단 「버전 폴더」 를 TRAIN_SRC 로 잡으면 안 된다. 거기도 config.py 가 있어서
    #   v 가 무엇이든 같은 폴더를 돌려준다. jsd200 검증이 bl200 스크립트를 복사한다.
    if os.path.isfile(os.path.join(TRAIN_SRC, 'config.py')):
        if re.match(r'^j.*_\d+k_', os.path.basename(TRAIN_SRC)):
            raise SystemExit(
                'TRAIN_SRC 가 버전 폴더다 (%s). '
                '그룹 폴더(runs/<그룹>) 나 평평한 2_train/steps/ 을 가리켜야 한다.' % TRAIN_SRC)
        return TRAIN_SRC
    h = _first([
        os.path.join(TRAIN_SRC, v),                               # 정리본/그룹 아래
        os.path.join(ROOT, 'runs', '*', v),                       # 지금 자리
        os.path.join(ROOT, 'runs', '*', '_archive_versions', v),  # 보관함
        os.path.join(ROOT, 'runs', '*', '*', v),                  # 한 겹 더 들어간 경우
        os.path.join(ROOT, 'scripts_%s' % v),                     # 옛 배치
    ])
    if h:
        return h
    c = os.path.join(ROOT, 'scripts')
    return c if os.path.isdir(c) else None


def script_dir(v):
    """예전 이름 그대로: ROOT 기준 상대경로."""
    h = script_dir_abs(v)
    return os.path.relpath(h, ROOT) if h else None


def dmr_csv(v):
    """DMR_confirmed_<버전>.csv 가 실제로 있는 자리. 없으면 None."""
    n = 'DMR_confirmed_%s.csv' % v
    R = os.path.join(DATA, 'results')
    return _first([
        os.path.join(R, v, 'step02_DMR', n),                # 지금 자리
        os.path.join(R, '*', v, 'step02_DMR', n),           # results/<무엇>/<버전>
        os.path.join(R, '*', '*', v, 'step02_DMR', n),      # _archive/old_versions/<버전>
        os.path.join(R, '*', '*', '*', v, 'step02_DMR', n), # 한 겹 더
    ], want_dir=False)

