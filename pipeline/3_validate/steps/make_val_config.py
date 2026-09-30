#!/usr/bin/env python
# coding: utf-8
"""외부검증용 작업 폴더와 config.py 를 만든다 (코호트 무관).

make_snu_config.py 를 일반화한 것. 코호트 정의를 COHORT_DIR(기본 ~/cohorts)의 <이름>.conf 에서
읽으므로 SNU 말고 다른 데이터(EGAD 등)도 같은 코드로 돌릴 수 있다.

버전의 config.py 를 그대로 읽어 뒤에 오버라이드 블록만 붙인다.
깊이·믹싱 방식·엔트로피 계산은 손대지 않는다. 학습 때와 같아야 하기 때문이다.
바꾸는 것은 입력 BAM 과 출력 경로뿐이다.

경로는 {코호트}_val/{검체}/... 로 잡는다. 코호트가 snu 면 기존 snu_val/... 과
같아져서 이미 만든 결과를 그대로 쓴다.

사용: python make_val_config.py <코호트> <검체> <version> [반복수]
"""
# 아무도 import 하지 않는다 (run_one_sample.sh 가 인자를 주고 부른다).
#   사고가 가드 없는 스크립트를 importlib 로 부른 데서 났다.
if __name__ != '__main__':
    raise ImportError(__file__ + ' 은(는) 스크립트다 — import 하지 않는다')

import io
import os
import shutil
_LAYOUT = [None]   # 코호트 리드 방식 자동 감지 (2026-08-11)
import sys

# 폴더가 runs/<그룹>/<버전> 으로 옮겨졌고, e·f·g 세트도
# 검증에 써야 해서 고정 딕셔너리를 없앴다. 찾는 규칙은 paths.py 한 곳에 있다.
from paths import (script_dir, script_dir_abs, VERSION_CONF, WORK_ROOT_CONF,
                   COHORT_DIR_CONF)


# 본 실행과 이름을 다르게 복사한다. 러너들이 pgrep 으로 프로세스를 세어
# 동시 실행을 제한하는데, 이름이 같으면 검증이 본 실행의 자리를 빼앗는다.
# STEPS 를 config.conf 의 STEP_MAP 에서 만든다. 두 정본이면 어긋난다.
#   STEP_MAP 형식:  <검체폴더에 놓일 이름>:<원본 파일명>;...
#   못 읽으면 아래 기본값으로 간다.
_STEPS_DEFAULT = {'01_Preprocessing_for_ML_RawData_Mixing.py': 'val03_1.py',
         '02_Preprocessing_for_ML_RawData_Mixing.py': 'val03_2.py',
         '03_Preprocessing_for_ML_RawData_Mixing.py': 'val03_3.py',
         '04_coverage_check.py': 'val04_0.py',
         '05_Feature_Matrix.py': 'val04_1.py',
         # 리드 엔트로피·MHL 을 검증에서도 만들 수 있게 추가.
         # 없는 버전에서는 위 복사 루프가 조용히 넘어간다.
         # readent·mhl 은 (마) 결정으로 안 쓴다. 복사하지 않는다.
         #   되살리려면 그 파일을 2_train/steps/ 에 먼저 두어야 한다 — 이 저장소에는 없다.
         # 'step04_1b_readfeat.py': 'val04_1b.py',
         # JSD 피처. 이 목록에 없으면
         #   검증에서 JSD-2 를 못 쟀다. 파일은 판 폴더가 아니라 한 단계 위에 있다.
         '06_jsdfeat.py': 'val04_1c.py',
         # [실행차단 교정] LLR 이 이 목록에 없어 검증이 llr 을 아예 안 만들었다.
         #   사전등록 1차 조합 mean|entropy|jsd|llr|pdr 이 통째로 계산될 수 없는 상태였다.
         #   옛 판은 별도 드라이버로 따로 깔아 돌렸다. 그래서 안 보였다.
         '07_llrfeat.py': 'val04_1d.py'}

try:
    from paths import STEP_MAP_CONF as _SM
except Exception:
    _SM = ''
if _SM:
    STEPS = {}
    for _pair in _SM.split(';'):
        if ':' not in _pair:
            continue
        _left, _right = _pair.split(':', 1)
        STEPS[_right.strip()] = _left.strip() + '.py'
    if not STEPS:
        sys.exit('STEP_MAP 을 읽었지만 항목이 없습니다. config.conf 를 보십시오.')
else:
    STEPS = _STEPS_DEFAULT

OVERRIDE = '''

# =====================================================================
# 외부검증 오버라이드: make_val_config.py 자동 생성. 손대지 말 것.
#
# 학습에 쓰지 않은 정상 검체를 배경으로, 암 검체 1종을 학습 때와 똑같은
# 비율·깊이·믹싱 방식으로 섞는다. 조건이 같아야 결과를 비교할 수 있다.
# =====================================================================
COHORT = '%(cohort)s'
SAMPLE = '%(sample)s'
DMR_SOURCE_VERSION = VERSION      # step02 는 돌리지 않는다
NUM_FILES = %(nrep)d              # 비율당 반복
%(ratios)s

# 중복 추출로 채우지 않는다.
#   학습은 암이 세포주라 풀이 깊어 이 경로를 탄 적이 없다(계획표·로그로 확인).
#   검증에서 중복으로 채우면 그 비율의 복제본 10판이 서로 독립이 아니게 되고,
#   판 사이 흩어짐이 작아져 성능이 실제보다 좋게 보인다.
#   대신 '부족한 비율만' 건너뛴다. 검체 전체를 버리지 않는다.
ALLOW_REPLACEMENT = False
SKIP_SHORT_RATIOS = True

FULL_NORMAL_BAM = '%(normal)s'
# 음성 배경을 공여자별로 가른다. 반복 i 가 배경 i 를 쓴다.
#   비어 있으면 예전처럼 병합본 하나에서 전부 뽑는다.
NORMAL_BAMS = %(normals)s
FULL_GBM_BAM = '%(cancer)s'

# 원본 config 의 guard_not_original() 은 출력 폴더 이름이 VERSION 이기를
# 요구한다. 그 보호장치를 끄지 않고 구조를 규칙에 맞춘다.
_VAL = f'{_R}/{COHORT}_val/{SAMPLE}'
SD = f'{_VAL}/sample_data/panels/{VERSION.split("_")[0]}/{VERSION}'
RS = f'{_VAL}/results/{VERSION}'

DMR_DIR = f'{RS}/step02_DMR'
# 본 실행에서 확정한 목록을 그대로 읽는다. 검증용으로 다시 뽑지 않는다
def _find_dmr(*_vs):
    """DMR 파일을 찾는다. 이름 후보를 순서대로, 자리 후보를 순서대로.
    2026-08-14: 폴더를 보관함으로 옮겼다고 검증이 끊기던 것을 막는다.
    SHARE_FROM 으로 패널을 공유하는 버전은 자기 이름의 파일이 없으므로
    VERSION 다음에 _S 도 본다."""
    import glob as _g
    _seen = []
    for _v in [_x for _x in _vs if _x]:
        _n = f'DMR_confirmed_{_v}.csv'
        _seen.append(_n)
        for _p in (f'{_R}/results/{_v}/step02_DMR/{_n}',
                   f'{_R}/results/*/{_v}/step02_DMR/{_n}',
                   f'{_R}/results/*/*/{_v}/step02_DMR/{_n}',
                   f'{_R}/results/*/*/*/{_v}/step02_DMR/{_n}'):
            _h = sorted(_g.glob(_p))
            if _h:
                return _h[0]
    raise SystemExit('[DMR 없음] results 아래 어디에도 없습니다: '
                     + ' , '.join(_seen))


# j 판 패널은 results/j_panel_dmr/<PANEL>/ 에 있다.
#   _find_dmr 은 step02_DMR 폴더만 찾으므로 j 판에서는 못 찾는다.
#   PANEL 이 있고 그 자리에 파일이 있으면 그것을 쓰고, 없으면 옛 방식으로 간다.
import os as _os2
_pj = ''
if 'PANEL' in dir():
    # [실행차단 교정] 자리가 'j_panel_dmr' 하나로 박혀 있었다.
    #   새 판(j15*)의 패널은 'j_panel_dmr15_panels' 에 있는데, 옛 자리에 같은 이름의
    #   옛 파일이 남아 있어서 검증이 조용히 옛 패널을 읽었다.
    #   실측: 옛 bl200 과 새 bl200 은 200칸 중 35칸만 겹친다 (열에 여덟이 다르다).
    #   그대로 뒀으면 학습은 새 패널, 검증은 옛 패널이 된다.
    #   자리를 VERSION 으로 가른다. 추측하지 않는다.
    # 'j15' 를 못 박으면 GEN 을 바꿨을 때 「옛 자리」(j_panel_dmr)로 내려간다.
    #   그 자리에 같은 이름의 옛 파일이 남아 있으면 조용히 옛 패널을 읽는다.
    # 이 줄은 「생성되는 config.py 안에서」 돈다. 그 파일은
    #   meth_config 를 _cfg 로 import 하므로 이름이 맞아야 한다.
    #   토큰 치환은 문자열 안을 못 건드려서 영어화 때 여기만 남았다.
    _g = _cfg.GEN
    _cands = ([_R + '/results/dmr/j_panel_dmr' + _g + '_panels/']
              if VERSION.startswith('j' + _g) else
              [_R + '/results/dmr/j_panel_dmr/'])
    for _c in _cands:
        _t = _c + PANEL + '/DMR_confirmed_' + PANEL + '.csv'
        if _os2.path.exists(_t):
            _pj = _t
            break
    if not _pj:
        raise SystemExit('중단: 패널을 못 찾았습니다 — ' + VERSION
                         + ' · ' + ' · '.join(_cands))
if _pj and _os2.path.exists(_pj):
    DMR_DIR = _os2.path.dirname(_pj)
    DMR_CSV = _pj
else:
    DMR_CSV = _find_dmr(VERSION, _S)

# 학습 판의 SHARE_FROM 이 그대로 복사돼 들어온다. 아래에서 MIX_SD 를
#   덮어쓰므로 지금도 맞게 돌지만, 변수가 남아 있으면 읽는 사람이 헷갈리고
#   누가 아래 줄을 지우면 검증이 조용히 학습 판의 조각 풀을 읽는다. 못 박는다.
SHARE_FROM = None
_S = VERSION
MIX_SD = f'{SD}/step03_Preprocessing_for_ML'
MUT_SORTED = f'{MIX_SD}/mut_sorted_{VERSION}.bam'
WT_SORTED = f'{MIX_SD}/wt_sorted_{VERSION}.bam'
MUT_NAMES = f'{MIX_SD}/mut_names_{VERSION}.txt'
WT_NAMES = f'{MIX_SD}/wt_names_{VERSION}.txt'

MIX_OUT = f'{RS}/step03_Preprocessing_for_ML'
BISMARK_OUT = f'{SD}/step04_ML_classifier'
ML_OUT = f'{RS}/step04_ML_classifier'
READ_ENT_OUT = f'{RS}/step04_read_entropy'

# val03_2 가 만들지 못한 비율을 적어 둔다. 뒷단계(val03_3 · val04_1 · val04_1b ·
# val04_2)는 모두 MUT_RATIOS / MUT_READ_COUNTS 를 돌기 때문에, 여기서 한 번
# 걸러내면 각 스크립트를 따로 고칠 필요가 없다.
import os as _os
_skip_f = MIX_OUT + '/skipped_ratios_' + VERSION + '.txt'
if _os.path.exists(_skip_f):
    _skipped = set(float(_x) for _x in open(_skip_f).read().split())
    MUT_RATIOS = [_r for _r in MUT_RATIOS if _r not in _skipped]
    MUT_READ_COUNTS = [int(_r * TARGET_DEPTH_READS) for _r in MUT_RATIOS]
    # 안내문은 화면(stdout)이 아니라 오류 쪽으로 보낸다.
    # run_val.sh 가 python -c 로 숫자만 받아가는데, 여기 찍히면 값에 섞인다.
    import sys as _sys
    print('  [config] 풀 부족으로 건너뛴 비율 '
          + ', '.join(str(_x) for _x in sorted(_skipped)), file=_sys.stderr)
'''


def load_conf(cohort):
    """COHORT_DIR 의 <이름>.conf 를 읽는다 (기본은 홈의 cohorts). KEY=값 한 줄씩, # 은 주석."""
    # 조회 자리가 run_validate.sh 와 갈려 있었다. 같은 순서로 맞춘다 —
    #   ~/cohorts 먼저(기관 자료), 없으면 저장소의 cohorts/ (예시·동봉본).
    # COHORT_CONF 로 못 박을 수 있게 한다. 작성자 기계에서는
    #   ~/cohorts 가 먼저라 저장소의 예시가 「검사된 적 없는 파일」 로 남는다.
    _env = os.environ.get('COHORT_CONF')
    if _env:
        if not os.path.exists(_env):
            sys.exit('COHORT_CONF 가 가리키는 파일이 없습니다: %s' % _env)
        p, _cands = _env, [_env]
    else:
        _cands = [os.path.join(COHORT_DIR_CONF, '%s.conf' % cohort),
                  os.path.join(os.path.dirname(os.path.abspath(__file__)),
                               'cohorts', '%s.conf' % cohort)]
        p = next((c for c in _cands if os.path.exists(c)), None)
    if p is None:
        sys.exit('코호트 정의 없음: ' + ' 또는 '.join(_cands))
    # 어느 파일을 읽었는지 「절대경로로」 남긴다. 둘 중 어디였나가 나중에 문제가 된다.
    print('  코호트 정의  %s' % os.path.abspath(p))
    d = {}
    for ln in io.open(p, encoding='utf-8'):
        ln = ln.strip()
        if not ln or ln.startswith('#') or '=' not in ln:
            continue
        k, v = ln.split('=', 1)
        d[k.strip()] = v.strip()
    for k in ('NORMAL', 'CANCER', 'SAMPLES'):
        if k not in d:
            sys.exit('%s 에 %s 가 없습니다' % (p, k))
    return d


def main():
    if len(sys.argv) < 4:
        sys.exit(__doc__)
    cohort, sample, v = sys.argv[1], sys.argv[2], sys.argv[3]
    conf = load_conf(cohort)
    # 깊이마다 배경 목록과 반복 수가 다르다. 판 이름에서 깊이를 읽는다.
    #   얕은 공여자는 깊은 판을 혼자 못 채운다. 섞으면 '사람' 과 '병합본' 이 한 표에 섞인다.
    _dep = next((x for x in ('5k', '15k', '50k', '100k') if ('_%s_' % x) in v), '')
    _dk = {'5k': '5K', '15k': '15K', '50k': '50K', '100k': '100K'}.get(_dep, '')
    _nrep_key = 'NREP_' + _dk if _dk and ('NREP_' + _dk) in conf else 'NREP'
    _norm_key = 'NORMALS_' + _dk if _dk and ('NORMALS_' + _dk) in conf else 'NORMALS'
    nrep = int(sys.argv[4]) if len(sys.argv) > 4 else int(conf.get(_nrep_key, 10))
    if script_dir(v) is None:
        sys.exit('알 수 없는 버전: %s' % v)
    if sample not in conf['SAMPLES'].split():
        sys.exit('%s 코호트에 없는 검체: %s' % (cohort, sample))

    # WORK_ROOT 를 따른다 (config.conf 가 「검체 작업 폴더의 부모」 로
    #   안내하는 값이다). 여기만 ~ 를 박아 두면 run_validate.sh 의 채점 루프가
    #   다른 자리를 보고 전 검체를 조용히 건너뛴다.
    home = os.environ.get('WORK_ROOT') or WORK_ROOT_CONF or os.path.expanduser('~')
    # ~/Methylation_rev1 로 되짚지 않는다 — paths 가 이미 절대경로를 안다.
    #   되짚으면 TRAIN_SRC 를 딴 데로 옮겼을 때 조용히 옛 코드를 복사한다.
    src = script_dir_abs(v)
    # ── [판이름 ↔ VERSION] ────────────────────────────────
    #   평평 배치(steps/ 에 config.py 하나)에서는 어떤 v 를 줘도 같은 config.py 를
    #   복사한다. 그 config.py 의 VERSION·깊이는 config.conf 에서 온다.
    #   작업 폴더는 v 로 만들고 산출 자리는 VERSION 으로 만들기 때문에, 둘이
    #   어긋나면 「50k 폴더에서 5k 깊이로 섞어 5k 자리에 쓰는」 일이 조용히 일어난다.
    #   meth_config.py 의 깊이 가드는 conf 안의 VERSION↔깊이만 보므로 이걸 못 잡는다.
    _flat = os.path.isfile(os.path.join(src, 'config.py'))
    if _flat and v != VERSION_CONF:
        sys.exit('중단: 판이름과 설정의 VERSION 이 다릅니다.'
                 '  인자 v=%s · config.conf VERSION=%s' % (v, VERSION_CONF)
                 + '  평평 배치에서는 config.py 가 하나라서 둘이 같아야 합니다.'
                   '  config.conf 의 VERSION 을 바꾸거나 VERSION=%s 로 주십시오.' % v)
    work = os.path.join(home, '%s_val' % cohort, v, sample)
    os.makedirs(work, exist_ok=True)

    for f, newname in STEPS.items():
        s, d = os.path.join(src, f), os.path.join(work, newname)
        if os.path.islink(d) or os.path.exists(d):
            os.remove(d)
        # 없는 파일을 조용히 넘기면 그 피처만 빠진 채 끝까지 돌고,
        #   채점이 결측을 평균으로 메워 AUC 0.500 이 된다.
        #   피처를 빼는 것은 STEPS 항목을 지워서 한다. 여기 있는데 파일이 없으면 항상 설정 오류다.
        if not os.path.exists(s):
            sys.exit('복사할 파일이 없다: %s — STEPS 에는 있는데 TRAIN_SRC 에 없다. '
                     '이름을 바꿨거나 TRAIN_SRC 가 딴 데를 가리킨다.' % s)
        if True:
            # 복사한다. 심볼릭 링크면 파이썬이 링크를 따라간 원본 폴더를
            # sys.path[0] 으로 잡아 원본 config.py 를 읽어버린다.
            shutil.copy2(s, d)
            if _LAYOUT[0] is None:
                from val_single import detect_layout
                try:
                    _cn = conf.get('CANCER', '').replace('{sample}', sample)
                except NameError:
                    _cn = ''
                _LAYOUT[0] = str(conf.get('LAYOUT') or
                                 detect_layout(conf.get('NORMAL', ''), _cn))
                print('  리드 방식: %s' % _LAYOUT[0])
            if _LAYOUT[0].startswith('single'):
                from val_single import to_single
                to_single(d)

    base = io.open(os.path.join(src, 'config.py'), encoding='utf-8').read()
    _rt = str(conf.get('RATIOS', '')).strip()
    if _rt:
        _rt = ('MUT_RATIOS = [' + _rt + ']\n'
               'MUT_READ_COUNTS = [int(r * TARGET_DEPTH_READS) for r in MUT_RATIOS]')
        print('  비율 재정의: ' + _rt.split('\n')[0])
    else:
        _rt = '# RATIOS 미지정 - 학습과 동일'
    io.open(os.path.join(work, 'config.py'), 'w', encoding='utf-8').write(
        base + OVERRIDE % {'cohort': cohort, 'sample': sample, 'nrep': nrep,
                           'ratios': _rt,
                           'normal': conf['NORMAL'],
                           'normals': repr([x for x in conf.get(_norm_key, '').split() if x]),
                           'cancer': conf['CANCER'].replace('{sample}', sample)})

    sys.path.insert(0, work)
    import importlib
    if 'config' in sys.modules:
        del sys.modules['config']
    C = importlib.import_module('config')
    print('  %s / %s / %s' % (cohort, sample, v))
    print('    깊이 %d · 비율 %d종 · 반복 %d  -> 검체 %d개'
          % (C.TARGET_DEPTH_READS, len(C.MUT_RATIOS), C.NUM_FILES,
             len(C.MUT_RATIOS) * C.NUM_FILES))
    print('    정상 %s' % os.path.basename(C.FULL_NORMAL_BAM))
    print('    암   %s' % os.path.basename(C.FULL_GBM_BAM))
    print('    결과 %s' % C.RS)
    for n, p in (('정상 BAM', C.FULL_NORMAL_BAM), ('암 BAM', C.FULL_GBM_BAM),
                 ('DMR', C.DMR_CSV)):
        if not os.path.exists(p):
            sys.exit('    X 없음: %s  %s' % (n, p))
    print('    OK 입력 3종 확인 · %s' % work)
    return 0


if __name__ == '__main__':
    sys.exit(main())
