# -*- coding: utf-8 -*-
"""4단계 공통 (make_tables.py, primary_metric.py 가 같이 씀)

person(검체 이름) → 사람 이름 (한 사람이 표 여럿 쥐지 않게)
  기본은 pat_1_..., ind_3_... 규칙. 다른 자료면 PERSON_RE 로, 괄호 하나 = 사람 이름
  예: PERSON_RE='^([A-Z0-9]+)-'
find_allcombo(접두, width, seen) → all_combos 채점 파일 목록
  없으면 옛 이름 _모든조합.csv 로 찾음 (_모든조합_포화전.csv 중간 산출은 안 걸림)
  seen[판_깊이] = (이름세대, width(첫 파일)) 를 남김 (gen_warn 이 읽음)
"""
import glob
import os
import re

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
        print('  ! PERSON_RE 가 검체 %d개 중 하나도 못 맞췄습니다. 사람 단위 묶기가 검체 단위와'
              ' 같아집니다. 환경변수 PERSON_RE 를 자료의 이름 규칙에 맞추세요.' % _pmiss[0])


def find_allcombo(prefix, width, seen):
    _k = prefix.split('/')[-1][len('y_prob_all_'):]
    fs = glob.glob(prefix + '_all_combos.csv')
    if fs:
        seen[_k] = ('all_combos', width(fs[0]))
        return fs
    fs = glob.glob(prefix + '_모든조합.csv')
    if fs:
        seen[_k] = ('모든조합', width(fs[0]))
        print('   [옛 이름으로 읽음] %s_모든조합.csv (%d개)' % (prefix.split('/')[-1], len(fs)))
    return fs
