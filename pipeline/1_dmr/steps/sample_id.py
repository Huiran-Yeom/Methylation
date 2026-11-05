# -*- coding: utf-8 -*-
"""검체 이름 규칙 (1단계 01, 02, 03, 11 과 size_cv 01 이 같이 씀)

sid(파일 경로) → 검체 이름
  꼬리표(.bismark.cov.gz .cov.gz .sort.bam .bam)를 떼고 첫 _ 앞까지 (cov, BAM 꼬리표 달라도 짝 맞게)
  규칙이 안 맞으면 SID_RE (첫 괄호가 검체 이름, 안 걸리면 기본 규칙). 예:
    SID_RE='(.+?)__'   H1876_0d__SRR...  -> H1876_0d
                       H1876__SRR...     -> H1876

grp(검체 이름, cancer) → 단위 (반복본을 세포주 하나로 묶음)
  기본은 SNU<n>-<반복> → SNU<n>. 다른 자료는 GROUP_RE (첫 괄호가 세포주 이름). 예:
    GROUP_RE='(.+?)_\\d+d$'        H1876_0d -> H1876, H1876 -> H1876
  GROUP_RE 는 암 쪽에만 씀. 정상은 한 사람당 한 파일 전제
  안 맞는 이름은 그대로 한 단위 (안 주면 반복본이 서로 다른 세포주로 세어짐)
"""
import os
import re

_SRE = os.environ.get('SID_RE')
_srx = re.compile(_SRE) if _SRE else None
_GRE = os.environ.get('GROUP_RE')
_grx = re.compile(_GRE) if _GRE else None
SUFFIXES = ('.bismark.cov.gz', '.cov.gz', '.sort.bam', '.bam')


def sid(p):
    b = os.path.basename(p)
    for suf in SUFFIXES:
        if b.endswith(suf):
            b = b[:-len(suf)]
            break
    if _srx is not None:
        m = _srx.match(b)
        if m and m.groups():
            return m.group(1)
    return b.split('_')[0]


def grp(s, cancer):
    if cancer and _grx is not None:
        m = _grx.match(s)
        return m.group(1) if (m and m.groups()) else s
    return s.split('-')[0] if s.upper().startswith('SNU') and '-' in s else s
