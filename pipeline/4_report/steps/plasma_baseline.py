# -*- coding: utf-8 -*-
"""plasma_metric.py 기준선 두 개를 사람마다 계산 (패널 없이 원자료에서)

  reads  사람별 총 리드 수. 기본은 BAM 1차 정렬 리드 수 (보조, 2차 정렬 제외)
         --reads-tsv 주면 그 값(sample, reads) 사용 (예: 정렬 전 원 리드 수)
  gmeth  전역 평균 메틸화 = bismark cov 의 sum(meth) / sum(meth + unmeth)

  기준선도 패널만큼 가르면 암 신호와 플랫폼, 깊이 차이를 구별 못 함
  입력은 검증 CANCER 슬롯과 같은 사람 자료여야 함

사용법:  COHORT=<혈장코호트> python plasma_baseline.py --cov '/환자/{sample}.bismark.cov.gz' '/대조/{sample}.bismark.cov.gz' \\
            [--bam '/path/{sample}.bam' | --reads-tsv reads.tsv] [--conf x.conf] -o baseline.csv
"""
if __name__ != '__main__':
    raise ImportError(__file__ + ' 은(는) 스크립트다 — import 하지 않는다')

import argparse
import io
import os
import sys

import pandas as pd

ap = argparse.ArgumentParser()
ap.add_argument('--cov', required=True, nargs='+')   # 여러 개면 앞에서부터 있는 것 사용
g = ap.add_mutually_exclusive_group()
g.add_argument('--bam')
g.add_argument('--reads-tsv')
ap.add_argument('--conf')
ap.add_argument('-o', '--out', required=True)
a = ap.parse_args()

COHORT = os.environ.get('COHORT') or sys.exit('COHORT 를 주십시오')
_h = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
# 코호트 설정 찾는 순서: --conf → COHORT_CONF → COHORT_DIR(기본 ~/cohorts) → 저장소 예시
_env = os.environ.get('COHORT_CONF')
cands = [a.conf] if a.conf else ([_env] if _env else []) + [os.path.join(os.environ.get('COHORT_DIR') or os.path.expanduser('~/cohorts'), '%s.conf' % COHORT),
                                 os.path.join(_h, '3_validate', 'steps', 'cohorts', '%s.conf' % COHORT)]
p = next((c for c in cands if c and os.path.exists(c)), None) or sys.exit('코호트 정의 없음: ' + ' 또는 '.join(cands))
conf = {}
for ln in io.open(p, encoding='utf-8'):
    ln = ln.strip()
    if ln and not ln.startswith('#') and '=' in ln:
        k, v = ln.split('=', 1)
        conf[k.strip()] = v.strip()
SAMPLES = conf['SAMPLES'].split()
bam = a.bam or conf.get('CANCER')


def gmeth(s):
    p = next((c.format(sample=s) for c in a.cov if os.path.exists(c.format(sample=s))), None)
    if p is None:
        sys.exit('%s 의 cov 가 어느 패턴에도 없다' % s)
    me = un = 0
    for ch in pd.read_csv(p, sep='\t', header=None, usecols=[4, 5], chunksize=2_000_000):
        me += ch[4].sum(); un += ch[5].sum()
    return me / (me + un)


def bam_reads(s):
    import pysam
    with pysam.AlignmentFile(bam.format(sample=s), 'rb') as f:
        return sum(1 for r in f.fetch(until_eof=True) if not (r.is_secondary or r.is_supplementary))


if a.reads_tsv:
    rt = pd.read_csv(a.reads_tsv, sep='\t').set_index('sample')['reads']
    missing = [s for s in SAMPLES if s not in rt.index]
    if missing:
        sys.exit('--reads-tsv 에 없는 사람: %s' % missing[:10])
rows = []
for s in SAMPLES:
    r = int(rt[s]) if a.reads_tsv else bam_reads(s)
    rows.append((s, r, gmeth(s)))
    print('%s\treads %d\tgmeth %.4f' % rows[-1])
pd.DataFrame(rows, columns=['sample', 'reads', 'gmeth']).to_csv(a.out, index=False)
print('씀', os.path.abspath(a.out), len(rows))
