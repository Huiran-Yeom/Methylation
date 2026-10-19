# -*- coding: utf-8 -*-
"""plasma_metric.py 의 기준선 두 개를 사람마다 잰다 (패널 없이 원자료에서)

  reads  사람별 총 리드 수. 기본은 BAM 의 1차 정렬 리드 수(보조·2차 정렬 제외).
         --reads-tsv 를 주면 그 값(sample, reads)을 쓴다. 예: 처리 로그의 「reads read」(정렬 전 원 리드 수).
  gmeth  전역 평균 메틸화 = bismark cov 의 Σmeth / Σ(meth + unmeth).

  기준선이 패널만큼 환자와 정상을 가르면, 패널 점수는 「암 신호」와 「플랫폼 · 깊이 ·
  전역 메틸화 차이」를 가르지 못한다. 그래서 plasma_metric.py 는 이 둘을 같은 순열로 잰다.
  입력은 검증의 CANCER 슬롯과 같은 사람의 자료여야 한다.

사용:  COHORT=<혈장코호트> python plasma_baseline.py --cov '/환자/{sample}.bismark.cov.gz' '/정상/{sample}.bismark.cov.gz' \\
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
ap.add_argument('--cov', required=True, nargs='+')   # 여러 개면 앞에서부터 있는 것을 쓴다
g = ap.add_mutually_exclusive_group()
g.add_argument('--bam')
g.add_argument('--reads-tsv')
ap.add_argument('--conf')
ap.add_argument('-o', '--out', required=True)
a = ap.parse_args()

COHORT = os.environ.get('COHORT') or sys.exit('COHORT 를 주십시오')
_h = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
cands = [a.conf] if a.conf else [os.path.expanduser('~/cohorts/%s.conf' % COHORT),
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
