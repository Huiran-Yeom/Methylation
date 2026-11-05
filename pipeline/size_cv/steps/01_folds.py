# -*- coding: utf-8 -*-
"""size_cv 01: 접기 나누기
  정상 cov, 암 cov 를 읽어 검체 id, 단위, 층을 정하고 단위째로 접기에 나눔
  층마다 단위를 바이트 순 정렬 → 1..NFOLD 에 번갈아 넣음 (시드 없음)
  접기별 학습 cov 링크, 빼 둔 정상 / 암 목록을 만듦. 원자료는 읽기만 하고 링크만 만듦
사용: python 01_folds.py --out <CV_ROOT> --nfold 5 [--strata <검체,층 csv>] [--dry]
  정상 cov = NORMAL_SET, 암 cov = CANCER_COV_DIR, 암 BAM = CANCER_BAM_DIR, 정상 BAM = NORMAL_BAM_DIR (옛 이름 GBM_* 도 받음)
  검체 id = SID_RE, 단위 = GROUP_RE (1단계 01 과 같은 규칙, 정상은 사람 = 단위)
  --strata 없으면 정상 / 암 두 층
"""
if __name__ != '__main__':
    raise ImportError(__file__ + ' 은(는) 스크립트다 — import 하지 않는다')
import glob
import os
import sys

import pandas as pd

a = sys.argv
OUT = a[a.index('--out') + 1]
NFOLD = int(a[a.index('--nfold') + 1]) if '--nfold' in a else 5
STRATA = a[a.index('--strata') + 1] if '--strata' in a else ''
DRY = '--dry' in a
ENV = os.environ
_P = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path[:0] = [os.path.join(_P, '0_setup'), os.path.join(_P, '1_dmr', 'steps')]
import meth_config as _cfg   # noqa: E402  (암 쪽 자료 이름, CANCER_LABEL)
import sample_id as _sid     # noqa: E402  (검체 이름, 단위 규칙)


ND, CD = ENV.get('NORMAL_SET', ''), _cfg.cancer_env('CANCER_COV_DIR', 'GBM_COV_DIR')
NB, CB = ENV.get('NORMAL_BAM_DIR', ''), _cfg.cancer_env('CANCER_BAM_DIR', 'GBM_BAM_DIR')
for k, v in (('NORMAL_SET', ND), ('CANCER_COV_DIR', CD), ('NORMAL_BAM_DIR', NB), ('CANCER_BAM_DIR', CB)):
    if not v or not os.path.isdir(v):
        sys.exit('멈춤: %s 폴더가 없음 — %s' % (k, v))
CANCER = _cfg.CANCER_LABEL
sid = _sid.sid


def unit(s, cancer):
    # 1단계 01 의 grp() 와 같은 규칙 (sample_id.py). GROUP_RE 는 암 쪽에만
    return _sid.grp(s, cancer)


def pairs(cov_dir, bam_dir, what):
    covs = sorted(glob.glob(cov_dir + '/*.cov.gz'))
    bams = {sid(b) for b in glob.glob(bam_dir + '/*.bam')}
    miss = [c for c in covs if sid(c) not in bams]
    if miss:
        sys.exit('멈춤: %s 짝 BAM 없는 검체 %d개' % (what, len(miss)))
    return covs


rows, src = [], {}
for c in pairs(ND, NB, '정상'):
    s = sid(c)
    rows.append((s, unit(s, False), 'Normal'))
    src[s] = c
for c in pairs(CD, CB, '암'):
    s = sid(c)
    rows.append((s, unit(s, True), CANCER))
    src[s] = c
F = pd.DataFrame(rows, columns=['sample', 'unit', 'stratum'])
if not F['sample'].is_unique:
    sys.exit('멈춤: 검체 id 가 겹침 (SID_RE 확인)')
if STRATA:
    st = pd.read_csv(STRATA).set_index('sample')['stratum']
    cm = F['stratum'] != 'Normal'
    miss = sorted(set(F.loc[cm, 'sample']) - set(st.index))
    if miss:
        sys.exit('멈춤: --strata 에 없는 암 검체 %d개' % len(miss))
    F.loc[cm, 'stratum'] = F.loc[cm, 'sample'].map(st).astype(str)
bad = F.groupby('unit')['stratum'].nunique()
if (bad > 1).any():
    sys.exit('멈춤: 한 단위가 층 둘 이상에 걸침 %d개' % int((bad > 1).sum()))

fold = {}
for t in sorted(F['stratum'].unique(), key=lambda x: x.encode('utf-8')):
    us = sorted(F.loc[F['stratum'] == t, 'unit'].unique(), key=lambda u: u.encode('utf-8'))
    for i, u in enumerate(us):
        fold[(t, u)] = i % NFOLD + 1
F['fold'] = [fold[(t, u)] for t, u in zip(F['stratum'], F['unit'])]
print('검체 · 단위:', {t: (int((F.stratum == t).sum()), int(F.loc[F.stratum == t, 'unit'].nunique()))
                     for t in sorted(F['stratum'].unique())})
print('접기별 단위 수:', {k: {t: int(F[(F.fold == k) & (F.stratum == t)]['unit'].nunique()) for t in sorted(F['stratum'].unique())}
                       for k in range(1, NFOLD + 1)})
for k in range(1, NFOLD + 1):
    h = F[F.fold == k]
    if not (h.stratum == 'Normal').any() or not (h.stratum != 'Normal').any():
        sys.exit('멈춤: 접기 %d 에 빼 둘 정상 또는 암 단위가 없음 (단위 수가 접기 수보다 적음)' % k)
if DRY:
    sys.exit(0)
fc = os.path.join(OUT, 'folds.csv')
if os.path.exists(fc) and not pd.read_csv(fc, dtype=str).equals(F.astype(str)):
    sys.exit('멈춤: %s 가 이미 있고 이번 배정과 다름. 다른 CV_ROOT 를 쓰거나 그 폴더를 지우고 다시' % fc)


def ln(target, rel):
    dst = os.path.join(OUT, rel)
    os.makedirs(os.path.dirname(dst), exist_ok=True)
    if os.path.lexists(dst):
        if os.path.realpath(dst) != os.path.realpath(target):
            sys.exit('멈춤: 이미 다른 곳을 가리키는 링크 — %s' % rel)
        return
    os.symlink(os.path.realpath(target), dst)


os.makedirs(OUT, exist_ok=True)
F.to_csv(os.path.join(OUT, 'folds.csv'), index=False)
for k in range(1, NFOLD + 1):
    tr, ho = F[F.fold != k], F[F.fold == k]
    for s, t in zip(tr['sample'], tr['stratum']):
        ln(src[s], 'fold%d/train_%s_cov/%s' % (k, 'normal' if t == 'Normal' else 'cancer', os.path.basename(src[s])))
    d = os.path.join(OUT, 'fold%d' % k)
    hn = ho[ho.stratum == 'Normal']
    hc = ho[ho.stratum != 'Normal']
    open(d + '/heldout_normal.txt', 'w').write('\n'.join(hn['sample']) + '\n')
    open(d + '/heldout_cancer_units.txt', 'w').write('\n'.join(sorted(hc['unit'].unique())) + '\n')
    open(d + '/heldout_cancer_samples.txt', 'w').write('\n'.join(hc['sample']) + '\n')
    # 학습 링크 = 이번 배정의 학습 검체 그대로인지 (옛 링크가 남으면 빼 둔 검체가 학습에 섞임)
    for w in ('normal', 'cancer'):
        want = {os.path.basename(src[s]) for s, t in zip(tr['sample'], tr['stratum']) if (t == 'Normal') == (w == 'normal')}
        have = set(os.listdir('%s/train_%s_cov' % (d, w)))
        if have != want:
            sys.exit('멈춤: fold%d/train_%s_cov 링크가 배정과 다름 (남는 %d · 빠진 %d)' % (k, w, len(have - want), len(want - have)))
print('저장: %s (접기 %d)' % (OUT, NFOLD))
