#!/usr/bin/env python
# coding: utf-8
"""
step04_1: bedGraph → Mean / Entropy feature 행렬

기준 : Methylation/scripts/step04_ML_classifier.py (최종 코드) 의
       `0) bed 파일 전처리`, `1) Mean`, `2) Entropy`

**최종 코드 구조를 그대로 따른다.** F드라이브 250226 의 가중치 엔트로피
(entropy_site / entropy_read / entropy_siteAread)는 최종 코드에 없으므로 넣지 않는다.
2026-09-25: 리드 단위 엔트로피(readent)·MHL 은 (마) 결정으로 쓰지 않는다.
  계산하던 step04_1b_readfeat.py 는 이 저장소에 없다(안 씀). 여기서는 안 만든다.

고친 것 (엔트로피 관련만)
  1. **엔트로피 계산식**: `SITE_ENTROPY_METHOD`
       'value_counts' : 최종 코드 방식. 서로 다른 실수값을 각각 하나의 범주로 센다.
                        📊 문제: 값 종류 수는 커버리지의 대리 지표다. GBM 리드가 섞이면
                        값 종류가 늘어 엔트로피가 올라간다(mut_0 5,555종 vs mut_5 6,441종).
                        즉 "외부 리드가 추가됐는지" 감지기로 작동할 수 있다.
       'histogram'    : 베타값을 10구간으로 묶어 계산(F:\\학부연구\\250207 방식.
                        F 드라이브는 옛 작업 기계로 이 저장소에 없다).
                        값 종류 증가에 둔감해 위 문제에 면역이다.
                        ★ 2026-08-20 부터 config 는 'histogram' 이다. 위 인공물을
                        피하려고 바꿨다. 08-13 에 적힌 'value_counts 다' 는 지금 틀리다.
                        지금 나오는 entropy 열은 히스토그램 방식이다.
  2. **좌표계**: bedGraph 는 0-based, step02 의 DMR 좌표는 cov(1-based) 기준이다.
     읽을 때 +1 해서 좌표계를 맞춘다. 최종 코드는 이 보정이 없어 블록 경계에서
     약 1% 가 다른 블록에 배정됐다.
  3. **파일 정렬**: 파일명 숫자순. 최종 코드는 크기 내림차순이라 Normal_i 가
     sampled_reads_i 와 무관해졌다. 숫자순이면 일치한다.
  4. **결측 블록 처리**: `MIN_SAMPLE_COVERAGE`
     최종 코드의 `dropna(how='any')` 는 샘플 하나라도 비면 블록을 버린다.
     depth 가 658배였을 때는 문제가 없었지만 총 5,000 reads 에서는 위험하다.
     📊 250226 실측(5샘플): 5개 전부 커버한 블록 849/1,139 = 74.5%
     반드시 `04_coverage_check.py` 를 먼저 돌려 실제 잔존 블록 수를 확인할 것.
       1.0  = 최종 코드와 동일
       0.8  = ★ 지금 설정. 20%까지 결측을 허용하고 아래에서 열평균으로 채운다.
       0.95 = 95% 이상 샘플에서 커버된 블록 유지

원본 불변: bedGraph 는 읽기만 하고 출력은 v3 폴더로만 간다.
"""

# 2026-09-26: 09-23 사고에서 가드 없는 이 부류가 importlib 로 「확인」 되다 실제로 돌았다.
#   82개 파일이 덮였다(내용은 결정적이라 같았지만 그건 운이었다).
if __name__ != '__main__':
    raise ImportError(__file__ + ' 은(는) 스크립트다 — import 하지 않는다')

import os
import sys

import numpy as np
import pandas as pd

import config as C

INPUT_FORMAT = getattr(C, 'INPUT_FORMAT', 'bedgraph')
FILE_ORDER = getattr(C, 'FILE_ORDER', 'numeric')
MIN_SAMPLE_COVERAGE = getattr(C, 'MIN_SAMPLE_COVERAGE', 1.0)

# 엔트로피 방식과 출력 이름을 명령줄로 바꾼다 (2026-08-14).
# value_counts 는 값 종류 수 = 커버리지 대리지표라 0.1% 에서 AUC 0.50 이었다.
# histogram(10구간) 은 종류 수에 둔감하므로 따로 만들어 비교한다.
_a = sys.argv
if '--entropy' in _a:
    C.SITE_ENTROPY_METHOD = _a[_a.index('--entropy') + 1]
ENT_NAME = _a[_a.index('--suffix') + 1] if '--suffix' in _a else 'entropy'

if '--ratios' in sys.argv:
    C.MUT_READ_COUNTS = [int(float(x) * C.TARGET_DEPTH_READS) for x in
                         sys.argv[sys.argv.index('--ratios') + 1].split(',')]
    print(f'[CLI] 대상 비율 -> {C.MUT_READ_COUNTS}')

KINDS = ('mean', 'entropy')


def load_dmr_set():
    df = pd.read_csv(C.DMR_CSV)
    chrom = df['Feature'].str.split('_').str.get(0).astype(str)
    blk = df['Feature'].str.split('_').str.get(1).astype(float).astype(int)
    return set(zip(chrom, blk))


def list_files(folder):
    import re

    suffix = '.bismark.cov.gz' if INPUT_FORMAT == 'cov' else '.bedGraph.gz'
    items = [
        (fn, os.path.getsize(os.path.join(folder, fn)))
        for fn in os.listdir(folder)
        if fn.endswith(suffix) and os.path.getsize(os.path.join(folder, fn)) > 0
    ]
    if FILE_ORDER == 'numeric':
        def num(fn):
            m = re.search(r'\d+', fn)
            return int(m.group()) if m else float('inf')

        items.sort(key=lambda x: num(x[0]))
    else:
        items = sorted(set(items), key=lambda x: x[1], reverse=True)
    return [fn for fn, _ in items]


def read_beta(path):
    """(chr, pos, beta[0~1]). pos 는 항상 1-based 로 맞춘다."""
    if INPUT_FORMAT == 'cov':
        df = pd.read_csv(
            path, sep='\t', header=None, usecols=[0, 1, 3],
            names=['chr', 'pos', 'beta'],
            dtype={'chr': str, 'pos': np.int64, 'beta': np.float64},
        )
    else:
        # 첫 줄은 `track type=bedGraph` (2026-07-30 실측 확인). 최종 코드의 df[1:] 와 동일
        df = pd.read_csv(
            path, sep='\t', header=None, skiprows=1, usecols=[0, 1, 3],
            names=['chr', 'pos', 'beta'],
            dtype={'chr': str, 'pos': np.int64, 'beta': np.float64},
        )
        df['pos'] = df['pos'] + 1  # 0-based -> 1-based (DMR 좌표계와 일치)
    df['beta'] = df['beta'] / 100.0
    return df


def _ent_value_counts(x):
    x = pd.Series(x).dropna()
    if x.empty:
        return np.nan
    p = x.value_counts(normalize=True)
    return float(-np.sum(p * np.log2(p)))


def _ent_histogram(x, bins=10):
    x = pd.Series(x).dropna()
    if x.empty:
        return np.nan
    h = np.histogram(x, bins=bins, range=(0.0, 1.0))[0]   # 08-20: 눈금을 0~1 로 고정
    h = h[h > 0]
    if h.size == 0:
        return 0.0
    p = h / h.sum()
    return float(-np.sum(p * np.log2(p)))


ENT = _ent_value_counts if C.SITE_ENTROPY_METHOD == 'value_counts' else _ent_histogram


def sample_features(path, dmr):
    """한 샘플의 블록별 (mean, entropy). 최종 코드의 groupby 계산과 동일."""
    df = read_beta(path)
    df['blk'] = df['pos'] // C.BLOCK_SIZE
    keys = list(zip(df['chr'].to_numpy(), df['blk'].to_numpy()))
    df = df[[k in dmr for k in keys]]
    if df.empty:
        return {}, {}
    g = df.groupby(['chr', 'blk'])['beta']
    mean = {f'{c}_{b}': v for (c, b), v in g.mean().items()}
    ent = {f'{c}_{b}': v for (c, b), v in g.apply(lambda x: ENT(x.to_numpy())).items()}
    return mean, ent


def collect(folder, prefix, dmr):
    files = list_files(folder)
    means, ents, mapping = {}, {}, []
    for i, fn in enumerate(files, 1):
        sid = f'{prefix}{i}'
        m, e = sample_features(os.path.join(folder, fn), dmr)
        means[sid], ents[sid] = m, e
        mapping.append({'sample_id': sid, 'file': fn})
        if i % 200 == 0:
            print(f'    {prefix} {i}/{len(files)}')
    return means, ents, mapping


def to_frame(rows_normal, rows_gbm, label):
    """행=샘플, 열=블록. MIN_SAMPLE_COVERAGE 로 결측 블록을 정리한다."""
    df = pd.DataFrame.from_dict({**rows_normal, **rows_gbm}, orient='index')
    n_total = len(df)
    before = df.shape[1]
    if MIN_SAMPLE_COVERAGE >= 1.0:
        df = df.dropna(axis=1, how='any')          # 최종 코드와 동일
    else:
        need = int(np.ceil(MIN_SAMPLE_COVERAGE * n_total))
        keep = df.notna().sum(axis=0) >= need
        df = df.loc[:, keep]
        # 남은 결측은 그 블록의 샘플 평균으로 채운다(모델이 NaN 을 못 받으므로)
        # ★ 누출 주의: 이 평균은 정상·암을 모두 합쳐 낸 값이고, 학습/검증을
        #    나누기 전에 채운다. 옳게 하려면 2_train/steps/08 에서 나눈 뒤 학습셋 평균으로
        #    채워야 한다. 📊 2026-08-13 실측: mean 행렬은 607열 중 2열(0.3%)만
        #    해당돼 영향이 작다. readent 는 크므로 반드시 확인할 것.
        df = df.fillna(df.mean(axis=0))
    print(f'    {label}: 블록 {before} -> {df.shape[1]} (커버율 기준 {MIN_SAMPLE_COVERAGE})')
    df['Type'] = [0 if s.startswith('Normal') else 1 for s in df.index]
    return df


def main():
    C.guard_not_original()
    C.banner('step04_1  Mean / Entropy feature 행렬')
    C.ensure_dirs()

    dmr = load_dmr_set()
    print(f'DMR {len(dmr)}개 · 입력 {INPUT_FORMAT} · 정렬 {FILE_ORDER}')
    print(f'엔트로피 {C.SITE_ENTROPY_METHOD} · 최소 커버율 {MIN_SAMPLE_COVERAGE}')
    if C.SITE_ENTROPY_METHOD == 'value_counts':
        print('  ! 경고: value_counts 는 값 종류 수(= 커버리지 대리지표)에 민감합니다.')
    print()

    normal_dir = f'{C.BISMARK_OUT}/mut_0_reads'
    if not os.path.isdir(normal_dir):
        raise SystemExit(f'정상 폴더가 없습니다. 2_train/steps/03 먼저.\n  {normal_dir}')
    print('정상(mut_0_reads) 처리 중...')
    n_mean, n_ent, n_map = collect(normal_dir, 'Normal', dmr)
    print(f'  {len(n_mean)}개 샘플\n')

    all_map = list(n_map)
    for p in C.MUT_READ_COUNTS:
        if p == 0:
            continue
        gbm_dir = f'{C.BISMARK_OUT}/mut_{p}_reads'
        if not os.path.isdir(gbm_dir):
            print(f'건너뜀: {gbm_dir} 없음')
            continue
        print(f'GBM(mut_{p}_reads) 처리 중...')
        g_mean, g_ent, g_map = collect(gbm_dir, 'GBM', dmr)
        all_map += [{**r, 'mut_reads': p} for r in g_map]

        for kind, nr, gr in (('mean', n_mean, g_mean), (ENT_NAME, n_ent, g_ent)):
            df = to_frame(nr, gr, kind)
            out = f'{C.ML_OUT}/{p}_{kind}_{C.VERSION}.csv'
            df.to_csv(out, index=True, index_label='sample')
            print(f'    저장 {os.path.basename(out)}  ({df.shape[0]}행 x {df.shape[1]-1}블록)')
        print()

    pd.DataFrame(all_map).to_csv(f'{C.ML_OUT}/sample_file_mapping_{C.VERSION}.csv', index=False)
    print(f'샘플↔파일 매핑 저장: sample_file_mapping_{C.VERSION}.csv')
    print('step04_1 완료.')


if __name__ == '__main__':
    sys.exit(main())
