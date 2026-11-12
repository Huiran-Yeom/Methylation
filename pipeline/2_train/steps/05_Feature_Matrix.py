#!/usr/bin/env python
# coding: utf-8
"""
2_train/steps/05: bedGraph → 블록별 mean, entropy 행렬 (비율마다 csv)

entropy 는 SITE_ENTROPY_METHOD (기본 histogram, 베타값 10구간)
bedGraph 는 0-based 라 +1 해서 DMR 좌표(1-based)에 맞춤, 파일은 이름 숫자순
학습은 MIN_SAMPLE_COVERAGE 미만 커버 블록 버림 (남은 빈칸은 08, 09 에서 채움)

옵션: --entropy <방식>  --suffix <출력 이름>  --ratios 0.001,0.01
"""

# 스크립트 전용 - import 하면 본체가 돌아 자료를 덮어씀
if __name__ != '__main__':
    raise ImportError(__file__ + ' 은(는) 스크립트다. import 하지 않는다')

import gzip
import os
import sys

import numpy as np
import pandas as pd

import config as C

INPUT_FORMAT = getattr(C, 'INPUT_FORMAT', 'bedgraph')
FILE_ORDER = getattr(C, 'FILE_ORDER', 'numeric')
MIN_SAMPLE_COVERAGE = getattr(C, 'MIN_SAMPLE_COVERAGE', 1.0)
COLSET_FROM = getattr(C, 'COLSET_FROM', 'all')     # all | normal - 어느 행 커버율로 블록 고를지

# 엔트로피 방식, 출력 이름 명령줄로 변경 (value_counts 는 커버리지 따라가서 비교용)
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
    """(chr, pos, beta[0~1]) 읽음. pos 는 1-based 로 맞춤"""
    if INPUT_FORMAT == 'cov':
        df = pd.read_csv(
            path, sep='\t', header=None, usecols=[0, 1, 3],
            names=['chr', 'pos', 'beta'],
            dtype={'chr': str, 'pos': np.int64, 'beta': np.float64},
        )
    else:
        # 첫 줄이 track 일 때만 건너뜀 (03 이 다시 쓴 파일엔 없음, 무조건 건너뛰면 첫 CpG 빠짐)
        with gzip.open(path, 'rt') as _h:
            _skip = 1 if _h.readline().startswith('track') else 0
        df = pd.read_csv(
            path, sep='\t', header=None, skiprows=_skip, usecols=[0, 1, 3],
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
    h = np.histogram(x, bins=bins, range=(0.0, 1.0))[0]   # 구간 0~1 고정
    h = h[h > 0]
    if h.size == 0:
        return 0.0
    p = h / h.sum()
    return float(-np.sum(p * np.log2(p)))


ENT = _ent_value_counts if C.SITE_ENTROPY_METHOD == 'value_counts' else _ent_histogram


def sample_features(path, dmr):
    """한 샘플의 블록별 (mean, entropy)"""
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


# 검증 실행인지 - 파일 이름(val*) 또는 METH_VALIDATE=1 로 판단
_VAL = os.path.basename(__file__).startswith('val') or os.environ.get('METH_VALIDATE') == '1'


def to_frame(rows_normal, rows_gbm, label):
    """행=샘플, 열=블록. 학습이면 MIN_SAMPLE_COVERAGE 로 블록 거름"""
    df = pd.DataFrame.from_dict({**rows_normal, **rows_gbm}, orient='index')
    n_total = len(df)
    before = df.shape[1]
    if _VAL:
        # 검증은 블록 안 버리고 안 채움 - 빈칸은 score_one_sample 이 학습 평균으로 채움
        pass
    else:
        # 커버율 재는 행. normal 이면 정상 행만 봐서 비율마다 같은 블록이 남음
        ref = df[df.index.astype(str).str.startswith('Normal')] if COLSET_FROM == 'normal' else df
        if len(ref) == 0:
            raise SystemExit(f'COLSET_FROM={COLSET_FROM} 인데 정상 행이 0개입니다 ({label}). 행 이름이 Normal 로 시작해야 합니다')
        if MIN_SAMPLE_COVERAGE >= 1.0:
            keep = ref.notna().all(axis=0)              # dropna(how='any') 와 같음
        else:
            keep = ref.notna().sum(axis=0) >= int(np.ceil(MIN_SAMPLE_COVERAGE * len(ref)))
        df = df.loc[:, keep].copy()                    # 아래 df['Type'] 대입이 원본에 안 가게 복사
        del ref
        # 남은 빈칸은 여기서 안 채움 (08, 09 의 Imputer 가 학습 fold 평균으로 채움)
        if COLSET_FROM == 'normal':
            # 정상은 덮이고 암 행에서 많이 비는 블록 수 (빈칸 자체가 신호가 될 수 있음)
            tum = df[~df.index.astype(str).str.startswith('Normal')]
            thin = int((tum.notna().mean(axis=0) < MIN_SAMPLE_COVERAGE).sum()) if len(tum) else 0
            print(f'    {label}: 암 행 커버율 {MIN_SAMPLE_COVERAGE} 미만인 남은 블록 {thin}개')
    print(f'    {label}: 블록 {before} -> {df.shape[1]} (커버율 기준 {MIN_SAMPLE_COVERAGE} · {COLSET_FROM} 행)')
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
        # 'GBM' 은 행 이름 접두사 (GBM1, GBM2 ...). 07 이 이 이름으로 읽어서 바꾸면 안 됨
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
