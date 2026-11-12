#!/usr/bin/env python
# coding: utf-8
"""
2_train/steps/04: 블록 커버리지 확인 (05 전에 실행)

깊이가 낮으면 샘플 하나만 비어도 그 블록이 빠지니 몇 개 남는지 미리 셈
블록 x 샘플 커버 여부 세고, 요구 커버율별 남는 블록과 샘플 수 늘릴 때 감소 곡선 출력
--limit N   앞 N개 샘플만 확인

출력: {ML_OUT}/coverage_check_{VERSION}.csv       (블록별 커버 샘플 수)
      {ML_OUT}/coverage_curve_{VERSION}.csv       (샘플 수별 전부 커버 블록)
      {ML_OUT}/coverage_thresholds_{VERSION}.csv  (요구 커버율별 남는 블록)
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

# 확인할 샘플 수. None 이면 전부
SAMPLE_LIMIT = None
if '--limit' in sys.argv:
    SAMPLE_LIMIT = int(sys.argv[sys.argv.index('--limit') + 1])
# 정상만 섞은 mut_0 기준으로 봄
FOLDERS = ['mut_0_reads']


def load_dmr_set():
    df = pd.read_csv(C.DMR_CSV)
    chrom = df['Feature'].str.split('_').str.get(0).astype(str)
    blk = df['Feature'].str.split('_').str.get(1).astype(float).astype(int)
    return set(zip(chrom, blk))


def covered_blocks(path, dmr, fmt):
    """이 파일이 커버한 DMR 블록 집합"""
    opener = gzip.open(path, 'rt')
    out = set()
    with opener as f:
        for line in f:
            if line.startswith('track'):      # track 줄 있으면 건너뜀
                continue
            parts = line.rstrip('\n').split('\t')
            if len(parts) < 4:
                continue
            chrom = parts[0]
            pos = int(parts[1]) + (1 if fmt == 'bedgraph' else 0)
            key = (chrom, pos // C.BLOCK_SIZE)
            if key in dmr:
                out.add(key)
    return out


def list_files(folder, fmt):
    import re

    suffix = '.bismark.cov.gz' if fmt == 'cov' else '.bedGraph.gz'
    files = [f for f in os.listdir(folder) if f.endswith(suffix)]

    def num(fn):
        m = re.search(r'\d+', fn)
        return int(m.group()) if m else float('inf')

    return sorted(files, key=num)


def main():
    C.guard_not_original()
    C.banner('2_train/steps/04  블록 커버리지 진단')
    C.ensure_dirs()

    fmt = getattr(C, 'INPUT_FORMAT', 'cov')
    dmr = load_dmr_set()
    print(f'DMR {len(dmr)}개 블록 · 입력 형식 {fmt}\n')

    counts = {}          # 블록 -> 커버한 샘플 수
    per_sample = []      # 샘플별 커버 블록 집합 (감쇠 곡선용)

    for folder_name in FOLDERS:
        folder = f'{C.BISMARK_OUT}/{folder_name}'
        if not os.path.isdir(folder):
            print(f'건너뜀: {folder} 없음')
            continue
        files = list_files(folder, fmt)
        n = min(len(files), SAMPLE_LIMIT or len(files))
        print(f'{folder_name}: {n}개 파일 검사')
        for i, fn in enumerate(files[:n], 1):
            blocks = covered_blocks(os.path.join(folder, fn), dmr, fmt)
            per_sample.append(blocks)
            for b in blocks:
                counts[b] = counts.get(b, 0) + 1
            if i % 100 == 0:
                print(f'  {i}/{n}')

    if not per_sample:
        raise SystemExit('읽은 파일이 없습니다. 2_train/steps/03 출력 경로를 확인하세요.')

    N = len(per_sample)
    # 인덱스로 먼저 정렬 후 안정 정렬 - 동률 블록 순서를 실행마다 같게
    cov = (pd.Series(counts, dtype=int).sort_index()
             .sort_values(ascending=False, kind='mergesort'))
    print(f'\n샘플 {N}개 · 한 번이라도 커버된 DMR 블록 {len(cov)}개 / {len(dmr)}개')
    print(f'샘플당 커버 블록 중간값 {int(np.median([len(s) for s in per_sample]))}개')

    # 요구 커버율별 남는 블록
    print(f"\n{'요구 커버율':>10} {'요구 샘플수':>10} {'잔존 블록':>9}")
    print('-' * 34)
    rows = []
    for frac in (1.0, 0.99, 0.95, 0.90, 0.80, 0.50):
        need = int(np.ceil(frac * N))
        keep = int((cov >= need).sum())
        rows.append({'min_fraction': frac, 'min_samples': need, 'blocks_kept': keep})
        mark = '  <- dropna(how="any") 와 동일' if frac == 1.0 else ''
        print(f'{frac*100:>9.0f}% {need:>10} {keep:>9}{mark}')

    # 감소 곡선: 샘플 수 늘릴 때 전부 커버 블록 수
    print('\n샘플 수를 늘릴 때 "전부 커버" 블록의 변화 (포화 여부 확인)')
    print(f"{'샘플수':>7} {'전부커버 블록':>13} {'직전 대비':>9}")
    print('-' * 33)
    marks = [k for k in (1, 2, 3, 5, 10, 20, 50, 100, 200, 500, 1000, 2000) if k <= N]
    if N not in marks:
        marks.append(N)
    inter = None
    prev = None
    curve = []
    for k in marks:
        inter = set.intersection(*per_sample[:k]) if k > 1 else set(per_sample[0])
        cnt = len(inter)
        delta = '' if prev is None else f'{cnt-prev:+d}'
        print(f'{k:>7} {cnt:>13} {delta:>9}')
        curve.append({'n_samples': k, 'blocks_all_covered': cnt})
        prev = cnt

    # 인덱스가 (염색체, 블록번호) 2단이라 열로 풀어서 저장
    out1 = f'{C.ML_OUT}/coverage_check_{C.VERSION}.csv'
    _cov = cov.rename('n_samples_covered').reset_index()
    _cov.columns = ['chr', 'block', 'n_samples_covered']
    _cov.insert(0, 'Feature', _cov['chr'].astype(str) + '_' + _cov['block'].astype(str))
    _cov.to_csv(out1, index=False)
    out2 = f'{C.ML_OUT}/coverage_curve_{C.VERSION}.csv'
    pd.DataFrame(curve).to_csv(out2, index=False)
    pd.DataFrame(rows).to_csv(f'{C.ML_OUT}/coverage_thresholds_{C.VERSION}.csv', index=False)
    print(f'\n저장 {out1}\n     {out2}')

    # 판단
    keep_any = int((cov >= N).sum())
    print('\n' + '=' * 66)
    if keep_any == 0:
        print('!! dropna(how="any") 로는 블록이 0개 남습니다.')
        print('   config.MIN_SAMPLE_COVERAGE 를 1.0 미만으로 낮추거나 depth 를 올려야 합니다.')
    elif keep_any < 30:
        print(f'!  dropna(how="any") 로 {keep_any}개만 남습니다. feature 가 너무 적습니다.')
        print('   MIN_SAMPLE_COVERAGE 를 0.95 정도로 낮추는 것을 검토하세요.')
    else:
        print(f'   dropna(how="any") 로 {keep_any}개 유지됩니다. 원본 설정을 그대로 써도 됩니다.')
    if len(curve) >= 3 and curve[-1]['blocks_all_covered'] == curve[-2]['blocks_all_covered']:
        print('   감쇠 곡선이 평평해졌습니다. 샘플을 더 늘려도 크게 줄지 않습니다.')
    else:
        print('   감쇠 곡선이 아직 줄고 있습니다. 전체 샘플로 다시 확인하세요.')
    print('=' * 66)
    print('2_train/steps/04 완료.')


if __name__ == '__main__':
    sys.exit(main())
