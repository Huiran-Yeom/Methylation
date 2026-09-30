#!/usr/bin/env python
# coding: utf-8
"""
2_train/steps/04: 블록 커버리지 진단.  **step04_1 을 돌리기 전에 반드시 먼저 실행할 것**

왜 필요한가
  원본 step04 는 feature 행렬을 만든 뒤 `dropna(axis=0, how='any')` 로 결측 블록을
  버린다. 즉 **샘플 하나라도 커버가 없으면 그 블록은 전체에서 탈락**한다.
  원본은 depth 가 658배였기 때문에 모든 블록이 모든 샘플에서 커버돼 문제가 없었다.

  📊 실측 (옛 작업 기계의 5_base_Normal.csv: 총 5,000 reads, 5샘플)
      ※ 그 파일은 이 저장소에 없다. 숫자만 인용한 것이다.
      자리 15,418개 · 블록 1,139개 · 샘플당 커버 자리 약 76%
      5개 샘플 전부 커버한 블록 849개 = 74.5%   (샘플-블록 커버 확률 q≈0.885)
      요구 샘플 수를 늘리면  ≥1:1139  ≥2:1078  ≥3:1021  ≥4:955  ≥5:849

  샘플이 2,000개면 `how='any'` 로 살아남는 블록이 몇 개일지 5샘플 데이터로는 알 수 없다.
  블록마다 커버리지 성향이 다르기 때문이다. **그래서 직접 세어봐야 한다.**

이 스크립트가 하는 일
  1. 커버리지 파일을 읽어 블록 x 샘플 커버 여부를 센다
  2. "N개 샘플 이상에서 커버된 블록" 곡선을 출력한다
  3. MIN_SAMPLE_COVERAGE 후보값마다 남는 블록 수를 보여준다
  4. 샘플 수를 늘려가며 잔존 블록이 어떻게 줄어드는지(감쇠 곡선) 측정한다
     → 실제로 포화하는지, 계속 줄어드는지가 여기서 드러난다

출력: {ML_OUT}/coverage_check_{VERSION}.csv  (블록별 커버 샘플 수)
      {ML_OUT}/coverage_curve_{VERSION}.csv  (요구 샘플 수별 잔존 블록)
"""

# 2026-09-26: 09-23 사고에서 가드 없는 이 부류가 importlib 로 「확인」 되다 실제로 돌았다.
#   82개 파일이 덮였다(내용은 결정적이라 같았지만 그건 운이었다).
if __name__ != '__main__':
    raise ImportError(__file__ + ' 은(는) 스크립트다 — import 하지 않는다')

import gzip
import os
import sys

import numpy as np
import pandas as pd

import config as C

# 진단용 표본 수. None 이면 config.NUM_FILES 전부. 먼저 100~200 으로 감을 잡는 것을 권장.
SAMPLE_LIMIT = None
if '--limit' in sys.argv:
    SAMPLE_LIMIT = int(sys.argv[sys.argv.index('--limit') + 1])
# 어느 조건을 볼지. 정상 대조군이 기준이므로 mut_0 이 기본.
FOLDERS = ['mut_0_reads']


def load_dmr_set():
    df = pd.read_csv(C.DMR_CSV)
    chrom = df['Feature'].str.split('_').str.get(0).astype(str)
    blk = df['Feature'].str.split('_').str.get(1).astype(float).astype(int)
    return set(zip(chrom, blk))


def covered_blocks(path, dmr, fmt):
    """이 파일이 커버한 DMR 블록 집합."""
    opener = gzip.open(path, 'rt')
    out = set()
    with opener as f:
        if fmt == 'bedgraph':
            next(f, None)  # track 헤더
        for line in f:
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
    # 2026-09-29: 동률 블록의 순서가 glob 순서에 달려 파일 md5 가 실행마다 달랐다
    #   (내용은 같고 행 순서만 다름). 인덱스로 먼저 정렬한 뒤 안정 정렬하면 결정적이다.
    cov = (pd.Series(counts, dtype=int).sort_index()
             .sort_values(ascending=False, kind='mergesort'))
    print(f'\n샘플 {N}개 · 한 번이라도 커버된 DMR 블록 {len(cov)}개 / {len(dmr)}개')
    print(f'샘플당 커버 블록 중간값 {int(np.median([len(s) for s in per_sample]))}개')

    # --- 요구 커버율별 잔존 블록 ---
    print(f"\n{'요구 커버율':>10} {'요구 샘플수':>10} {'잔존 블록':>9}")
    print('-' * 34)
    rows = []
    for frac in (1.0, 0.99, 0.95, 0.90, 0.80, 0.50):
        need = int(np.ceil(frac * N))
        keep = int((cov >= need).sum())
        rows.append({'min_fraction': frac, 'min_samples': need, 'blocks_kept': keep})
        mark = '  <- dropna(how="any") 와 동일' if frac == 1.0 else ''
        print(f'{frac*100:>9.0f}% {need:>10} {keep:>9}{mark}')

    # --- 감쇠 곡선: 샘플 수를 늘려가며 '전부 커버' 블록이 어떻게 줄어드나 ---
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

    # 인덱스가 (염색체, 블록번호) 2단이므로 index_label 도 두 개여야 한다.
    # 하나만 주면 헤더 2개 / 데이터 3개가 되어 읽을 때 염색체가 인덱스로 먹힌다.
    out1 = f'{C.ML_OUT}/coverage_check_{C.VERSION}.csv'
    _cov = cov.rename('n_samples_covered').reset_index()
    _cov.columns = ['chr', 'block', 'n_samples_covered']
    _cov.insert(0, 'Feature', _cov['chr'].astype(str) + '_' + _cov['block'].astype(str))
    _cov.to_csv(out1, index=False)
    out2 = f'{C.ML_OUT}/coverage_curve_{C.VERSION}.csv'
    pd.DataFrame(curve).to_csv(out2, index=False)
    pd.DataFrame(rows).to_csv(f'{C.ML_OUT}/coverage_thresholds_{C.VERSION}.csv', index=False)
    print(f'\n저장 {out1}\n     {out2}')

    # --- 판단 ---
    keep_any = int((cov >= N).sum())
    print('\n' + '=' * 66)
    if keep_any == 0:
        print('🔴 dropna(how="any") 로는 블록이 0개 남습니다.')
        print('   config.MIN_SAMPLE_COVERAGE 를 1.0 미만으로 낮추거나 depth 를 올려야 합니다.')
    elif keep_any < 30:
        print(f'🟡 dropna(how="any") 로 {keep_any}개만 남습니다. feature 가 너무 적습니다.')
        print('   MIN_SAMPLE_COVERAGE 를 0.95 정도로 낮추는 것을 검토하세요.')
    else:
        print(f'🟢 dropna(how="any") 로 {keep_any}개 유지됩니다. 원본 설정을 그대로 써도 됩니다.')
    if len(curve) >= 3 and curve[-1]['blocks_all_covered'] == curve[-2]['blocks_all_covered']:
        print('   감쇠 곡선이 평평해졌습니다. 샘플을 더 늘려도 크게 줄지 않습니다.')
    else:
        print('   감쇠 곡선이 아직 줄고 있습니다. 전체 샘플로 다시 확인하세요.')
    print('=' * 66)
    print('2_train/steps/04 완료.')


if __name__ == '__main__':
    sys.exit(main())
