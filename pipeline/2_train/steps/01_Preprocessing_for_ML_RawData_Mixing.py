#!/usr/bin/env python
# coding: utf-8
"""
2_train/steps/01: DMR 영역 리드만 추출해 가벼운 BAM 을 만든다.

원본 : Methylation/scripts/01_Preprocessing_for_ML_RawData_Mixing.py
변경 : 경로를 config.py 로 옮긴 것뿐. **계산 로직은 원본과 동일**하다.
       - DMR 블록마다 암/정상 리드를 세고
       - 양쪽 모두 MIN_READS_PER_BLOCK(5) 이상인 블록만 통과시키고
       - 통과한 블록의 리드를 새 BAM 에 쓰고 samtools sort + index

출력은 전부 v3 폴더로 간다. 원본 v1/v2 는 건드리지 않는다.
"""

# 2026-09-26: 09-23 사고에서 가드 없는 이 부류가 importlib 로 「확인」 되다 실제로 돌았다.
#   82개 파일이 덮였다(내용은 결정적이라 같았지만 그건 운이었다).
if __name__ != '__main__':
    raise ImportError(__file__ + ' 은(는) 스크립트다 — import 하지 않는다')

import os
import subprocess
import sys

import pandas as pd

import config as C


def load_dmr():
    df = pd.read_csv(C.DMR_CSV)
    g = pd.DataFrame(
        {
            'chr': df['Feature'].str.split('_').str.get(0).astype(str),
            'blk': df['Feature'].str.split('_').str.get(1).astype(float).astype(int),
        }
    ).dropna()
    return g


def main():
    import pysam  # 지연 임포트: 서버에만 있으면 됨

    C.guard_not_original()
    C.banner('2_train/steps/01  DMR 영역 리드 추출')
    C.ensure_dirs()

    if not os.path.exists(C.DMR_CSV):
        raise SystemExit(f'DMR 목록이 없습니다. 1_dmr/steps/08_export_panels.py 를 먼저 실행하세요.\n  {C.DMR_CSV}')
    g = load_dmr()
    print(f'DMR 블록 {len(g)}개')

    n = C.BLOCK_SIZE
    mut_combined = f'{C.MIX_SD}/mut_combined_{C.VERSION}.bam'
    wt_combined = f'{C.MIX_SD}/wt_combined_{C.VERSION}.bam'

    mut_bam = pysam.AlignmentFile(C.FULL_GBM_BAM, 'rb')
    wt_bam = pysam.AlignmentFile(C.FULL_NORMAL_BAM, 'rb')

    block_counts = []
    kept = 0
    with pysam.AlignmentFile(mut_combined, 'wb', header=mut_bam.header) as mo, \
         pysam.AlignmentFile(wt_combined, 'wb', header=wt_bam.header) as wo:
        for _, row in g.iterrows():
            chrom = row['chr']
            start = int(row['blk']) * n
            end = start + n
            try:
                mut_n = sum(1 for _ in mut_bam.fetch(chrom, start, end))
                wt_n = sum(1 for _ in wt_bam.fetch(chrom, start, end))
            except ValueError:
                # 참조 이름이 BAM 에 없는 경우
                continue
            # 원본과 동일한 조건: 양쪽 모두 5개 이상
            if mut_n >= C.MIN_READS_PER_BLOCK and wt_n >= C.MIN_READS_PER_BLOCK:
                for r in mut_bam.fetch(chrom, start, end):
                    mo.write(r)
                for r in wt_bam.fetch(chrom, start, end):
                    wo.write(r)
                block_counts.append(
                    {
                        'chr': chrom,
                        'start_range': start,
                        'end_range': end,
                        'mut_reads': mut_n,
                        'wt_reads': wt_n,
                        'total_reads': mut_n + wt_n,
                    }
                )
                kept += 1

    mut_bam.close()
    wt_bam.close()
    print(f'통과 블록 {kept} / {len(g)}  (양쪽 모두 {C.MIN_READS_PER_BLOCK} reads 이상)')

    bc = pd.DataFrame(block_counts)
    bc_path = f'{C.MIX_SD}/block_read_counts_{C.VERSION}.txt'
    bc.to_csv(bc_path, sep='\t', index=False)
    print(f'블록별 리드 수 저장: {bc_path}')
    if len(bc):
        print(
            f"  정상 블록당 중간값 {bc['wt_reads'].median():.0f} reads / "
            f"암 {bc['mut_reads'].median():.0f} reads"
        )

    print('\nsamtools sort + index...')
    # ---- 2026-09-14: 검증 음성 배경을 **공여자별로** 따로 뽑는다 ----
    #   왜: 음성 반복본 10개가 병합 BAM 하나에서 나왔다. Normal1~10 이 사람 10명이
    #       아니라 같은 통에서 10번 뜬 것이라 **음성 쪽 사람 간 변동이 0** 이었고,
    #       특이도가 실제보다 좋게 보였다.
    #   블록 집합은 위에서 병합본으로 이미 정했다. 여기서는 그 블록만 배경별로 떠낸다.
    #   그래야 배경마다 칸이 달라지지 않는다.
    #   C.NORMAL_BAMS 가 없으면(=학습) 통째로 건너뛴다.
    _bgs = getattr(C, 'NORMAL_BAMS', None)
    if _bgs:
        print('')
        print('배경별 정상 추출: %d개' % len(_bgs))
        _blocks = [(r['chr'], int(r['blk']) * n) for _, r in g.iterrows()]
        for _k, _bp in enumerate(_bgs, 1):
            _out = C.WT_SORTED.replace('.bam', '_bg%02d.bam' % _k)
            if os.path.exists(_out + '.bai'):
                print('  bg%02d 이미 있음' % _k)
                continue
            _tmp = _out + '.tmp.bam'
            _b = pysam.AlignmentFile(_bp, 'rb')
            _cnt = 0
            with pysam.AlignmentFile(_tmp, 'wb', header=_b.header) as _o:
                for _c, _s in _blocks:
                    try:
                        for _r in _b.fetch(_c, _s, _s + n):
                            _o.write(_r)
                            _cnt += 1
                    except ValueError:
                        continue
            _b.close()
            subprocess.run(['samtools', 'sort', '-o', _out, _tmp], check=True)
            subprocess.run(['samtools', 'index', _out], check=True)
            os.remove(_tmp)
            print('  bg%02d  리드 %s  %s' % (_k, format(_cnt, ','), os.path.basename(_bp)))

    for combined, sorted_path in ((mut_combined, C.MUT_SORTED), (wt_combined, C.WT_SORTED)):
        subprocess.run(['samtools', 'sort', '-o', sorted_path, combined], check=True)
        subprocess.run(['samtools', 'index', sorted_path], check=True)
        print(f'  {sorted_path}')

    # 2_train/steps/02 가 쓰는 쿼리 이름 목록: 짝이 정확히 2개인 것만 남긴다
    print('\n리드쌍 이름 목록 생성 (짝 2개인 것만)...')
    for bam_path, out_txt, label in (
        (C.MUT_SORTED, C.MUT_NAMES, '암'),
        (C.WT_SORTED, C.WT_NAMES, '정상'),
    ):
        counts = {}
        with pysam.AlignmentFile(bam_path, 'rb') as b:
            for r in b.fetch(until_eof=True):
                counts[r.query_name] = counts.get(r.query_name, 0) + 1
        paired = [k for k, v in counts.items() if v in (1, 2)]
        with open(out_txt, 'w') as f:
            f.write('\n'.join(paired) + '\n')
        dropped = len(counts) - len(paired)
        print(
            f'  {label}: {len(paired):,}쌍  (제외 {dropped:,}개, '
            f'{dropped/max(1,len(counts))*100:.1f}% — 짝이 2개가 아님)'
        )

    print('\n2_train/steps/01 완료.')


if __name__ == '__main__':
    sys.exit(main())
