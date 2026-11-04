#!/usr/bin/env python
# coding: utf-8
"""
2_train/steps/01: DMR 블록 리드만 뽑아 작은 BAM 만듦

DMR 블록마다 암, 정상 리드 세고 양쪽 모두 MIN_READS_PER_BLOCK 이상인 블록만 남김
남은 블록 리드를 새 BAM 에 쓰고 samtools sort + index
"""

# 스크립트 전용 - import 하면 본체가 돌아 자료를 덮어씀
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
    import pysam  # 서버에만 있으면 돼서 여기서 import

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
    seen_mut, seen_wt = set(), set()
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
                # BAM 에 없는 염색체
                continue
            # 암, 정상 모두 MIN_READS_PER_BLOCK 이상
            if mut_n >= C.MIN_READS_PER_BLOCK and wt_n >= C.MIN_READS_PER_BLOCK:
                # 두 블록에 걸친 리드는 한 번만 씀 (두 번 쓰면 아래 이름당 레코드 수 거르기에서 빠짐)
                for r in mut_bam.fetch(chrom, start, end):
                    key = (r.query_name, r.flag, r.reference_id, r.reference_start)
                    if key not in seen_mut:
                        seen_mut.add(key)
                        mo.write(r)
                for r in wt_bam.fetch(chrom, start, end):
                    key = (r.query_name, r.flag, r.reference_id, r.reference_start)
                    if key not in seen_wt:
                        seen_wt.add(key)
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
    # 검증용 정상 배경은 사람별로 따로 뽑음
    # 병합 BAM 하나에서 뽑으면 사람 간 차이가 없어져 특이도가 좋게 나옴
    # 블록은 위에서 정한 것 그대로 사용. NORMAL_BAMS 없으면(학습) 건너뜀
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
            _seen = set()          # 두 블록에 걸친 리드는 한 번만 씀
            with pysam.AlignmentFile(_tmp, 'wb', header=_b.header) as _o:
                for _c, _s in _blocks:
                    try:
                        for _r in _b.fetch(_c, _s, _s + n):
                            _key = (_r.query_name, _r.flag, _r.reference_id, _r.reference_start)
                            if _key in _seen:
                                continue
                            _seen.add(_key)
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

    # 02 가 쓰는 리드 이름 목록 (이름당 레코드 1~2개인 것만)
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
