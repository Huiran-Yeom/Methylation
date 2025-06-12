#!/usr/bin/env python
# coding: utf-8

# In[ ]:


import os
import random
import pysam
import shutil
import tarfile
from collections import defaultdict

# 랜덤 시드 설정 (재현성 확보용)
random.seed(42)

# 설정값
target_depth = 5000
mut_ratios = [0, 0.001, 1, 0.005, 0.01, 0.02, 0.025, 0.03, 0.05, 0.1]
num_files = 1000
BLOCK_SIZE = 100

mut_bam_path = "/ssd_data/Methylation/sample_data/step03_Preprocessing_for_ML/mut_sorted_v1.bam"
wt_bam_path = "/ssd_data/Methylation/sample_data/step03_Preprocessing_for_ML/wt_sorted_v1.bam"
base_output_dir = "/ssd_data/Methylation/results/step03_Preprocessing_for_ML/v1/"
mut_name_file = "/ssd_data/Methylation/sample_data/step03_Preprocessing_for_ML/mut_names_v1.txt"
wt_name_file = "/ssd_data/Methylation/sample_data/step03_Preprocessing_for_ML/wt_names_v1.txt"

# 참조 길이 및 읽기 길이 계산
with pysam.AlignmentFile(mut_bam_path, "rb") as mut_bam:
    reference_lengths = dict(zip(mut_bam.references, mut_bam.lengths))
    sample_read = next(mut_bam.fetch(until_eof=True))
    read_length = sample_read.query_length

# 블록 구하기
def get_actual_blocks(bam_path):
    blocks = set()
    with pysam.AlignmentFile(bam_path, "rb") as bam:
        for read in bam.fetch(until_eof=True):
            if read.is_unmapped:
                continue
            chrom = read.reference_name
            start = (read.reference_start // BLOCK_SIZE) * BLOCK_SIZE
            blocks.add((chrom, start))
    return blocks

mut_blocks = get_actual_blocks(mut_bam_path)
wt_blocks = get_actual_blocks(wt_bam_path)
all_blocks = mut_blocks.union(wt_blocks)

total_blocks = len(all_blocks)
baseline_total_pairs = total_blocks * target_depth // 2  # 페어 기준

print(f"실제 사용 블록 수: {total_blocks}")
print(f"총 필요한 read 쌍 수: {baseline_total_pairs}")

# 쿼리 이름 저장
def save_query_names(bam_path, output_txt):
    with pysam.AlignmentFile(bam_path, "rb") as bam, open(output_txt, "w") as f:
        names = set()
        for read in bam.fetch(until_eof=True):
            if read.is_paired and read.query_name not in names:
                f.write(read.query_name + "\n")
                names.add(read.query_name)

save_query_names(mut_bam_path, mut_name_file)
save_query_names(wt_bam_path, wt_name_file)

# 쿼리 이름 로딩
with open(mut_name_file) as f:
    mut_names = [line.strip() for line in f]

with open(wt_name_file) as f:
    wt_names = [line.strip() for line in f]

# read 객체 메모리 로딩
def load_reads_by_name(bam_path):
    reads = defaultdict(list)
    with pysam.AlignmentFile(bam_path, "rb") as bam:
        for read in bam.fetch(until_eof=True):
            reads[read.query_name].append(read)
    return reads

print("read 메모리 로딩 중...")
mut_reads_dict = load_reads_by_name(mut_bam_path)
wt_reads_dict = load_reads_by_name(wt_bam_path)
print("메모리 로딩 완료.")

# 블록 카운트 함수
def count_blocks(read, label, block_counts):
    if read.is_unmapped:
        return
    chrom = read.reference_name
    start = (read.reference_start // BLOCK_SIZE) * BLOCK_SIZE
    block_counts[(chrom, start)][label] += 1

# 선택된 read 쌍 쓰기
def write_selected_reads(reads_dict, names_set, outbam, label, block_counts):
    for name in names_set:
        read_pair = reads_dict.get(name, [])
        if len(read_pair) != 2:
            continue
        for read in read_pair:
            outbam.write(read)
            count_blocks(read, label, block_counts)

# 블록 필터링
def filter_blocks_by_min_reads(block_counts, min_reads=5):
    return {block: counts for block, counts in block_counts.items() if (counts['mut'] + counts['wt']) >= min_reads}

# 샘플 BAM 생성
def create_sampled_bam(mut_ratio):
    folder_name = f"mut_{int(mut_ratio * 5000)}_reads"
    output_dir = os.path.join(base_output_dir, folder_name)
    os.makedirs(output_dir, exist_ok=True)
    txt_dir = os.path.join(output_dir, "txt")
    os.makedirs(txt_dir, exist_ok=True)
    mut_pair_count = int(baseline_total_pairs * mut_ratio)
    wt_pair_count = baseline_total_pairs - mut_pair_count
    print(f"Mut ratio: {mut_ratio:.3f} -> Mut pairs: {mut_pair_count}, WT pairs: {wt_pair_count}")
    for i in range(num_files):
        sampled_mut_set = (random.sample(mut_names, mut_pair_count)
                           if mut_pair_count <= len(mut_names)
                           else random.choices(mut_names, k=mut_pair_count))
        sampled_wt_set = (random.sample(wt_names, wt_pair_count)
                          if wt_pair_count <= len(wt_names)
                          else random.choices(wt_names, k=wt_pair_count))
        output_bam_path = os.path.join(output_dir, f"sampled_reads_{i+1}.bam")
        header = {
            'HD': {'VN': '1.0', 'SO': 'coordinate'},
            'SQ': [{'SN': k, 'LN': v} for k, v in reference_lengths.items()]
        }
        block_counts = defaultdict(lambda: {"mut": 0, "wt": 0})
        with pysam.AlignmentFile(output_bam_path, "wb", header=header) as outbam:
            write_selected_reads(wt_reads_dict, sampled_wt_set, outbam, "wt", block_counts)
            write_selected_reads(mut_reads_dict, sampled_mut_set, outbam, "mut", block_counts)
        # 블록 통계 저장
        filtered_block_counts = filter_blocks_by_min_reads(block_counts)
        stats_path = os.path.join(txt_dir, f"block_read_counts_{i+1}.txt")
        with open(stats_path, "w") as f:
            f.write("chr\tstart\tend\tmutant_count\twt_count\n")
            for (chrom, start) in sorted(filtered_block_counts.keys()):
                end = start + BLOCK_SIZE
                counts = filtered_block_counts[(chrom, start)]
                f.write(f"{chrom}\t{start}\t{end}\t{counts['mut']}\t{counts['wt']}\n")
        if (i + 1) % 10 == 0 or i == num_files - 1:
            print(f"[{mut_ratio:.3%}] 샘플 {i+1} 생성 완료: {output_bam_path}")
    # txt 폴더 압축 및 삭제
    if os.path.exists(txt_dir):
        archive_path = os.path.join(output_dir, "txt.tar.gz")
        with tarfile.open(archive_path, "w:gz") as tar:
            tar.add(txt_dir, arcname="txt")
        shutil.rmtree(txt_dir)
        print(f"[{mut_ratio:.3%}] txt 폴더 압축 완료 및 삭제됨: {archive_path}")

# 실행
if __name__ == '__main__':
    for mut_ratio in mut_ratios:
        create_sampled_bam(mut_ratio)



# In[ ]:




