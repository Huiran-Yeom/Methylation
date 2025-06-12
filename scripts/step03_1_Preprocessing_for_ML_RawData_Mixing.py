#!/usr/bin/env python
# coding: utf-8

# In[ ]:


import pysam
import pandas as pd
import subprocess

# CSV 파일 읽기
n = 100
th = 4
df2 = pd.read_csv(rf'/ssd_data/Methylation/results/step02_DMR/{n}기준_{th}이상_feature_importances(inner)_v1.csv')

# Feature에서 chr와 start 정보 추출
g_list = df2[['Feature']].copy()
g_list['chr'] = df2['Feature'].str.split('_').str.get(0)
g_list[f'start_div_{n}'] = df2['Feature'].str.split('_').str.get(1)
g_list = g_list[['chr', f'start_div_{n}']].dropna()
g_list['chr'] = g_list['chr'].astype(str)
g_list[f'start_div_{n}'] = g_list[f'start_div_{n}'].astype(int)

# BAM 파일 경로
mut_bam_file_path = "/ssd_data/Methylation/sample_data/step03_Preprocessing_for_ML/merged_glio_cell_all_sorted.bam"
wt_bam_file_path = "/ssd_data/Methylation/sample_data/step03_Preprocessing_for_ML/merged_normal_cfDNA_v2_sorted.bam"

# BAM 파일 열기
mut_bam = pysam.AlignmentFile(mut_bam_file_path, "rb")
wt_bam = pysam.AlignmentFile(wt_bam_file_path, "rb")

# 출력 파일 경로
mut_combined_output_path = "/ssd_data/Methylation/sample_data/step03_Preprocessing_for_ML/mut_combined_v1.bam"
wt_combined_output_path = "/ssd_data/Methylation/sample_data/step03_Preprocessing_for_ML/wt_combined_v1.bam"
mut_sorted_path = "/ssd_data/Methylation/sample_data/step03_Preprocessing_for_ML/mut_sorted_v1.bam"
wt_sorted_path = "/ssd_data/Methylation/sample_data/step03_Preprocessing_for_ML/wt_sorted_v1.bam"
mut_index_path = mut_sorted_path + ".bai"
wt_index_path = wt_sorted_path + ".bai"
block_read_count_output_path = "/ssd_data/Methylation/sample_data/step03_Preprocessing_for_ML/block_read_counts_v1.txt"
mut_block_read_count_output_path = "/ssd_data/Methylation/sample_data/step03_Preprocessing_for_ML/mut_block_read_counts_v1.txt"
wt_block_read_count_output_path = "/ssd_data/Methylation/sample_data/step03_Preprocessing_for_ML/wt_block_read_counts_v1.txt"

# 블록별 리드 수를 저장할 리스트
block_read_counts = []
mut_block_read_counts = []
wt_block_read_counts = []

# Mut와 WT 출력 파일 생성 (블록에 포함되는 리드 추출)
with pysam.AlignmentFile(mut_combined_output_path, "wb", header=mut_bam.header) as mut_out_bam, \
     pysam.AlignmentFile(wt_combined_output_path, "wb", header=wt_bam.header) as wt_out_bam:
    for idx, row in g_list.iterrows():
        chr_name = row['chr']
        start_range = int(row[f'start_div_{n}']) * n
        end_range = (int(row[f'start_div_{n}']) + 1) * n
        print(f"Processing: {chr_name}, {start_range}, {end_range}")  # 디버깅용 출력
        mut_count = sum(1 for _ in mut_bam.fetch(chr_name, start_range, end_range))
        wt_count = sum(1 for _ in wt_bam.fetch(chr_name, start_range, end_range))
        # Mut와 WT 리드가 모두 5개 이상인 경우만 저장
        if mut_count >= 5 and wt_count >= 5:
            mut_block_read_counts.append({
                'chr': chr_name,
                'start_range': start_range,
                'end_range': end_range,
                'mut_reads': mut_count
            })
            wt_block_read_counts.append({
                'chr': chr_name,
                'start_range': start_range,
                'end_range': end_range,
                'wt_reads': wt_count
            })
            for read in mut_bam.fetch(chr_name, start_range, end_range):
                mut_out_bam.write(read)
            for read in wt_bam.fetch(chr_name, start_range, end_range):
                wt_out_bam.write(read)
            block_read_counts.append({
                'chr': chr_name,
                'start_range': start_range,
                'end_range': end_range,
                'mut_reads': mut_count,
                'wt_reads': wt_count,
                'total_reads': mut_count + wt_count
            })
            print(f"Processed block {chr_name}:{start_range}-{end_range} - Mut: {mut_count}, WT: {wt_count}")

# BAM 파일 닫기
mut_bam.close()
wt_bam.close()

# 블록별 리드 수를 텍스트 파일로 저장
with open(block_read_count_output_path, "w") as count_file:
    count_file.write("chr\tstart_range\tend_range\tmut_reads\twt_reads\ttotal_reads\n")
    for block in block_read_counts:
        count_file.write(
            f"{block['chr']}\t{block['start_range']}\t{block['end_range']}\t"
            f"{block['mut_reads']}\t{block['wt_reads']}\t{block['total_reads']}\n"
        )

# Mut와 WT 블록별 리드 수를 별도로 저장
with open(mut_block_read_count_output_path, "w") as mut_count_file:
    mut_count_file.write("chr\tstart_range\tend_range\tmut_reads\n")
    for block in mut_block_read_counts:
        mut_count_file.write(
            f"{block['chr']}\t{block['start_range']}\t{block['end_range']}\t{block['mut_reads']}\n"
        )

with open(wt_block_read_count_output_path, "w") as wt_count_file:
    wt_count_file.write("chr\tstart_range\tend_range\twt_reads\n")
    for block in wt_block_read_counts:
        wt_count_file.write(
            f"{block['chr']}\t{block['start_range']}\t{block['end_range']}\t{block['wt_reads']}\n"
        )

print("Read count files saved successfully.")
print("Starting BAM sort and index...")

# Samtools sort + index (subprocess 사용)
subprocess.run(["samtools", "sort", "-o", mut_sorted_path, mut_combined_output_path], check=True)
subprocess.run(["samtools", "index", mut_sorted_path], check=True)
subprocess.run(["samtools", "sort", "-o", wt_sorted_path, wt_combined_output_path], check=True)
subprocess.run(["samtools", "index", wt_sorted_path], check=True)

print(f"Mut reads saved and indexed: {mut_sorted_path}")
print(f"WT reads saved and indexed: {wt_sorted_path}")


