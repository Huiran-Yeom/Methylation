#!/usr/bin/env python
# coding: utf-8

# ### 원하는 비율 조합 파일 생성 
# - total reads: 5000 reads
# - GBM percent: 0.1%, 0.5%, 1%, 2%, 2.5%, 3%, 5%, 10%, 100%

# In[ ]:


import os
import subprocess

mut_ratios = [0, 5, 5000, 25, 50, 100, 125, 150, 250, 500]  # mutation 비율 리스트
num_files = 1000  # 각 비율당 샘플 개수
OUTPUT_BASE = "/ssd_data/Methylation/sample_data/step04_ML_classifier/v1/"  # 출력 베이스 경로
GENOME_FOLDER = "./data/Bisulfite_Genome/"  # bismark genome preparation 디렉토리

def process_bam(mut_reads, i):
    input_bam = f"/ssd_data/Methylation/results/step03_Preprocessing_for_ML/v1/mut_{mut_reads}_reads/sampled_reads_{i}.bam"
    output_dir = f"{OUTPUT_BASE}/mut_{mut_reads}_reads"
    if not os.path.exists(input_bam):
        print(f"❌ {input_bam} does not exist. Skipping.")
        return
    os.makedirs(output_dir, exist_ok=True)
    # Bismark methylation extractor 명령어
    bismark_cmd = [
        "bismark_methylation_extractor", 
        "--paired-end", 
        "--ignore_3prime", "1",
        "--genome_folder", GENOME_FOLDER,
        "--bedGraph", 
        "--gzip", 
        "--output_dir", output_dir, 
        input_bam
    ]
    print(f"🔄 Running: {' '.join(bismark_cmd)}")
    try:
        subprocess.run(bismark_cmd, check=True)
        print(f"✅ Finished mut_{mut_reads}_reads, sample {i}")
    except subprocess.CalledProcessError as e:
        print(f"❌ Error in mut_{mut_reads}_reads, sample {i}: {e}")

# 실행 루프
for mut_reads in mut_ratios:
    for i in range(1, num_files + 1):
        process_bam(mut_reads, i)



# In[ ]:




