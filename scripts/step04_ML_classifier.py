#!/usr/bin/env python
# coding: utf-8

# ### 0) bed 파일 전처리

# In[ ]:


import pandas as pd
import os

type = '.bedGraph.gz'
n = 100

def filename(dir, bedfile, type):
    for filename in os.listdir(dir):
        if filename.endswith(type):
            file_path = os.path.join(dir, filename)
            file_size = os.path.getsize(file_path)
            if file_size > 0:
                bedfile.append((filename, file_size))
    return bedfile

##### 정상 샘플 #####
directory1 = f'/ssd_data/Methylation/sample_data/step04_ML_classifier/v1/mut_0_reads/'
bed_files1 = []
bed_files1 = filename(directory1, bed_files1, type)
bed_files1 = list(set(bed_files1))
bed_files1 = [file[0] for file in sorted(bed_files1, key=lambda x: x[1], reverse=True)]
print(len(bed_files1))
df_1 = pd.read_csv(rf'{directory1}{bed_files1[0]}', sep='\t', header=None, names=['chr', 'start', 'end', 'Normal1'], dtype={'chr': str})
df_1 = df_1[1:]
df_1 = df_1[df_1.columns[:4]]
df_1['Normal1'] = df_1['Normal1']/100
grouped1_2 = df_1

for i in range(2, len(bed_files1) + 1):
    df1_1 = pd.read_csv(rf'{directory1}{bed_files1[i - 1]}', sep='\t', header=None, 
            names=['chr', 'start', 'end', f'Normal{i}'], dtype={'chr': str})
    df1_1 = df1_1[1:]
    df1_1 = df1_1[df1_1.columns[:4]]
    df1_1[f'Normal{i}'] = df1_1[f'Normal{i}']/100
    grouped1_2 = pd.merge(grouped1_2, df1_1, on=['chr', f'start', 'end'], how='outer')


print(grouped1_2)
grouped1_2.to_csv(f'/ssd_data/Methylation/results/step04_ML_classifier/v1/base_Normal_v1.csv', index=False)


# In[ ]:


##### 암 샘플 #####
p_list = [5, 5000, 25, 50, 100, 125, 150, 250, 500]  # 원하는 p 값 리스트

for p in p_list:
    directory2 = f'/ssd_data/Methylation/sample_data/step04_ML_classifier/v1/mut_{p}_reads/'
    bed_files2 = []
    bed_files2 = filename(directory2, bed_files2, type)
    bed_files2 = list(set(bed_files2))
    bed_files2 = [file[0] for file in sorted(bed_files2, key=lambda x: x[1], reverse=True)]
    print(len(bed_files2))
    df_2 = pd.read_csv(rf'{directory2}{bed_files2[0]}', sep='\t', header=None, names=['chr', 'start', 'end', 'GBM1'], dtype={'chr': str})
    df_2 = df_2[1:]
    df_2 = df_2[df_2.columns[:4]]
    df_2['GBM1'] = df_2['GBM1']/100
    grouped2_2 = df_2
    for i in range(2, len(bed_files2) + 1):
        df2_1 = pd.read_csv(rf'{directory2}{bed_files2[i - 1]}', sep='\t', header=None, 
                            names=['chr', 'start', 'end', f'GBM{i}'], dtype={'chr': str})
        df2_1 = df2_1[1:]
        df2_1 = df2_1[df2_1.columns[:4]]
        df2_1[f'GBM{i}'] = df2_1[f'GBM{i}']/100
        grouped2_2 = pd.merge(grouped2_2, df2_1, on=['chr', f'start', 'end'], how='outer')
    print(grouped2_2)
    grouped2_2.to_csv(f'/ssd_data/Methylation/results/step04_ML_classifier/v1/{p}_base_GBM_v1.csv', index=False)



# ### 1) Mean

# In[ ]:


import pandas as pd

th = 4
df2 = pd.read_csv(rf'/ssd_data/Methylation/results/step02_DMR/{n}기준_{th}이상_feature_importances(inner)_v1.csv')

# Feature에서 chr와 start 정보 추출
g_list = df2[['Feature']].copy()
g_list['chr'] = df2['Feature'].str.split('_').str.get(0)
g_list[f'start_div_{n}'] = df2['Feature'].str.split('_').str.get(1)
g_list = g_list[['chr', f'start_div_{n}']].dropna()
g_list['chr'] = g_list['chr'].astype(str)
g_list[f'start_div_{n}'] = g_list[f'start_div_{n}'].astype(int)

for p in p_list:
    ##### 정상 샘플 #####
    grouped1_2 = pd.read_csv(f'/ssd_data/Methylation/results/step04_ML_classifier/v1/base_Normal_v1.csv')
    grouped1_2[f'start_div_{n}']=grouped1_2['start']//n
    grouped1_2 = grouped1_2.merge(g_list, on=['chr',f'start_div_{n}'], how='inner')
    # 'chr'과 'start_div_100'을 기준으로 그룹화하여 'Normal1'부터 'Normal1000'까지의 평균을 계산
    normal_cols = [col for col in grouped1_2.columns if 'Normal' in col]
    mean1 = grouped1_2.groupby(['chr', f'start_div_{n}'])[normal_cols].mean().reset_index()
    # print(mean1.head())
    ##### 암 샘플 #####
    grouped2_2 = pd.read_csv(f'/ssd_data/Methylation/results/step04_ML_classifier/v1/{p}_base_GBM_v1.csv')
    grouped2_2[f'start_div_{n}']=grouped2_2['start']//n
    grouped2_2 = grouped2_2.merge(g_list, on=['chr',f'start_div_{n}'], how='inner')
    # 'chr'과 'start_div_100'을 기준으로 그룹화하여 'GBM1'부터 'GBM1000'까지의 평균을 계산
    gbm_cols = [col for col in grouped2_2.columns if 'GBM' in col]
    mean2 = grouped2_2.groupby(['chr', f'start_div_{n}'])[gbm_cols].mean().reset_index()
    # print(mean2.head())
    ##### 병합 #####
    grouped_mean = pd.merge(mean1, mean2, on=['chr', f'start_div_{n}'], how='inner')
    grouped_mean = grouped_mean.dropna(axis=0, how='any')
    # print(grouped_mean.head())
    df_mean = grouped_mean
    df_mean['chr_start'] = df_mean['chr'].astype(str) + '_' + df_mean[f'start_div_{n}'].astype(str)
    df_mean = df_mean.drop(['chr',f'start_div_{n}'], axis=1)
    df_mean = df_mean.set_index('chr_start').transpose()
    df_mean['Type'] = df_mean.index.map(lambda x: 0 if 'Normal' in str(x) else (1 if 'GBM' in str(x) else np.nan))
    print(df_mean.head())
    df_mean.to_csv(f'/ssd_data/Methylation/results/step04_ML_classifier/v1/{p}_mean_v1.csv', index=False)



# ### 2) Entropy

# In[ ]:


import pandas as pd
import numpy as np

# 엔트로피 계산 함수
def calculate_entropy(series):
    series = series.dropna()
    if series.empty:
        return np.nan
    value_counts = series.value_counts(normalize=True)
    entropy = -np.sum(value_counts * np.log2(value_counts))
    return entropy


for p in p_list:
    ##### 정상 샘플 #####
    grouped1_2 = pd.read_csv(f'/ssd_data/Methylation/results/step04_ML_classifier/v1/base_Normal_v1.csv')
    grouped1_2[f'start_div_{n}']=grouped1_2['start']//n
    grouped1_2 = grouped1_2.merge(g_list, on=['chr',f'start_div_{n}'], how='inner')
    grouped1 = grouped1_2.groupby(['chr', f'start_div_{n}'])
    entropy_resultss1 = {}
    for i in range(1, len(bed_files1) + 1):
        col_name = f'Normal{i}'
        entropy_resultss1[col_name] = grouped1[col_name].apply(lambda x: calculate_entropy(x))
    entropy_df1 = pd.DataFrame(entropy_resultss1).reset_index()
    # print(entropy_df1.head())
    ##### 암 샘플 #####
    grouped2_2 = pd.read_csv(f'/ssd_data/Methylation/results/step04_ML_classifier/v1/{p}_base_GBM_v1.csv')
    grouped2_2[f'start_div_{n}']=grouped2_2['start']//n
    grouped2_2 = grouped2_2.merge(g_list, on=['chr',f'start_div_{n}'], how='inner')
    grouped2 = grouped2_2.groupby(['chr', f'start_div_{n}'])
    entropy_resultss2 = {}
    for i in range(1, len(bed_files2) + 1):
        col_name = f'GBM{i}'
        entropy_resultss2[col_name] = grouped2[col_name].apply(lambda x: calculate_entropy(x))
    entropy_df2 = pd.DataFrame(entropy_resultss2).reset_index()
    # print(entropy_df2.head())
    ##### 병합 #####
    grouped_entropy = pd.merge(entropy_df1, entropy_df2, on=['chr', f'start_div_{n}'], how='inner')
    grouped_entropy  = grouped_entropy.dropna(axis=0, how='any')
    # print(grouped_entropy.head())
    df_entropy = grouped_entropy 
    df_entropy['chr_start'] = df_entropy['chr'].astype(str) + '_' + df_entropy[f'start_div_{n}'].astype(str)
    df_entropy = df_entropy.drop(['chr',f'start_div_{n}'], axis=1)
    df_entropy = df_entropy.set_index('chr_start').transpose()
    df_entropy['Type'] = df_entropy.index.map(lambda x: 0 if 'Normal' in str(x) else (1 if 'GBM' in str(x) else np.nan))
    print(df_entropy.head())
    df_entropy.to_csv(f'/ssd_data/Methylation/results/step04_ML_classifier/v1/{p}_entropy_v1.csv', index=False)



# ### 3) 모델링

# In[ ]:


import pandas as pd
import itertools
from sklearn.model_selection import train_test_split
from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import RandomForestClassifier
from xgboost import XGBClassifier
from sklearn.svm import SVC
from sklearn.metrics import accuracy_score

# p_list = [5, 5000, 25, 50, 100, 125, 150, 250, 500]  # 원하는 p 값 리스트
m_list = ['mean', 'entropy']

# 조합 + 단일 feature set 포함
all_m_cases = [(m,) for m in m_list] + list(itertools.combinations(m_list, 2))

for p in p_list:
    sheet_dict = {}  # 시트별로 저장할 결과 딕셔너리
    for m_case in all_m_cases:
        try:
            if len(m_case) == 1:
                m1 = m_case[0]
                file_path = f'/ssd_data/Methylation/results/step04_ML_classifier/v1/{p}_{m1}_v1.csv'
                df = pd.read_csv(file_path)
                df = df.rename(columns=lambda x: f"{m1}_{x}")
                df['Type'] = df[f'{m1}_Type']
                df = df.drop([f'{m1}_Type'], axis=1)
                combo_name = m1
            else:
                m1, m2 = m_case
                file_path1 = f'/ssd_data/Methylation/results/step04_ML_classifier/v1/{p}_{m1}_v1.csv'
                df1 = pd.read_csv(file_path1)   
                df1 = df1.rename(columns=lambda x: f"{m1}_{x}")
                file_path2 = f'/ssd_data/Methylation/results/step04_ML_classifier/v1/{p}_{m2}_v1.csv'
                df2 = pd.read_csv(file_path2)   
                df2 = df2.rename(columns=lambda x: f"{m2}_{x}")
                df = pd.concat([df1, df2], axis=1)
                df = df.drop_duplicates()
                df['Type'] = df[f'{m1}_Type']
                df = df.drop([f'{m1}_Type', f'{m2}_Type'], axis=1)
                df.to_csv(f'/ssd_data/Methylation/results/step04_ML_classifier/v1/{p}_{m1}_{m2}_v1.csv', index=False)
                combo_name = f'{m1}|{m2}'
            # 데이터 분리
            train_df, val_df = train_test_split(df, test_size=0.2, random_state=42)
            X_train = train_df.drop('Type', axis=1)
            y_train = train_df['Type']
            X_val = val_df.drop('Type', axis=1)
            y_val = val_df['Type']
            # 모델 리스트
            models = {
                'Logistic Regression': LogisticRegression(random_state=42),
                'Random Forest': RandomForestClassifier(random_state=42),
                'XGBoost': XGBClassifier(random_state=42, use_label_encoder=False, eval_metric='logloss'),
                'SVM': SVC(random_state=42)
            }
            # 결과 저장용 리스트
            results_list = []
            for model_name, model in models.items():
                model.fit(X_train, y_train)
                y_pred = model.predict(X_val)
                accuracy = accuracy_score(y_val, y_pred)
                results_list.append({
                    'Model': model_name,
                    'Feature Set': combo_name,
                    'Accuracy': accuracy
                })
            # 시트 저장용
            sheet_df = pd.DataFrame(results_list)
            sheet_dict[combo_name] = sheet_df
        except Exception as e:
            print(f"❌ p={p}, {m_case} 처리 중 오류 발생: {e}")
    # ExcelWriter로 시트별 저장
    output_file = f'/ssd_data/Methylation/results/step04_ML_classifier/v1/Acc_{p}_v1.xlsx'
    with pd.ExcelWriter(output_file, engine='openpyxl') as writer:
        for sheet_name, df in sheet_dict.items():
            # 시트 이름은 31자 이하여야 하므로 자르기
            safe_sheet_name = sheet_name[:31]
            df.to_excel(writer, sheet_name=safe_sheet_name, index=False)
    print(f"✅ p={p}에 대한 모든 결과가 시트별로 저장되었습니다 → {output_file}")



