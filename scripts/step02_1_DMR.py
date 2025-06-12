#!/usr/bin/env python
# coding: utf-8

# ### cov file load

# In[ ]:


import os

# 디렉토리 경로
directory1 = "/ssd_data/Methylation/sample_data/step02_DMR/Normal/"
directory2 = "/ssd_data/Methylation/sample_data/step02_DMR/GBM/"

# 파일 패턴에 맞는 파일 중 상위 5개 추출 함수
def get_top5_largest_files_by_pattern(directory, suffix_pattern):
    matched_files = []
    for filename in os.listdir(directory):
        if filename.endswith(suffix_pattern):
            filepath = os.path.join(directory, filename)
            if os.path.isfile(filepath):
                matched_files.append((filename, os.path.getsize(filepath)))
    # 크기 기준 내림차순 정렬 후 상위 5개 선택
    top5 = sorted(matched_files, key=lambda x: x[1], reverse=True)[:5]
    return [fname for fname, _ in top5]  # 파일 이름만 추출하여 반환

# 상위 5개 파일 이름만 리스트에 저장
normal_samples = get_top5_largest_files_by_pattern(directory1, '_1_val_1_bismark_bt2_pe.sort.bismark.cov.gz')
cancer_samples = get_top5_largest_files_by_pattern(directory2, '_1_val_1_bismark_bt2_pe.sort.bismark.cov.gz')

# 결과 출력
print("정상 샘플 상위 5개 파일:")
for fname in normal_samples:
    print(fname)

print("\n암 샘플 상위 5개 파일:")
for fname in cancer_samples:
    print(fname)

# 결과 저장
with open('/ssd_data/Methylation/results/step02_DMR/top5_normal_samples.txt', 'w') as f:
    for fname in normal_samples:
        f.write(fname + '\n')

with open('/ssd_data/Methylation/results/step02_DMR/top5_cancer_samples.txt', 'w') as f:
    for fname in cancer_samples:
        f.write(fname + '\n')

print("\n저장 완료: top5_normal_samples.txt, top5_cancer_samples.txt")


# ### 데이터 처리(100bp & 4 count 이상)

# In[ ]:


def group_bp(df, col, i, n):
    df[f'{col}_count{i}'] = df[f'{col}{i}'].notna().astype(int)  # Normal 열의 결측값 확인
    df[f'start_div_{n}'] = df['start'] // n
    df[f'end_div_{n}'] = df['end'] // n
    grouped = df.groupby(['chr', f'start_div_{n}'])[f'{col}_count{i}'].sum().reset_index()
    return grouped



# In[ ]:


import pandas as pd

n = 100
th = 4

df_1 = pd.read_csv(rf'{directory1}{normal_samples[0]}', sep='\t', header=None, 
                   names=['chr', 'start', 'end', 'Normal1','tn','cn'], dtype={'chr': str})
# df_1 = df_1[(df_1['tn']+df_1['cn']) > th]
df_1 = df_1[df_1.columns[:4]]
df_1['Normal_count1'] = df_1.iloc[:, 3:].notna().sum(axis=1)
df_1[f'start_div_{n}'] = df_1['start'] // n
df_1[f'end_div_{n}'] = df_1['end'] // n

grouped1 = df_1.groupby(['chr', f'start_div_{n}'])['Normal_count1'].sum().reset_index()

for i in range(2, len(normal_samples) + 1):
    df1_1 = pd.read_csv(rf'{directory1}{normal_samples[i - 1]}', sep='\t', header=None, 
                         names=['chr', 'start', 'end', f'Normal{i}','tn','cn'], dtype={'chr': str})
    # df1_1 = df1_1[(df1_1['tn']+df1_1['cn']) > th]
    df1_1 = df1_1[df1_1.columns[:4]]
    grouped1_1 = group_bp(df1_1, f'Normal', i, n)
    grouped1 = pd.merge(grouped1, grouped1_1, on=['chr', f'start_div_{n}'], how='inner')

grouped1


# In[ ]:


df_2 = pd.read_csv(rf'{directory2}{cancer_samples[0]}', sep='\t', header=None, 
                   names=['chr', 'start', 'end', 'GBM1','tn','cn'], dtype={'chr': str})
# df_2 = df_2[(df_2['tn']+df_2['cn']) > th]
df_2 = df_2[df_2.columns[:4]]
df_2['GBM_count1'] = df_2.iloc[:, 3:].notna().sum(axis=1)
df_2[f'start_div_{n}'] = df_2['start'] // n
df_2[f'end_div_{n}'] = df_2['end'] // n
grouped2 = df_2.groupby(['chr', f'start_div_{n}'])['GBM_count1'].sum().reset_index()

for i in range(2, len(cancer_samples) + 1):
    df2_1 = pd.read_csv(rf'{directory2}{cancer_samples[i - 1]}', sep='\t', header=None, 
                        names=['chr', 'start', 'end', f'GBM{i}','tn','cn'], dtype={'chr': str})
    # df2_1 = df2_1[(df2_1['tn']+df2_1['cn']) > th]
    df2_1 = df2_1[df2_1.columns[:4]]
    grouped2_1 = group_bp(df2_1, 'GBM', i, n)
    grouped2 = pd.merge(grouped2, grouped2_1, on=['chr', f'start_div_{n}'], how='inner')

grouped2


# In[ ]:


df = pd.merge(grouped1, grouped2, on=['chr', f'start_div_{n}'], how='inner')
# df = df.drop(f'end_div_{n}', axis=1)
df


# In[ ]:


df.to_csv(f'/ssd_data/Methylation/results/step02_DMR/sample_div_{n}_v1.csv', index=False)


# In[ ]:


th=4

df2 = df[(df.iloc[:, 2:] >= th).all(axis=1)]
df2


# In[ ]:


df2.to_csv(f'/ssd_data/Methylation/results/step02_DMR/{n}기준_{th}이상_v1.csv', index=False)


# In[ ]:


def mean_bp(df, col, i, n, g_list):
    df[f'start_div_{n}'] = df['start'] // n
    grouped = df.groupby(['chr', f'start_div_{n}'])[f'{col}{i}'].mean().reset_index()
    grouped_df = pd.merge(g_list, grouped, on=['chr', f'start_div_{n}'], how='inner')
    return grouped_df



# In[ ]:


import pandas as pd

n=100
th=4

df2 = pd.read_csv(f'/ssd_data/Methylation/results/step02_DMR/{n}기준_{th}이상_v1.csv')
g_list = df2[['chr','start_div_100']]

df_1 = pd.read_csv(rf'{directory1}{normal_samples[0]}', sep='\t', header=None, names=['chr', 'start', 'end', 'Normal1','tn','cn'], dtype={'chr': str})
df_1 = df_1[df_1.columns[:4]]
df_1[f'start_div_{n}'] = df_1['start'] // n
grouped1_2 = df_1.groupby(['chr', f'start_div_{n}'])['Normal1'].mean().reset_index()
grouped1_2 = pd.merge(g_list, grouped1_2, on=['chr', f'start_div_{n}'], how='inner')

for i in range(2, len(normal_samples)+1):
    df1_1 = pd.read_csv(rf'{directory1}{normal_samples[i - 1]}', sep='\t', header=None, 
                        names=['chr', 'start', 'end', f'Normal{i}','tn','cn'], dtype={'chr': str})
    df1_1 = df1_1[df1_1.columns[:4]]
    grouped1_3 = mean_bp(df1_1, 'Normal', i, n, g_list)
    grouped1_2 = pd.merge(grouped1_2, grouped1_3, on=['chr', f'start_div_{n}'], how='inner')

grouped1_2


# In[ ]:


import pandas as pd

n=100
th=4

df2 = pd.read_csv(f'/ssd_data/Methylation/results/step02_DMR/{n}기준_{th}이상_v1.csv')
g_list = df2[['chr','start_div_100']]

df_2 = pd.read_csv(rf'{directory2}{cancer_samples[0]}', sep='\t', header=None, names=['chr', 'start', 'end', 'GBM1','tn','cn'], dtype={'chr': str})
df_2 = df_2[df_2.columns[:4]]
df_2[f'start_div_{n}'] = df_2['start'] // n
grouped2_2 = df_2.groupby(['chr', f'start_div_{n}'])['GBM1'].mean().reset_index()
grouped2_2 = pd.merge(g_list, grouped2_2, on=['chr', f'start_div_{n}'], how='inner')

for i in range(2, len(cancer_samples)+1):
    df2_1 = pd.read_csv(rf'{directory2}{cancer_samples[i - 1]}', sep='\t', header=None, names=['chr', 'start', 'end', f'GBM{i}','tn','cn'], dtype={'chr': str})
    df2_1 = df2_1[df2_1.columns[:4]]
    grouped2_3 = mean_bp(df2_1, 'GBM', i, n, g_list)
    grouped2_2 = pd.merge(grouped2_2, grouped2_3, on=['chr', f'start_div_{n}'], how='inner')

grouped2_2


# In[ ]:


df3 = pd.merge(grouped1_2, grouped2_2, on=['chr', f'start_div_{n}'], how='inner')
df3.iloc[:, 2:] = df3.iloc[:, 2:] / 100
df3['chr_start'] = df3['chr'].astype(str) + '_' + df3[f'start_div_{n}'].astype(str)
df3 = df3.drop(['chr',f'start_div_{n}'], axis=1)
df3 = df3.set_index('chr_start').transpose()
df3['Type'] = df3.index.to_series().apply(lambda x: 0 if 'Normal' in str(x) else 1)

df3


# In[ ]:


df3.to_csv(f'/ssd_data/Methylation/results/step02_DMR/{n}기준_{th}이상_평균_v1.csv', index=False)


# ### 처리한 data load

# In[ ]:


import pandas as pd

n=100
th=4

df3 =pd.read_csv(f'/ssd_data/Methylation/results/step02_DMR/{n}기준_{th}이상_평균_v1.csv')


# ### RandomForest

# In[ ]:


import numpy as np
import pandas as pd
from sklearn.ensemble import RandomForestClassifier

X = df3[df3.columns[:-1]]
# X = X.drop(columns=['Unnamed: 0'])
y = df3['Type']

rf = RandomForestClassifier(random_state=42)
rf.fit(X, y)

importances = rf.feature_importances_

indices = np.argsort(importances)[::-1]
RandomForest_df = pd.DataFrame({
    'Feature': X.columns[indices],
    'RandomForest': importances[indices].round(4)
})

RandomForest_df = RandomForest_df[RandomForest_df['RandomForest'] != 0]
RandomForest_df = RandomForest_df.head(1000)
RandomForest_df


# ### XGB

# In[ ]:


import xgboost as xgb

X = df3[df3.columns[:-1]]
# X = X.drop(columns=['Unnamed: 0'])
y = df3['Type']

# XGBoost 모델 훈련
model = xgb.XGBClassifier(random_state=42)
model.fit(X, y)

# Feature importances 얻기
importances = model.feature_importances_

# 인덱스를 중요도에 따라 정렬
indices = np.argsort(importances)[::-1]
XGB_df = pd.DataFrame({
    'Feature': X.columns[indices],
    'XGBoost': importances[indices].round(4)
})

XGB_df = XGB_df[XGB_df['XGBoost'] != 0]
XGB_df = XGB_df[:1000]
XGB_df


# ### LogisticRegression

# In[ ]:


from sklearn.linear_model import LogisticRegression

# Logistic Regression 모델 훈련
model = LogisticRegression(max_iter=1000)
model.fit(X, y)

# Feature importances (계수) 얻기
importances = model.coef_[0]
indices = np.argsort(np.abs(importances))[::-1]  # 절댓값을 기준으로 정렬

# Feature 이름과 중요도 출력
print("Feature importances:")
LogisticRegression_df = pd.DataFrame({
    'Feature': X.columns[indices],
    'Logistic_Regression': importances[indices].round(4)
})

LogisticRegression_df = LogisticRegression_df[LogisticRegression_df['Logistic_Regression'] != 0]
LogisticRegression_df = LogisticRegression_df[:1000]
LogisticRegression_df


# ### T-Test

# In[ ]:


from scipy import stats

t_test_resultss = {}
X = df3[df3.columns[:-1]]
# X = X.drop(columns=['Unnamed: 0'])
y = df3['Type']
type_0 = X[y == 0]
type_1 = X[y == 1]

for column in X.columns:
    t_stat, p_val = stats.ttest_ind(type_0[column], type_1[column], equal_var=False)
    t_test_resultss[column] = p_val

t_test_df = pd.DataFrame({
    'Feature': list(t_test_resultss.keys()),
    'p_value': list(t_test_resultss.values())
}).sort_values(by='p_value')

## 조건
T_test__df = t_test_df[t_test_df['p_value'] < 0.05]
T_test__df = T_test__df.head(1000)
T_test__df


# ### LogisticRegression & T-test 공통 중요 피처 확인

# In[ ]:


common_features = pd.merge(LogisticRegression_df,T_test__df,on='Feature',how='inner')
# feature_importances_df = pd.merge(feature_importances_df,feature_importances_df3,on='Feature',how='inner')
# feature_importances_df = feature_importances_df.sort_values(by='Feature')
common_features


print("Logistic Regression Importance for All Features:\n", LogisticRegression_df)
print("\nT-test p-values for All Features:\n", T_test__df)
print("\nCommon Features (Importance > 0 & p-value < 0.05) with Importance and p-value:\n", common_features)


# In[ ]:


common_features.to_csv(f'/ssd_data/Methylation/results/step02_DMR/{n}기준_{th}이상_feature_importances(inner)_v1.csv', index=False)


# In[ ]:




