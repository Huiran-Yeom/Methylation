# -*- coding: utf-8 -*-
"""1단계 단계들이 같이 쓰는 작은 함수

gen(cfg)                       패널 세대: 환경변수 GEN, 없으면 config.conf
need_cancer(n, cancer)         암 검체(단위)가 0개면 멈춤 (CANCER_LABEL 과 samples.csv type 어긋남)
unit_matrix(B, M, cancer)      검체 β 행렬 → 단위(세포주, 사람) 평균. (X, y, 단위 이름)
fill_col_mean(X)               빈 칸을 열 평균으로 (열이 다 비면 0)
"""
import os
import sys

import numpy as np
import pandas as pd


def gen(cfg):
    # 산출 폴더 이름에 GEN(패널 세대) 붙임. 못 읽으면 중단
    g = os.environ.get('GEN', '')
    if not g:
        try:
            g = cfg.GEN
        except Exception:
            raise SystemExit('GEN(패널 세대)을 못 읽었습니다 — config.conf 의 GEN 을 채우거나 '
                             '환경변수로 주십시오. 비워 두면 옛 세대 폴더를 가리킵니다.')
    return g


def need_cancer(n, cancer, msg=None):
    # CANCER_LABEL 과 samples.csv type 이 안 맞으면 암 0개 → 중단
    if int(n) == 0:
        sys.exit(msg or ('samples.csv 에 type=%s 인 검체가 없습니다. CANCER_LABEL 과 1단계가 쓴 값이'
                         ' 어긋났습니다 (현재 CANCER_LABEL=%s).' % (cancer, cancer)))


def unit_matrix(B, M, cancer, msg=None):
    # 세포주 단위로 평균 (같은 세포주 반복 파일을 독립으로 세지 않게). 정상은 각자 한 단위
    X = B.T.reindex(M['sample']).values.astype(np.float32)
    y = (M['type'] == cancer).astype(int).values
    need_cancer(y.sum(), cancer, msg)
    cols = [str(i) for i in range(X.shape[1])]
    D = pd.DataFrame(X, columns=cols); D['_g'] = M['group'].values; D['_y'] = y
    A = D.groupby('_g', sort=True).mean()
    return A[cols].values.astype(np.float32), (A['_y'].values > 0.5).astype(int), np.array(A.index)


def fill_col_mean(X):
    mu = np.nanmean(X, axis=0); mu = np.where(np.isfinite(mu), mu, 0.0)
    return np.where(np.isfinite(X), X, mu)
