#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
폐 jsdb200 의 entropy 가 「얼마나」 갈리는지 잰다. 추기 33항.

왜
  entropy 단독 CV 가 0.1% 와 0.5% 에서 0.6011 / 0.6012 다. 종양 비율 5배에
  +0.0001 이라 종양 기원이 아니다. 그러면 「왜 갈리나」보다 「얼마나 갈리나」가
  먼저다. 변동폭 대비 치우침이 아주 작은데 일관되면 그것은 생물이 아니라
  기술 차이의 서명이고, 원인을 못 찾아도 그 문장은 보고할 수 있다.

  거울상 선례: 변동계수 1% 축의 상관 0 은 「없음」이 아니라 「못 쟀음」이다.
  여기는 반대쪽이다 — 변동이 거의 없는 축의 AUC 0.60 은 「셈」이 아니라
  「미세하고 일관됨」일 수 있다.

내는 것
  칸별 변동계수 · 칸별 표준화 차이(d) · 칸 하나로 낸 AUC 의 분포
  그리고 칸 평균 하나로 낸 AUC — 다변량 0.60 이 어디서 오는지 가른다
"""
import sys

import numpy as np
import pandas as pd

ARMS = [('0.1%', '5'), ('0.5%', '25')]
TARGETS = [
    ('폐', 'j16jsdb200', '<METH_ROOT>/lung/results/panels'),
    ('폐', 'j16bl200', '<METH_ROOT>/lung/results/panels'),
    ('GBM', 'j15bl200', '<METH_ROOT>/results/판'),
]


def auc(pos, neg):
    """Mann-Whitney U 로 낸 AUC. 동점은 0.5."""
    pos, neg = np.asarray(pos, float), np.asarray(neg, float)
    m = ~np.isnan(pos)
    pos = pos[m]
    neg = neg[~np.isnan(neg)]
    if pos.size == 0 or neg.size == 0:
        return np.nan
    allv = np.concatenate([pos, neg])
    r = pd.Series(allv).rank().to_numpy()
    rp = r[:pos.size].sum()
    return (rp - pos.size * (pos.size + 1) / 2.0) / (pos.size * neg.size)


def main():
    for tag, panel, base in TARGETS:
        d = panel + '_5k_uni'
        for lab, mr in ARMS:
            f = '%s/%s/%s/step04_ML_classifier/%s_entropy_%s.csv' % (
                base, panel, d, mr, d)
            try:
                df = pd.read_csv(f, index_col=0)
            except Exception as e:
                print('%s %s %s 못 읽음: %s' % (tag, panel, lab, e))
                continue
            # 주의: Type 은 라벨이다. 08_ML_classifier 도 drop 한다.
            #   특징으로 넣으면 AUC 가 1.0 이 된다.
            if 'Type' not in df.columns:
                raise SystemExit('Type 열이 없다 — 파일 형식이 바뀌었다: ' + f)
            df = df.drop(columns=['Type'])
            isneg = df.index.to_series().str.startswith('Normal').to_numpy()
            X = df.to_numpy(dtype=float)
            print('=' * 72)
            print('%s %s · %s · 칸 %d · 행 %d (음성 %d / 양성 %d)'
                  % (tag, panel, lab, X.shape[1], X.shape[0],
                     isneg.sum(), (~isneg).sum()))
            print('=' * 72)

            mu = np.nanmean(X, axis=0)
            sd = np.nanstd(X, axis=0, ddof=1)
            with np.errstate(divide='ignore', invalid='ignore'):
                cv = np.where(np.abs(mu) > 1e-12, sd / np.abs(mu), np.nan)
            n = X[isneg], X[~isneg]
            mn, mp = np.nanmean(n[0], axis=0), np.nanmean(n[1], axis=0)
            sn, sp = np.nanstd(n[0], axis=0, ddof=1), np.nanstd(n[1], axis=0, ddof=1)
            pooled = np.sqrt((sn ** 2 + sp ** 2) / 2.0)
            with np.errstate(divide='ignore', invalid='ignore'):
                dd = np.where(pooled > 1e-12, (mp - mn) / pooled, np.nan)
            aucs = np.array([auc(n[1][:, j], n[0][:, j]) for j in range(X.shape[1])])

            def q(a, name, fmt='%.4f'):
                a = a[np.isfinite(a)]
                if a.size == 0:
                    print('  %-16s 값 없음' % name)
                    return
                print(('  %-16s 중앙 ' + fmt + ' · 사분 ' + fmt + '~' + fmt
                       + ' · 폭 ' + fmt + '~' + fmt)
                      % (name, np.median(a), np.percentile(a, 25),
                         np.percentile(a, 75), a.min(), a.max()))

            q(cv * 100, '변동계수 %', '%.2f')
            q(dd, '표준화 차이 d')
            q(aucs, '칸 하나 AUC')
            sign = dd[np.isfinite(dd)]
            print('  d 의 부호        양수 %d / %d (%.0f%%)'
                  % (int((sign > 0).sum()), sign.size,
                     100.0 * (sign > 0).mean()))
            good = aucs[np.isfinite(aucs)]
            print('  칸 하나 AUC 0.60 넘는 칸  %d / %d'
                  % (int((good > 0.60).sum()), good.size))

            # 칸 평균 하나로 낸 AUC — 다변량 0.60 이 「분산된 미세 치우침」인지
            rowmean = np.nanmean(X, axis=1)
            a_mean = auc(rowmean[~isneg], rowmean[isneg])
            print('  칸 평균 한 축 AUC  %.4f' % a_mean)
            # 꼬리 기여: 그 축에서 양끝 10행을 빼고 다시
            o = np.argsort(rowmean)
            keep = np.ones(rowmean.size, bool)
            keep[o[:10]] = False
            keep[o[-10:]] = False
            a_trim = auc(rowmean[keep & ~isneg], rowmean[keep & isneg])
            print('  양끝 10행 빼면     %.4f' % a_trim)
            print()
    return 0


if __name__ == '__main__':
    sys.exit(main())
