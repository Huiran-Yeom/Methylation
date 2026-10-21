# GBM (본 연구) — 학습 교차검증 요약 (최종 코드 · 0.5% · 실제 패널)

학습 CV AUC(`cv_auc`) · 특징 5종(mean · entropy · jsd · llr · pdr)의 31조합 · 조합마다 모델 넷(LR · RF · SVM · XGBoost) 중 최고.
전체 값(모든 비율 · 무작위 패널 rand · randb 포함)은 같은 폴더의 `cv_all.csv`.

| 패널 | 깊이 | 1차 조합 mean\|entropy\|jsd\|llr\|pdr | llr 뺀 mean\|entropy\|jsd\|pdr | 31조합 중 최고 |
|---|---|---|---|---|
| bl200 | 5k | 0.698 | 0.630 | 0.727 (mean) |
| bl200 | 50k | 0.992 | 1.000 | 1.000 (mean) |
| jsd200 | 5k | 0.653 | 0.572 | 0.659 (mean\|llr) |
| jsd200 | 50k | 0.999 | 1.000 | 1.000 (mean) |
| jsdb200 | 5k | 0.720 | 0.566 | 0.720 (entropy\|jsd\|llr) |
| jsdb200 | 50k | 0.976 | 0.991 | 0.995 (mean\|entropy\|jsd) |

**주의** — 학습 CV 는 패널을 고른 같은 자료로 잰 값이라 성능 근거가 아니다. 판정은 외부검증 1차 지표로만 한다.
