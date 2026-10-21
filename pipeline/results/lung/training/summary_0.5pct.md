# 폐 SCLC (뒷받침) — 학습 교차검증 요약 (최종 코드 · 0.5% · 실제 패널)

학습 CV AUC(`cv_auc`) · 특징 5종(mean · entropy · jsd · llr · pdr)의 31조합 · 조합마다 모델 넷(LR · RF · SVM · XGBoost) 중 최고.
전체 값(모든 비율 · 무작위 패널 rand · randb 포함)은 같은 폴더의 `cv_all.csv`.

| 패널 | 깊이 | 1차 조합 mean\|entropy\|jsd\|llr\|pdr | llr 뺀 mean\|entropy\|jsd\|pdr | 31조합 중 최고 |
|---|---|---|---|---|
| bl200 | 5k | 0.762 | 0.774 | 0.786 (entropy\|jsd\|llr\|pdr) |
| bl200 | 50k | 0.998 | 1.000 | 1.000 (entropy) |
| jsd200 | 5k | 0.884 | 0.747 | 0.914 (mean\|jsd\|llr\|pdr) |
| jsd200 | 50k | 1.000 | 1.000 | 1.000 (entropy\|jsd) |
| jsdb200 | 5k | 0.815 | 0.650 | 0.831 (entropy\|jsd\|llr) |
| jsdb200 | 50k | 1.000 | 1.000 | 1.000 (entropy\|jsd) |
| jsdb200 | 100k | 1.000 | 1.000 | 1.000 (entropy) |

**주의** — 학습 CV 는 패널을 고른 같은 자료로 잰 값이라 성능 근거가 아니다. 판정은 외부검증 1차 지표로만 한다.
