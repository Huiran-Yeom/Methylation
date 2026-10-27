# LUAD 다음 판 1안 · 10런 — 학습 교차검증 요약 (0.5% · 실제 패널)

학습 CV AUC(`cv_auc`) · 특징 5종(mean · entropy · jsd · llr · pdr)의 31조합 · 조합마다 모델 넷(LR · RF · SVM · XGBoost) 중 최고.
전체 값(모든 비율 · 무작위 패널 rand · randb 포함)은 같은 폴더의 `cv_all.csv` · 보기용 `cv_all.md`.

| 패널 | 깊이 | 1차 조합 mean\|entropy\|jsd\|llr\|pdr | llr 뺀 mean\|entropy\|jsd\|pdr | 31조합 중 최고 |
|---|---|---|---|---|
| bl200 | 5k | 0.746 | 0.758 | 0.768 (entropy) |
| bl200 | 50k | 0.980 | 1.000 | 1.000 (mean\|entropy\|jsd\|pdr) |
| bl200 | 100k | 1.000 | 1.000 | 1.000 (mean 외 같은 값 15개) |
| jsd200 | 5k | 0.718 | 0.627 | 0.744 (mean\|llr\|pdr) |
| jsd200 | 50k | 0.957 | 0.882 | 0.963 (mean\|jsd\|llr) |
| jsd200 | 100k | 0.996 | 0.997 | 0.999 (mean) |

**주의** — 학습 CV 는 패널을 고른 같은 자료로 잰 값이라 성능 근거가 아니다. 판정은 외부검증 1차 지표로만 한다.
