# GBM (본 연구) · 결함 수정 2판 — 학습 교차검증 요약 (0.5% · 실제 패널)

학습 CV AUC(`cv_auc`) · 특징 5종(mean · entropy · jsd · llr · pdr)의 31조합 · 조합마다 모델 넷(LR · RF · SVM · XGBoost) 중 최고.
전체 값(모든 비율 · 무작위 패널 rand · randb 포함)은 같은 폴더의 `cv_all.csv` · 보기용 `cv_all.md`.

| 패널 | 깊이 | 1차 조합 mean\|entropy\|jsd\|llr\|pdr | llr 뺀 mean\|entropy\|jsd\|pdr | 31조합 중 최고 |
|---|---|---|---|---|
| bl200 | 5k | 0.656 | 0.672 | 0.675 (mean\|entropy\|llr) |
| bl200 | 50k | 0.978 | 0.998 | 1.000 (mean\|pdr) |
| jsd200 | 5k | 0.719 | 0.601 | 0.752 (jsd\|llr\|pdr) |
| jsd200 | 50k | 0.996 | 1.000 | 1.000 (mean\|pdr 외 같은 값 2개) |
| jsdb200 | 5k | 0.670 | 0.529 | 0.670 (mean\|entropy\|jsd\|llr\|pdr) |
| jsdb200 | 50k | 0.982 | 0.985 | 0.999 (mean\|jsd) |

**주의** — 학습 CV 는 패널을 고른 같은 자료로 잰 값이라 성능 근거가 아니다. 판정은 외부검증 1차 지표로만 한다.
