# GBM 학습 교차검증 요약 — 0.5% · 실제 패널

학습 CV AUC (`cv_<비율>_<판>.csv` 의 `cv_auc`) · 조합마다 모델 넷 중 최고(09 단계의 BEST 와 같은 규칙).
1차 = 이번 판 원래 실행 · 민감도 = 결함 수정판(LLR 참조 누수 · 묶음 평균 채우기 고침). 전체 값은 같은 폴더의 `cv_all_*.csv`(무작위 패널 rand · randb 포함).

## 1차

| 패널 | 깊이 | 1차 조합 mean\|entropy\|jsd\|llr\|pdr | llr 뺀 mean\|entropy\|jsd\|pdr | 조합 수 | 전체 조합 중 최고 |
|---|---|---|---|---|---|
| bl200 | 5k | 0.686 | 0.520 | 31 | 0.695 (jsd\|llr) |
| bl200 | 15k | 0.798 | 0.691 | 31 | 0.815 (mean\|llr) |
| bl200 | 50k | 0.982 | 0.982 | 31 | 0.988 (mean\|jsd\|pdr) |
| bl200 | 100k | 1.000 | 1.000 | 31 | 1.000 (entropy\|jsd\|llr) |
| jsd200 | 5k | 0.797 | 0.674 | 127 | 0.805 (mean\|readent\|mhl\|llr) |
| jsd200 | 15k | 0.902 | 0.754 | 127 | 0.917 (mean\|entropy\|readent\|jsd\|llr\|pdr) |
| jsd200 | 50k | 1.000 | 0.989 | 127 | 1.000 (readent\|mhl\|jsd\|llr\|pdr) |
| jsdb200 | 5k | 0.713 | 0.505 | 127 | 0.722 (mean\|readent\|mhl\|jsd\|llr\|pdr) |
| jsdb200 | 15k | 0.871 | 0.548 | 127 | 0.873 (readent\|jsd\|llr) |
| jsdb200 | 50k | 0.985 | 0.910 | 127 | 0.986 (mean\|entropy\|mhl\|jsd\|llr) |

## 민감도(결함 수정판)

| 패널 | 깊이 | 1차 조합 mean\|entropy\|jsd\|llr\|pdr | llr 뺀 mean\|entropy\|jsd\|pdr | 조합 수 | 전체 조합 중 최고 |
|---|---|---|---|---|---|
| bl200 | 5k | 0.698 | 0.630 | 31 | 0.727 (mean) |
| bl200 | 50k | 0.992 | 1.000 | 31 | 1.000 (mean) |
| jsd200 | 5k | 0.653 | 0.572 | 31 | 0.660 (mean\|jsd\|llr) |
| jsd200 | 50k | 0.999 | 1.000 | 31 | 1.000 (mean) |
| jsdb200 | 5k | 0.720 | 0.565 | 31 | 0.720 (entropy\|jsd\|llr) |
| jsdb200 | 50k | 0.976 | 0.991 | 31 | 0.995 (mean\|entropy\|jsd) |

**주의** — 학습 CV 는 패널을 고른 같은 자료로 잰다(사전등록 29-10 · 학습 CV 를 성능 근거로 쓰지 않는다).
1차 판에는 LLR 참조 누수가 있다(민감도 판에서 고침). GBM 1차 판에는 옛 특징(readent · mhl)이 든 조합도 있다. 판정은 외부검증 1차 지표로만 한다.
