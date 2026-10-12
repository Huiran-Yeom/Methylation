# 결과 성능 표 — cfDNA 메틸화 분류기 (2026-10-12)

이 폴더는 `pipeline/` 코드로 돌린 **학습 · 외부검증의 성능 표**만 담는다. 원자료 · 검체별 점수 · 공여자 정보는 없다.

- **본 연구는 GBM 이다.** 폐(SCLC 학습 → NSCLC 위주 검증)는 같은 파이프라인이 다른 암에서도 같은 양상인지 보는 뒷받침이다.
- 외부검증은 모두 **반복 10판**이다. 검체마다 비율 칸당 10번 섞어 채점했다.
  - GBM: 178검체 × 10
  - 폐 In silico: 18명 × 10
  - 폐 Patient: 48명 × 10

```
results/
  GBM/
    training/     summary_0.5pct.md · cv_all_primary.csv · cv_all_sensitivity_fixed.csv
    validation/   patient178_tables_primary.txt · patient178_primary_metric_sensitivity_fixed.txt
  lung/
    training/     summary_0.5pct.md · cv_all_primary.csv · cv_all_sensitivity_fixed.csv
    validation/   insilico18_tables_primary.txt · patient48_tables_primary.txt · *_primary_metric_sensitivity_fixed.txt
```

## 1차 지표 (사전등록) — 외부검증 · 깊이 5k · 종양 0.5% · 검체 단위 AUC · 조합 mean|entropy|jsd|llr|pdr

| 코호트 | 판 | bl200 | jsd200 | jsdb200 | 판정 |
|---|---|---|---|---|---|
| **GBM** 환자 178검체(100명) | 1차 | 0.523 | 0.542 | 0.525 | 0.70 미달 → 되돌림 |
| **GBM** | 민감도(결함 수정판) | 0.592 | 0.573 | 0.575 | 같음 |
| 폐 In silico (NSCLC 종양 조직 18명) | 1차 | 0.512 | 0.502 | 0.520 | 0.70 미달 → 되돌림 |
| 폐 In silico | 민감도 | 0.545 | 0.540 | 0.478 | 같음 |
| 폐 Patient (환자 혈장 48명) | 1차 | 0.517 | 0.507 | 0.514 | 0.70 미달 → 되돌림 |
| 폐 Patient | 민감도 | 0.493 | 0.481 | 0.520 | 같음 |

- 사전등록 문턱은 0.70 이다. 세 판 모두 못 넘어 대표 패널을 정하지 않는다(되돌림).
- 민감도 판은 결함 둘을 고쳐 다시 돌린 것이다(LLR 참조 누수 · 묶음 평균 채우기 · 코드 `e8d5fe5d`). 1차 결과를 덮어쓰지 않고 나란히 보고한다.

## 파일 읽는 법
- `validation/*_tables_primary.txt`
  - 패널 · 깊이(5k · 50k · 100k)마다 9조합 × 7비율(0.1~5%)의 검체 단위 · 사람 단위 AUC 표가 있다.
  - 무작위 패널(rand · randb) 표와 1차 지표 블록(짝지은 부트스트랩 구간 · 특이도 · llr 뺀 조합)이 뒤에 붙어 있다.
- `validation/*_sensitivity_fixed.txt` — 결함 수정판의 1차 지표 블록.
- `training/summary_0.5pct.md` — 학습 교차검증 AUC 요약(0.5% · 조합마다 모델 넷 중 최고).
- `training/cv_all_*.csv` — 학습 교차검증 전체. 열은 `panel · depth · mut_reads · fraction_pct · feature_set · model · cv_auc · n` 이다.

## 해석의 한계
- **학습 CV 는 성능 근거가 아니다.** 패널을 고른 같은 자료로 쟀고, 1차 판에는 LLR 참조 누수가 있다. 판정은 외부검증으로만 한다.
- 학습 안(처음 보는 세포주 · 정상인)에서는 일반화된다. 학습–검증 괴리는 세포주 → 환자 단계에서 생긴다.
- 폐는 학습이 SCLC 세포주이고 검증은 NSCLC 위주라 조직형이 다르다. 폐 패널은 혈액이 아닌 폐 조직 DNA 를 잡는다.
- 임상 진단 성능 · 패널 우열 · 암 특이 신호는 이 결과로 주장하지 않는다.
