# 5_analysis — 사전등록 추기의 숫자를 낸 분석 스크립트

`0_setup`~`4_report` 의 파이프라인 밖에서, 파이프라인 산출물을 읽어 사전등록 추기
(29~45)와 보고서 검증에 쓴 숫자를 낸 스크립트다. 파이프라인을 다시 돌리지 않는다.

## 실행본과의 차이 — 경로와 검체 ID 만 가렸다

서버에서 실제로 돈 파일에서 아래 문자열만 자리표시자로 바꿨다. 동작이 달라진 곳은 없다.
그래서 **그대로는 돌지 않는다** — 자리표시자를 자기 환경 값으로 채워야 한다.

| 자리표시자 | 뜻 |
|---|---|
| `<METH_ROOT>` | 파이프라인 산출물 뿌리 (`0_setup` 의 `METH_ROOT`) |
| `<METH_ROOT_HELDOUT>` | LLR 참조 누수를 고친 재학습판의 뿌리 |
| `<WORK>` | 작업 폴더 (중간 산출물 · 표 · 병기 표) |
| `<정상A>` `<정상B>` | 학습 정상 15명 중 같은 사람인지 의심한 두 런 |
| `<검체>` | 재현 시험에 쓴 Patient 검체 하나 |
| `<정상 7명>` | Patient 검증의 정상 배경 7명 |
| `<사내1>`~`<사내4>` | 학습 정상 중 연구실 공여자 4명 (이름이라 가렸다) |

검체 ID 를 가린 이유는 `3_validate/steps/cohorts/sample_cohort.conf.example` 과 같다.

## 스크립트

머리말(docstring)에 적힌 사전등록 추기 번호를 옮겼다. 번호가 없는 것은 머리말을 본다.

| 파일 | 하는 일 | 추기 |
|---|---|---|
| `llr_leak_scale.py` | LLR 자기참조 누수 크기를 척도로 나눈다 | 29 |
| `llr_refsize_8.py` | 참조 크기 15 → 8 의 순수 효과 | 머리말 |
| `entropy_dispersion.py` | 폐 jsdb200 의 entropy 가 얼마나 갈리나 | 33 |
| `onefeat_vs_many.py` | 칸 평균 한 칸 대 197칸 | 머리말 |
| `null_5050.py` | 신호가 없을 때의 AUC — 귀무 바닥 | 머리말 |
| `matched_null.py` | 관측 CV 와 귀무 바닥을 같은 클래스 크기에서 | 머리말 |
| `val_floor.py` | 검증 AUC 의 재어진 우연 기준선 | 38 |
| `heldout_compare.py` | LLR 누수 수정 전후 학습 CV 짝비교 | 27 |
| `floor_split.py` | 귀무 바닥을 두 조각(모델 선택 · CV)으로 | 머리말 |
| `matched_null_heldout.py` | 수정판의 귀무 바닥 | 머리말 |
| `delta_null.py` | 수정 전후 Δ 의 귀무 분포 | 42 |
| `srr_same_person.py` · `srr_homozygote_test.py` | 학습 정상 두 런이 같은 사람인지 | 머리말 |
| `rep7_same.py` | 폐 Patient 반복 7→10 재현 시험 | 44 |
| `lung10_compare.py` | 폐 Patient 7판 대 10판 | 44 |
| `lung_stage.py` · `lung_stage_depth.py` · `lung_stage_grid.py` | 폐 Patient 병기별 보조 분석 | 45 |
| `verify_html.py` | 보고서 HTML 숫자를 서버 표와 전수 대조 | — |
| `build_lung_report.py` | 폐 외부검증 보고서를 서버 표에서 다시 만든다 | 44 |

## 넣지 않은 것

- 실제 코호트 정의 (`~/cohorts/*.conf`) — 검체 목록과 통제접근 데이터셋 목록이 들어 있다
- 자료 · 표 · 중간 산출물 — 데이터는 저장소에 넣지 않는다
- 서버 작업 감시용 스크립트
