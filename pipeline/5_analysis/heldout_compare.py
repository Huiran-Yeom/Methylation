#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
LLR 참조 누수 수정 전후 학습 CV 짝비교.

이 파일은 수정 후 산출물(<METH_ROOT_HELDOUT> 아래 cv_*.csv)이 「0개」일 때
작성했다 — 숫자를 보기 전에 집계 규칙을 고정하기 위해서다.
사전등록 추기 27항의 보고 규칙을 그대로 구현한다.

규칙 (추기 27항)
  - 수정 전과 후를 「같은 표」에 나란히 싣는다.
  - 주 지표는 수정판으로 한다.
  - 어느 판이 올라가도 그 하나만으로 판 순위를 매기지 않는다.
    → 이 스크립트는 판 사이 순위·최댓값·「최고 판」을 「출력하지 않는다」.

미리 정한 핵심 대조
  llr 을 쓰는 조합과 쓰지 않는 조합을 갈라서 본다.
  누수가 원인이면 llr 조합만 내려가고 llr 없는 조합은 거의 그대로여야 한다.
  둘이 같이 내려가면 원인이 누수가 아니라 분할 자체(사람 수 15→7/8)다.
  이 해석 규칙도 숫자를 보기 전에 적었다.
"""
import csv
import os
import sys

NL = chr(10)

BEFORE = '<METH_ROOT>/results/판'
AFTER = '<METH_ROOT_HELDOUT>/results/판'
PANELS = ['j15bl200', 'j15jsd200', 'j15jsdb200']
DEPTH = '5k_uni'
MUT = '25'                                    # 5,000 리드의 0.5%
PRIMARY = 'mean|entropy|jsd|llr|pdr'          # 사전등록 1차 조합
NO_LLR = 'mean|entropy|jsd|pdr'               # 거기서 llr 만 뺀 것 (추기 30-2)
ONLY_LLR = 'llr'

# 검증 1차 지표 (5k · 0.5% · 검체 단위 AUC). 수정 「전」 값이고 동결돼 있다.
VALID = {'j15bl200': 0.5227, 'j15jsd200': 0.5420, 'j15jsdb200': 0.5254}


def cv_path(root, panel):
    d = '%s_%s' % (panel, DEPTH)
    return os.path.join(root, panel, d, 'step04_ML_classifier',
                        'cv_%s_%s.csv' % (MUT, d))


def read_cv(path):
    """(feature_set, model) -> cv_auc. 비율이 0.5 인 줄만 받는다."""
    out = {}
    with open(path, newline='', encoding='utf-8') as fh:
        for row in csv.DictReader(fh):
            if abs(float(row['fraction_pct']) - 0.5) > 1e-9:
                continue
            out[(row['feature_set'], row['model'])] = float(row['cv_auc'])
    if not out:
        raise SystemExit('0.5%% 줄이 없다: %s' % path)
    return out


def best_over_models(table, combo):
    """09_save_best_j 가 고르는 방식 — 그 조합의 모델 최댓값."""
    vals = [v for (c, _m), v in table.items() if c == combo]
    return max(vals) if vals else None


def main():
    missing = [cv_path(r, p) for p in PANELS for r in (BEFORE, AFTER)
               if not os.path.exists(cv_path(r, p))]
    if missing:
        raise SystemExit('아직 없다 (%d개):%s%s'
                         % (len(missing), NL, NL.join(missing)))

    rows = []
    for panel in PANELS:
        b = read_cv(cv_path(BEFORE, panel))
        a = read_cv(cv_path(AFTER, panel))

        shared = sorted(set(b) & set(a))
        only_b = sorted(set(b) - set(a))
        only_a = sorted(set(a) - set(b))

        pb = best_over_models(b, PRIMARY)
        pa = best_over_models(a, PRIMARY)
        rows.append((panel, pb, pa, len(shared), len(only_b), len(only_a), b, a))

    print('=' * 74)
    print('1차 지표 조합 — 깊이 5k · 비율 0.5%% · %s' % PRIMARY)
    print('모델 4개 중 최댓값 (09_save_best_j 와 같은 방식)')
    print('=' * 74)
    print('%-12s %10s %10s %10s   %s' % ('판', '수정 전', '수정 후', 'Δ', '칸'))
    for panel, pb, pa, ns, nb, na, _b, _a in rows:
        if pb is None or pa is None:
            print('%-12s %10s %10s %10s   공통 %d · 전만 %d · 후만 %d'
                  % (panel, '없음' if pb is None else '%.4f' % pb,
                     '없음' if pa is None else '%.4f' % pa, '—', ns, nb, na))
        else:
            print('%-12s %10.4f %10.4f %+10.4f   공통 %d · 전만 %d · 후만 %d'
                  % (panel, pb, pa, pa - pb, ns, nb, na))

    print()
    print('=' * 74)
    print('추기 30-4 사전 예측 대조 — 예측은 11:17(산출물 0개)에 적었다')
    print('=' * 74)
    print('%-12s %8s %8s %8s %8s %8s %8s'
          % ('판', '수정전', '예측', '실측', '검증', '예측간극', '실측간극'))
    for panel, _pb, _pa, _ns, _nb, _na, b, a in rows:
        pred = best_over_models(b, NO_LLR)        # llr 뺀 조합의 수정 「전」 CV
        before = best_over_models(b, PRIMARY)
        after = best_over_models(a, PRIMARY)
        v = VALID.get(panel)
        if None in (pred, before, after) or v is None:
            print('%-12s  칸이 없다 — 조합 집합을 확인하라' % panel)
            continue
        print('%-12s %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f'
              % (panel, before, pred, after, v, abs(pred - v), abs(after - v)))
    print()
    print('예측 = 수정 전 `%s` 의 CV. 누수를 빼면 여기로 내려간다고 적었다.' % NO_LLR)
    print('실측간극이 거의 0 이면 학습–검증 괴리가 LLR 누수로 설명된다.')
    print('남으면 다른 원인이 더 있다 — 「누수가 전부였다」로 적지 않는다.')
    print('실측이 예측보다 낮으면 참조 축소(29-5)가 겹친 것이다.')

    print()
    print('=' * 74)
    print('`llr` 단독 · `llr` 뺀 조합 · 1차 조합 — 전후')
    print('=' * 74)
    print('%-12s %-24s %9s %9s %9s' % ('판', '조합', '수정전', '수정후', 'Δ'))
    for panel, _pb, _pa, _ns, _nb, _na, b, a in rows:
        for c in (ONLY_LLR, NO_LLR, PRIMARY):
            vb, va = best_over_models(b, c), best_over_models(a, c)
            if vb is None or va is None:
                print('%-12s %-24s %9s %9s %9s'
                      % (panel, c, '없음' if vb is None else '%.4f' % vb,
                         '없음' if va is None else '%.4f' % va, '—'))
                continue
            print('%-12s %-24s %9.4f %9.4f %+9.4f' % (panel, c, vb, va, va - vb))
        print()

    print('=' * 74)
    print('미리 정한 대조 — llr 쓰는 조합 대 안 쓰는 조합')
    print('조합별로 모델 최댓값을 내고, 그 Δ 들의 중앙값과 부호를 센다')
    print('=' * 74)
    for panel, _pb, _pa, _ns, _nb, _na, b, a in rows:
        combos = sorted({c for c, _m in b} & {c for c, _m in a})
        with_llr, without_llr = [], []
        for c in combos:
            d = best_over_models(a, c) - best_over_models(b, c)
            if 'llr' in c.split('|'):
                with_llr.append((c, d))
            else:
                without_llr.append((c, d))
        print('[%s]' % panel)
        for label, group in (('llr 있음', with_llr), ('llr 없음', without_llr)):
            if not group:
                print('   %-9s 조합 0개' % label)
                continue
            ds = sorted(d for _c, d in group)
            mid = ds[len(ds) // 2] if len(ds) % 2 else (ds[len(ds) // 2 - 1] + ds[len(ds) // 2]) / 2
            down = sum(1 for d in ds if d < 0)
            print('   %-9s 조합 %2d · 중앙 Δ %+.4f · 내려간 것 %d/%d · 폭 %+.4f~%+.4f'
                  % (label, len(ds), mid, down, len(ds), ds[0], ds[-1]))
        print()

    print('=' * 74)
    print('조합 전부 — 수정 전 / 수정 후 / Δ (모델 최댓값)')
    print('=' * 74)
    for panel, _pb, _pa, _ns, _nb, _na, b, a in rows:
        combos = sorted({c for c, _m in b} & {c for c, _m in a})
        print('[%s] 공통 조합 %d개' % (panel, len(combos)))
        for c in combos:
            vb, va = best_over_models(b, c), best_over_models(a, c)
            mark = ' *' if c == PRIMARY else '  '
            print('  %s%-34s %7.4f %7.4f %+8.4f' % (mark, c, vb, va, va - vb))
        print()

    print('판 사이 비교·순위·「최고 판」은 일부러 내지 않는다 (추기 27항).')
    print('이 값은 학습 CV 다. 검증 성능이 아니고, 검증 재채점은 입력이')
    print('없어져서 못 한다 (추기 28항).')


if __name__ == '__main__':
    sys.exit(main())
