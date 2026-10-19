# -*- coding: utf-8 -*-
"""학습 혼합군을 A / B 로 나눈다 — 교차 채점(4_report/steps/cross_score.py)의 준비

  학습 CV 는 같은 세포주 · 같은 정상인의 리드를 접기 양쪽에 나눠 쓰므로, 처음 보는 세포주 ·
  정상인으로의 일반화를 재지 못한다. 혼합군을 사람 · 세포주 단위로 둘로 나눠 A 로 학습한 모델이
  B 를, B 로 학습한 모델이 A 를 맞히는지 보면 그것을 학습 자료 안에서 잴 수 있다.

  단위: 정상은 사람, 암은 group(세포주) — 11_pool.py · 05_llrref.py 와 같다.
  가르기: meth_config.split_ref_mix 로 참조군 / 혼합군을 가르고(LLR 참조와 같은 가르기),
          혼합군을 같은 함수로 다시 A / B 로 가른다(바이트 순 번갈아 · 시드 없음).
  --layout-tsv (열: sample, layout)를 주면 정상은 layout(SE · PE 등)마다 따로 번갈아 나눠,
          두 무리의 라이브러리 구성을 맞춘다(층마다 시작하는 쪽을 번갈아 둔다). 주지 않으면 맞추지 않는다.
  LLR_HELDOUT=1 이어야 한다. 0 이면 LLR 참조표에 A · B 가 다 들어가 교차가 성립하지 않는다.
  --use-* 에는 원래 풀을 만든 목록을 준다. 주지 않으면 samples.csv 전부가 기준이 된다.
  주의: A · B 는 «분류기가 처음 보는» 무리일 뿐이다. 패널 선정(1단계)과 jsd · pdr 의 정상 기준 분포는
        A · B 를 모두 보고 만들어졌다. 그래서 «파이프라인 전체가 처음 보는 세포주 · 정상인» 은 아니다.

  출력: <out>/use_normal_A.txt · use_cancer_A.txt · use_normal_B.txt · use_cancer_B.txt
        (검체 이름 한 줄에 하나 · --use-* 목록의 순서)

사용:  python split_ab.py --out <폴더> [--cand <j_candidates 폴더>]
           [--use-normal <목록>] [--use-cancer <목록>] [--layout-tsv <tsv>]
다음:  트리 둘을 만든다. 1단계 산출(results/dmr 아래 j_pool* 를 뺀 것)을 심링크로 걸고,
       트리마다 11_pool.py --use use_<normal|cancer>_<A|B>.txt 로 풀을 만든 뒤 2단계를 돌린다.
       그다음 cross_score.py 로 두 트리를 서로 채점한다(pipeline/README.md).
"""
# 이 파일은 스크립트다. 다른 코드가 import 하면 본체가 그대로 돈다.
if __name__ != '__main__':
    raise ImportError(__file__ + ' 은(는) 스크립트다 — import 하지 않는다')

import os as _os, sys as _sys
_h = _os.path.dirname(_os.path.abspath(__file__))
for _ in range(5):
    _c = _os.path.join(_h, '0_setup')
    if _os.path.isdir(_c):
        _sys.path.insert(0, _c); break
    _h = _os.path.dirname(_h)
else:
    _e = _os.environ.get('METH_CONF_DIR')
    if _e and _os.path.isdir(_e):
        _sys.path.insert(0, _e)
import meth_config as _cfg

import argparse
import os
import sys

import pandas as pd

ap = argparse.ArgumentParser()
ap.add_argument('--out', required=True)
ap.add_argument('--cand', default='%s/results/dmr/j_candidates%s' % (_cfg.METH_ROOT, _cfg.GEN))
ap.add_argument('--use-normal')
ap.add_argument('--use-cancer')
ap.add_argument('--layout-tsv')
ap.add_argument('--force', action='store_true', help='이미 있는 목록을 덮어쓴다')
a = ap.parse_args()

if _cfg.LLR_HELDOUT != 1:
    sys.exit('LLR_HELDOUT=%s — 1 이어야 교차가 성립한다(0 이면 LLR 참조표에 A · B 가 다 들어간다)' % _cfg.LLR_HELDOUT)
S = pd.read_csv(os.path.join(a.cand, 'samples.csv'))
typ = dict(zip(S['sample'].astype(str), S['type'].astype(str)))
grp = dict(zip(S['sample'].astype(str), S['group'].astype(str)))
cancer = sorted(set(typ.values()) - {'Normal'})
if len(cancer) != 1:
    sys.exit('samples.csv 의 암 종류가 하나가 아니다: %s' % cancer)
lay = None
if a.layout_tsv:
    L = pd.read_csv(a.layout_tsv, sep='\t')
    lay = dict(zip(L['sample'].astype(str), L['layout'].astype(str)))
    lack = [s for s in typ if typ[s] == 'Normal' and s not in lay]
    if lack:
        sys.exit('--layout-tsv 에 없는 정상: %s' % lack[:10])
os.makedirs(a.out, exist_ok=True)
_have = [x for x in ('use_normal_A.txt', 'use_normal_B.txt', 'use_cancer_A.txt', 'use_cancer_B.txt')
         if os.path.exists(os.path.join(a.out, x))]
if _have and not a.force:
    sys.exit('이미 있다: %s — 그 목록으로 풀을 만들었다면 출처가 어긋난다. 덮으려면 --force' % _have)

for kind, usef, name in (('Normal', a.use_normal, 'normal'), (cancer[0], a.use_cancer, 'cancer')):
    use = open(usef).read().split() if usef else sorted(s for s in typ if typ[s] == kind)
    _bad = [s for s in use if typ.get(s) != kind]
    if _bad:
        sys.exit('%s 목록에 samples.csv 의 %s 가 아닌 이름이 있다: %s' % (usef, kind, _bad[:10]))
    unit = (lambda s: s) if kind == 'Normal' else (lambda s: grp[s])
    allu = [unit(s) for s in typ if typ[s] == kind]
    ref, mix = _cfg.split_ref_mix(allu)
    if lay is not None and kind == 'Normal':
        ga, gb = [], []
        for i, st in enumerate(sorted(set(lay[u] for u in mix))):
            x, y = _cfg.split_ref_mix([u for u in mix if lay[u] == st])
            if i % 2:                     # 홀수 층이 모두 A 를 하나씩 더 갖지 않게 시작 쪽을 번갈아 둔다
                x, y = y, x
            ga += x; gb += y
    else:
        ga, gb = _cfg.split_ref_mix(mix)
    if set(ga) & set(gb) or set(ref) & (set(ga) | set(gb)):
        sys.exit('나누기가 겹친다 — 멈춤')
    for g, us in (('A', set(ga)), ('B', set(gb))):
        pick = [s for s in use if unit(s) in us]
        with open(os.path.join(a.out, 'use_%s_%s.txt' % (name, g)), 'w') as h:
            h.write('\n'.join(pick) + '\n')
        extra = ''
        if lay is not None and kind == 'Normal':
            extra = ' · ' + ' '.join('%s %d' % (st, sum(lay[u] == st for u in us)) for st in sorted(set(lay.values())))
        print('%-8s %s  단위 %2d · 검체 %2d   (참조군 %d · 혼합군 %d)%s'
              % (kind, g, len(us), len(pick), len(ref), len(mix), extra))
print('씀', os.path.abspath(a.out))
