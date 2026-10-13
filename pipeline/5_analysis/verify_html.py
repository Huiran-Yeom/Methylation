#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""업로드된 HTML 보고서의 숫자를 서버 표와 전수 대조한다.

서버 표(silico.txt · patient.txt)는 피처·모델마다 두 줄을 낸다 — 「검체」와 「사람」.
HTML 은 본문에 「AUC는 사람 단위로 보고합니다」라고 적혀 있으므로 「사람」 줄과
맞아야 한다. 양쪽을 (패널, 깊이, 조합) -> 비율 7개로 만들어 칸마다 비교한다.

맞으면 조용하고, 다르면 그 칸을 전부 찍는다. 「대부분 맞다」로 넘기지 않는다.
"""
import io
import json
import re
import sys

RATIOS = ['0.1%', '0.5%', '1%', '2%', '2.5%', '3%', '5%']


def parse_server(path, tag):
    """{(판, 깊이, 조합): {'검체': [...], '사람': [...]}}"""
    out = {}
    panel = depth = None
    cur = None
    for ln in io.open(path, encoding='utf-8', errors='replace'):
        m = re.match(r'\s*#####\s+(\S+)', ln)
        if m:
            panel = m.group(1)
            continue
        m = re.match(r'\s*==\s+(\S+)\s+리드\s+·\s+(\S+)', ln)
        if m:
            depth, panel = m.group(1), m.group(2)
            continue
        # 피처 줄: 이름 모델 수7 '검체'
        m = re.match(r'^(\S+)\s+(\S+)\s+((?:-?\d+\.\d+\s+){7})(검체|사람)\s*$', ln)
        if m:
            cur = (tag, panel, depth, m.group(1))
            out.setdefault(cur, {})[m.group(4)] = [float(x) for x in m.group(3).split()]
            continue
        # 이어지는 줄: 수7 '검체'/'사람'
        m = re.match(r'^\s+((?:-?\d+\.\d+\s+){7})(검체|사람)\s*$', ln)
        if m and cur:
            out[cur][m.group(2)] = [float(x) for x in m.group(1).split()]
    return out


def main():
    tabs = json.load(io.open(sys.argv[1], encoding='utf-8'))
    srv = {}
    for p in sys.argv[2:]:
        import os
        srv.update(parse_server(p, os.path.basename(p)))
    print('서버 표에서 읽은 (판,깊이,조합) 칸: %d' % len(srv))

    # HTML 의 비율 표만 고른다 (머리행이 조합/모델/비율)
    cand = [t for t in tabs if t and t[0][:2] == ['조합', '모델']]
    print('HTML 비율 표: %d개' % len(cand))

    total = same = diff = notfound = 0
    bad = []
    for t in cand:
        for row in t[1:]:
            if len(row) < 9:
                continue
            combo = row[0].replace(' †', '').replace('†', '').strip()
            vals = []
            for x in row[2:9]:
                try:
                    vals.append(float(x))
                except ValueError:
                    vals.append(None)
            # 서버에서 같은 조합의 「사람」 줄을 찾는다 (판·깊이는 모름 → 전부 훑어 일치 찾기)
            hit = None
            for (tg, pn, dp, cb), v in srv.items():
                if cb != combo or '사람' not in v:
                    continue
                if all(a is None or abs(a - b) < 5e-4
                       for a, b in zip(vals, v['사람'])):
                    hit = (tg, pn, dp)
                    break
            total += 1
            if hit:
                same += 1
            else:
                # 조합 이름은 있는데 값이 안 맞는 경우와 아예 없는 경우를 가른다
                have = [k for k in srv if k[3] == combo and '사람' in srv[k]]
                if have:
                    diff += 1
                    bad.append((combo, vals, have[:3]))
                else:
                    notfound += 1
                    bad.append((combo, vals, []))
    print()
    print('대조한 줄 %d · 서버와 「사람」 줄이 일치 %d · 불일치 %d · 조합 없음 %d'
          % (total, same, diff, notfound))
    if bad:
        print()
        print('=== 안 맞은 줄')
        for combo, vals, have in bad[:40]:
            print('  %-18s %s' % (combo, ' '.join('%.3f' % v if v is not None else '—'
                                                  for v in vals)))
            if have:
                print('      서버에 같은 조합은 있다: %s' % (have,))
    return 0


if __name__ == '__main__':
    sys.exit(main())
