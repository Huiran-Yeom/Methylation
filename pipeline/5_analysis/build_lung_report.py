# -*- coding: utf-8 -*-
"""폐 외부검증 보고서를 서버 표에서 다시 만든다. 2026-10-07

틀    업로드된 「폐암 외부검증 결과 (4).html」. 머리·CSS·글·스크립트는 그대로 쓰고
      표와 그림은 «전부» 서버 표에서 새로 채운다. 손으로 넣는 숫자가 없다.
입력  입력/silico.txt · patient.txt   (4_report/make_report.sh 출력)
      입력/gbm.txt                    (GBM 외부검증 표 · <WORK>/표/02_GBM/stored_5k_v2.txt
                                      · 코호트 patient178_j51 · 5k·50k·100k 가 다 들어 있다.
                                      fix4_50k.txt 는 다른 코호트라 쓰지 않는다)
      입력/dev_cv.json                (models/<조합>/pooled/BEST.joblib 의 cv_auc)
      입력/template.html
바꾸는 것
      0장 1차 지표(새) · 3~6장 머리 경계선 · 세 판(jsdb200) · 반복 설명 · 7장 한계(새)
검사  python build_lung_report.py 입력 출력.html --regress
      → 기존 bl200·jsd200 표를 다시 그려 틀과 바이트 단위로 같은지 본다.
        다르면 rc=1. 렌더러가 틀을 그대로 따르는지의 증명이다.
"""
import datetime
import hashlib
import io
import json
import math
import os
import re
import sys

COMBOS = ['ME', 'ME+JSD', 'LLR', 'ME+LLR', 'ME+JSD+LLR', 'PDR', 'ME+PDR',
          'ME+JSD+PDR', 'ME+JSD+LLR+PDR']
COMBO_KEY = {'ME': 'mean|entropy', 'ME+JSD': 'mean|entropy|jsd', 'LLR': 'llr',
             'ME+LLR': 'mean|entropy|llr', 'ME+JSD+LLR': 'mean|entropy|jsd|llr',
             'PDR': 'pdr', 'ME+PDR': 'mean|entropy|pdr',
             'ME+JSD+PDR': 'mean|entropy|jsd|pdr',
             'ME+JSD+LLR+PDR': 'mean|entropy|jsd|llr|pdr'}
PRIMARY = 'ME+JSD+LLR+PDR'
RATIOS = ['0.1%', '0.5%', '1%', '2%', '2.5%', '3%', '5%']
DEPTHS = ['5k', '50k', '100k']
PANELS = ['bl200', 'jsd200', 'jsdb200']
NULLS = ['rand1200', 'rand2200', 'rand3200', 'randb1200', 'randb2200', 'randb3200']
TAG = {'bl200': 'q&lt;0.05 &amp; |dβ|≥30 → 부트스트랩 순위',
       'jsd200': 'JSD 순위 상위',
       'jsdb200': 'JSD 부트스트랩 빈도 상위 · 귀무 문턱 통과'}
MODEL = {'Ran': 'RF', 'Log': 'LR', 'XGB': 'XGB', 'SVM': 'SVM'}
DAG = (' <span class="dag" title="동점이 많아 순위를 주장할 수 없는 줄">†</span>')
NL = '\n'


# ── 서버 표 읽기 ──────────────────────────────────────────────
def parse(path):
    """-> rows{(판,깊이,조합): {model, 검체[7], 사람[7], tie}}, heads{(판,깊이): {...}}"""
    rows, heads = {}, {}
    panel = depth = cur = None
    for ln in io.open(path, encoding='utf-8'):
        m = re.match(r'==\s+(\S+)\s+리드\s+·\s+(\S+)\s+·\s+\S+\s+·\s+채점\s+(\d+)/(\d+)'
                     r'\s+·\s+비율당 양성/음성 각\s+(\d+)', ln)
        if m:
            depth, panel = m.group(1), m.group(2)
            heads[(panel, depth)] = {'scored': int(m.group(3)), 'total': int(m.group(4)),
                                     'per': int(m.group(5))}
            cur = None
            continue
        if ln.startswith('=='):
            cur = None
            continue
        m = re.match(r'^(\S+)\s+(\S+)\s+((?:-?\d+\.\d+\s+){7})검체\s*$', ln)
        if m:
            cur = (panel, depth, m.group(1))
            rows[cur] = {'model': MODEL.get(m.group(2), m.group(2)), 'tie': False,
                         '검체': [float(x) for x in m.group(3).split()]}
            continue
        m = re.match(r'^\s+((?:-?\d+\.\d+\s+){7})사람\s*$', ln)
        if m and cur:
            rows[cur]['사람'] = [float(x) for x in m.group(1).split()]
            continue
        if cur and '동점이 많아' in ln:
            rows[cur]['tie'] = True
    if not rows:
        sys.exit('표를 하나도 못 읽었다: ' + path)
    return rows, heads


def md5(path):
    return hashlib.md5(open(path, 'rb').read()).hexdigest()


# ── 칸 그리기 (틀의 마크업을 그대로 따른다) ─────────────────────
def alpha(v):
    return min(max((v - 0.5) / 0.5, 0.0), 1.0)


THEAD = ('<thead><tr><th scope="col">조합</th><th scope="col">모델</th>'
         + ''.join('<th scope="col" class="num">%s</th>' % r for r in RATIOS)
         + '</tr></thead>')


def ratio_table(rows, gen, panel, depth, unit='사람'):
    body = []
    for c in COMBOS:
        r = rows.get(('j%s%s' % (gen, panel), depth, c))
        if r is None:
            continue
        vals = r[unit]
        first = next((i for i, v in enumerate(vals) if v >= 0.7), None)
        tds = []
        for i, v in enumerate(vals):
            cls = 'num' + (' lo' if v < 0.5 else '') + (' hit' if i == first else '')
            tds.append('<td class="%s" style="--a:%.2f">%.3f</td>' % (cls, alpha(v), v))
        body.append('<tr><th scope="row">%s%s</th><td class="mdl">%s</td>%s</tr>'
                    % (c, DAG if r['tie'] else '', r['model'], ''.join(tds)))
    return '<table>%s<tbody>%s</tbody></table>' % (THEAD, ''.join(body))


def figure(rows, heads, gen, panel, depth):
    h = heads.get(('j%s%s' % (gen, panel), depth))
    if h is None:
        return ('<figure class="panel"><figcaption>%s <span class="tag">%s</span></figcaption>'
                '<p class="nline">이 깊이의 채점 파일이 없다</p></figure>' % (panel, TAG[panel]))
    return ('<figure class="panel"><figcaption>%s <span class="tag">%s</span></figcaption>'
            '<div class="tbl">%s</div><p class="nline">채점 %d/%d · 사람 단위 집계 '
            '(한 사람의 반복 검체를 묶음)</p></figure>'
            % (panel, TAG[panel], ratio_table(rows, gen, panel, depth),
               h['scored'], h['total']))


def trio(rows, heads, gen, depth, panels):
    return ('<div class="trio">%s</div>'
            % ''.join(figure(rows, heads, gen, p, depth) for p in panels))


# ── 0장 1차 지표 ────────────────────────────────────────────
def parse_primary(path):
    """표 끝의 「=== 1차 지표」 블록 (4_report 의 primary_metric.py 출력).
    음성은 모든 비율의 음성 행을 합친다 — 사전등록 계산법. 비율별 표 칸과 계산이 달라
    같은 칸이라도 값이 조금 다르다(GBM jsdb200 0.5254 대 0.526). 0장은 이 블록만 쓴다."""
    s = io.open(path, encoding='utf-8').read()
    i = s.find('=== 1차 지표')
    if i < 0:
        sys.exit('1차 지표 블록이 없다: ' + path)
    blk = s[i:]
    out = {}
    for m in re.finditer(r'^\s+(bl200|jsd200|jsdb200)\s+AUC\s+([\d.]+)\s+'
                         r'\(양성 (\d+)행 · 음성 (\d+)행\)', blk, re.M):
        out.setdefault(m.group(1), (float(m.group(2)), int(m.group(3)), int(m.group(4))))
    m = re.search(r'교집합 검체 (\d+)개 · 사람 (\d+)명', blk)
    out['_n'] = (int(m.group(1)), int(m.group(2))) if m else None
    out['_rollback'] = '되돌림' in blk and '발동' in blk[:blk.find('쌍별') if '쌍별' in blk else None]
    if not all(p in out for p in PANELS):
        sys.exit('1차 지표 블록에 세 판이 다 없다: %s %s' % (path, sorted(out)))
    return out


def section_primary(data):
    trs = []
    for name, pm, rows, gen in data:
        tds = []
        for p in PANELS:
            v, npos, nneg = pm[p]
            r = rows.get(('j%s%s' % (gen, p), '5k', PRIMARY))
            tds.append('<td class="num" style="--a:%.2f">%.4f<span class="m2">%s</span></td>'
                       % (alpha(v), v, r['model'] if r else ''))
        npos, nneg = pm[PANELS[0]][1:]
        below = all(pm[p][0] < 0.70 for p in PANELS)
        trs.append('<tr><th scope="row">%s</th>%s<td class="num">%s / %s</td><td>%s</td></tr>'
                   % (name, ''.join(tds), format(npos, ','), format(nneg, ','),
                      ('<b>미달</b> · 되돌림 발동' if pm['_rollback'] else '<b>미달</b>')
                      if below else '<b>넘은 판 있음</b>'))
    return NL.join([
        '<section id="primary">',
        '  <h2><span class="kind val">0</span>1차 지표 — 사전등록한 단 하나의 칸</h2>',
        '  <p class="sub">이 보고서에서 판정에 쓰는 값은 이 표뿐입니다. '
        '<b>검체 단위 · 깊이 5k · 암 비율 0.5% · 조합 ME+JSD+LLR+PDR</b>. '
        '양성은 0.5% 검체, 음성은 <b>모든 비율의 음성 검체를 합친 것</b>입니다(사전등록 계산법). '
        '모델은 학습 CV 로 고른 것(숫자 옆 회색)입니다. 세 판을 모두 싣습니다.</p>',
        '  <figure class="panel"><div class="tbl"><table><thead><tr><th scope="col">자료</th>'
        + ''.join('<th scope="col" class="num">%s</th>' % p for p in PANELS)
        + '<th scope="col" class="num">양성 / 음성 행</th>'
          '<th scope="col">보조 바닥 0.70</th></tr></thead><tbody>'
        + ''.join(trs) + '</tbody></table></div></figure>',
        '  <ul class="plain">',
        '    <li><b>판정.</b> 세 자료 모두 세 판이 0.70 에 못 미쳐, 파이프라인이 사전등록의 '
        '되돌림 발동을 출력했습니다. 사전등록의 주 판정인 <b>귀무 문턱 판정은 계산하지 '
        '않았습니다</b> — 0.70 은 보조 바닥입니다.</li>',
        '    <li><b>순위를 매기지 않습니다.</b> 되돌림이 발동하면 순위를 매기지 않는다고 결과 '
        '전에 정했습니다. 판·조합·깊이 사이에서 「가장 좋은 구성」을 고르지 않습니다. '
        '3~6장의 표에서 큰 값을 골라 읽지 마십시오.</li>',
        '    <li>3~6장의 비율별 표는 음성을 같은 비율의 음성만 쓰므로, 같은 칸이라도 이 표와 '
        '값이 조금 다릅니다.</li>',
        '  </ul>',
        '</section>'])


# ── 3장 학습 CV ─────────────────────────────────────────────
def section_dev_tables(dev):
    figs = []
    for p in PANELS:
        body = []
        for c in COMBOS:
            tds = []
            for d in DEPTHS:
                x = dev.get('j16%s|%s|%s' % (p, d, c))
                tds.append('<td class="num na">—</td>' if x is None else
                           '<td class="num" style="--a:%.2f">%.3f<span class="m2">%s</span></td>'
                           % (alpha(x[0]), x[0], x[1]))
            body.append('<tr><th scope="row">%s</th>%s</tr>' % (c, ''.join(tds)))
        figs.append('<figure class="panel"><figcaption>%s <span class="tag">%s</span></figcaption>'
                    '<div class="tbl"><table><thead><tr><th scope="col">조합</th>%s</tr></thead>'
                    '<tbody>%s</tbody></table></div></figure>'
                    % (p, TAG[p], ''.join('<th scope="col" class="num">%s</th>' % d
                                          for d in DEPTHS), ''.join(body)))
    return '<div class="trio">%s</div>' % ''.join(figs)


# ── 6장 비교 표·그림 ────────────────────────────────────────
def gbm_table(sets):
    head1 = ''.join('<th scope="col" class="grp" colspan="3">%s</th>' % n for n, _, _ in sets)
    head2 = ''.join('<th scope="col" class="num">%s</th>' % p for _ in sets for p in PANELS)
    body = []
    for c in COMBOS:
        tds = []
        for _, rows, gen in sets:
            for p in PANELS:
                r = rows.get(('j%s%s' % (gen, p), '50k', c))
                if r is None:
                    tds.append('<td class="num na">—</td>')
                else:
                    v = r['사람'][RATIOS.index('1%')]
                    tds.append('<td class="num" style="--a:%.2f">%.3f</td>' % (alpha(v), v))
        body.append('<tr><th scope="row">%s</th>%s</tr>' % (c, ''.join(tds)))
    return ('<figure class="panel"><figcaption>50k 리드 · 암 비율 1%% · 사람 단위 AUC</figcaption>'
            '<div class="tbl"><table><thead><tr><th></th>%s</tr><tr><th scope="col">조합</th>%s'
            '</tr></thead><tbody>%s</tbody></table></div></figure>' % (head1, head2, ''.join(body)))


def xpos(r):
    lo, hi = math.log10(0.1), math.log10(5.0)
    return 40 + (math.log10(float(r.rstrip('%'))) - lo) / (hi - lo) * 470


def ypos(v):
    # 틀의 격자선(0.5→190.7 · 0.6→155.3 · 0.7→120.0 · 1.0→14.0)에서 역산한 식.
    #   190.7 로 반올림한 상수를 쓰면 점 21개 중 2개가 0.1px 어긋난다.
    return 14.0 + (1.0 - v) * 1060.0 / 3.0


def chart_series(sets):
    """(이름, 값7, 클래스, 사각점여부) 세 줄을 그린다. 틀의 마크업과 좌표식을 따른다."""
    out, last = [], []
    for i, (name, vals) in enumerate(sets, start=1):
        pts = [(xpos(r), ypos(v)) for r, v in zip(RATIOS, vals)]
        out.append('<polyline class="line s%d%s" points="%s"/>'
                   % (i, ' dash' if i == 3 else '', ' '.join('%.1f,%.1f' % p for p in pts)))
        for (x, y), r, v in zip(pts, RATIOS, vals):
            if i == 3:
                dot = ('<rect class="dot s3" x="%.1f" y="%.1f" width="7" height="7"/>'
                       % (x - 3.5, y - 3.5))
            else:
                dot = ('<circle class="dot s%d" cx="%.1f" cy="%.1f" r="%s"/>'
                       % (i, x, y, '4' if i == 1 else '3.5'))
            out.append('<g class="pt" data-tip="%s · %s · AUC %.3f">%s'
                       '<circle class="hit" cx="%.1f" cy="%.1f" r="11"/></g>'
                       % (name, r, v, dot, x, y))
        last.append([pts[-1][1] + 3.5, i, name, vals[-1]])
    last.sort()
    for k in range(1, len(last)):                       # 겹치지 않게 12px 띄운다
        if last[k][0] - last[k - 1][0] < 12:
            last[k][0] = last[k - 1][0] + 12
    for y, i, name, v in last:
        out.append('<text class="dlab s%d" x="512" y="%.1f">%s %.3f</text>' % (i, y, name, v))
    return ''.join(out)


# ── 조립 ───────────────────────────────────────────────────
def replace_once(s, old, new, what):
    n = s.count(old)
    if n != 1:
        sys.exit('틀에서 「%s」 자리가 %d곳이다 (1곳이어야 한다)' % (what, n))
    return s.replace(old, new)


def between(s, start, end, what):
    i = s.find(start)
    j = s.find(end, i + len(start))
    if i < 0 or j < 0:
        sys.exit('틀에서 「%s」 구간을 못 찾았다' % what)
    return i, j


def reps_of(heads, gen):
    h = heads.get(('j%sbl200' % gen, '5k'))
    return h['per'] // h['scored'] if h and h['scored'] else None


BOUNDARY = '<p class="boundary">%s</p>'
CSS = ('.boundary{margin:2px 0 6px;padding:8px 12px;border-left:3px solid var(--warn);'
       'background:var(--warn-bg);color:var(--ink);font-size:13px;max-width:none;'
       'border-radius:4px} ul.plain{margin:10px 0 0;padding-left:18px;display:grid;gap:6px;'
       'max-width:80ch} #limits ol{margin:0;padding-left:20px;display:grid;gap:8px;max-width:80ch}')


def build(indir, regress):
    P = lambda f: os.path.join(indir, f)
    tpl = io.open(P('template.html'), encoding='utf-8').read()
    sil, silh = parse(P('silico.txt'))
    pat, path_ = parse(P('patient.txt'))
    g5, _ = parse(P('gbm.txt'))
    g50 = g5
    dev = json.load(io.open(P('dev_cv.json'), encoding='utf-8'))
    s = tpl

    if regress:
        # † 는 뺀다 — 틀의 † 기준은 기록이 없고, 여기서는 make_tables.py 의 규칙을 쓴다.
        #   그 밖의 모든 글자(값·모델·색·lo·hit·채점 문구)가 바이트 단위로 같아야 한다.
        strip = lambda x: re.sub(r' <span class="dag"[^>]*>†</span>', '', x)
        bad = 0
        for sec, rows, heads in (('silico', sil, silh), ('patient', pat, path_)):
            i, j = between(s, '<section id="%s">' % sec, '</section>', sec)
            old = strip(s[i:j])
            for d in DEPTHS:
                for p in ('bl200', 'jsd200'):
                    new = strip(figure(rows, heads, '16', p, d))
                    if new not in old:
                        bad += 1
                        print('  다름 %s %s %s' % (sec, d, p))
        print('회귀 검사(† 제외): bl200·jsd200 그림 12개 중 틀과 다른 것 %d' % bad)
        return None if bad else 'ok'

    rp, rs = reps_of(path_, '16'), reps_of(silh, '16')

    # 머리·목차·사양 글
    k = s.rfind('</style>')       # 틀에 style 이 둘이다 — 마지막(본문 CSS) 끝에 붙인다
    s = s[:k] + CSS + s[k:]
    s = replace_once(s, '<b>패널</b> bl200 · jsd200 (각 200블록)',
                     '<b>패널</b> bl200 · jsd200 · jsdb200 (각 200블록)', '머리 패널')
    s = replace_once(s, '<nav class="toc" aria-label="목차">',
                     '<nav class="toc" aria-label="목차"><a href="#primary">0. 1차 지표</a>', '목차')
    s = replace_once(s, '<a href="#gbm">6. GBM 비교</a></nav>',
                     '<a href="#gbm">6. GBM 비교</a><a href="#limits">7. 한계</a></nav>', '목차 끝')
    s = replace_once(s, '<dd>bl200 · jsd200 — 각 100bp 블록 200개</dd>',
                     '<dd>bl200 · jsd200 · jsdb200 — 각 100bp 블록 200개</dd>', '사양 패널')
    s = replace_once(s, 'bl200(통계 → 부트스트랩 순위)과 jsd200(JSD 순위)을 만들었습니다.',
                     'bl200(통계 → 부트스트랩 순위) · jsd200(JSD 순위) · jsdb200(JSD 부트스트랩 '
                     '빈도 · 귀무 문턱 통과)을 만들었습니다.', '방법 패널')
    s = replace_once(s, 'AUC는 <b>사람 단위</b>로 보고합니다 (한 사람의 반복 검체 점수를 묶어 집계).',
                     '3~6장 AUC는 <b>사람 단위</b>입니다 (한 사람의 반복 검체 점수를 묶어 집계). '
                     '0장 1차 지표만 사전등록대로 <b>검체 단위</b>입니다.', '방법 단위')
    s = replace_once(s, '<td>비율당 7회 — 공여자가 7명이라 그 이상은 독립이 아님</td>',
                     '<td>비율당 %d회 · 반복마다 배경 7명 <b>전원</b>에서 비례로 뽑는다 '
                     '(반복 사이 차이는 사람 차이가 아니라 뽑기 잡음)</td>' % rp, 'Patient 반복')
    s = replace_once(s, '반복 비율당 10회 · 배경 10명 전원 비례 혼합',
                     '반복 비율당 %d회 · 배경 10명 전원 비례 혼합' % rs, 'In silico 반복') \
        if '반복 비율당 10회 · 배경 10명 전원 비례 혼합' in s else s

    old_key = '흐린 글자 = 0.5 미만 · † = 동점 많아 순위 불가'
    if s.count(old_key) != 2:
        sys.exit('틀의 † 범례가 2곳이 아니다: %d' % s.count(old_key))
    s = s.replace(old_key, '흐린 글자 = 0.5 미만 · † = 검체 점수 동점으로 AUC 가 0.02 넘게 '
                           '정해지지 않는 칸이 있는 줄 (make_tables.py 규칙)')

    # 0장
    s = replace_once(s, '<section id="data">',
                     section_primary([('GBM · 조직 178', parse_primary(P('gbm.txt')), g5, '15'),
                                      ('폐 In silico · 조직 18', parse_primary(P('silico.txt')), sil, '16'),
                                      ('폐 Patient · cfDNA 48', parse_primary(P('patient.txt')), pat, '16')])
                     + NL + '<section id="data">', '0장 자리')

    # 3장
    i, j = between(s, '<section id="dev">', '</section>', 'dev')
    sec = s[i:j]
    k = sec.find('<div class="trio">')
    sec = (sec[:k].replace('</h2>', '</h2>' + NL + '  ' + BOUNDARY % (
        '<b>보조 지표 · 학습 교차검증.</b> 9개 비율(10%·100% 포함)을 <b>합쳐</b> 학습한 값이라 '
        '특정 농도의 값이 아닙니다. 저장된 모델 파일(<code>BEST.joblib</code>) 스스로 '
        '「5% 이상이 1.000 이라 평균이 올라간다. 보고에 그대로 쓰지 않는다」고 적어 두었습니다. '
        '4개 모델 중 최고를 고른 값이고, 패널을 고른 세포주로 다시 잰 값입니다. '
        '검증 값과 나란히 읽지 마십시오.'), 1)
           + section_dev_tables(dev) + NL)
    s = s[:i] + sec + s[j:]

    # 4·5장
    note = ('<b>보조 지표.</b> 사람 단위 · 깊이·비율·조합이 0장과 달라 1차 지표와 <b>비교할 수 '
            '없습니다</b>. 표 하나가 9조합 × 7비율 = 63칸이고 세 판 · 세 깊이를 합치면 567칸입니다 — '
            '큰 값을 골라 읽으면 고르는 것만으로 값이 오릅니다. 분모는 %s 입니다.')
    for sec_id, rows, heads, who in (('silico', sil, silh, 'In silico 18명'),
                                     ('patient', pat, path_, 'Patient 48명')):
        i, j = between(s, '<section id="%s">' % sec_id, '</section>', sec_id)
        sec = s[i:j]
        k = sec.find('<h3>')
        parts = [sec[:k].replace('</h2>', '</h2>' + NL + '  ' + BOUNDARY % (note % who), 1)]
        for d in DEPTHS:
            parts.append('  <h3>%s 리드</h3>' % d + NL + trio(rows, heads, '16', d, PANELS))
        s = s[:i] + NL.join(parts) + NL + s[j:]

    # 6장
    i, j = between(s, '<section id="gbm">', '</section>', 'gbm')
    sec = s[i:j]
    sec = sec.replace('</h2>', '</h2>' + NL + '  ' + BOUNDARY % (
        '<b>보조 지표.</b> 아래 그림의 구성(jsd200 · LLR · 50k)은 <b>예시</b>입니다 — '
        '판·조합 사이 순위로 고른 것으로 읽지 마십시오(0장).'), 1)
    a, b = sec.find('<svg'), sec.find('</svg>')
    svg = sec[a:b]
    keep = svg[:svg.find('<polyline')]
    series = [('GBM', g50.get(('j15jsd200', '50k', 'LLR'))),
              ('폐 In silico', sil.get(('j16jsd200', '50k', 'LLR'))),
              ('폐 Patient', pat.get(('j16jsd200', '50k', 'LLR')))]
    if any(r is None for _, r in series):
        sys.exit('6장 그림 줄이 없다: ' + str([n for n, r in series if r is None]))
    sec = sec[:a] + keep + chart_series([(n, r['사람']) for n, r in series]) + sec[b:]
    a = sec.find('<figure class="panel"><figcaption>50k 리드 · 암 비율 1%')
    b = sec.find('</figure>', a) + len('</figure>')
    sec = sec[:a] + gbm_table([('GBM · 조직 178검체', g50, '15'),
                               ('폐 In silico · 조직 18검체', sil, '16'),
                               ('폐 Patient · cfDNA 48검체', pat, '16')]) + sec[b:]
    s = s[:i] + sec + s[j:]

    # 7장 한계 + 출처
    stamp = ' · '.join('%s %s' % (f, md5(P(f))[:12]) for f in
                       ('silico.txt', 'patient.txt', 'gbm.txt', 'dev_cv.json'))
    limits = NL.join([
        '<section id="limits">',
        '  <h2><span class="kind val">7</span>한계</h2>',
        '  <ol>',
        '    <li><b>음성이 독립이 아닙니다.</b> Patient 의 음성은 전부 같은 정상 7명'
        '(<정상 7명>)에서 나왔고, In silico 배경 10명에 그 7명이 들어 있습니다. '
        '반복을 늘려도 독립 음성은 7명입니다.</li>',
        '    <li><b>반복은 배경 전원에서 비례로 뽑습니다.</b> 반복 사이 차이는 사람 간 차이가 '
        '아니라 뽑기 잡음이라, 특이도가 사람 간 변동을 과소평가합니다.</li>',
        '    <li><b>학습과 검증의 비율 축이 다릅니다.</b> 학습은 9개(0.1~100%), 검증은 7개'
        '(0.1~5%)입니다. 3장 값에는 10%·100% 가 들어 있습니다.</li>',
        '    <li><b>LLR 특징에 자기참조 누수가 있습니다.</b> 학습에 쓴 정상·세포주가 LLR 참조표에 '
        '그대로 들어갑니다. GBM 에서 크기를 쟀고(사전등록 추기 30~43), 폐는 같은 코드라 구조가 '
        '같습니다 — 크기는 재지 않았습니다. 학습 CV 를 떠받치는 것이 이 특징입니다.</li>',
        '    <li><b>패널을 고른 자료로 학습 CV 를 쟀습니다.</b> 3장 값에 해당하고, 4·5장 검증 값에는 '
        '해당하지 않습니다.</li>',
        '    <li><b>사전등록의 주 판정(귀무 문턱)은 계산하지 않았습니다.</b> 0장의 판정은 '
        '보조 바닥 0.70 으로 한 것입니다.</li>',
        '  </ol>',
        '  <p class="sub" style="margin-top:14px">생성 %s · build_lung_report.py · 입력 md5 %s</p>'
        % (datetime.datetime.now().strftime('%Y-%m-%d %H:%M'), stamp),
        '</section>'])
    i, _ = between(s, '<section id="gbm">', '</section>', 'gbm 끝')
    j = s.find('</section>', i) + len('</section>')
    s = s[:j] + NL + limits + s[j:]
    return s


def main():
    if len(sys.argv) < 3:
        sys.exit('쓰임: python build_lung_report.py <입력폴더> <출력.html> [--regress]')
    regress = '--regress' in sys.argv
    out = build(sys.argv[1], regress)
    if regress:
        sys.exit(0 if out else 1)
    io.open(sys.argv[2], 'w', encoding='utf-8', newline='\n').write(out)
    print('씀', os.path.abspath(sys.argv[2]), len(out))


if __name__ == '__main__':
    main()
