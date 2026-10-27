# -*- coding: utf-8 -*-
"""학습 교차검증 전체 표(cv_all.csv)를 GitHub 에서 읽기 좋은 Markdown 으로 옮긴다.

쓰는 법: python cv_md.py <cv_all.csv> [> cv_all.md]
  - 패널 · 깊이마다 표 하나: 행 = 31조합, 열 = 혼합 비율, 칸 = 그 조합의 모델 넷 중 최고 cv_auc(괄호 안 모델 첫 글자).
  - CSV 가 원본이다. 이 표는 보기용이고 값은 CSV 에서 그대로 옮긴다(넷째 자리).
  - 끝에서 표의 모든 칸을 CSV 에서 다시 계산해 맞대고, 다르면 멈춘다(assert).
"""
import io, sys
import pandas as pd

ABBR = {'Logistic Regression': 'L', 'Random Forest': 'R', 'SVM': 'S', 'XGBoost': 'X'}
FEATS = ['mean', 'entropy', 'jsd', 'llr', 'pdr']


def order(fs):
    t = fs.split('|')
    return (len(t), [FEATS.index(x) for x in t])


def build(d):
    best = d.loc[d.groupby(['panel', 'depth', 'fraction_pct', 'feature_set'])['cv_auc'].idxmax()]
    out = ['# 학습 교차검증 — 조합마다 모델 넷 중 최고 cv_auc', '',
           '원본은 같은 폴더의 `cv_all.csv`(모든 모델 · 모든 칸)다. 괄호는 최고 모델: L 로지스틱 회귀 · R 랜덤 포레스트 · S SVM · X XGBoost(값이 같으면 CSV 에서 먼저 나온 모델).',
           '학습 CV 는 패널을 고른 같은 자료로 잰 값이라 성능 근거가 아니다. 판정은 외부검증 1차 지표로만 한다.', '']
    cells = {}
    dk = lambda s: int(s.rstrip('k'))
    for p in sorted(best['panel'].unique()):
        for dep in sorted(best.loc[best.panel == p, 'depth'].unique(), key=dk):
            b = best[(best.panel == p) & (best.depth == dep)]
            fr = sorted(b['fraction_pct'].unique())
            out.append('## %s · %s' % (p, dep)); out.append('')
            out.append('| 조합 | ' + ' | '.join('%g%%' % f for f in fr) + ' |')
            out.append('|---|' + '---:|' * len(fr))
            for fs in sorted(b['feature_set'].unique(), key=order):
                row = []
                for f in fr:
                    r = b[(b.feature_set == fs) & (b.fraction_pct == f)]
                    if r.empty:
                        row.append('—'); continue
                    v, m = float(r['cv_auc'].iloc[0]), ABBR.get(r['model'].iloc[0], r['model'].iloc[0])
                    cells[(p, dep, f, fs)] = v
                    row.append('%.4f (%s)' % (v, m))
                out.append('| %s | %s |' % (fs.replace('|', '\\|'), ' | '.join(row)))
            out.append('')
    return '\n'.join(out) + '\n', cells


def parse(md):
    """만든 Markdown 을 다시 읽어 (판, 깊이, 비율, 조합) → 표시값(문자열) 을 돌려준다(R23 Codex 5 · 출력식까지 검사)."""
    got, key, fr = {}, None, None
    for l in md.splitlines():
        if l.startswith('## '):
            p, dep = l[3:].split(' · '); key = (p, dep); fr = None
        elif l.startswith('| 조합 |'):
            fr = [float(x.strip().rstrip('%')) for x in l.strip('|').split('|')[1:]]
        elif fr and l.startswith('| ') and not l.startswith('|---'):
            c = [x.strip() for x in l.replace('\\|', '\x00').strip('|').split('|')]
            fs = c[0].replace('\x00', '|')
            assert len(c) == len(fr) + 1, '칸 수가 머리와 다르다: %r' % l
            for f, x in zip(fr, c[1:]):
                if x != '—':
                    got[key + (f, fs)] = x.split(' ')[0]
    return got


if __name__ == '__main__':
    d = pd.read_csv(sys.argv[1])
    md, cells = build(d)
    chk = d.groupby(['panel', 'depth', 'fraction_pct', 'feature_set'])['cv_auc'].max()
    got = parse(md)
    assert len(got) == len(chk), '칸 수가 다르다: 표 %d · CSV %d' % (len(got), len(chk))
    for k, v in chk.items():
        assert got.get(k) == '%.4f' % v, '값이 다르다: %s · CSV %.4f · 표 %s' % (k, v, got.get(k))
    sys.stdout.buffer.write(md.encode('utf-8'))
