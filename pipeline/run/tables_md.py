# -*- coding: utf-8 -*-
"""4단계 표(make_report.sh 텍스트 출력)를 Markdown 으로 옮김

사용법: python tables_md.py <..._tables.txt> [> ..._tables.md]
  - 판별 표('== Nk 리드' 블록)는 Markdown 표로, 1차 지표 절('=== 1차 지표' 부터)은 원문 그대로 고정폭 블록에
  - 숫자는 안 바꿈. 끝에서 원문과 결과의 숫자 열을 맞대 보고 다르면 중단 (assert)
"""
import io, re, sys

NUM = re.compile(r'[-+]?\d+(?:[.,]\d+)*%?')
VAL = re.compile(r'\d\.\d{3}|·')          # 판별 표 칸: 소수 셋째 자리 값, 또는 값 없음 표시
UNITS = ('검체', '사람', '↑동점폭')


def convert(text):
    out, lines = [], text.splitlines()
    i, n, missing = 0, len(lines), []

    def flush_missing():
        if missing:
            out.append('- 채점 파일 없음: ' + ' · '.join(missing))
            out.append('')
            missing.clear()

    while i < n:
        ln = lines[i]
        if ln.startswith('=== 1차 지표'):                 # 1차 지표 절은 원문 그대로
            flush_missing()
            out.append('## ' + ln.strip('= ').strip())
            out.append('')
            out.append('```text')
            out.extend(lines[i + 1:])
            out.append('```')
            break
        if ln.startswith('==='):
            flush_missing()
            out.append('## ' + ln.strip('= ').strip()); out.append('')
        elif ln.startswith('#####'):
            flush_missing()
            out.append('### ' + ln.lstrip('#').strip()); out.append('')
        elif re.match(r'== \S+ 리드: ', ln):
            missing.append(ln[3:].split(':')[0].strip())
        elif ln.startswith('== '):
            flush_missing()
            out.append('#### ' + ln[3:].strip()); out.append('')
            i += 1
            head = None
            while i < n and not lines[i].startswith('==') and not lines[i].startswith('#'):   # 다음 블록이나 판 머리에서 끝
                l = lines[i]
                if not l.strip():
                    i += 1; continue
                if l.startswith('피처'):
                    head = l.split()[2:]
                    if out and out[-1] != '':
                        out.append('')                      # 앞 설명 줄과 표 사이 빈 줄
                    out.append('| 피처 | 모델 | 단위 | ' + ' | '.join(head) + ' |')
                    out.append('|---|---|---|' + '---:|' * len(head))
                elif head is None:                           # 표 앞 설명 줄
                    out.append('- ' + l.strip())
                elif re.match(r'^\S', l) or re.match(r'^\s+[-\d.]', l):   # 값 줄 (새 피처, 또는 같은 피처의 다음 단위)
                    f = l.split() if re.match(r'^\S', l) else [''] * 2 + l.split()
                    vals, unit = f[2:2 + len(head)], ' '.join(f[2 + len(head):])
                    # 칸 수, 값 꼴, 단위 하나라도 어긋나면 열 밀린 것이라 중단
                    assert len(vals) == len(head) and all(VAL.fullmatch(v) for v in vals) and unit in UNITS,                         '%d행 값 줄 꼴이 다르다: %r' % (i + 1, l)
                    out.append('| %s | %s | %s | %s |' % (f[0], f[1], unit, ' | '.join(vals)))
                else:                                        # 표 안 주석 줄은 그 자리에 한 행으로 (순서 유지)
                    out.append('| | | ↳ %s |%s' % (l.strip().replace('|', '\\|'), ' |' * len(head)))
                i += 1
            out.append('')
            continue
        else:
            if ln.strip():
                flush_missing(); out.append(ln)
        i += 1
    flush_missing()
    return '\n'.join(out) + '\n'


def numbers(s):
    return NUM.findall(s)


if __name__ == '__main__':
    src = io.open(sys.argv[1], encoding='utf-8').read()
    md = convert(src)
    a, b = numbers(src), numbers(md)
    assert a == b, '숫자 열이 다르다: 원문 %d개 · 결과 %d개 · 첫 차이 %s' % (
        len(a), len(b), next(((k, x, y) for k, (x, y) in enumerate(zip(a, b)) if x != y), '길이'))
    sys.stdout.buffer.write(md.encode("utf-8"))
