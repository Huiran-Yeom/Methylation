# LAYOUT=single 코호트용: 검증 사본(val03_*.py)의 짝 요구를 없앤다 (2026-08-11)
# 학습 스크립트는 건드리지 않는다.
import io, os

# detect_layout 정의가 두 번 있었다. 둘 다 같은 일을 하고 「뒤의 것이」
#   앞의 것을 덮고 있었다. 앞의 것은 죽은 코드다. 지워도 동작은 그대로다.
#   (지운 쪽 독스트링이 짧아, 아래 남긴 쪽이 더 설명이 낫다.)

def to_single(path):
    b = os.path.basename(path)
    s = io.open(path, encoding='utf-8').read()
    hits = []
    def rep(old, new, tag):
        nonlocal s
        if old in s:
            s = s.replace(old, new); hits.append(tag)
    if b == 'val03_1.py':
        rep("paired = [k for k, v in counts.items() if v == 2]",
            "paired = list(counts)   # 단일 리드: 짝 요구 없음", "짝요구제거")
    elif b == 'val03_2.py':
        rep("            reads[r.query_name].append(r)",
            "            reads['%s#%d' % (r.query_name, len(reads))] = [r]", "1리드=1단위")
        rep("        if len(rs) != 2:\n            drop += 1; continue\n", "", "len2필터제거")
        rep("        r = rs[0] if rs[0].reference_start <= rs[1].reference_start else rs[1]",
            "        r = rs[0]", "대표리드")
        rep("target = C.TARGET_DEPTH_READS // 2",
            "target = C.TARGET_DEPTH_READS   # 단일: 1리드=1단위", "목표깊이")
        rep("total*2:,} reads", "total:,} reads", "보고1")
        rep("'total_reads': total*2", "'total_reads': total", "보고2")
        # ── 균일 무작위판 (3·7·4) 전용
        rep("    paired = {k: v for k, v in reads.items() if len(v) == 2}",
            "    paired = dict(reads)   # 단일 리드: 짝 요구 없음", "균일:len2필터제거")
        rep("        base = C.TARGET_DEPTH_READS // 2",
            "        base = C.TARGET_DEPTH_READS   # 단일: 1리드=1단위", "균일:목표깊이")
        for _k in ('mut', 'wt'):
            rep("        assert written['%s'] == r['%s_pairs'] * 2, (written['%s'], r['%s_pairs'] * 2)" % (_k,_k,_k,_k),
                "        assert written['%s'] == r['%s_pairs'], (written['%s'], r['%s_pairs'])" % (_k,_k,_k,_k),
                "균일:쌍검사_" + _k)
    elif b == 'val03_3.py':
        rep("        '--paired-end',\n", "        '-s',\n", "단일모드")
    if hits:
        io.open(path, 'w', encoding='utf-8').write(s)
        print('  [single] %-11s %s' % (b, ' · '.join(hits)))
    elif b in ('val03_1.py', 'val03_2.py', 'val03_3.py'):
        print('  [single] %-11s !! 바꾼 것 없음. 코드가 다를 수 있음' % b)


def detect_layout(*bams):
    """종양·배경 중 하나라도 단일 리드면 'single'. 한 BAM 안에 두 방식이
    섞이면 bismark_methylation_extractor 가 처리할 수 없기 때문이다."""
    import pysam
    for p in bams:
        if not p or not os.path.exists(p):
            continue
        paired = False
        with pysam.AlignmentFile(p, 'rb') as b:
            for i, r in enumerate(b.fetch(until_eof=True)):
                if r.is_paired:
                    paired = True; break
                if i > 2000:
                    break
        if not paired:
            return 'single'
    return 'paired'
