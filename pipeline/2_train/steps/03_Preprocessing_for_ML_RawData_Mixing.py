#!/usr/bin/env python
# coding: utf-8
"""2_train/steps/03 (PE·SE 혼합판): 2026-08-20

m 판은 학습 자료가 전부 페어드엔드라 `--paired-end` 하나로 끝났다.
j 판은 정상 20명 중 7명(GSM 5 · SRR 2)이 단일말단이라 그 방법을 못 쓴다.

  --single-end 로 통일하면 안 되는 이유
    페어드엔드 짝이 겹친 구간(cfDNA 166bp 에 2x100 리드면 약 34bp)이 두 번 세어진다.
    커버리지가 비균일하게 부풀어, 우리가 계속 싸운 아티팩트가 하나 더 생긴다.

  그래서 섞인 BAM 을 flag 로 나눠 두 번 돌리고 cov 를 자리별로 합친다.
"""
# 2026-09-26: 09-23 사고에서 가드 없는 이 부류가 importlib 로 「확인」 되다 실제로 돌았다.
#   82개 파일이 덮였다(내용은 결정적이라 같았지만 그건 운이었다).
if __name__ != '__main__':
    raise ImportError(__file__ + ' 은(는) 스크립트다 — import 하지 않는다')

import gzip, os, subprocess, sys, time
import config as C

BATCH_SIZE = int(os.environ.get('BATCH_SIZE', 40))
SLICE_I    = int(os.environ.get('SLICE_I', 1))
SLICE_N    = int(os.environ.get('SLICE_N', 1))
BISMARK_MULTICORE = int(os.environ.get('BISMARK_MULTICORE', 1))
EXTRA_ARGS = []


def sh(cmd):
    subprocess.run(cmd, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.PIPE)


def split_pe_se(bam, tmp):
    """이름순으로 정렬한 뒤 짝이 온전한 것만 PE 로 보낸다.

    `samtools view -f 1` 만으로는 (가) 짝이 나란히 놓였는지 (나) 짝이 있기는 한지
    둘 다 보장되지 않는다. bismark PE 모드는 두 줄을 순서대로 읽어 이름이 다르면
    즉시 멈추므로, 여기서 직접 세어 나눈다.
    """
    _b = os.path.basename(bam)
    if _b.endswith('.bam'):
        _b = _b[:-4]
    ns = os.path.join(tmp, _b + '_ns.bam')
    sh(['samtools', 'sort', '-n', '-@', '4', '-m', '1G', '-o', ns, bam])

    pe_sam = os.path.join(tmp, _b + '_pe.sam')
    se_sam = os.path.join(tmp, _b + '_se.sam')
    pe = os.path.join(tmp, _b + '_pe.bam')
    se = os.path.join(tmp, _b + '_se.bam')

    # 페어드 관련 비트: paired · proper · mate_unmapped · mate_reverse · read1 · read2
    PAIRBITS = 0x1 | 0x2 | 0x8 | 0x20 | 0x40 | 0x80

    n_pe = n_se = 0
    p = subprocess.Popen(['samtools', 'view', '-h', ns],
                         stdout=subprocess.PIPE, universal_newlines=True)
    fp = open(pe_sam, 'w')
    fs = open(se_sam, 'w')

    def to_single(cols):
        """페어드 비트를 끄고 mate 정보를 지운다."""
        f = int(cols[1]) & ~PAIRBITS
        cols[1] = str(f)
        cols[6] = '*'   # RNEXT
        cols[7] = '0'   # PNEXT
        cols[8] = '0'   # TLEN
        return cols

    def flush(buf):
        nonlocal n_pe, n_se
        if not buf:
            return
        ok = False
        if len(buf) == 2:
            f0, f1 = int(buf[0][1]), int(buf[1][1])
            if (f0 & 0x1) and (f1 & 0x1) and ((f0 & 0x40) != (f1 & 0x40)):
                ok = True
        if ok:
            r1, r2 = (buf[0], buf[1]) if (int(buf[0][1]) & 0x40) else (buf[1], buf[0])
            fp.write('\t'.join(r1) + '\n')
            fp.write('\t'.join(r2) + '\n')
            n_pe += 2
        else:
            for c in buf:
                fs.write('\t'.join(to_single(list(c))) + '\n')
                n_se += 1

    buf = []
    last = None
    for line in p.stdout:
        line = line.rstrip('\n')
        if line.startswith('@'):
            fp.write(line + '\n')
            fs.write(line + '\n')
            continue
        cols = line.split('\t')
        if cols[0] != last:
            flush(buf)
            buf = []
        last = cols[0]
        buf.append(cols)
    flush(buf)
    p.stdout.close()
    p.wait()
    fp.close()
    fs.close()

    sh(['samtools', 'view', '-b', '-o', pe, pe_sam])
    sh(['samtools', 'view', '-b', '-o', se, se_sam])
    for f in (ns, pe_sam, se_sam):
        try:
            os.remove(f)
        except OSError:
            pass
    return pe, se


def build_cmd(bams, out_dir, paired, use_multicore, use_genome):
    cmd = ['bismark_methylation_extractor',
           '--paired-end' if paired else '--single-end',
           '--ignore_3prime', '1', *EXTRA_ARGS,
           '--bedGraph', '--gzip', '--output_dir', out_dir]
    if use_genome:    cmd += ['--genome_folder', C.GENOME_FOLDER]
    if use_multicore: cmd += ['--multicore', str(BISMARK_MULTICORE)]
    return cmd + bams


def read_cov(path):
    d = {}
    if not path or not os.path.exists(path):
        return d
    with gzip.open(path, 'rt') as f:
        for ln in f:
            c, s, e, _pct, m, u = ln.rstrip('\n').split('\t')
            k = (c, s, e)
            a = d.get(k)
            if a is None: d[k] = [int(m), int(u)]
            else:         a[0] += int(m); a[1] += int(u)
    return d


def merge_cov(pe_cov, se_cov, out_cov, out_bg):
    """두 cov 를 자리별로 더한다. 메틸·비메틸 개수를 합치고 비율을 다시 계산한다."""
    d = read_cov(pe_cov)
    for k, v in read_cov(se_cov).items():
        a = d.get(k)
        if a is None: d[k] = list(v)
        else:         a[0] += v[0]; a[1] += v[1]
    keys = sorted(d, key=lambda k: (k[0], int(k[1])))
    with gzip.open(out_cov, 'wt') as fc, gzip.open(out_bg, 'wt') as fb:
        for k in keys:
            m, u = d[k]
            pct = 100.0 * m / (m + u) if (m + u) else 0.0
            fc.write('%s\t%s\t%s\t%.6f\t%d\t%d\n' % (k[0], k[1], k[2], pct, m, u))
            fb.write('%s\t%d\t%s\t%.6f\n' % (k[0], int(k[1]) - 1, k[2], pct))
    return len(keys)


def main():
    C.guard_not_original()
    C.banner('2_train/steps/03  Bismark 메틸화 추출 (PE·SE 혼합)')
    C.ensure_dirs()
    use_genome = os.path.isdir(C.GENOME_FOLDER)
    use_multicore = BISMARK_MULTICORE > 1
    t0, done, skip = time.time(), 0, 0

    for p in C.MUT_READ_COUNTS:
        in_dir  = '%s/mut_%d_reads' % (C.MIX_OUT, p)
        out_dir = '%s/mut_%d_reads' % (C.BISMARK_OUT, p)
        if not os.path.isdir(in_dir):
            print('건너뜀: %s 없음' % in_dir); continue
        tmp = out_dir + '/_split'
        os.makedirs(out_dir, exist_ok=True); os.makedirs(tmp, exist_ok=True)

        todo = []
        for i in range(1, C.NUM_FILES + 1):
            if SLICE_N > 1 and (i - 1) % SLICE_N != (SLICE_I - 1): continue
            bam  = '%s/sampled_reads_%d.bam' % (in_dir, i)
            fin  = '%s/sampled_reads_%d.bismark.cov.gz' % (out_dir, i)
            if not os.path.exists(bam): continue
            if os.path.exists(fin) and os.path.getsize(fin) > 0: skip += 1; continue
            todo.append((i, bam))
        print('mut_%d_reads: 처리할 것 %d개' % (p, len(todo)))
        if not todo: continue

        parts = {}
        for i, bam in todo:
            parts[i] = split_pe_se(bam, tmp)
        pes = [v[0] for v in parts.values() if v[0]]
        ses = [v[1] for v in parts.values() if v[1]]
        print('  나눔: 페어드엔드 %d개 · 단일말단 %d개' % (len(pes), len(ses)))

        for lst, paired in ((pes, True), (ses, False)):
            for s in range(0, len(lst), BATCH_SIZE):
                b = lst[s:s + BATCH_SIZE]
                try:
                    sh(build_cmd(b, tmp, paired, use_multicore, use_genome))
                except subprocess.CalledProcessError as e:
                    msg = (e.stderr or b'').decode('utf-8', 'replace').strip().splitlines()
                    print('  ! %s 배치 %d 실패: %s'
                          % ('PE' if paired else 'SE', s, msg[-1] if msg else e))

        for i, _ in todo:
            pe, se = parts[i]
            cov = lambda f: (f[:-4] + '.bismark.cov.gz') if f else None
            n = merge_cov(cov(pe) and '%s/%s' % (tmp, os.path.basename(cov(pe))),
                          cov(se) and '%s/%s' % (tmp, os.path.basename(cov(se))),
                          '%s/sampled_reads_%d.bismark.cov.gz' % (out_dir, i),
                          '%s/sampled_reads_%d.bedGraph.gz' % (out_dir, i))
            if n: done += 1
        el = time.time() - t0
        print('  합침 %d개 · 누적 %d개 · %.1f분' % (len(todo), done, el / 60))

    print('\n처리 %d개 · 건너뜀 %d개 · %.1f분' % (done, skip, (time.time() - t0) / 60))
    print('2_train/steps/03 완료.')


if __name__ == '__main__':
    main()
