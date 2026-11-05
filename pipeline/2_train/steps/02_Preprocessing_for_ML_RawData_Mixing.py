#!/usr/bin/env python
# coding: utf-8
"""
2_train/steps/02: 비율마다 암, 정상 리드 섞어 BAM 만듦 (in-silico mixing)

DEPTH_MODE='total' 이면 샘플당 총 리드 TARGET_DEPTH_READS, 'per_block' 이면 블록 수 곱함
풀 모자라면 중단 (중복 뽑기 안 함), 명목/실제 비율은 계획표 CSV 로 저장
"""

# 스크립트 전용 - import 하면 본체가 돌아 자료를 덮어씀
if __name__ != '__main__':
    raise ImportError(__file__ + ' 은(는) 스크립트다 — import 하지 않는다')

import os
import random
import shutil
import sys
import tarfile
from collections import defaultdict

import config as C

# 명령줄 옵션
#   --plan-only        BAM 안 만들고 계획표만
#   --files 10         비율당 10개만 (시험용)
#   --ratios 0,0.001   특정 비율만
PLAN_ONLY = '--plan-only' in sys.argv


# 계획 계산
def avg_reads_per_frag(bam_path, names):
    """조각(이름)당 평균 리드 수 (페어드엔드 2, 단일말단 1)"""
    import pysam
    s = set(names)
    n = 0
    with pysam.AlignmentFile(bam_path, 'rb') as b:
        for r in b.fetch(until_eof=True):
            if r.query_name in s:
                n += 1
    return n / max(1, len(s))


def build_plan(n_mut_pairs, n_wt_pairs, total_blocks, avg_mut=2.0, avg_wt=2.0):
    # 목표는 리드 수 기준 (단일말단 섞이면 조각 x 2 가 리드 수가 아님)
    if C.DEPTH_MODE == 'total':
        T = C.TARGET_DEPTH_READS
    elif C.DEPTH_MODE == 'per_block':
        T = total_blocks * C.TARGET_DEPTH_READS
    else:
        raise ValueError(f"DEPTH_MODE 는 'total' 또는 'per_block' — 받은 값 {C.DEPTH_MODE}")
    base = int(round(T / max(1e-9, avg_wt)))   # 표시용 (정상 기준 조각 수)

    rows, problems = [], []
    for ratio in C.MUT_RATIOS:
        mut_reads_t = T * ratio
        mut_pairs = int(round(mut_reads_t / max(1e-9, avg_mut)))
        wt_pairs = int(round((T - mut_reads_t) / max(1e-9, avg_wt)))
        mut_reads = mut_pairs * avg_mut
        total_reads = mut_reads + wt_pairs * avg_wt
        actual = mut_reads / total_reads * 100 if total_reads else 0.0
        nominal = ratio * 100
        sm = max(0, mut_pairs - n_mut_pairs)
        sw = max(0, wt_pairs - n_wt_pairs)
        if sm or sw:
            problems.append(
                f'  {nominal:g}% : '
                + (f'암 {sm:,}쌍 부족(요청 {mut_pairs:,}/보유 {n_mut_pairs:,}) ' if sm else '')
                + (f'정상 {sw:,}쌍 부족(요청 {wt_pairs:,}/보유 {n_wt_pairs:,})' if sw else '')
            )
        rows.append(
            {
                'folder': C.folder_name(ratio),
                'ratio': ratio,
                'nominal_pct': nominal,
                'mut_pairs': mut_pairs,
                'wt_pairs': wt_pairs,
                'mut_reads': int(round(mut_reads)),
                'total_reads': total_reads,
                'actual_pct': actual,
                'rel_err_pct': None if nominal == 0 else (actual - nominal) / nominal * 100,
                'pool_short_mut': sm,
                'pool_short_wt': sw,
            }
        )
    return rows, base, problems


def print_plan(rows, base, n_mut, n_wt, total_blocks):
    print(f"DEPTH_MODE = '{C.DEPTH_MODE}'   TARGET_DEPTH_READS = {C.TARGET_DEPTH_READS:,}")
    if C.DEPTH_MODE == 'per_block':
        print(f'  블록 {total_blocks:,}개 x {C.TARGET_DEPTH_READS:,} / 2 = 샘플당 {base:,}쌍')
    else:
        print(f'  샘플당 {base:,}쌍 = {base*2:,} reads   (250430 Method 문서 스펙)')
    print(f'  풀 — 암 {n_mut:,}쌍 / 정상 {n_wt:,}쌍')
    print('-' * 74)
    print(f"{'폴더':>16} {'명목%':>8} {'암 쌍':>9} {'암 리드':>9} {'총 리드':>10} {'실제%':>8} {'오차':>7}")
    for r in rows:
        err = '-' if r['rel_err_pct'] is None else (
            '정확' if abs(r['rel_err_pct']) < 1e-9 else f"{r['rel_err_pct']:+.0f}%"
        )
        flag = '  <- 풀 부족' if (r['pool_short_mut'] or r['pool_short_wt']) else ''
        print(
            f"{r['folder']:>16} {r['nominal_pct']:>8.3f} {r['mut_pairs']:>9,} "
            f"{r['mut_reads']:>9,} {r['total_reads']:>10,} {r['actual_pct']:>8.3f} {err:>7}{flag}"
        )
    print()


# 실행
def load_pairs_exactly_two(bam_path):
    import pysam

    reads = defaultdict(list)
    with pysam.AlignmentFile(bam_path, 'rb') as bam:
        for r in bam.fetch(until_eof=True):
            reads[r.query_name].append(r)
    total = len(reads)
    paired = {k: v for k, v in reads.items() if len(v) in (1, 2)}
    print(
        f'  {os.path.basename(bam_path)}: 이름 {total:,} 중 짝 2개 {len(paired):,} '
        f'(제외 {total-len(paired):,})'
    )
    return paired


def load_names_exactly_two(bam_path):
    """load_pairs_exactly_two 와 같은 이름, 같은 순서. 리드 객체 대신 이름 → 리드 수만 들고 있음
    배경(NORMAL_BAMS) 있으면 병합 정상 풀은 이름 수 세는 데만 써서 이걸로 충분
    병합 풀이 크면 리드 객체를 다 들 때 메모리가 수 G 넘음"""
    import pysam

    cnt = defaultdict(int)
    with pysam.AlignmentFile(bam_path, 'rb') as bam:
        for r in bam.fetch(until_eof=True):
            cnt[r.query_name] += 1
    total = len(cnt)
    paired = {k: v for k, v in cnt.items() if v in (1, 2)}
    print(
        f'  {os.path.basename(bam_path)}: 이름 {total:,} 중 짝 2개 {len(paired):,} '
        f'(제외 {total-len(paired):,})'
    )
    return paired


def count_blocks(paths):
    import pysam

    blocks = set()
    for p in paths:
        with pysam.AlignmentFile(p, 'rb') as bam:
            for r in bam.fetch(until_eof=True):
                if not r.is_unmapped:
                    blocks.add(
                        (r.reference_name, (r.reference_start // C.BLOCK_SIZE) * C.BLOCK_SIZE)
                    )
    return len(blocks)


def sample_names(names, k, label, rng=random):
    if k == 0:
        return []
    if k <= len(names):
        return rng.sample(names, k)
    if not C.ALLOW_REPLACEMENT:
        raise RuntimeError(
            f'{label} 풀 부족: {k:,}쌍 요청 / {len(names):,}쌍 보유 ({k/len(names):.1f}배).\n'
            f'  ALLOW_REPLACEMENT=False 이므로 중단합니다.\n'
            f'  해결: DEPTH_MODE/TARGET_DEPTH_READS 를 낮추거나 여러 샘플을 pooling 하세요.'
        )
    print(f'  ! 경고: {label} 중복 허용 추출 ({k/len(names):.1f}배)')
    return rng.choices(names, k=k)


def _cli():
    """명령줄 인자로 config 값 임시 변경 (시험 실행용)"""
    a = sys.argv
    if '--files' in a:
        C.NUM_FILES = int(a[a.index('--files') + 1])
        print(f'[CLI] NUM_FILES = {C.NUM_FILES}')
    if '--ratios' in a:
        C.MUT_RATIOS = [float(x) for x in a[a.index('--ratios') + 1].split(',')]
        C.MUT_READ_COUNTS = [int(r * C.TARGET_DEPTH_READS) for r in C.MUT_RATIOS]
        print(f'[CLI] MUT_RATIOS = {C.MUT_RATIOS}')


def main():
    """비율마다 암, 정상 조각 섞어 BAM 씀

    입력 (전부 config)
      MUT_SORTED, WT_SORTED   암, 정상 풀 BAM
      NORMAL_BAMS             배경(사람)별 BAM. 없으면 병합본 사용
      MUT_RATIOS, NUM_FILES   비율 목록, 비율당 반복 수
      SEED                    rng = Random(SEED*10^10 + 비율*5000*10^4 + i)

    출력
      MIX_OUT/mixing_plan_<모드>_<판>.csv     계획표
      MIX_OUT/mut_<n>_reads/sampled_reads_<i>.bam
      MIX_OUT/skipped_ratios_<판>.txt         풀 모자라 건너뛴 비율
      MIX_OUT/mut_<n>_reads/background_map_<판>.csv   반복별 배경 사용 기록

    검증도 이 파일을 복사해 씀 - rng 식 바꾸면 학습과 검증 섞기가 달라짐
    """
    import pysam  # 서버에만 있으면 돼서 여기서 import

    _cli()
    C.guard_not_original()
    C.banner('2_train/steps/02  In-silico mixing (depth 수정판)')
    C.ensure_dirs()
    random.seed(C.SEED)  # 전역 시드 (보조)

    for p in (C.MUT_SORTED, C.WT_SORTED):
        if not os.path.exists(p):
            raise SystemExit(f'입력 BAM 이 없습니다. 2_train/steps/01 을 먼저 실행하세요.\n  {p}')

    # 입력 읽기
    print('BAM 로딩...')
    mut_reads = load_pairs_exactly_two(C.MUT_SORTED)
    # 배경 있으면 병합 정상 풀에서 안 뽑으니 이름만 읽음 (메모리 절약)
    _have_bg = bool(getattr(C, 'NORMAL_BAMS', []) or [])
    wt_reads = (load_names_exactly_two if _have_bg else load_pairs_exactly_two)(C.WT_SORTED)
    # 배경(사람)별 정상 조각. 없으면(학습) 빈 목록 → 병합본 사용
    _bgnames = []
    _bgreads = []   # 배경별 리드 사전
    _bgavg = []     # 배경별 조각당 리드
    for _k in range(1, len(getattr(C, 'NORMAL_BAMS', []) or []) + 1):
        _f = C.WT_SORTED.replace('.bam', '_bg%02d.bam' % _k)
        if os.path.exists(_f):
            _d = load_pairs_exactly_two(_f)
            _bgnames.append(list(_d))
            _bgreads.append(_d)   # 리드도 배경에서 가져옴
            _bgavg.append(sum(len(v) for v in _d.values()) / max(len(_d), 1))
        else:
            # 배경 BAM 없으면 중단 - 그냥 넘어가면 병합본으로 조용히 돌아도 티가 안 남
            sys.exit('중단: 배경 %d 의 BAM 이 없습니다 — %s'
                     '   2_train/steps/01 이 배경별 BAM 을 만들었는지 보십시오.' % (_k, _f))
    if _bgnames:
        print('배경 %d개 · 조각 수 %s'
              % (len(_bgnames), ' · '.join(format(len(x), ',') for x in _bgnames)))
    mut_names, wt_names = list(mut_reads), list(wt_reads)
    total_blocks = count_blocks([C.MUT_SORTED, C.WT_SORTED])
    print(f'  실제 사용 블록 수 {total_blocks:,}\n')

    am = avg_reads_per_frag(C.MUT_SORTED, mut_names)
    aw = avg_reads_per_frag(C.WT_SORTED, wt_names)
    print('  조각당 평균 리드: 암 %.2f · 정상 %.2f' % (am, aw))
    # 계획 세우기
    rows, base, problems = build_plan(len(mut_names), len(wt_names), total_blocks, am, aw)
    print_plan(rows, base, len(mut_names), len(wt_names), total_blocks)

    plan_csv = f'{C.MIX_OUT}/mixing_plan_{C.DEPTH_MODE}_{C.VERSION}.csv'
    cols = list(rows[0].keys())
    with open(plan_csv, 'w', encoding='utf-8') as f:
        f.write(','.join(cols) + '\n')
        for r in rows:
            f.write(','.join('' if r[c] is None else str(r[c]) for c in cols) + '\n')
    print(f'계획표 저장: {plan_csv}')

    # 계획 점검 (풀 모자란 비율)
    if problems:
        print('\n' + '!' * 74)
        print('풀이 부족한 비율이 있습니다:')
        for p in problems:
            print(p)
        print('!' * 74)

    if PLAN_ONLY:
        print('\nPLAN_ONLY=True: BAM 을 만들지 않았습니다.')
        return
    if problems and getattr(C, 'SKIP_SHORT_RATIOS', False):
        # 모자란 비율만 빼고 진행 (중복 채우기 안 함, 전체 중단도 안 함)
        _short = [r['ratio'] for r in rows if r['pool_short_mut'] or r['pool_short_wt']]
        rows = [r for r in rows if not (r['pool_short_mut'] or r['pool_short_wt'])]
        os.makedirs(C.MIX_OUT, exist_ok=True)
        with open(f'{C.MIX_OUT}/skipped_ratios_{C.VERSION}.txt', 'w') as _h:
            _h.write('\n'.join(str(x) for x in _short) + '\n')
        C.MUT_RATIOS = [r['ratio'] for r in rows]
        C.MUT_READ_COUNTS = [int(x * C.TARGET_DEPTH_READS) for x in C.MUT_RATIOS]
        print('\n건너뛴 비율: ' + ', '.join(f'{x*100:g}%' for x in _short))
        print('  남은 비율 ' + str(len(rows)) + '개 — 목록은 skipped_ratios 파일에')
        problems = []
    if problems and not C.ALLOW_REPLACEMENT:
        raise SystemExit('풀 부족으로 중단했습니다.')

    # 암 헤더로 정상, 배경 리드도 쓰니까 참조 서열 목록이 같아야 함
    C._cfg.check_same_refs([C.MUT_SORTED, C.WT_SORTED]
                         + [f for f in (C.WT_SORTED.replace('.bam', '_bg%02d.bam' % k) for k in range(1, 100))
                            if os.path.exists(f)], '섞기 입력')
    # 출력 헤더 (리드쌍 연달아 쓰니 SO:unsorted, GO:query)
    with pysam.AlignmentFile(C.MUT_SORTED, 'rb') as b:
        hdr = b.header.to_dict()
    hdr['HD'] = {'VN': '1.6', 'SO': 'unsorted', 'GO': 'query'}
    header = pysam.AlignmentHeader.from_dict(hdr)

    # 본체: 비율 x 반복마다 뽑아 BAM 씀 (rng 식 바꾸면 학습, 검증 섞기가 달라짐)
    for r in rows:
        out_dir = f"{C.MIX_OUT}/{r['folder']}"
        txt_dir = f'{out_dir}/txt'
        os.makedirs(txt_dir, exist_ok=True)
        print(
            f"\n[{r['nominal_pct']:g}%] {r['folder']} — 암 {r['mut_pairs']:,}쌍 / "
            f"정상 {r['wt_pairs']:,}쌍 (실제 {r['actual_pct']:.3f}%)"
        )
        # 이어하기: BAM 다 있고 txt.tar.gz 까지 있으면 그 비율 건너뜀
        _n = sum(1 for _i in range(1, C.NUM_FILES + 1)
                 if os.path.exists(f'{out_dir}/sampled_reads_{_i}.bam')
                 and os.path.getsize(f'{out_dir}/sampled_reads_{_i}.bam') > 0)
        if _n >= C.NUM_FILES and os.path.exists(f'{out_dir}/txt.tar.gz'):
            print(f'    건너뜀: 이미 {_n}개 완성')
            shutil.rmtree(txt_dir, ignore_errors=True)
            continue
        _bgmap = []
        for i in range(1, C.NUM_FILES + 1):
            # 파일마다 독립 시드 - 일부만 다시 돌려도 같은 BAM 나옴
            # 비율 x5000 은 0.0002 미만 비율이면 0 과 겹침 (바꾸면 기존 결과와 달라짐)
            # 검증은 검체 이름도 시드에 넣음 (안 넣으면 검체마다 정상이 똑같아짐), 학습은 _smp=''
            _smp = (os.path.basename(os.getcwd())
                    if os.path.basename(__file__).startswith('val03_2') else '')
            # 학습 파일을 _val 경로에서 돌리면 중단 (학습 산출물이 바뀜)
            assert _smp or '_val' not in os.getcwd(), \
                '학습인데 _val 경로다. 학습 산출물이 바뀐다'
            _sbase = C.SEED * 10**10 + int(r['ratio'] * 5000) * 10**4 + i
            rng = random.Random(_sbase if not _smp else '%d|%s' % (_sbase, _smp))
            # 반복 하나를 배경 K명 전원에서 균등하게 뽑아 만듦 (한 사람으로는 깊이를 못 채움)
            # 몫 못 채운 사람 부족분은 여유 있는 사람에게 나눔, 그래도 모자라면 fallback 기록
            # 반복마다 같은 K명이라 반복 간 차이에 사람 간 차이는 안 들어감 (그래서 특이도가 사람 간 변동을 덜 반영함)
            _wt_reads_t = r['total_reads'] - r['mut_reads']
            if _bgnames:
                _K = len(_bgnames)
                _take, _short = [0] * _K, 0.0
                for _j in range(_K):
                    _awj = _bgavg[_j] if _bgavg[_j] > 0 else aw
                    _need = int(round((_wt_reads_t / float(_K)) / max(1e-9, _awj)))
                    _take[_j] = min(_need, len(_bgnames[_j]))
                    _short += (_need - _take[_j]) * _awj
                for _ in range(5):          # 부족분 재배분, 안 줄면 멈춤
                    if _short <= 0.5:
                        break
                    _sur = [(len(_bgnames[_j]) - _take[_j], _j) for _j in range(_K)]
                    _tot = sum(_s for _s, _ in _sur if _s > 0)
                    if _tot <= 0:
                        break
                    _moved = 0.0
                    for _s, _j in _sur:
                        if _s <= 0:
                            continue
                        _awj = _bgavg[_j] if _bgavg[_j] > 0 else aw
                        _add = min(_s, int(round(_short * (_s / float(_tot))
                                                 / max(1e-9, _awj))))
                        if _add <= 0:
                            continue
                        _take[_j] += _add
                        _moved += _add * _awj
                    if _moved <= 0:
                        break
                    _short -= _moved
                # 사람마다 몫만큼 뽑아 이어 붙임 (리드 이름은 사람끼리 안 겹침)
                s_wt, _usereads = [], {}
                for _j in range(_K):
                    if _take[_j] <= 0:
                        continue
                    s_wt.extend(sample_names(_bgnames[_j], _take[_j],
                                             '정상bg%02d' % (_j + 1), rng))
                    _usereads.update(_bgreads[_j])
                # 계획한 정상 조각 수 (아래 리드 수 검사에서 사용)
                _wp = sum(_take)
                _bgid = 0
                _bgmap.append({'ratio': r['ratio'], 'rep': i, 'bg': 0,
                               'bg_file': '혼합%d명' % sum(1 for _t in _take if _t > 0),
                               'fallback': int(_short > 0.5),
                               'take': '|'.join(str(_t) for _t in _take)})
                if _short > 0.5:
                    print('  ! 배경 부족: 비율 %g 반복 %d 에서 리드 %d개 모자람'
                          % (r['ratio'], i, int(_short)))
            else:
                # 학습: 배경별 BAM 없음, 병합본 사용
                _usereads = wt_reads
                _wp = int(round(_wt_reads_t / max(1e-9, aw)))
                s_wt = sample_names(wt_names, _wp, '정상', rng)
                _bgmap.append({'ratio': r['ratio'], 'rep': i, 'bg': 0,
                               'bg_file': '병합', 'fallback': 0, 'take': ''})
            # 배경 몫은 그 배경의 조각당 리드로 계산함 (병합본 평균 쓰면 최대 2배 어긋남)
            s_mut = sample_names(mut_names, r['mut_pairs'], '암', rng)
            blk = defaultdict(lambda: {'mut': 0, 'wt': 0})
            written = {'mut': 0, 'wt': 0}
            out_bam = f'{out_dir}/sampled_reads_{i}.bam'
            with pysam.AlignmentFile(out_bam, 'wb', header=header) as out:
                for label, names, src in (('wt', s_wt, _usereads),
                                          ('mut', s_mut, mut_reads)):
                    for name in names:
                        for read in src[name]:
                            out.write(read)
                            written[label] += 1
                            if not read.is_unmapped:
                                key = (
                                    read.reference_name,
                                    (read.reference_start // C.BLOCK_SIZE) * C.BLOCK_SIZE,
                                )
                                blk[key][label] += 1
            # 쓴 리드 수가 계획 근처인지 확인 (여유 10%)
            for _lab, _exp in (('wt', r['total_reads'] - r['mut_reads']),
                               ('mut', r['mut_reads'])):
                if _exp > 0:
                    # 최소 여유 6리드 - 조각 몇 개뿐인 낮은 비율에서 단일말단 하나로 걸리지 않게
                    _n = r.get('mut_pairs', 0) if _lab == 'mut' else _wp
                    _tol = max(0.10 * _exp, 6.0, 2.0 * (float(_n) ** 0.5))
                    assert abs(written[_lab] - _exp) <= _tol, (_lab, written[_lab], _exp, _tol)
                else:
                    assert written[_lab] == 0, (_lab, written[_lab])

            with open(f'{txt_dir}/block_read_counts_{i}.txt', 'w') as f:
                f.write('chr\tstart\tend\tmutant_count\twt_count\n')
                for chrom, start in sorted(
                    k for k, v in blk.items() if v['mut'] + v['wt'] >= C.MIN_READS_PER_BLOCK
                ):
                    c = blk[(chrom, start)]
                    f.write(f"{chrom}\t{start}\t{start+C.BLOCK_SIZE}\t{c['mut']}\t{c['wt']}\n")

            if i % 100 == 0 or i == C.NUM_FILES:
                print(f'    {i}/{C.NUM_FILES}')

        # 반복별 배경 사용 기록 저장
        if _bgmap:
            import pandas as _pd          # 여기서만 씀
            _bm = _pd.DataFrame(_bgmap)
            _bmp = f'{out_dir}/background_map_{C.VERSION}.csv'
            _bm.to_csv(_bmp, index=False)
            _fb = int(_bm['fallback'].sum())
            print('    background_map 저장 %s · 되돌린 반복 %d개' % (os.path.basename(_bmp), _fb))

        with tarfile.open(f'{out_dir}/txt.tar.gz', 'w:gz') as tar:
            tar.add(txt_dir, arcname='txt')
        shutil.rmtree(txt_dir)

    print('\n2_train/steps/02 완료.')


if __name__ == '__main__':
    sys.exit(main())
