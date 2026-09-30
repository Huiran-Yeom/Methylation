#!/usr/bin/env python
# coding: utf-8
"""
2_train/steps/02: In-silico mixing.  depth 정의를 250430 Method 문서 스펙으로 고친 판

원본 : Methylation/scripts/02_Preprocessing_for_ML_RawData_Mixing.py

무엇이 달라졌나
  1. depth  원본은 `baseline_total_pairs = total_blocks * target_depth // 2` 라서
     블록 수(약 658)가 곱해져 샘플당 1,645,000쌍(3,289,864 reads)이 만들어졌다.
     250430 Method 문서는 "각 샘플은 총 2,500쌍(5,000개)" 이라고 명시한다.
     `DEPTH_MODE='total'` 이 문서 스펙이고, 10개 비율 전부 문서 표와 일치함을 확인했다.
     (`per_block` 으로 두면 원본 동작을 재현한다)
  2. 풀 부족 시 멈춘다. 원본은 조용히 `random.choices`(중복 허용)로 넘어가
     같은 리드쌍을 평균 21번 복제했다. 복제는 통계적 정보량을 늘리지 않는다.
  3. 짝이 2개인 이름만 미리 걸러 뽑는다. 원본은 뽑은 뒤 건너뛰어 요청보다 적게
     기록됐다. 0.1%(2쌍)에서 한 쌍만 빠지면 비율이 50% 틀어진다.
  4. BAM 헤더의 정렬 표기를 사실대로 적는다(원본은 SO:coordinate 라 적고 정렬 안 함).
     리드쌍을 연달아 쓰므로 SO:unsorted / GO:query 가 맞다.
  5. 명목 비율과 실제 비율을 계획표 CSV 로 남긴다. Methods 에 그대로 쓸 수 있다.

바뀌지 않은 것
  - 소수점 내림(int), MUT_RATIOS 값, 폴더 이름 규칙 f'mut_{int(r*5000)}_reads'
  - random.seed(SEED), 비율당 NUM_FILES 개 생성, 블록 통계 txt 저장
"""

# 이 파일은 스크립트다. 다른 코드가 import 하면 본체가 그대로 돌아
# 실제 자료를 덮어쓴다. 실행은 `python <파일>` 로만 한다.
if __name__ != '__main__':
    raise ImportError(__file__ + ' 은(는) 스크립트다 — import 하지 않는다')

import os
import random
import shutil
import sys
import tarfile
from collections import defaultdict

import config as C

# True 면 BAM 을 만들지 않고 계획표만 출력한다. 처음 한 번은 True 로 확인 권장.
#   명령줄로도 지정 가능:  --plan-only        계획표만
#                          --files 10         비율당 10개만 (시험용)
#                          --ratios 0,0.001   특정 비율만
PLAN_ONLY = '--plan-only' in sys.argv


# =====================================================================
# 계획 계산
# =====================================================================
def avg_reads_per_frag(bam_path, names):
    """조각(이름) 하나가 실제로 리드 몇 개인지. 페어드엔드 2 · 단일말단 1."""
    import pysam
    s = set(names)
    n = 0
    with pysam.AlignmentFile(bam_path, 'rb') as b:
        for r in b.fetch(until_eof=True):
            if r.query_name in s:
                n += 1
    return n / max(1, len(s))


def build_plan(n_mut_pairs, n_wt_pairs, total_blocks, avg_mut=2.0, avg_wt=2.0):
    # 목표를 '조각 수' 가 아니라 '리드 수' 로 잡는다.
    #   단일말단이 섞이면 조각당 리드가 1~2 로 흔들려 '조각 x 2 = 리드' 가 깨진다.
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


# =====================================================================
# 실행
# =====================================================================
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
    """명령줄 인자로 config 를 임시 조정한다 (파일 수정 없이 시험 실행)."""
    a = sys.argv
    if '--files' in a:
        C.NUM_FILES = int(a[a.index('--files') + 1])
        print(f'[CLI] NUM_FILES = {C.NUM_FILES}')
    if '--ratios' in a:
        C.MUT_RATIOS = [float(x) for x in a[a.index('--ratios') + 1].split(',')]
        C.MUT_READ_COUNTS = [int(r * C.TARGET_DEPTH_READS) for r in C.MUT_RATIOS]
        print(f'[CLI] MUT_RATIOS = {C.MUT_RATIOS}')


def main():
    """비율마다 암·정상 조각을 섞어 BAM 을 쓴다. 272줄 · 구역 여섯.

    받는 것 (전부 config)
      MUT_SORTED · WT_SORTED   암·정상 풀 BAM (조각 이름으로 묶여 있다)
      NORMAL_BAMS              배경(공여자)별 BAM. 비면 병합본 하나로 간다
      MUT_RATIOS · NUM_FILES   비율 목록 · 비율당 복제본 수
      SEED                     rng = Random(SEED*10^10 + 비율*10^4 + i)

    내는 것
      MIX_OUT/mixing_plan_<모드>_<판>.csv   계획표 (비율마다 조각 수)
      BISMARK_OUT/mut_<n>_reads/sampled_reads_<i>.bam
      MIX_OUT/skipped_ratios_<판>.txt       풀이 모자라 건너뛴 비율
      MIX_OUT/background_map_<판>.csv              어느 배경을 몇 조각 썼나

    구역
      1) 준비        설정·가드·배너·폴더·시드
      2) 입력 읽기   BAM 로딩 · 배경별 조각 · 블록 수 · 조각당 리드
      3) 계획        build_plan -> 계획표 저장
      4) 계획 점검   모자란 비율: 건너뛰기 / 중단 / PLAN_ONLY
      5) 출력 헤더   정렬 상태를 unsorted 로 바꿔 쓴다
      6) 본체        비율 x 복제본 마다 뽑아 BAM 쓰기  <- 182줄. 여기가 덩어리다

    ※ 검증(run_one_sample.sh)이 이 파일을 검체 폴더로 복사해 돌린다.
      구역 6 의 rng 식을 바꾸면 학습과 검증이 서로 다른 섞기를 한다.
    """
    import pysam  # 지연 임포트: 서버에만 있으면 됨

    _cli()
    C.guard_not_original()
    C.banner('2_train/steps/02  In-silico mixing (depth 수정판)')
    C.ensure_dirs()
    random.seed(C.SEED)  # (전역 시드도 유지. 보조)

    for p in (C.MUT_SORTED, C.WT_SORTED):
        if not os.path.exists(p):
            raise SystemExit(f'입력 BAM 이 없습니다. 2_train/steps/01 을 먼저 실행하세요.\n  {p}')

    # ── 2) 입력 읽기 - BAM · 배경별 조각 · 블록 수 ────────────────────────────────────────
    print('BAM 로딩...')
    mut_reads = load_pairs_exactly_two(C.MUT_SORTED)
    wt_reads = load_pairs_exactly_two(C.WT_SORTED)
    # 배경(공여자)별 정상 조각. 반복 i 가 배경 i 를 쓴다.
    #   없으면(=학습) 빈 목록이라 아래에서 병합본으로 간다.
    _bgnames = []
    _bgreads = []   # 2026-09-14 배경별 리드 사전
    _bgavg = []     # 2026-09-14 배경별 조각당 리드 (EGAD 1.99 · G3 1.00)
    for _k in range(1, len(getattr(C, 'NORMAL_BAMS', []) or []) + 1):
        _f = C.WT_SORTED.replace('.bam', '_bg%02d.bam' % _k)
        if os.path.exists(_f):
            _d = load_pairs_exactly_two(_f)
            _bgnames.append(list(_d))
            _bgreads.append(_d)   # 2026-09-14 리드까지 배경에서 가져온다
            _bgavg.append(sum(len(v) for v in _d.values()) / max(len(_d), 1))
        else:
            # [관문] 선언한 배경 BAM 이 없으면 멈춘다.
            #   print 만 하면 그 반복이 조용히 병합본으로 돌아가고 fallback 도
            #   0 으로 남는다(_bgid 가 0 이라). 178검체가 전부 합본 음성으로 돌아도
            #   표에서는 정상으로 보인다. 설정 잘못이지 자료 한계가 아니다.
            #   조각이 모자라 되돌리는 것은 정상 경로이고 그건 fallback 로 남는다.
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
    # ── 3) 계획 세우기 ────────────────────────────────────────
    rows, base, problems = build_plan(len(mut_names), len(wt_names), total_blocks, am, aw)
    print_plan(rows, base, len(mut_names), len(wt_names), total_blocks)

    plan_csv = f'{C.MIX_OUT}/mixing_plan_{C.DEPTH_MODE}_{C.VERSION}.csv'
    cols = list(rows[0].keys())
    with open(plan_csv, 'w', encoding='utf-8') as f:
        f.write(','.join(cols) + '\n')
        for r in rows:
            f.write(','.join('' if r[c] is None else str(r[c]) for c in cols) + '\n')
    print(f'계획표 저장: {plan_csv}')

    # ── 4) 계획 점검 - 모자란 비율을 어떻게 할까 ────────────────────────────────────────
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
        # 부족한 비율만 빼고 나머지는 정상 진행한다.
        #   중복으로 채우면 복제본이 독립이 아니게 되고,
        #   검체 전체를 중단하면 멀쩡한 낮은 비율까지 잃는다.
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

    # ── 5) 출력 헤더 ────────────────────────────────────────
    with pysam.AlignmentFile(C.MUT_SORTED, 'rb') as b:
        hdr = b.header.to_dict()
    hdr['HD'] = {'VN': '1.6', 'SO': 'unsorted', 'GO': 'query'}
    header = pysam.AlignmentHeader.from_dict(hdr)

    # ── 6) 본체: 비율 × 복제본 마다 뽑아 BAM 을 쓴다 ────────────────────
    #        182줄. 함수로 빼낼 첫 후보다. 빼낼 때 rng 식을 건드리면
    #        학습과 검증이 서로 다른 섞기를 한다. 재현 시험으로 확인할 것.
    for r in rows:
        out_dir = f"{C.MIX_OUT}/{r['folder']}"
        txt_dir = f'{out_dir}/txt'
        os.makedirs(txt_dir, exist_ok=True)
        print(
            f"\n[{r['nominal_pct']:g}%] {r['folder']} — 암 {r['mut_pairs']:,}쌍 / "
            f"정상 {r['wt_pairs']:,}쌍 (실제 {r['actual_pct']:.3f}%)"
        )
        # 이어하기: 이미 완성된 비율은 통째로 건너뛴다 (일치판 135줄과 같은 취지).
        # txt.tar.gz 까지 있어야 완성으로 본다. BAM 만 있고 tar 가 없으면 다시 만든다.
        _n = sum(1 for _i in range(1, C.NUM_FILES + 1)
                 if os.path.exists(f'{out_dir}/sampled_reads_{_i}.bam')
                 and os.path.getsize(f'{out_dir}/sampled_reads_{_i}.bam') > 0)
        if _n >= C.NUM_FILES and os.path.exists(f'{out_dir}/txt.tar.gz'):
            print(f'    건너뜀: 이미 {_n}개 완성')
            shutil.rmtree(txt_dir, ignore_errors=True)
            continue
        _bgmap = []
        for i in range(1, C.NUM_FILES + 1):
            # 파일별 독립 시드: 부분 실행·재시작·순서 변경에도 같은 BAM 이 나온다.
            # 전역 시드 1회 방식은 앞 파일의 소비량에 결과가 달려 재현성이 깨진다.
            # 파일마다 독립 시드: 원본은 파일 맨 위 random.seed(42) 하나뿐이라
            # 중단·재시작·병렬 실행에서 결과가 달라졌다.
            # ×5000 은 원본 폴더 이름 규칙(mut_{ratio*5000}_reads)에서 온 것이고,
            # 지금 쓰는 10개 비율에서는 0·5·25·50·100·125·150·250·500·5000 으로
            # 전부 구분된다. 0.0001 처럼 촘촘한 비율을 추가하면 int() 에서
            # 0 으로 뭉개져 0% 와 시드가 겹친다. 그때는 v5/v6 처럼 ×10**6 으로
            # 바꿔야 하며, 바꾸면 기존 파일과 재현이 깨지므로 전체 재실행이 필요하다.
            # [시드] 검체 식별자를 넣는다.
            #   옛 시드는 (SEED, 비율, 반복) 뿐이라 검체가 축에 없었다.
            #   그래서 178검체의 음성이 리드 수준에서 100% 같았다(교집합 2,757/2,757).
            #   점수에서도 확인된다. 음성 178검체의 검체평균 SD 가 llr 에서 0.0001 이다.
            #   학습(2_train/steps/02)은 검체 개념이 없으므로 _smp 가 빈 문자열이고 옛 시드 그대로다.
            #   → 학습 산출물은 비트 단위로 안 바뀐다. 검증만 바뀐다.
            _smp = (os.path.basename(os.getcwd())
                    if os.path.basename(__file__).startswith('val03_2') else '')
            # [대칭] 검증 쪽에는 「검증인데 _smp 없음」 assert 가 있다.
            #   반대(학습인데 _smp 있음)는 막혀 있지 않았다. 학습을 경로에 _val/ 이 낀
            #   자리에서 한 번이라도 돌리면 학습 산출물이 조용히 바뀐다. 제일 비싼 사고다.
            # [정정] 앞의 조건은 구조상 불가능했다 — _smp 는 파일명이
            #   val03_2* 일 때만 채워지므로 학습에서는 늘 '' 이고 검증에서는 늘 찬다.
            #   그래서 이 가드는 학습을 지키지 못하고 검증만 통째로 죽였다.
            #   (09-22 14:51 추가 · 09-23 에 LODO 첫 검체가 여기서 멈춰 드러났다.
            #    바로 다음 줄의 검증 시드 경로도 같이 도달 불가였다.)
            #   막으려던 진짜 사고는 「학습 파일을 _val 경로에서 돌리는 것」 이다.
            assert _smp or '_val' not in os.getcwd(), \
                '학습인데 _val 경로다. 학습 산출물이 바뀐다'
            _sbase = C.SEED * 10**10 + int(r['ratio'] * 5000) * 10**4 + i
            rng = random.Random(_sbase if not _smp else '%d|%s' % (_sbase, _smp))
            # 반복 i 는 배경 (i-1)%K 를 쓴다 — 음성이 서로 다른 사람에서 나온다.
            #   조각이 모자라면 병합본으로 되돌리고 그 사실을 기록한다.
            #   (깊은 판에서는 공여자 하나가 목표 깊이를 못 채울 수 있다.)
            # [혼합] 반복 하나를 전원에서 비율대로 뽑아 만든다.
            #   옛 방식 `_j = (i-1) % K` 는 한 반복 = 한 사람이었다. 그래서 목표 깊이를
            #   혼자 못 채우는 사람은 통째로 빠졌다. 50k 를 혼자 채울 GSE 공여자가
            #   하나도 없어(최대 bg07 49,599 < 50,000) EGAD 3명만 남았다.
            #   새 방식: 사람마다 균등 몫(리드 기준 총정상리드/K)을 뽑고, 몫을 못 채우는
            #   사람은 있는 만큼만 내고 부족분을 여유 있는 사람에게 비례 배분한다.
            #
            #   [대가 · 사전등록에 적었다] 반복 10개가 전부 같은 K명을 포함하므로
            #   반복 사이 차이가 「사람 간 차이」가 아니라 「뽑기 잡음」이 된다.
            #   즉 특이도가 사람 간 변동을 과소평가한다. 알고 쓰는 것이다.
            #
            #   사람별로 실제 몇 조각을 냈는지 background_map에 take 로 남긴다 —
            #   안 남기면 나중에 "정말 K명이었나" 를 셀 수 없다.
            _wt_reads_t = r['total_reads'] - r['mut_reads']
            if _bgnames:
                _K = len(_bgnames)
                _take, _short = [0] * _K, 0.0
                for _j in range(_K):
                    _awj = _bgavg[_j] if _bgavg[_j] > 0 else aw
                    _need = int(round((_wt_reads_t / float(_K)) / max(1e-9, _awj)))
                    _take[_j] = min(_need, len(_bgnames[_j]))
                    _short += (_need - _take[_j]) * _awj
                for _ in range(5):          # 부족분 재배분. 잔여가 안 줄면 멈춘다.
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
                # 사람마다 그 몫만큼 뽑아 이어 붙인다.
                #   리드 이름은 공여자끼리 0.00% 겹친다(5명 쌍 전수 확인): 합쳐도 안 덮어쓴다.
                s_wt, _usereads = [], {}
                for _j in range(_K):
                    if _take[_j] <= 0:
                        continue
                    s_wt.extend(sample_names(_bgnames[_j], _take[_j],
                                             '정상bg%02d' % (_j + 1), rng))
                    _usereads.update(_bgreads[_j])
                # 아래 리드 수 검사가 쓰는 「계획한 정상 조각 수」.
                #   혼합에서는 사람별 몫의 합이다. 이 줄이 없으면 411행에서 NameError.
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
                # 학습 경로: 배경별 BAM 이 없다(NORMAL_BAMS 미설정). 병합본 그대로.
                _usereads = wt_reads
                _wp = int(round(_wt_reads_t / max(1e-9, aw)))
                s_wt = sample_names(wt_names, _wp, '정상', rng)
                _bgmap.append({'ratio': r['ratio'], 'rep': i, 'bg': 0,
                               'bg_file': '병합', 'fallback': 0, 'take': ''})
            # 계획은 병합본 평균(약 1.7)으로 세웠는데 실제 배경은
            #   EGAD 1.99(PE 99%) · G3 1.00(SE 100%) 이다. 그대로 쓰면 리드가
            #   최대 2배 어긋난다. 그 배경의 조각당 리드로 다시 계산한다.
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
            # 조각당 리드가 1~2 로 섞이므로 '조각 x 2' 로 검사하면 안 된다.
            #   실제로 쓴 리드가 계획한 리드 수 근처인지로 본다.
            #   조각 3,416개면 리드 수 표준편차가 약 24 (5,000 의 0.5%). 여유 10%.
            for _lab, _exp in (('wt', r['total_reads'] - r['mut_reads']),
                               ('mut', r['mut_reads'])):
                if _exp > 0:
                    # 절대 하한 6리드: 0.1% 판은 암 조각이 3개뿐이라 그중 하나가
                    # 단일말단이면 6 대신 5가 된다. 리드 1개 차이인데 비율로는 17%다.
                    # 이 검사의 목적은 '두 배로 썼다' 같은 큰 사고를 잡는 것이다.
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

        # 어느 반복본이 어느 공여자 배경에서 나왔는지 남긴다.
        #   이게 없으면 나중에 음성을 사람 단위로 묶을 수 없다.
        if _bgmap:
            import pandas as _pd          # 2026-09-14: 이 파일엔 pandas 가 없다. 지역 임포트.
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
