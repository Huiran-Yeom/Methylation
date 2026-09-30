#!/bin/bash
# ── 자체 점검: 한 번 돌리면 표 하나가 나옵니다 ──────────────────────────
#
#   목적: 47개 파일을 읽지 않고 「제대로 되어 있나」 를 확인한다.
#   읽기만 합니다. 아무것도 만들지 않고 아무것도 고치지 않습니다.
#
#   사용
#     bash selfcheck.sh
#
#   PASS 가 전부면 「구조가 성립한다」 는 뜻입니다. 「결과가 맞다」 는 뜻은 아닙니다 —
#   그건 README 의 「Reproducibility」 표(md5 대조)가 답합니다.
set -u

_HERE=$(cd "$(dirname "${BASH_SOURCE[0]:-$0}")" && pwd)
cd "$_HERE" || exit 1

PASS=0; FAIL=0
ok()   { PASS=$((PASS+1)); printf '  PASS  %s\n' "$1"; }
bad()  { FAIL=$((FAIL+1)); printf '  FAIL  %s\n' "$1"; [ $# -gt 1 ] && printf '        %s\n' "$2"; }

echo "════ 1. 있어야 할 것이 있나 ════"
for f in 0_setup/config.conf 0_setup/meth_config.py \
         1_dmr/run_dmr.sh 2_train/run_train.sh \
         3_validate/run_validate.sh 4_report/make_report.sh \
         README.md requirements.txt; do
  [ -f "$f" ] && ok "$f" || bad "$f 가 없습니다"
done

echo
echo "════ 2. 단계 파일 수 ════"
#   2026-09-29: 「비었나」만 보면 11개 중 10개를 지워도 통과한다. 기대 개수를 적어 둔다.
#   단계를 더하거나 빼면 이 숫자도 같이 고친다. 그게 이 검사의 요점이다.
while IFS=: read -r d want; do
  n=$(find "$d" -maxdepth 1 -type f \( -name '*.py' -o -name '*.sh' \) 2>/dev/null | wc -l)
  if [ "$n" -eq "$want" ]; then ok "$d  $n개"
  else bad "$d  $n개 (기대 $want개)" "지웠거나 더했다면 selfcheck.sh 의 기대 개수도 고치십시오"; fi
done <<LIST
1_dmr/steps:16
2_train/steps:10
3_validate/steps:5
4_report/steps:2
LIST

echo
echo "════ 3. 문법 (전부) ════"
_py=$(command -v python || command -v python3 || true)
if [ -n "$_py" ]; then
  n=0; e=0
  while IFS= read -r f; do
    n=$((n+1))
    "$_py" -c "import ast,io,sys;ast.parse(io.open(sys.argv[1],encoding='utf-8').read())" "$f" 2>/dev/null \
      || { e=$((e+1)); echo "        문법 오류: $f"; }
  done < <(find . -name '*.py')
  [ "$e" -eq 0 ] && ok "파이썬 $n개 전부 통과" || bad "파이썬 $e/$n 개 문법 오류"
else
  bad "python 을 못 찾아 파이썬 문법을 못 봤습니다"
fi
n=0; e=0
while IFS= read -r f; do
  n=$((n+1)); bash -n "$f" 2>/dev/null || { e=$((e+1)); echo "        문법 오류: $f"; }
done < <(find . -name '*.sh')
[ "$e" -eq 0 ] && ok "셸 $n개 전부 통과" || bad "셸 $e/$n 개 문법 오류"

echo
echo "════ 4. 설정이 한 곳에서 읽히나 ════"
if [ -n "$_py" ]; then
  _out=$("$_py" - <<'PY' 2>&1
import sys, os
sys.path.insert(0, '0_setup')
try:
    import meth_config as C
except SystemExit as e:
    print('SETUPFAIL', e); raise SystemExit(0)
except Exception as e:
    print('SETUPFAIL', type(e).__name__, e); raise SystemExit(0)
keys = [k for k in open('0_setup/config.conf', encoding='utf-8').read().split('\n')
        if k.strip() and not k.strip().startswith('#') and '=' in k]
miss = [k.split('=')[0].strip() for k in keys if not hasattr(C, k.split('=')[0].strip())]
print('KEYS', len(keys))
print('EXPOSED_MISSING', ','.join(miss) if miss else '-')
print('VERSION', C.VERSION, 'GEN', C.GEN, 'DEPTH', C.TARGET_DEPTH_READS)
print('VERNAME', C.version_name(C.PANELS[0], '5k'))
PY
)
  case "$_out" in
    *SETUPFAIL*) bad "설정을 못 읽습니다" "$_out" ;;
    *) ok "설정 읽힘: $(echo "$_out" | grep '^VERSION')"
       _m=$(echo "$_out" | awk '/^EXPOSED_MISSING/{print $2}')
       [ "$_m" = "-" ] && ok "conf 키 전부 노출됨 ($(echo "$_out" | awk '/^KEYS/{print $2}')개)" \
                       || bad "노출 안 된 conf 키: $_m" ;;
  esac
fi

echo
echo "════ 5. 이름 규칙 ════"
_ko=$(find . -name '*[가-힣]*' 2>/dev/null | wc -l)
[ "$_ko" -eq 0 ] && ok "한글 파일·폴더 이름 0개" || bad "한글 이름 $_ko개"
_id=$(grep -rlP '^\s*[A-Za-z_]*[가-힣][A-Za-z0-9_가-힣]*\s*=' --include='*.py' . 2>/dev/null | wc -l)
[ "$_id" -eq 0 ] && ok "한글 파이썬 식별자 0개" || bad "한글 식별자가 있는 파일 $_id개"

echo
echo "════ 6. 서로를 부르는 이름이 실재하나 ════"
#   2026-09-29: $S/ · $_HERE/ 만 보던 것을 넓혔다. run_train.sh 와 make_report.sh 는
#   `cd "$S" && python -u <파일>` 꼴로 「맨 이름」 을 넘기고, run_one_sample.sh 는
#   "$TRAIN_SRC/<파일>" 을 쓴다. 셋 다 옛 규칙에 하나도 걸리지 않아, 단계 파일을
#   지워도 이 묶음이 통과했다. 파이썬으로 본다(한글 경로에서 grep 이 새기도 한다).
if [ -n "$_py" ]; then
  _out=$("$_py" - <<'HD6'
import io, os, re
REF = re.compile(r'''(?:\$(?:S|_HERE|TRAIN_SRC|ANALYSIS_DIR)/|python3?\s+-u\s+|python3?\s+|bash\s+|import\s+)["']?([A-Za-z0-9_][A-Za-z0-9_./-]*\.(?:py|sh))''')
have = {}
for r, ds, fs in os.walk('.'):
    ds[:] = [d for d in ds if d not in ('.git', '__pycache__')]
    for f in fs:
        have.setdefault(f, []).append(r.replace(os.sep, '/'))
# 검체 작업폴더에 놓이는 이름은 트리에 없다 (검증이 복사해 만든다)
WORK = re.compile(r'^(val0|config\.py$|_n\.py$)')
miss = []
for r, ds, fs in os.walk('.'):
    ds[:] = [d for d in ds if d not in ('.git', '__pycache__')]
    for f in sorted(fs):
        if not f.endswith(('.sh', '.py')): continue
        p = os.path.join(r, f).replace(os.sep, '/')
        for i, L in enumerate(io.open(p, encoding='utf-8', errors='replace'), 1):
            if L.lstrip().startswith('#'): continue
            for m in REF.finditer(L.split('#')[0]):
                n = os.path.basename(m.group(1))
                if n in have or WORK.match(n) or n == os.path.basename(p): continue
                miss.append('%s:%d -> %s' % (p, i, n))
miss = sorted(set(miss))
print('REFMISS %d' % len(miss))
for x in miss: print('   ', x)
HD6
)
  _n=$(echo "$_out" | awk '/^REFMISS/{print $2}')
  if [ "${_n:-1}" -eq 0 ]; then ok "드라이버·단계가 부르는 파일 전부 실재"
  else bad "없는 파일을 부르는 자리 ${_n}곳" "$(echo "$_out" | sed 1d)"; fi
fi

echo
echo "════ 7. 공개하면 안 되는 것이 없나 ════"
_pii=$(grep -rlE '\b(pat_[0-9]{3}_N[0-9]+|ind_[0-9]+_N[0-9]+)' . 2>/dev/null | wc -l)
[ "$_pii" -eq 0 ] && ok "검체 식별자 0건" || bad "검체 식별자가 있는 파일 $_pii개"
_conf=$(find . -path '*/cohorts/*' -name '*.conf' ! -name '*.example' 2>/dev/null | wc -l)
[ "$_conf" -eq 0 ] && ok "실제 코호트 정의 0개 (예시만)" || bad "실제 코호트 정의 $_conf개 — 빼십시오"

echo
echo "════ 8. 단일인용 heredoc 이 쓰는 변수가 export 되나 ════"
#   <<'X' 는 쓸 때 치환되지 않는다. 그 안의 스크립트는 「따로」 도니까 환경에서 찾는다.
#   밖에서 export 안 했고 내부가 set -u 면 그 자리에서 죽는다 — 출력이 /dev/null 이면 조용히.
#   2026-09-29: 정리 중 score_parallel.sh 의 _APPLY 가 정확히 이 꼴로 깨져 있었다.
if [ -n "$_py" ]; then
  _out=$("$_py" - <<'HD8'
import io, os, re
REF=re.compile(r'\$\{?([A-Za-z_][A-Za-z0-9_]*)\}?')
ENV={'PATH','HOME','PWD','USER','SHELL','LANG','IFS','TMPDIR','RANDOM','LINENO',
     'PYTHONPATH','OMP_NUM_THREADS','MKL_NUM_THREADS','OPENBLAS_NUM_THREADS'}
try:  # conf 의 키는 셸 머리글의 `set -a` 가 export 한다
    ENV |= set(L.split('=')[0].strip() for L in
               io.open('0_setup/config.conf', encoding='utf-8')
               if '=' in L and not L.strip().startswith('#'))
except Exception: pass
bad=[]
for r,ds,fs in os.walk('.'):
    ds[:]=[d for d in ds if d!='.git']
    for f in sorted(fs):
        if not f.endswith('.sh'): continue
        p=os.path.join(r,f)
        lines=io.open(p,encoding='utf-8',errors='replace').read().split(chr(10))
        exp=set(m.group(1) for L in lines
                for m in re.finditer(r'\bexport\s+([A-Za-z_][A-Za-z0-9_]*)', L))
        i=0
        while i<len(lines):
            if lines[i].lstrip().startswith('#'): i+=1; continue   # 주석의 <<'X' 는 heredoc 이 아니다
            _L=lines[i].split('#')[0]          # 줄 끝 주석의 <<'X' 도 heredoc 이 아니다
            # 2026-09-29: <<'T' 만 보던 것을 넓혔다. <<"T" 와 <<\T 도 치환을 막으므로
            #   같은 결함이 나고, 옛 규칙에는 걸리지 않았다.
            _BS = chr(92); _Q = chr(34); _A = chr(39)
            _NM = '([A-Za-z_]' + _BS + 'w*)'
            m=re.search('<<-?' + _BS + 's*(?:' + _A + _NM + _A + '|' + _Q + _NM + _Q
                        + '|' + _BS + _BS + _NM + ')', _L)
            if not m: i+=1; continue
            tag=m.group(1) or m.group(2) or m.group(3); body=[]; j=i+1
            while j<len(lines) and lines[j].strip()!=tag: body.append(lines[j]); j+=1
            if j>=len(lines): i+=1; continue      # 닫히는 줄이 없다 = 헛맞음
            inner=set()
            for L in body:
                inner |= set(mm.group(1) for mm in re.finditer(
                    r'(?:^|[;&|(]|\bexport\s+|\bfor\s+)\s*([A-Za-z_][A-Za-z0-9_]*)(?:=|\s+in\b)', L))
            for L in body:
                for mm in REF.finditer(L):
                    n=mm.group(1)
                    if n in ENV or n in inner or n in exp: continue
                    bad.append('%s:%d  <<%s 안의 $%s' % (p.replace(os.sep,'/'), i+1, tag, n))
            i=j+1
bad=sorted(set(bad))
print('HDBAD %d' % len(bad))
for b in bad: print('   ', b)
HD8
)
  _n=$(echo "$_out" | awk '/^HDBAD/{print $2}')
  if [ "${_n:-1}" -eq 0 ]; then ok "단일인용 heredoc 의 변수 전부 export 됨"
  else bad "export 안 된 변수 ${_n}건" "$(echo "$_out" | sed 1d)"; fi
fi
echo
echo "════ 9. 산출물 경로에 한글이 없나 ════"
#   코드 주석은 한글이다(의도). 그러나 「경로·파일명」 에 한글이 있으면 영어권 사용자가
#   만든 산출물을 읽을 수 없고, 저장소의 이름 규칙이 깨진다.
#   판정: 공백 없는 토큰에 한글 + (구분자 2개 이상 또는 아는 확장자). 따옴표 안팎을 다 본다.
if [ -n "$_py" ]; then
  _out=$("$_py" - <<'HD9'
import io, os, re
KO  = re.compile('[\uac00-\ud7a3]')
TOK = re.compile(r'[^\s\'"`,;()\[\]{}=|&<>]+')
EXT = ('.log','.lock','.sh','.py','.csv','.txt','.conf','.parquet','.joblib','.bam',
       '.tsv','.json','.bed','.gz','.bai','.md')
FRAG = re.compile(chr(34) + "([^" + chr(34) + chr(10) + "]*)" + chr(34) + "|" + chr(39) + "([^" + chr(39) + chr(10) + "]*)" + chr(39))
hits = []
for r, ds, fs in os.walk('.'):
    ds[:] = [d for d in ds if d not in ('.git', '__pycache__')]
    for f in sorted(fs):
        if not f.endswith(('.py', '.sh', '.conf', '.example')): continue
        p = os.path.join(r, f).replace(os.sep, '/')
        for i, L in enumerate(io.open(p, encoding='utf-8', errors='replace'), 1):
            if L.lstrip().startswith('#'): continue
            for m in TOK.finditer(L.split('#')[0]):
                t = m.group(0)
                if not KO.search(t): continue
                # 경로: 구분자가 둘 이상이거나 아는 확장자로 끝난다.
                #   「요청 3/보유 5」 같은 말은 구분자가 하나뿐이라 걸리지 않는다.
                if t.count('/') >= 2 or t.endswith(EXT):
                    hits.append('%s:%d  %s' % (p, i, t[:64]))
        # 2026-09-29: 위는 「한 토큰」 만 본다. 경로를 이어 붙이면 한글이 조각으로
        #   숨는다. `+ '_cellline'` 자리에 `'_세포주'` 가 있던 것을 못 잡았다.
        #   따옴표 안의 「공백 없는」 한글 조각 중 앞이 _ . / 인 것을 경로 조각으로 본다.
        #   (「암」·「정상」 같은 분류 라벨은 앞에 그런 글자가 없어 걸리지 않는다.)
        for i, L in enumerate(io.open(p, encoding='utf-8', errors='replace'), 1):
            if L.lstrip().startswith('#'): continue
            for m in FRAG.finditer(L.split('#')[0]):
                t = m.group(1) or m.group(2) or ''
                if not KO.search(t) or ' ' in t or chr(9) in t: continue
                if t[:1] in ('_', '.', '/') or t.endswith('/'):
                    hits.append('%s:%d  조각 %r' % (p, i, t[:48]))
print('KOPATH %d' % len(hits))
for h in hits: print('   ', h)
HD9
)
  _n=$(echo "$_out" | awk '/^KOPATH/{print $2}')
  if [ "${_n:-1}" -eq 0 ]; then ok "산출물 경로에 한글 0건"
  else bad "경로에 한글 ${_n}건" "$(echo "$_out" | sed 1d)"; fi
fi
echo
echo "════ 10. 단계 이름 목록 셋이 일치하나 ════"
#   검증은 학습 스크립트를 「복사해」 돌린다. 그 이름 짝(val03_1 <-> 01_...py)이
#   세 곳에 적혀 있다. config.conf 의 STEP_MAP · make_val_config.py 의 기본값 ·
#   run_one_sample.sh 의 기본값. 셋이 어긋나면 그 피처만 빠진 채 완주하고
#   채점이 결측을 학습 평균으로 메워 AUC 0.500 이 난다. 에러는 안 난다.
if [ -n "$_py" ]; then
  _out=$("$_py" - <<'HD10'
import io, re
NL2 = chr(10)
_PAIR = chr(39) + "([^" + chr(39) + "]+[.]py)" + chr(39) + "[ ]*:[ ]*" + chr(39) + "([^" + chr(39) + "]+)[.]py" + chr(39)
def parse_pairs(t):
    d = {}
    for chunk in re.split('[;' + chr(10) + ']', t):
        chunk = chunk.strip().strip('"').strip()
        if not chunk or ':' not in chunk: continue
        k, v = chunk.split(':', 1)
        k, v = k.strip(), v.strip().strip('"')
        if k and v.endswith('.py'): d[k] = v
    return d
conf = io.open('0_setup/config.conf', encoding='utf-8').read()
m = re.search('STEP_MAP=(.*?)' + chr(10) + '[A-Za-z_#]', conf + chr(10) + '#', re.S)
A = parse_pairs(m.group(1)) if m else {}
sh = io.open('3_validate/steps/run_one_sample.sh', encoding='utf-8').read()
m = re.search(r'_pairs="\$\{STEP_MAP:-(.*?)\}"', sh, re.S)
B = parse_pairs(m.group(1)) if m else {}
py = io.open('3_validate/steps/make_val_config.py', encoding='utf-8').read()
m = re.search('_STEPS_DEFAULT = ' + chr(123) + '(.*?)' + chr(125), py, re.S)
C = {}
if m:
    # 주석 줄을 걷어낸다 — 꺼 둔 항목(# 'step04_1b_readfeat.py': ...)을 읽으면
    #   있지도 않은 불일치를 만든다.
    _body = NL2.join(x for x in m.group(1).split(NL2) if not x.strip().startswith('#'))
    for k, v in re.findall(_PAIR, _body):
        C[v] = k                      # {val03_1: 01_....py} 로 방향을 맞춘다
bad = []
for name, d in (('config.conf STEP_MAP', A),
                ('run_one_sample.sh 기본값', B),
                ('make_val_config.py 기본값', C)):
    if not d: bad.append('%s 를 못 읽었다 (형식이 바뀌었나)' % name)
if not bad:
    keys = set(A) | set(B) | set(C)
    for k in sorted(keys):
        vals = {A.get(k), B.get(k), C.get(k)}
        if len(vals) > 1:
            bad.append('%s : conf=%s · sh=%s · py=%s' % (k, A.get(k), B.get(k), C.get(k)))
print('STEPBAD %d (항목 %d)' % (len(bad), len(A)))
for x in bad: print('   ', x)
HD10
)
  _n=$(echo "$_out" | awk '/^STEPBAD/{print $2}')
  _c=$(echo "$_out" | awk '/^STEPBAD/{print $3 $4}')
  if [ "${_n:-1}" -eq 0 ]; then ok "단계 이름 목록 셋 일치 $_c"
  else bad "단계 이름 목록이 어긋난다 ${_n}건" "$(echo "$_out" | sed 1d)"; fi
fi
echo
echo "════ 11. 저장소에 잡동사니가 없나 ════"
#   시험하다 남긴 임시 폴더·백업·편집기 파일이 그대로 올라가는 일이 있었다.
#   2026-09-29: 관문을 태우려고 만든 _t/ 가 트리에 남아 README 의 파일 수와 어긋났다.
_junk=$(find . -not -path './.git/*' -not -path '*__pycache__*' \
  \( -name '_*' -o -name '*.bak' -o -name '*.orig' -o -name '*~' -o -name '*.tmp' \
     -o -name '*.prev' -o -name '.DS_Store' -o -name '*.swp' \) \
  ! -name '_*.md' 2>/dev/null | sort)
if [ -z "$_junk" ]; then ok "임시·백업 파일 0개"
else bad "치워야 할 것 $(printf '%s\n' "$_junk" | wc -l)개" "$(printf '%s\n' "$_junk" | sed 's/^/        /')"; fi
echo
echo "════ 12. 경로가 맨 앞 「경로 블록」 에 모여 있나 ════"
#   각 단계 스크립트는 자기가 읽고 쓰는 자리를 파일 「맨 앞」 한 묶음에 적는다.
#   그래야 그 파일만 열면 무엇을 건드리는지 보이고, 배치가 다른 사람은 거기만 고친다.
#   블록 「밖」 에 자료 경로가 흩어져 있으면 그것을 못 찾는다. 그 상태를 잡는다.
if [ -n "$_py" ]; then
  _out=$("$_py" - <<'HD12'
import io, os
MARK_A, MARK_B = "이 단계가 쓰는 경로", chr(9562)
BAD = ("normal_cfDNA_public", "GBM_cell-line", "Bisulfite_Genome",
       "_input/normal_pub")
SKIP = ("config.conf", "RENAMES.txt", "README.md", "IO.md", "ENVIRONMENT.md",
        "selfcheck.sh", "meth_config.py")
hits, marked = [], 0
for r, ds, fs in os.walk("."):
    ds[:] = [d for d in ds if d not in (".git", "__pycache__")]
    for f in sorted(fs):
        if f in SKIP or not f.endswith((".py", ".sh")): continue
        p = os.path.join(r, f).replace(os.sep, "/")
        L = io.open(p, encoding="utf-8", errors="replace").read().split(chr(10))
        a = next((i for i, x in enumerate(L) if MARK_A in x), None)
        b = next((i for i, x in enumerate(L) if MARK_B in x), None) if a is not None else None
        if a is not None: marked += 1
        for i, x in enumerate(L):
            if x.lstrip().startswith("#"): continue
            body = x.split("#")[0]
            if a is not None and b is not None and a <= i <= b: continue   # 블록 「안」 은 제자리
            for t in BAD:
                if t in body:
                    hits.append("%s:%d  %s" % (p, i + 1, body.strip()[:58]))
print("OUTSIDE %d (블록 표시된 파일 %d개)" % (len(hits), marked))
for h in hits: print("   ", h)
HD12
)
  _n=$(echo "$_out" | awk '/^OUTSIDE/{print $2}')
  if [ "${_n:-1}" -eq 0 ]; then ok "자료 경로 전부 맨 앞 블록 안에 $(echo "$_out" | sed -n '1s/.*(\(.*\))/()/p')"
  else bad "블록 「밖」 에 흩어진 경로 ${_n}곳" "$(echo "$_out" | sed 1d)"; fi
fi
echo
echo "════════════════════════════════════════"
printf '  PASS %d · FAIL %d\n' "$PASS" "$FAIL"
if [ "$FAIL" -eq 0 ]; then
  echo "  구조가 성립합니다. 「결과가 맞다」 는 README 의 Reproducibility 표를 보십시오."
  exit 0
else
  echo "  위 FAIL 을 보십시오."
  exit 1
fi
