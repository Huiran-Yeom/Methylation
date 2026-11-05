#!/bin/bash
# 판마다 다른 설정
# profile <이름> 이 정하는 값:
#   DS, O(대기열, 로그 폴더), COHORTS, NV, NL, MEMW, REAL(실판), DEL_ROOTS(중간 산출 지워도 되는 METH_ROOT)
#   OUT(표 폴더), STAGE1(1단계 드라이버 이름, 안내용), EXC_*(미리 정해 둔 예외 칸, 있으면)
profile(){
  local name=$1
  unset EXC_C EXC_V EXC_S EXC_N
  case "$name" in
    r20_gbm|r20_lung)
      DS=${name#r20_}; O=$HOME/tmp/r20/run; NL=2; MEMW=26; STAGE1=r20_run.sh
      if [ "$DS" = gbm ]; then
        env_gbm /ssd_data/Methylation_r20 r20_gbm; COHORTS="patient178_j51"; NV=12
        export NORMAL_SET=/ssd_data/Methylation/sample_data/_input/normal_pub15
        DEL_ROOTS="/ssd_data/Methylation_r20"; OUT=$HOME/표/00_보고서/결과_결함수정2판/GBM
      else
        env_lung /ssd_data/Methylation_r20/lung r20_lung; COHORTS="Patient_Test10 Silico_Test"; NV=8
        DEL_ROOTS="/ssd_data/Methylation_r20/lung"; OUT=$HOME/표/00_보고서/결과_결함수정2판/lung
      fi
      REAL="bl200 jsd200 jsdb200" ;;
    r1_10|r1_5)
      local s=${name#r1_}
      env_luad "$s" || return 1
      DS=luad$s; O=$HOME/tmp/r1/run; COHORTS="Patient_Test10 Silico_Test"; NV=8; NL=1; MEMW=22; STAGE1=r1_run.sh
      REAL=${PANELS//,/ }; DEL_ROOTS="/ssd_data/Methylation_luad_a1 /ssd_data/Methylation_luad_a1u"
      OUT=$HOME/표/00_보고서/결과_다음판_1안/LUAD_${s}런 ;;
    r2_10|r2_5)
      local s=${name#r2_}
      env_luad2 "$s" m || return 1
      DS=luad2_$s; O=$HOME/tmp/r2/run; COHORTS="Patient_Test10"; NV=8; NL=1; MEMW=22; STAGE1=r2_run.sh
      REAL=${PANELS//,/ }; DEL_ROOTS="/ssd_data/Methylation_luad_a2 /ssd_data/Methylation_luad_a2u"
      OUT=$HOME/표/00_보고서/결과_다음판_2안/LUAD_${s}런 ;;
    *) echo "판 이름: r20_gbm r20_lung r1_10 r1_5 r2_10 r2_5"; return 1 ;;
  esac
  case "$name" in r1_10|r2_10)     # 예외 칸: 이 한 칸만 45/48 (3명 풀 부족) 허용
    EXC_C=Patient_Test10; EXC_V=j17bl200_100k_uni; EXC_S=45/48
    EXC_N="G301266_GSM9078998 G301266_GSM9079001 G301266_GSM9079004" ;;
  esac
}

person_re(){   # 4단계 사람 단위 묶기 정규식
  case "$1" in Silico_Test) echo '^(Person[0-9]+)__';; Patient_Test10) echo '^(G301266_GSM[0-9]+)$';; *) echo '';; esac
}
