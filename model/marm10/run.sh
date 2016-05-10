#!/bin/bash
# vim: fdm=marker fdl=0
#
# usage example: ./run.sh fm 8
#

# it is critical important to set this environment variable
# in order to keep everything in the corrent folder
export DATAPATH=`pwd`/

allowed_tasks=( \
  "fm" "essfwi" "enfwi" "noise" \
  "fm-sw" "essfwi-sw" "enfwi-sw" "noise-sw" \
  "fm-swintel" "essfwi-swintel" "enfwi-swintel" "noise-swintel" \
  )

# check the parameters#{{{
task=$1
task_is_valid=0
for t in ${allowed_tasks[@]}; do
  if [[ $task == $t ]]; then task_is_valid=1; fi
done
if [[ $task_is_valid -eq 0 ]]; then
  echo "usage: $0 task [nthreads]"
  echo "allowed_tasks: ${allowed_tasks[@]}"
  echo "nthreads only used in GNU or Intel compilers"
  exit
fi
#}}}
# set omp threads#{{{
if [[ -z $2 ]]; then  nthread=8; else  nthread=$2; fi
export OMP_NUM_THREADS=$nthread
#}}}

scons task=$1
