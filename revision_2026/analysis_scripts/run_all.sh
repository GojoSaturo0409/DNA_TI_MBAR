#!/bin/bash
PY=<ENV_DIR>/rocklinc/bin/python; D=<STATE_DIR>/rocklin
top(){ case $1_$2 in 64pp_pb) echo <PROJECT_ROOT>/TI/system.prmtop;; 64pp_docked) echo <PROJECT_ROOT>/TI/rahul_ti_unbiased/docked/system_docked_cleaned.prmtop;; cpd_pb) echo <PROJECT_ROOT>/TI/cpd_unbiased/topologies/cpd_pb.prmtop;; cpd_docked) echo <PROJECT_ROOT>/TI/cpd_unbiased/topologies/cpd_docked.prmtop;; esac; }
sel(){ case $1 in *bhd2) echo "resid 392:406";; *bhd3) echo "resid 475:485";; *lesion) echo "resid 602";; *partner) echo "resid 584:585";; esac; }
run(){ s=$1; l=$2; st=${l%%_*}; o=$D/all/${s}_${l}; mkdir -p $o; $PY $D/run_rocklin.py $(top $s $st) $D/seeds/${s}_${l}/production.rst7 "$(sel $l)" $o/apbs > $o/rocklin.log 2>&1; echo "$s $l rc $?"; }
LIST="64pp:pb_bhd2 64pp:docked_bhd2 64pp:pb_bhd3 64pp:docked_bhd3 64pp:pb_partner 64pp:docked_partner cpd:pb_lesion cpd:docked_lesion cpd:pb_bhd2 cpd:docked_bhd2 cpd:pb_bhd3 cpd:docked_bhd3 cpd:pb_partner cpd:docked_partner"
set -- $LIST
while [ $# -gt 0 ]; do run ${1%%:*} ${1##*:} & [ $# -gt 1 ] && run ${2%%:*} ${2##*:} & wait; shift; [ $# -gt 0 ] && shift; done
