#!/bin/bash
# CPD unbiased: minimization + equilibration, same order/settings as stage A v3 but WITHOUT TI (no icfe/soft-core).
#   P1 restrained min (heavy atoms 1-620, k=100) 25,000 cycles; P2 unrestrained min 25,000 cycles; T1 steepest-descent min 25,000 cycles (no TI)
#   H heating 10->300 K NVT 200 ps; V NVT 500 ps; N NPT 2.4 ns   (1 fs, ntc=1/ntf=1 as in stage A)
# Needs env: SYS (pb|docked). Reuses heat.inc/chk.inc from ~/ti_rep1 (heat_segmented handles the GPU box-change halt).
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=2
#SBATCH --time=12:00:00
#SBATCH --partition=u22
#SBATCH --nodes=1
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=rahul.singal@research.iiit.ac.in
set -euo pipefail
N=$(hostname -s | tr -dc '0-9'); [ "$((10#${N:-0}))" -gt 44 ] || { echo "refusing to run on $(hostname -s)"; exit 1; }
module load u22/amber/24
R=$HOME/ti_rep1; C=$HOME/cpd_unbiased; S=/scratch/$USER/cpd_unbiased/$SYS; mkdir -p "$C/$SYS" "$S"; cd "$S"
PRM=$C/topologies/cpd_$SYS.prmtop; CRD=$C/topologies/cpd_$SYS.inpcrd
source "$R/heat.inc"; source "$R/chk.inc"; fail(){ echo "EQUIL FAILED at $1"; exit 1; }
trap 'for f in p1 p2 t1 h v; do [ -f $f.out ] && cp -f $f.out "$C/$SYS/"; done; [ -f npt.out ] && gzip -c npt.out > "$C/$SYS/npt.out.gz"; echo "$(hostname -s):$S" > "$C/$SYS/where.txt"; true' EXIT
RALL=':1-620&!@H='
cat > p1.in <<X
Restrained minimization (heavy atoms RAD4+DNA, k=100)
&cntrl
  imin=1, ntmin=1, maxcyc=25000, ncyc=10000, ntb=1, cut=10.0,
  ntr=1, restraint_wt=100.0, restraintmask='$RALL', ntpr=1000, ntwr=5000,
/
X
pmemd.cuda -O -i p1.in -p "$PRM" -c "$CRD" -ref "$CRD" -o p1.out -r p1.rst7 || fail P1; chk_min p1.out || fail P1-check
cat > p2.in <<X
Unrestrained minimization
&cntrl
  imin=1, ntmin=1, maxcyc=25000, ncyc=10000, ntb=1, cut=10.0, ntr=0, ntpr=1000, ntwr=5000,
/
X
pmemd.cuda -O -i p2.in -p "$PRM" -c p1.rst7 -o p2.out -r p2.rst7 || fail P2; chk_min p2.out || fail P2-check
# T1: steepest-descent minimization (as stage A "T1" but without TI). REQUIRED: P1/P2 end with "REPEATED LINMIN FAILURE" and their restart
# files do not hold the printed minimum (P2 final -4.7e5 but next stage starts at +4e5 / -1.5e5); heating from p2.rst7 blew up (TEMP ~1e6 K).
cat > t1.in <<X
Unrestrained steepest-descent minimization (no TI)
&cntrl
  imin=1, ntmin=2, maxcyc=25000, ntb=1, cut=10.0, ntr=0, ntpr=1000, ntwr=5000,
/
X
pmemd.cuda -O -i t1.in -p "$PRM" -c p2.rst7 -o t1.out -r t1.rst7 || fail T1; chk_min t1.out || fail T1-check
cat > h.in <<X
Heating 10 -> 300 K, NVT
&cntrl
  imin=0, ntx=1, irest=0, nstlim=200000, dt=0.001, ntb=1, ntt=3, gamma_ln=2.0, tempi=10.0, temp0=300.0, ig=-1,
  ntc=1, ntf=1, cut=10.0, ntpr=1000, ntwx=10000, ntwr=5000, nmropt=1,
/
&wt type='TEMP0', istep1=0, istep2=200000, value1=10.0, value2=300.0, /
&wt type='END', /
X
pmemd.cuda -O -i h.in -p "$PRM" -c t1.rst7 -o h.out -r h.rst7 -x h.nc || fail H; chk_md h.out 270 310 || fail H-check
cat > v.in <<X
NVT equilibration 300 K
&cntrl
  imin=0, ntx=5, irest=1, nstlim=500000, dt=0.001, ntb=1, ntt=3, gamma_ln=2.0, temp0=300.0, ig=-1,
  ntc=1, ntf=1, cut=10.0, ntpr=1000, ntwx=10000, ntwr=5000,
/
X
pmemd.cuda -O -i v.in -p "$PRM" -c h.rst7 -o v.out -r v.rst7 -x v.nc || fail V; chk_md v.out 290 310 || fail V-check
# NPT input written explicitly (stage A settings: 1 fs, ntc=1/ntf=1, 2.4 ns, Langevin + Berendsen, no TI block)
cat > npt.in <<X
NPT equilibration 300 K, 1 bar, 2.4 ns
&cntrl
  imin=0,
  ntx=5,
  irest=1,
  nstlim=2400000,
  dt=0.001,
  ntb=2,
  pres0=1.0,
  ntp=1,
  taup=2.0,
  temp0=300.0,
  ntt=3,
  gamma_ln=2.0,
  ig=-1,
  ntc=1,
  ntf=1,
  cut=10.0,
  ntpr=1000,
  ntwx=20000,
  ntwr=5000,
/
X
grep -q "nstlim=2400000," npt.in && grep -q "ntx=5," npt.in && grep -q "irest=1," npt.in || fail "N-input-check"
NSEG="$(for i in $(seq 10); do echo -n '500 '; done)$(for i in $(seq 10); do echo -n '2000 '; done)$(for i in $(seq 9); do echo -n '10000 '; done)85000 $(for i in $(seq 22); do echo -n '100000 '; done)"
HEAT_TEMPLATE="$S/npt.in" HEAT_NSTLIM=2400000 HEAT_START=v.rst7 heat_segmented npt 0 $NSEG || fail N
chk_md npt_s$(ls npt_s*.out | sed -E 's/npt_s([0-9]+)\.out/\1/' | sort -n | tail -1).out 290 310 0.97 1.12 || echo "WARNING: NPT end check failed (non-blocking by request; inspect npt.out.gz before trusting this system)"
cp -f npt.rst7 "$C/$SYS/equil.rst7"; echo "EQUIL DONE: $C/$SYS/equil.rst7"
