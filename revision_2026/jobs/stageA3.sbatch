#!/bin/bash
# Stage A v3 - minimization/equilibration order requested by the user, with the paper's (Raghunathan & Krishnan, JCIM 2025) minimization settings:
#   P1  restrained minimization, 25,000 cycles, heavy atoms of RAD4+DNA restrained (k = 100 kcal/mol/A^2)   [paper]
#   P2  unrestrained minimization, 25,000 cycles                                                           [paper]
#   T1  unrestrained minimization at lambda = 0.5 (TI soft-core on; steepest descent as soft-core requires), 25,000 cycles
#   H   heating 10 -> 300 K, NVT, 200 ps, lambda = 0.5
#   V   NVT equilibration, 300 K, 500 ps, lambda = 0.5
#   N   NPT equilibration, 300 K, 1 bar, 2.4 ns, lambda = 0.5  (Langevin + Berendsen, PME)  -> production.rst7
#   G   gate: TI windows 1 and 12, 5,000 steps each from production.rst7 (must stay at ~300 K)
# All MD at lambda = 0.5 uses soft-core, which AMBER only allows with ntc=1/ntf=1 -> dt = 1 fs.
# Durations of H/V/N and dt are NOT given in the paper's main text (SI Table S1) - values here are assumptions to verify.
# Needs env: LEG MASK STATE (pb|docked). Each stage is sanity-checked (temperature/density/energy), not only the pmemd exit code.
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=2
#SBATCH --time=12:00:00
#SBATCH --partition=u22
#SBATCH --nodes=1
set -euo pipefail
N=$(hostname -s | tr -dc '0-9'); [ "$((10#${N:-0}))" -gt 44 ] || { echo "refusing to run on $(hostname -s): node must be > gnode044"; exit 1; }
module load u22/amber/24
command -v pmemd.cuda >/dev/null || { echo "pmemd.cuda not found after module load"; exit 1; }
R=$HOME/ti_rep1; L=$R/$LEG; S=/scratch/$USER/ti_rep1/$LEG; IN=$R/amber_inputs
mkdir -p "$L" "$S"; cd "$S"; hostname > "$L/stageA_node.txt"
if [ "$STATE" = pb ]; then TOP=system; else TOP=system_docked_cleaned; fi
PRM=$R/topologies/$TOP.prmtop; CRD=$R/topologies/$TOP.inpcrd
RMASK=':1-620&!@H='          # RAD4 (565 res) + DNA (55 res, incl. lesion T64); ions (621-643) and waters are free
# Heavy files (trajectories, per-segment outputs/restarts) stay in this node's scratch ($S); home gets only small stage outputs,
# a gzipped npt.out, and a note of the node. production.rst7 (needed by windows on other nodes) is copied separately below.
trap 'mkdir -p "$L/prep"; echo "$(hostname -s):$S" > "$L/prep/where.txt"; for f in "$S"/p1.out "$S"/p2.out "$S"/t1.out "$S"/h.out "$S"/v.out "$S"/gate_w*.out "$S"/p1.in "$S"/p2.in "$S"/t1.in "$S"/h.in "$S"/v.in "$S"/npt.in; do [ -f "$f" ] && [ "$(stat -c %s "$f")" -lt 5000000 ] && cp -f "$f" "$L/prep"/; done; [ -f "$S/npt.out" ] && gzip -c "$S/npt.out" > "$L/prep/npt.out.gz"; true' EXIT
source "$R/heat.inc"; source "$R/chk.inc"
fail() { echo "STAGE A FAILED at $1"; exit 1; }
T0=$SECONDS; stamp() { echo "STAGE $1 done at $((SECONDS-T0)) s"; }

# TI block shared by all lambda=0.5 stages
TI="  icfe=1, clambda=0.5, ifsc=1,
  timask1='$MASK', scmask1='$MASK', timask2='', scmask2='',
  gti_lam_sch=1, gti_ele_sc=1, gti_vdw_sc=1, gti_cut_sc=2, gti_add_sc=1,"

# ---------- P1: restrained minimization ----------
cat > p1.in <<EOF
Restrained minimization (heavy atoms RAD4+DNA, k=100)
&cntrl
  imin=1, ntmin=1, maxcyc=25000, ncyc=10000,
  ntb=1, cut=10.0,
  ntr=1, restraint_wt=100.0, restraintmask='$RMASK',
  ntpr=1000, ntwr=5000,
/
EOF
pmemd.cuda -O -i p1.in -p "$PRM" -c "$CRD" -ref "$CRD" -o p1.out -r p1.rst7 -inf p1.mdinfo || fail P1
chk_min p1.out || fail P1-check; stamp P1

# ---------- P2: unrestrained minimization ----------
cat > p2.in <<EOF
Unrestrained minimization
&cntrl
  imin=1, ntmin=1, maxcyc=25000, ncyc=10000,
  ntb=1, cut=10.0, ntr=0,
  ntpr=1000, ntwr=5000,
/
EOF
pmemd.cuda -O -i p2.in -p "$PRM" -c p1.rst7 -o p2.out -r p2.rst7 -inf p2.mdinfo || fail P2
chk_min p2.out || fail P2-check; stamp P2

# ---------- T1: unrestrained minimization at lambda = 0.5 ----------
cat > t1.in <<EOF
Unrestrained minimization at lambda=0.5 (soft-core: steepest descent only)
&cntrl
  imin=1, ntmin=2, maxcyc=25000,
  ntb=1, cut=10.0, ntr=0,
  ntpr=1000, ntwr=5000,
$TI
/
EOF
pmemd.cuda -O -i t1.in -p "$PRM" -c p2.rst7 -o t1.out -r t1.rst7 -inf t1.mdinfo || fail T1
chk_min t1.out || fail T1-check; stamp T1

# ---------- H: heating 10 -> 300 K, NVT, 200 ps, lambda = 0.5 ----------
cat > h.in <<EOF
Heating 10 -> 300 K, NVT, lambda=0.5
&cntrl
  imin=0, ntx=1, irest=0, nstlim=200000, dt=0.001,
  ntb=1, ntt=3, gamma_ln=2.0, tempi=10.0, temp0=300.0, ig=-1,
  ntc=1, ntf=1, cut=10.0,
  ntpr=1000, ntwx=10000, ntwr=5000,
  nmropt=1,
$TI
/
&wt type='TEMP0', istep1=0, istep2=200000, value1=10.0, value2=300.0, /
&wt type='END', /
EOF
pmemd.cuda -O -i h.in -p "$PRM" -c t1.rst7 -o h.out -r h.rst7 -x h.nc -inf h.mdinfo || fail H
chk_md h.out 285 310 || fail H-check; stamp H

# ---------- V: NVT equilibration, 300 K, 500 ps, lambda = 0.5 ----------
cat > v.in <<EOF
NVT equilibration 300 K, lambda=0.5
&cntrl
  imin=0, ntx=5, irest=1, nstlim=500000, dt=0.001,
  ntb=1, ntt=3, gamma_ln=2.0, temp0=300.0, ig=-1,
  ntc=1, ntf=1, cut=10.0,
  ntpr=1000, ntwx=10000, ntwr=5000,
$TI
/
EOF
pmemd.cuda -O -i v.in -p "$PRM" -c h.rst7 -o v.out -r v.rst7 -x v.nc -inf v.mdinfo || fail V
chk_md v.out 290 310 || fail V-check; stamp V

# ---------- N: NPT equilibration, 300 K, 1 bar, 2.4 ns, lambda = 0.5 (restart segments: box changes) ----------
sed -e "s/:392-406/$MASK/g" "$IN/post_heat.in" > npt.in
grep -q "ntx=5," npt.in && grep -q "irest=1," npt.in && grep -q "nstlim=2400000," npt.in && grep -q "ntc=1," npt.in && grep -q "ntf=1," npt.in || fail "N-input-check"
NSEG="$(for i in 1 2 3 4 5 6 7 8 9 10; do echo -n '500 '; done)$(for i in 1 2 3 4 5 6 7 8 9 10; do echo -n '2000 '; done)$(for i in 1 2 3 4 5 6 7 8 9; do echo -n '10000 '; done)85000 $(for i in $(seq 1 22); do echo -n '100000 '; done)"
HEAT_TEMPLATE="$S/npt.in" HEAT_NSTLIM=2400000 HEAT_START=v.rst7 heat_segmented npt 0 $NSEG || fail N
chk_md npt_s$(ls npt_s*.out | sed -E 's/npt_s([0-9]+)\.out/\1/' | sort -n | tail -1).out 290 310 0.97 1.12 || fail N-check
cp -f npt.rst7 "$L/production.rst7"; stamp N

# ---------- G: gate - TI windows 1 and 12, 5000 steps each ----------
LAMS=(0.00922 0.04794 0.11505 0.20634 0.31608 0.43738 0.56262 0.68392 0.79366 0.88495 0.95206 0.99078)
for i in 1 12; do
  LAM=${LAMS[$((i-1))]}; W=gate_w$i
  sed -e "s/:392-406/$MASK/g" -e "s/nstlim=10000000,/nstlim=5000,/" -e "s/ntpr=1000,/ntpr=500,/" -e "s/clambda=0.00922,/clambda=$LAM,/" "$IN/ti.in" > $W.in
  grep -q "clambda=$LAM," $W.in || fail "gate lambda substitution w$i"
  pmemd.cuda -O -i $W.in -p "$PRM" -c "$L/production.rst7" -o $W.out -r $W.rst7 -inf $W.mdinfo || fail "gate w$i"
  chk_md $W.out 285 315 || fail "gate-check w$i"
done
echo "stage A done: $L/production.rst7  (total $((SECONDS-T0)) s)"
