#!/bin/bash
# CPD unbiased production: 50 ns, NPT 300 K / 1 bar, 2 fs with SHAKE (ntc=2/ntf=2). Needs env: SYS (pb|docked) REP (1-4).
# Replicas start from the same equilibrated restart; ig=-1 (random Langevin seed) makes them diverge.
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=2
#SBATCH --time=24:00:00
#SBATCH --partition=u22
#SBATCH --nodes=1
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=rahul.singal@research.iiit.ac.in
set -euo pipefail
N=$(hostname -s | tr -dc '0-9'); [ "$((10#${N:-0}))" -gt 44 ] || { echo "refusing to run on $(hostname -s)"; exit 1; }
module load u22/amber/24
source "$HOME/ti_rep1/chk.inc"
C=$HOME/cpd_unbiased; S=/scratch/$USER/cpd_unbiased/$SYS/rep$REP; mkdir -p "$S" "$C/$SYS/rep$REP"; cd "$S"
PRM=$C/topologies/cpd_$SYS.prmtop
cat > prod.in <<X
50 ns unbiased production
&cntrl
  imin=0, ntx=5, irest=1, nstlim=25000000, dt=0.002, ntb=2, ntp=1, pres0=1.0, taup=2.0,
  ntt=3, gamma_ln=2.0, temp0=300.0, ig=-1, ntc=2, ntf=2, cut=10.0,
  ntpr=5000, ntwx=25000, ntwr=500000,
/
X
# restart if a previous run of this replica was interrupted
if [ -f prod.rst7 ]; then echo "continuing from existing prod.rst7 is NOT automatic; rerun from equil.rst7"; fi
pmemd.cuda -O -i prod.in -p "$PRM" -c "$C/$SYS/equil.rst7" -o prod.out -r prod.rst7 -x prod.nc -inf prod.mdinfo
# sanity check (exit code alone misses blown-up runs): mean T 295-305 K, density 0.97-1.12
chk_md prod.out 295 305 0.97 1.12 || { echo "PROD CHECK FAILED $SYS rep$REP"; gzip -c prod.out > "$C/$SYS/rep$REP/prod.out.gz"; exit 1; }
gzip -c prod.out > "$C/$SYS/rep$REP/prod.out.gz"; echo "$(hostname -s):$S" > "$C/$SYS/rep$REP/where.txt"; echo "PROD DONE $SYS rep$REP"
