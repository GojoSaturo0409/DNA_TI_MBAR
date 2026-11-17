#!/bin/bash
# Two-stage TI production window (validated protocol, SESSION_NOTES sec 35-38).
#  windows 1-25  = decharge stage E (lambda 0->0.5, lambda_E.sch, scalpha=scbeta=0; exactly linear -> exact MBAR from DV/DL)
#  windows 26-50 = vdW stage V (separate lambda' 0->1, lambda_V.sch, scalpha=0.2)
# Env: SYS (64pp|cpd) LEG (pb_bhd2 ... docked_partner) WIN (1..50).
# Optional (tests only): NS0 NS1 NS2 (steps), NTWX, FR1 (first frame of s1 used for single points), FROFF (frame stride),
#   FE1/FEOFF (frames of s2 for the E-stage linearity check), TAG (separate output tree for pilots/tests).
# Segments: s0 = 100 ps + sanity check; s1 = 2.4 ns; s2 = 2.5 ns. Then single points (1-step minimisation, exact energies):
#  V windows: every FROFF-th frame from FR1 of s1 and all of s2, at all 25 V lambdas -> sp_energies.txt (exact MBAR input)
#  E windows: 3 frames at lambda-0.02, lambda, lambda+0.02 -> linearity proof for this leg
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=2
#SBATCH --time=24:00:00
#SBATCH --partition=u22
#SBATCH --nodes=1
set -uo pipefail
N=$(hostname -s | tr -dc '0-9'); [ "$((10#${N:-0}))" -gt 44 ] || { echo "refusing to run on $(hostname -s): node must be > gnode044"; exit 1; }
module load u22/amber/24
command -v pmemd.cuda >/dev/null && command -v cpptraj >/dev/null || { echo "pmemd.cuda/cpptraj not found"; exit 1; }
NS0=${NS0:-100000}; NS1=${NS1:-2400000}; NS2=${NS2:-2500000}; NTWX=${NTWX:-20000}; FR1=${FR1:-46}; FROFF=${FROFF:-4}; FE1=${FE1:-40}; FEOFF=${FEOFF:-40}; TAG=${TAG:-}
B=$HOME/prod2c; W=$(printf "w%02d" $WIN); O=$B/${TAG:+$TAG/}$SYS/$LEG/runs; S=/scratch/$USER/prod2c/${TAG:+$TAG/}$SYS/$LEG/$W
mkdir -p $O $S
[ -f $O/$W.done ] && { echo "$SYS $LEG $W already done"; exit 0; }
case $LEG in *bhd2) MASK=':392-406';; *bhd3) MASK=':475-485';; *lesion) MASK=':602';; *partner) MASK=':584-585';; *) echo "unknown leg $LEG"; exit 1;; esac
case $LEG in pb_*) TOP=system;; docked_*) TOP=system_docked_cleaned;; esac
if [ $SYS = 64pp ]; then R=$HOME/ti_rep1; else R=$HOME/ti_cpd; fi
SEED=$R/$LEG/production.rst7; PRM=$R/topologies/$TOP.prmtop
declare -A MD5=( [64pp/pb_bhd2]=551c6e9a [64pp/pb_bhd3]=10e0f1a3 [64pp/pb_lesion]=f9677a4c [64pp/pb_partner]=972ee576
 [64pp/docked_bhd2]=b73488dd [64pp/docked_bhd3]=e7aec56b [64pp/docked_lesion]=ca1eea56 [64pp/docked_partner]=dbe48acf
 [cpd/pb_bhd2]=c8ed5562 [cpd/pb_bhd3]=fbf0d14f [cpd/pb_lesion]=2715a522 [cpd/pb_partner]=84ec7375
 [cpd/docked_bhd2]=7f7fa50c [cpd/docked_bhd3]=44e6938e [cpd/docked_lesion]=03c704f0 [cpd/docked_partner]=0e8e7a88
 [64pp/system]=3c576a5d [64pp/system_docked_cleaned]=d961c524 [cpd/system]=ab6fbba2 [cpd/system_docked_cleaned]=23dfa291 )
[ "$(md5sum $SEED | cut -c1-8)" = "${MD5[$SYS/$LEG]}" ] || { echo "SEED md5 mismatch for $SYS/$LEG ($SEED)"; exit 1; }
[ "$(md5sum $PRM | cut -c1-8)" = "${MD5[$SYS/$TOP]}" ] || { echo "PRMTOP md5 mismatch for $SYS/$TOP ($PRM)"; exit 1; }
read ST LAM <<< "$(awk -F, -v w=$WIN 'NR>1 && $1==w {print $2, $3}' $B/lambdas.csv)"; [ -n "${LAM:-}" ] || { echo "no lambda for window $WIN"; exit 1; }
cd $S; cp -f $B/lambda_$ST.sch lambda.sch   # pmemd reads ./lambda.sch
mkin(){ sed -e "s/@LAM@/$LAM/" -e "s/@MASK@/$MASK/g" -e "s/@NSTLIM@/$1/" -e "s/@NTWX@/$NTWX/" $B/template_$ST.in > $2
  grep -q "clambda=$LAM," $2 && grep -q "timask1='$MASK'" $2 || { echo "input substitution failed ($2)"; exit 1; }; }
echo "$(hostname -s):$S" > $O/$W.where; echo "$SYS $LEG $W stage $ST lambda $LAM mask $MASK on $(hostname -s)"
run(){ # seg nstlim start_coords
  mkin $2 ${W}_$1.in
  pmemd.cuda -O -i ${W}_$1.in -p $PRM -c $3 -o ${W}_$1.out -r ${W}_$1.rst7 -x ${W}_$1.nc -inf ${W}_$1.mdinfo > ${W}_$1.stdout 2>&1
  gzip -c ${W}_$1.out > $O/${W}_$1.out.gz
  grep -q "Total wall time" ${W}_$1.out || { echo "$SYS $LEG $W $1 FAILED: $(grep -h -m1 -E 'Periodic box|ERROR|Error' ${W}_$1.out ${W}_$1.stdout)"; tail -3 ${W}_$1.stdout; exit 1; }
}
run s0 $NS0 $SEED
awk '/A V E R A G E S/{f=1} f{next} /TEMP\(K\)/ && !/Softcore/{for(i=1;i<NF;i++) if($i=="TEMP(K)"){t+=$(i+2);nt++}} /EELEC/{for(i=1;i<NF;i++) if($i=="EELEC"){e+=$(i+2);ne++}}
     /Density/{d=$3; if(d<dmin||dmin=="")dmin=d; if(d>dmax)dmax=d}
     END{T=t/nt; E=e/ne; ok=(T>295 && T<305 && E<-400000 && dmin>1.03 && dmax<1.09); printf "s0 check: <T> %.1f <EELEC> %.0f density %.4f-%.4f -> %s\n",T,E,dmin,dmax,(ok?"OK":"FAIL"); exit !ok}' ${W}_s0.out \
  | tee $O/${W}_s0.check || { echo "$SYS $LEG $W: s0 sanity check FAILED, stopping"; exit 1; }
run s1 $NS1 ${W}_s0.rst7
run s2 $NS2 ${W}_s1.rst7
# ---- single points (exact energies) ----
mkdir -p sp; cd sp; cp -f $B/lambda_$ST.sch lambda.sch
if [ $ST = V ]; then
  printf "parm $PRM\ntrajin ../${W}_s1.nc $FR1 last $FROFF\ntrajin ../${W}_s2.nc 1 last $FROFF\ntrajout fr.rst7 restart multi\nrun\n" > ext.in
  TARGETS=$(awk -F, 'NR>1 && $2=="V" {print $3}' $B/lambdas.csv)
else
  printf "parm $PRM\ntrajin ../${W}_s2.nc $FE1 last $FEOFF\ntrajout fr.rst7 restart multi\nrun\n" > ext.in
  # linearity check points stay inside the decharge stage [0,0.5]: delta = min(0.02, lam/2, (0.5-lam)/2)
  TARGETS=$(awk -v l=$LAM 'BEGIN{d=0.02; if(l/2<d)d=l/2; if((0.5-l)/2<d)d=(0.5-l)/2; printf "%.6f %.6f %.6f", l-d, l, l+d}')
fi
cpptraj -i ext.in > ext.log 2>&1 || { echo "cpptraj failed"; tail -5 ext.log; exit 1; }
FR=$(ls fr.rst7.* 2>/dev/null | sort -t. -k3 -n); [ -n "$FR" ] || { echo "no frames extracted"; exit 1; }
OUTF=$O/${W}_sp_energies.txt; echo "# $SYS $LEG $W stage $ST sampled_lambda $LAM ; columns: frame target_lambda E_total_kcal DV/DL" > $OUTF.tmp
for f in $FR; do for L in $TARGETS; do
  sed -e "s/@LAM@/$L/" -e "s/@MASK@/$MASK/g" $B/sp_$ST.in > sp.in
  pmemd.cuda -O -i sp.in -p $PRM -c $f -o sp.out > sp.stdout 2>&1 || { echo "single point failed ($f, $L)"; tail -3 sp.stdout; exit 1; }
  echo "${f##*.} $L $(awk -f $B/sp_energy.awk sp.out)" >> $OUTF.tmp
done; done
grep -q PARSE_ERROR $OUTF.tmp && { echo "single point parse error"; exit 1; }
mv $OUTF.tmp $OUTF; echo "single points: $(grep -vc '^#' $OUTF) energies ($(echo $FR | wc -w) frames x $(echo $TARGETS | wc -w) lambdas)"
touch $O/$W.done; echo "$SYS $LEG $W done ($(hostname -s))"
