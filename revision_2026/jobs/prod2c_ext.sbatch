#!/bin/bash
# Extension segment s3 of a finished two-stage production window (SESSION_NOTES sec 47, 2026-11-06).
# Continues from the backed-up wNN_s2.rst7 (~/prod2c/restart_backup/<sys>/<leg>/, md5-checked) with the SAME template,
# lambda schedule and mask as prod2c.sbatch (irest=1, ntx=5 -> velocities kept). NS3 steps (default 2.5 ns).
# V windows: every FROFF-th frame of s3 at all 29 V lambdas of lambdas.csv -> wNN_sp3_energies.txt (exact MBAR input).
# Output: ~/prod2c/ext3/<sys>/<leg>/runs/wNN_s3.out.gz, wNN_s3.check, [wNN_sp3_energies.txt], wNN_s3.done
# Env: SYS LEG WIN [NS3 NTWX FROFF TAG]
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=2
#SBATCH --time=12:00:00
#SBATCH --partition=u22
#SBATCH --nodes=1
set -uo pipefail
N=$(hostname -s | tr -dc '0-9'); [ "$((10#${N:-0}))" -gt 44 ] || { echo "refusing to run on $(hostname -s): node must be > gnode044"; exit 1; }
module load u22/amber/24
command -v pmemd.cuda >/dev/null && command -v cpptraj >/dev/null || { echo "pmemd.cuda/cpptraj not found"; exit 1; }
NS3=${NS3:-2500000}; NTWX=${NTWX:-20000}; FROFF=${FROFF:-4}; TAG=${TAG:-ext3}
B=$HOME/prod2c; W=$(printf "w%02d" $WIN); O=$B/$TAG/$SYS/$LEG/runs; S=/scratch/$USER/prod2c/$TAG/$SYS/$LEG/$W
mkdir -p $O $S
[ -f $O/${W}_s3.done ] && { echo "$SYS $LEG $W s3 already done"; exit 0; }
case $LEG in *bhd2) MASK=':392-406';; *bhd3) MASK=':475-485';; *lesion) MASK=':602';; *partner) MASK=':584-585';; *) echo "unknown leg $LEG"; exit 1;; esac
case $LEG in pb_*) TOP=system;; docked_*) TOP=system_docked_cleaned;; esac
if [ $SYS = 64pp ]; then R=$HOME/ti_rep1; else R=$HOME/ti_cpd; fi
PRM=$R/topologies/$TOP.prmtop
declare -A MD5=( [64pp/system]=3c576a5d [64pp/system_docked_cleaned]=d961c524 [cpd/system]=ab6fbba2 [cpd/system_docked_cleaned]=23dfa291 )
[ "$(md5sum $PRM | cut -c1-8)" = "${MD5[$SYS/$TOP]}" ] || { echo "PRMTOP md5 mismatch for $SYS/$TOP ($PRM)"; exit 1; }
RST=$B/restart_backup/$SYS/$LEG/${W}_s2.rst7
[ -s $RST ] && [ "$(md5sum < $RST | cut -c1-32)" = "$(cut -c1-32 $RST.md5)" ] || { echo "restart missing or md5 mismatch: $RST"; exit 1; }
read ST LAM <<< "$(awk -F, -v w=$WIN 'NR>1 && $1==w {print $2, $3}' $B/lambdas.csv)"; [ -n "${LAM:-}" ] || { echo "no lambda for window $WIN"; exit 1; }
cd $S; cp -f $B/lambda_$ST.sch lambda.sch
sed -e "s/@LAM@/$LAM/" -e "s/@MASK@/$MASK/g" -e "s/@NSTLIM@/$NS3/" -e "s/@NTWX@/$NTWX/" $B/template_$ST.in > ${W}_s3.in
grep -q "clambda=$LAM," ${W}_s3.in && grep -q "timask1='$MASK'" ${W}_s3.in || { echo "input substitution failed"; exit 1; }
echo "$(hostname -s):$S" > $O/${W}_s3.where; echo "$SYS $LEG $W s3 stage $ST lambda $LAM mask $MASK from $RST on $(hostname -s)"
pmemd.cuda -O -i ${W}_s3.in -p $PRM -c $RST -o ${W}_s3.out -r ${W}_s3.rst7 -x ${W}_s3.nc -inf ${W}_s3.mdinfo > ${W}_s3.stdout 2>&1
gzip -c ${W}_s3.out > $O/${W}_s3.out.gz
grep -q "Total wall time" ${W}_s3.out || { echo "$SYS $LEG $W s3 FAILED: $(grep -h -m1 -E 'Periodic box|ERROR|Error' ${W}_s3.out ${W}_s3.stdout)"; tail -3 ${W}_s3.stdout; exit 1; }
cp -f ${W}_s3.rst7 $O/   # 4 MB, allows a further extension
# same sanity numbers as the s0 check, over the whole s3 segment (record; continuity check on the analysis side)
awk '/A V E R A G E S/{f=1} f{next} /TEMP\(K\)/ && !/Softcore/{for(i=1;i<NF;i++) if($i=="TEMP(K)"){t+=$(i+2);nt++}} /EELEC/{for(i=1;i<NF;i++) if($i=="EELEC"){e+=$(i+2);ne++}}
     /Density/{d=$3; if(d<dmin||dmin=="")dmin=d; if(d>dmax)dmax=d}
     END{T=t/nt; E=e/ne; ok=(T>295 && T<305 && E<-400000 && dmin>1.03 && dmax<1.09); printf "s3 check: <T> %.1f <EELEC> %.0f density %.4f-%.4f -> %s\n",T,E,dmin,dmax,(ok?"OK":"FAIL"); exit !ok}' ${W}_s3.out \
  | tee $O/${W}_s3.check || { echo "$SYS $LEG $W: s3 sanity check FAILED"; exit 1; }
if [ $ST = V ]; then
  mkdir -p sp; cd sp; cp -f $B/lambda_V.sch lambda.sch
  printf "parm $PRM\ntrajin ../${W}_s3.nc 1 last $FROFF\ntrajout fr.rst7 restart multi\nrun\n" > ext.in
  TARGETS=$(awk -F, 'NR>1 && $2=="V" {print $3}' $B/lambdas.csv)
  cpptraj -i ext.in > ext.log 2>&1 || { echo "cpptraj failed"; tail -5 ext.log; exit 1; }
  FR=$(ls fr.rst7.* 2>/dev/null | sort -t. -k3 -n); [ -n "$FR" ] || { echo "no frames extracted"; exit 1; }
  OUTF=$O/${W}_sp3_energies.txt; echo "# $SYS $LEG $W stage V s3 sampled_lambda $LAM ; columns: frame target_lambda E_total_kcal DV/DL" > $OUTF.tmp
  for f in $FR; do for L in $TARGETS; do
    sed -e "s/@LAM@/$L/" -e "s/@MASK@/$MASK/g" $B/sp_V.in > sp.in
    pmemd.cuda -O -i sp.in -p $PRM -c $f -o sp.out > sp.stdout 2>&1 || { echo "single point failed ($f, $L)"; tail -3 sp.stdout; exit 1; }
    echo "${f##*.} $L $(awk -f $B/sp_energy.awk sp.out)" >> $OUTF.tmp
  done; done
  grep -q PARSE_ERROR $OUTF.tmp && { echo "single point parse error"; exit 1; }
  mv $OUTF.tmp $OUTF; echo "single points: $(grep -vc '^#' $OUTF) energies"
fi
touch $O/${W}_s3.done; echo "$SYS $LEG $W s3 done ($(hostname -s))"
