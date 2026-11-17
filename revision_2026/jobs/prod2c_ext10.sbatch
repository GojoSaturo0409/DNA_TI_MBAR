#!/bin/bash
# Bring one finished 5 ns (or 7.5 ns) two-stage window to exactly 10 ns (SESSION_NOTES sec 50, 2026-11-06; user chose
# "uniform 10 ns everywhere"). Same template / lambda schedule / mask as prod2c.sbatch (irest=1, ntx=5, velocities kept).
#  no s3 yet            -> s3 = 5.0 ns from the md5-checked s2 restart (restart_backup/[bhd2dense/]<sys>/<leg>/wNN_s2.rst7)
#  s3 of 2.5 ns exists  -> s4 = 2.5 ns from ext3/<sys>/<leg>/runs/wNN_s3.rst7
#  s3 of 5.0 ns exists  -> nothing to run (already 10 ns)
# V windows: every FROFF-th frame of the new segment at all 32 unique V lambdas of lambdas_bhd2dense.csv -> wNN_sp3/sp4_energies.txt
# Output: ~/prod2c/ext3/<sys>/<leg>/runs/wNN_<seg>.out.gz, .check, [sp], wNN_<seg>.done, wNN_10ns.done
# Env: SYS LEG WIN [NTWX FROFF]
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=2
#SBATCH --time=12:00:00
#SBATCH --partition=u22
#SBATCH --nodes=1
set -uo pipefail
N=$(hostname -s | tr -dc '0-9'); [ "$((10#${N:-0}))" -gt 44 ] || { echo "refusing to run on $(hostname -s): node must be > gnode044"; exit 1; }
module load u22/amber/24
command -v pmemd.cuda >/dev/null && command -v cpptraj >/dev/null || { echo "pmemd.cuda/cpptraj not found"; exit 1; }
NTWX=${NTWX:-20000}; FROFF=${FROFF:-4}
B=$HOME/prod2c; W=$(printf "w%02d" $WIN); O=$B/ext3/$SYS/$LEG/runs; mkdir -p $O
[ -f $O/${W}_10ns.done ] && { echo "$SYS $LEG $W already 10 ns"; exit 0; }
LAMS=$B/lambdas_bhd2dense.csv
if [ $WIN -ge 55 ]; then RST0=$B/restart_backup/bhd2dense/$SYS/$LEG/${W}_s2.rst7; [ -s $RST0 ] || RST0=$B/restart_backup/dense9/$SYS/$LEG/${W}_s2.rst7; else RST0=$B/restart_backup/$SYS/$LEG/${W}_s2.rst7; fi
# guard: a window that was run at 10 ns directly (s2 = 7.5 ns) needs no extension
for f in $B/$SYS/$LEG/runs/${W}_s2.out.gz $B/dense9/$SYS/$LEG/runs/${W}_s2.out.gz $B/bhd2dense/$SYS/$LEG/runs/${W}_s2.out.gz; do
  [ -s $f ] && [ "$(zcat $f | grep -m1 -o 'nstlim *= *[0-9]*' | tr -dc '0-9')" = 7500000 ] && { touch $O/${W}_10ns.done; echo "$SYS $LEG $W was run at 10 ns directly"; exit 0; }
done
if [ -f $O/${W}_s3.done ]; then
  N3=$(zcat $O/${W}_s3.out.gz | grep -m1 -o 'nstlim *= *[0-9]*' | tr -dc '0-9')
  case $N3 in
    5000000) touch $O/${W}_10ns.done; echo "$SYS $LEG $W: s3 is already 5 ns -> 10 ns"; exit 0;;
    2500000) SEG=s4; NS=2500000; RST=$O/${W}_s3.rst7; [ -s $RST ] || { echo "missing $RST"; exit 1; };;
    *) echo "unexpected s3 nstlim '$N3'"; exit 1;;
  esac
else
  SEG=s3; NS=5000000; RST=$RST0
  [ -s $RST ] && [ "$(md5sum < $RST | cut -c1-32)" = "$(cut -c1-32 $RST.md5)" ] || { echo "restart missing or md5 mismatch: $RST"; exit 1; }
fi
[ -f $O/${W}_$SEG.done ] && { touch $O/${W}_10ns.done; exit 0; }
case $LEG in *bhd2) MASK=':392-406';; *bhd3) MASK=':475-485';; *lesion) MASK=':602';; *partner) MASK=':584-585';; *) echo "unknown leg $LEG"; exit 1;; esac
case $LEG in pb_*) TOP=system;; docked_*) TOP=system_docked_cleaned;; esac
if [ $SYS = 64pp ]; then R=$HOME/ti_rep1; else R=$HOME/ti_cpd; fi
PRM=$R/topologies/$TOP.prmtop
declare -A MD5=( [64pp/system]=3c576a5d [64pp/system_docked_cleaned]=d961c524 [cpd/system]=ab6fbba2 [cpd/system_docked_cleaned]=23dfa291 )
[ "$(md5sum $PRM | cut -c1-8)" = "${MD5[$SYS/$TOP]}" ] || { echo "PRMTOP md5 mismatch for $SYS/$TOP ($PRM)"; exit 1; }
read ST LAM <<< "$(awk -F, -v w=$WIN 'NR>1 && $1==w {print $2, $3}' $LAMS)"; [ -n "${LAM:-}" ] || { echo "no lambda for window $WIN"; exit 1; }
S=/scratch/$USER/prod2c/ext3/$SYS/$LEG/${W}_$SEG; mkdir -p $S; cd $S; cp -f $B/lambda_$ST.sch lambda.sch
sed -e "s/@LAM@/$LAM/" -e "s/@MASK@/$MASK/g" -e "s/@NSTLIM@/$NS/" -e "s/@NTWX@/$NTWX/" $B/template_$ST.in > ${W}_$SEG.in
grep -q "clambda=$LAM," ${W}_$SEG.in && grep -q "timask1='$MASK'" ${W}_$SEG.in && grep -q "nstlim *= *$NS" ${W}_$SEG.in || { echo "input substitution failed"; exit 1; }
echo "$(hostname -s):$S" > $O/${W}_$SEG.where; echo "$SYS $LEG $W $SEG ($NS steps) stage $ST lambda $LAM mask $MASK from $RST on $(hostname -s)"
pmemd.cuda -O -i ${W}_$SEG.in -p $PRM -c $RST -o ${W}_$SEG.out -r ${W}_$SEG.rst7 -x ${W}_$SEG.nc -inf ${W}_$SEG.mdinfo > ${W}_$SEG.stdout 2>&1
gzip -c ${W}_$SEG.out > $O/${W}_$SEG.out.gz
grep -q "Total wall time" ${W}_$SEG.out || { echo "$SYS $LEG $W $SEG FAILED: $(grep -h -m1 -E 'Periodic box|ERROR|Error' ${W}_$SEG.out ${W}_$SEG.stdout)"; tail -3 ${W}_$SEG.stdout; exit 1; }
cp -f ${W}_$SEG.rst7 $O/
awk -v seg=$SEG '/A V E R A G E S/{f=1} f{next} /TEMP\(K\)/ && !/Softcore/{for(i=1;i<NF;i++) if($i=="TEMP(K)"){t+=$(i+2);nt++}} /EELEC/{for(i=1;i<NF;i++) if($i=="EELEC"){e+=$(i+2);ne++}}
     /Density/{d=$3; if(d<dmin||dmin=="")dmin=d; if(d>dmax)dmax=d}
     END{T=t/nt; E=e/ne; ok=(T>295 && T<305 && E<-400000 && dmin>1.03 && dmax<1.09); printf "%s check: <T> %.1f <EELEC> %.0f density %.4f-%.4f -> %s\n",seg,T,E,dmin,dmax,(ok?"OK":"FAIL"); exit !ok}' ${W}_$SEG.out \
  | tee $O/${W}_$SEG.check || { echo "$SYS $LEG $W: $SEG sanity check FAILED"; exit 1; }
if [ $ST = V ]; then
  mkdir -p sp; cd sp; cp -f $B/lambda_V.sch lambda.sch
  printf "parm $PRM\ntrajin ../${W}_$SEG.nc 1 last $FROFF\ntrajout fr.rst7 restart multi\nrun\n" > ext.in
  TARGETS=$(awk -F, 'NR>1 && $2=="V" && !s[$3]++ {print $3}' $LAMS)
  cpptraj -i ext.in > ext.log 2>&1 || { echo "cpptraj failed"; tail -5 ext.log; exit 1; }
  FR=$(ls fr.rst7.* 2>/dev/null | sort -t. -k3 -n); [ -n "$FR" ] || { echo "no frames extracted"; exit 1; }
  OUTF=$O/${W}_sp${SEG#s}_energies.txt; echo "# $SYS $LEG $W stage V $SEG sampled_lambda $LAM ; columns: frame target_lambda E_total_kcal DV/DL" > $OUTF.tmp
  for f in $FR; do for L in $TARGETS; do
    sed -e "s/@LAM@/$L/" -e "s/@MASK@/$MASK/g" $B/sp_V.in > sp.in
    pmemd.cuda -O -i sp.in -p $PRM -c $f -o sp.out > sp.stdout 2>&1 || { echo "single point failed ($f, $L)"; tail -3 sp.stdout; exit 1; }
    echo "${f##*.} $L $(awk -f $B/sp_energy.awk sp.out)" >> $OUTF.tmp
  done; done
  grep -q PARSE_ERROR $OUTF.tmp && { echo "single point parse error"; exit 1; }
  mv $OUTF.tmp $OUTF; echo "single points: $(grep -vc '^#' $OUTF) energies"
fi
touch $O/${W}_$SEG.done $O/${W}_10ns.done; echo "$SYS $LEG $W $SEG done -> 10 ns ($(hostname -s))"
