# Sourced by stageA.sbatch. Heating in restart segments, robust to the GPU "box dimensions changed too much" halt:
#  - first <cpu_steps> on a CPU AMBER build if one exists (falls back to GPU if the CPU run fails),
#  - then GPU segments; every segment restarts from the previous restart file (regenerates GPU grid cells);
#  - if a GPU segment fails, it is split in half and retried from the last good restart (down to MINSEG steps).
# usage: heat_segmented <prefix> <cpu_steps> <gpu_seg_steps...>   (needs PRM, R, min.rst7 in cwd; writes <prefix>.rst7 / <prefix>.out)
heat_segmented() {
  local pre=$1 cpu_steps=$2; shift 2
  local IN=$R/amber_inputs TPL="${HEAT_TEMPLATE:-$R/amber_inputs/heat.in}" NST="${HEAT_NSTLIM:-1200000}" prev="${HEAT_START:-min.rst7}" k=0 first=1 n f t0 CPU="" retries=0 MAXRETRY=40 MINSEG=50
  if   command -v pmemd.MPI   >/dev/null 2>&1; then CPU="mpirun --map-by :OVERSUBSCRIBE -np ${SLURM_CPUS_PER_TASK:-2} pmemd.MPI"
  elif command -v pmemd       >/dev/null 2>&1; then CPU="pmemd"
  elif command -v sander.MPI  >/dev/null 2>&1; then CPU="mpirun --map-by :OVERSUBSCRIBE -np ${SLURM_CPUS_PER_TASK:-2} sander.MPI"
  elif command -v sander      >/dev/null 2>&1; then CPU="sander"; fi
  echo "HEAT cpu_binary='${CPU:-none}' cpu_steps=$cpu_steps gpu_segments='$*'"
  mkseg() { # steps first? -> $f.in   (first=1: fresh start from min.rst7; first=0: continue velocities)
    local sedargs=(-e "s/nstlim=$NST,/nstlim=$1,/")
    [ "$2" = 0 ] && sedargs+=(-e "s/ntx=1,/ntx=5,/" -e "s/irest=0,/irest=1,/")
    sed "${sedargs[@]}" "$TPL" > "$f.in"
    grep -q "nstlim=$1," "$f.in" || { echo "HEAT nstlim substitution failed ($f)"; return 1; }
    [ "$2" = 0 ] && { grep -q "irest=1," "$f.in" && grep -q "ntx=5," "$f.in" || { echo "HEAT restart substitution failed ($f)"; return 1; }; }
    return 0; }
  if [ "$cpu_steps" -gt 0 ]; then
    if [ -n "$CPU" ]; then
      f=${pre}_s00; mkseg "$cpu_steps" 1 || return 1; t0=$SECONDS
      if timeout 900 $CPU -O -i $f.in -p "$PRM" -c $prev -o $f.out -r $f.rst7 -x $f.nc -inf $f.mdinfo; then
        echo "HEAT cpu segment: $cpu_steps steps in $((SECONDS-t0)) s"; prev=$f.rst7; first=0
      else
        echo "HEAT cpu segment FAILED (rc=$?) after $((SECONDS-t0)) s - falling back to GPU for these $cpu_steps steps"
        rm -f $f.rst7; set -- "$cpu_steps" "$@"
      fi
    else
      set -- "$cpu_steps" "$@"          # no CPU build: run the CPU steps as a short first GPU segment
    fi
  fi
  local Q=("$@") total=0
  while [ ${#Q[@]} -gt 0 ]; do
    n=${Q[0]}; Q=("${Q[@]:1}")
    k=$((k+1)); f=$(printf "%s_s%02d" "$pre" "$k"); mkseg "$n" "$first" || return 1; t0=$SECONDS
    if pmemd.cuda -O -i $f.in -p "$PRM" -c $prev -o $f.out -r $f.rst7 -x $f.nc -inf $f.mdinfo; then
      total=$((total+n)); echo "HEAT gpu segment $k: $n steps in $((SECONDS-t0)) s (steps done on GPU: $total)"; prev=$f.rst7; first=0
    else
      mv -f $f.out $f.failed 2>/dev/null; rm -f $f.rst7
      if [ "$n" -ge $((2*MINSEG)) ] && [ $retries -lt $MAXRETRY ]; then
        retries=$((retries+1)); h=$((n/2)); echo "HEAT gpu segment $k ($n steps) FAILED - splitting into $h + $((n-h)) and retrying (retry $retries)"
        Q=("$h" "$((n-h))" "${Q[@]}")
      else
        echo "HEAT gpu segment $k ($n steps) failed and cannot be split further (retries=$retries)"; return 1
      fi
    fi
  done
  cp -f "$prev" "$pre.rst7" && cat "${pre}"_s*.out > "$pre.out"
}
