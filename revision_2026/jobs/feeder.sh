#!/bin/bash
# prod2c feeder: submits queued windows (one sbatch per window) while the account has < LIMIT jobs in squeue.
# ~/prod2c/feeder/queue.txt lines: SYS LEG WIN NICE (in order). Submitted lines go to submitted.txt with the job id.
# FEED_DRY=1 only prints. Removes its own cron entry when the queue is empty. Run from cron under flock.
set -u
LIMIT=${LIMIT:-19}
D=$HOME/prod2c/feeder; Q=$D/queue.txt; EXC='gnode[001-008,011-012,017,020-023,025-029,031,033-036,038,040-044,061,066,068,069,076,077,090]'
cd $HOME/prod2c || exit 1; mkdir -p logs
log(){ echo "$(date '+%F %T') $*"; }
[ -s $Q ] || [ -s $D/submitted.txt ] || { log "queue and submitted list empty, removing cron entry"; crontab -l 2>/dev/null | grep -v 'prod2c/feeder/feeder.sh' | crontab -; exit 0; }
touch $Q
# ---- retry: a submitted window that left the queue without its .done marker is re-queued ONCE (then logged as given up) ----
if live=$(timeout 120 squeue -u $USER -h -o %i 2>/dev/null); then
  if [ -s $D/submitted.txt ]; then
    : > $D/submitted.keep
    while read -r S1 L1 W1 N1 J1; do
      [ -z "${J1:-}" ] && continue
      if echo "$live" | grep -qx "$J1"; then echo "$S1 $L1 $W1 $N1 $J1" >> $D/submitted.keep; continue; fi
      if [ -f $S1/$L1/runs/$(printf w%02d $W1).done ]; then echo "$S1 $L1 $W1 $N1 $J1" >> $D/completed.txt; continue; fi
      n=$(grep -c "^$S1 $L1 $W1\$" $D/retries.txt 2>/dev/null); n=${n:-0}
      if [ "$n" -lt 1 ]; then echo "$S1 $L1 $W1 $N1" >> $Q; echo "$S1 $L1 $W1" >> $D/retries.txt; log "RETRY $S1 $L1 w$(printf %02d $W1) (job $J1 ended without .done)"
      else echo "$S1 $L1 $W1 $N1 $J1" >> $D/failed.txt; log "GAVE UP $S1 $L1 w$(printf %02d $W1) (job $J1, second failure)"; fi
    done < $D/submitted.txt
    [ "${FEED_DRY:-0}" = 1 ] || mv $D/submitted.keep $D/submitted.txt
  fi
else log "squeue failed - skipping this round (no retries, no submissions)"; exit 0; fi
cur=$(squeue -u $USER -h | wc -l)
cp $Q $Q.new; : > $Q.keep
while read -r SYS LEG WIN NICE; do
  [ -z "${SYS:-}" ] && continue
  if [ $cur -ge $LIMIT ] || [ -f $SYS/$LEG/runs/$(printf w%02d $WIN).done ]; then
    [ -f $SYS/$LEG/runs/$(printf w%02d $WIN).done ] || echo "$SYS $LEG $WIN $NICE" >> $Q.keep; continue; fi
  n=$(printf %02d $WIN)
  cmd="sbatch --parsable --nice=$NICE --exclude=$EXC -J P3_${SYS}_${LEG}_w$n -o logs/P3_${SYS}_${LEG}_w${n}_%j.log --export=ALL,SYS=$SYS,LEG=$LEG,WIN=$WIN prod2c.sbatch"
  if [ "${FEED_DRY:-0}" = 1 ]; then log "[dry] cur=$cur $cmd" | cut -c1-200; echo "$SYS $LEG $WIN $NICE" >> $Q.keep; cur=$((cur+1)); continue; fi
  if id=$(eval $cmd 2>&1) && [[ $id =~ ^[0-9]+$ ]]; then log "queued $SYS $LEG w$n -> $id"; echo "$SYS $LEG $WIN $NICE $id" >> $D/submitted.txt; cur=$((cur+1))
  else log "sbatch failed for $SYS $LEG w$n: $id"; echo "$SYS $LEG $WIN $NICE" >> $Q.keep; fi
done < $Q.new
[ "${FEED_DRY:-0}" = 1 ] || mv $Q.keep $Q
rm -f $Q.new $Q.keep 2>/dev/null; true
