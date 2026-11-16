# Sanity checks on AMBER output (exit codes alone miss blown-up systems: pmemd.cuda can exit 0 at 10^5 K).
chk_md() { # out tmin tmax [dmin dmax]
  local o=$1 t m d
  # only instantaneous energy blocks: stop at the end-of-run "A V E R A G E S" / "R M S" blocks (their TEMP/Density are not instantaneous values)
  t=$(awk '/A V E R A G E S/{exit} /NSTEP =.*TEMP\(K\) *=/{print}' "$o" | tail -20 | sed -E 's/.*TEMP\(K\) *= *([^ ]+).*/\1/')
  [ -n "$t" ] || { echo "CHECK FAILED $o: no energy lines"; return 1; }
  case "$t" in *'*'*|*[Nn]a[Nn]*) echo "CHECK FAILED $o: overflow/NaN in TEMP"; return 1;; esac
  m=$(echo "$t" | awk '{s+=$1;n++} END{printf "%.1f", s/n}')
  awk -v m="$m" -v a="$2" -v b="$3" 'BEGIN{exit !(m>=a && m<=b)}' || { echo "CHECK FAILED $o: mean TEMP $m K outside [$2,$3]"; return 1; }
  if [ -n "${4:-}" ]; then
    d=$(awk '/A V E R A G E S/{exit} /Density/{print}' "$o" | tail -1 | sed -E 's/.*Density *= *([^ ]+).*/\1/')
    case "$d" in ''|*'*'*) echo "CHECK FAILED $o: bad density '$d'"; return 1;; esac
    awk -v d="$d" -v a="$4" -v b="$5" 'BEGIN{exit !(d>=a && d<=b)}' || { echo "CHECK FAILED $o: density $d outside [$4,$5]"; return 1; }
  fi
  echo "CHECK OK $o: mean TEMP $m K${d:+, density $d}"
}
chk_min() { # out  (final energy must be finite and below -2e5 kcal/mol for this 112k-atom system)
  local o=$1 e
  e=$(grep -A1 "NSTEP *ENERGY" "$o" | grep -vE "NSTEP|^--" | tail -1 | awk '{print $2}')
  case "$e" in ''|*'*'*|*[Nn]a[Nn]*) echo "CHECK FAILED $o: bad final minimization energy '$e'"; return 1;; esac
  awk -v e="$e" 'BEGIN{exit !(e < -2e5)}' || { echo "CHECK FAILED $o: final energy $e not below -2e5"; return 1; }
  echo "CHECK OK $o: final minimization energy $e"
}
