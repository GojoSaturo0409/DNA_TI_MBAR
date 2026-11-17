# total potential energy (sum of the 9 terms of the first energy block after RESULTS) and DV/DL from a pmemd 1-step minimisation output.
# If any term overflows the print field (******), the true energy is astronomically high -> print inf (weight 0 in MBAR).
/4\.  RESULTS/ {r=1}
r && /^ BOND/ {blk=1}
r && blk && n<3 { if ($0 ~ /= *\*+/) ovf=1; line=$0; while (match(line, /= *-?[0-9]+\.[0-9]+/)) { v=substr(line, RSTART+1, RLENGTH-1); s+=v; line=substr(line, RSTART+RLENGTH) } n++ }
r && /^ DV\/DL/ && !d { d=$3 }
END { if (n!=3) print "PARSE_ERROR"; else if (ovf) printf "inf %s\n", d; else printf "%.4f %s\n", s, d }
