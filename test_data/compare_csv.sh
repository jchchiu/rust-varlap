cd test_data/
awk -F',' '
NR==FNR {
  for (i=1; i<=NF; i++) a[FNR,i]=$i
  max_nf[FNR]=NF
  next
}
{
  max = (NF > max_nf[FNR] ? NF : max_nf[FNR])
  for (i=1; i<=max; i++) {
    if (a[FNR,i] != $i) {
      print "Row " FNR ", Col " i ": file1=\"" a[FNR,i] "\", file2=\"" $i "\""
    }
  }
}
' rust-varlap.output.csv variants_chrM_heavy_stress.varlap.csv > diff.txt