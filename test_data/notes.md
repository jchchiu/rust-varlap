to build release:
`cargo build -r`
then run at:
`./target/release/rust-varlap`
e.g.
`time ./target/release/rust-varlap -v /home/jch/genomic-data/HG00096_chr11/HG00096.chr11.vcf.gz -b /home/jch/genomic-data/HG00096_chr11/HG00096.chrom11.ILLUMINA.bwa.GBR.low_coverage.20120522.bam --varclass snv -c test_data/chr11-snv-rust-varlap.output.csv`

`time ./target/release/rust-varlap -v /home/jch/genomic-data/HG00096_chr11/HG00096.chr11.vcf.gz -b /home/jch/genomic-data/HG00096_chr11/HG00096.chrom11.ILLUMINA.bwa.GBR.low_coverage.20120522.bam --varclass indel -c test_data/chr11-indel-rust-varlap.output.csv`

varlap commands:
`source varlap_dev/bin/activate`
SNV
`time varlap --varclass SNV --format VCF -- /home/jch/genomic-data/HG00096_chr11/HG00096.chrom11.ILLUMINA.bwa.GBR.low_coverage.20120522.bam < /home/jch/genomic-data/HG00096_chr11/HG00096.chr11.vcf > chr11-snv.varlap.csv`
INDEL (with log)
`time varlap --log chr11-indel.log --varclass INDEL --format VCF -- /home/jch/genomic-data/HG00096_chr11/HG00096.chrom11.ILLUMINA.bwa.GBR.low_coverage.20120522.bam < /home/jch/genomic-data/HG00096_chr11/HG00096.chr11.vcf > chr11-indel.varlap.csv`


chr11 test data:
from:
http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/data/HG00096/alignment/
command used:
`bcftools mpileup --max-depth 100000 --no-BAQ -q 0 -Q 0 -f human_g1k_v37.fasta HG00096.chrom11.ILLUMINA.bwa.GBR.low_coverage.20120522.bam | bcftools call -mv -Oz -o HG00096.chr11.vcf.gz`
`tabix -p vcf HG00096.chr11.vcf.gz`

comparing results command
```
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
' /home/jch/git/rust-varlap/test_data/chr11-snv-rust-varlap.output.csv chr11-snv.varlap.csv > chr11-snv.diff.txt
```