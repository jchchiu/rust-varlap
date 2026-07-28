to build release:
`cargo build -r`
then run at:
`./target/release/rust-varlap`
e.g.
`time ./target/release/rust-varlap -v /home/jch/genomic-data/HG00096_chr11/HG00096.chr11.vcf.gz -b /home/jch/genomic-data/HG00096_chr11/HG00096.chrom11.ILLUMINA.bwa.GBR.low_coverage.20120522.bam --varclass snv -o test_data/chr11-snv-rust-varlap.output.csv`

`time ./target/release/rust-varlap -v /home/jch/genomic-data/HG00096_chr11/HG00096.chr11.vcf.gz -b /home/jch/genomic-data/HG00096_chr11/HG00096.chrom11.ILLUMINA.bwa.GBR.low_coverage.20120522.bam --varclass indel -o test_data/chr11-indel-rust-varlap.output.csv`

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

# HG00096 mapped exome test data:
### BAM/CRAM file from:
http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/data/HG00096/exome_alignment/

  FOR PYTHON:
  [ ]	HG00096.mapped.ILLUMINA.bwa.GBR.exome.20120522.bam	2015-04-30 15:31 	8.6G	 
  [ ]	HG00096.mapped.ILLUMINA.bwa.GBR.exome.20120522.bam.bai	2015-04-30 22:49 	6.5M
  `wget HG00096.mapped.ILLUMINA.bwa.GBR.exome.20120522.bam`
  `wget HG00096.mapped.ILLUMINA.bwa.GBR.exome.20120522.bam.bai`
  FOR RUST:
  [   ]	HG00096.mapped.ILLUMINA.bwa.GBR.exome.20120522.bam.cram	2015-05-05 01:18	2.1G	 
  [   ]	HG00096.mapped.ILLUMINA.bwa.GBR.exome.20120522.bam.cram.crai	2015-04-30 19:06	175K	 
  `wget HG00096.mapped.ILLUMINA.bwa.GBR.exome.20120522.bam.cram`
  `wget HG00096.mapped.ILLUMINA.bwa.GBR.exome.20120522.bam.cram.crai`

### FASTA REFERENCE FILE from the header of the cram:
  `wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/phase2_reference_assembly_sequence/hs37d5.fa.gz`
  Clear out any broken downloaded slices from your hidden cache folder
  `rm -rf ~/.cache/hts-ref`
  reindex
  `gunzip hs37d5.fa.gz`
  `samtools faidx hs37d5.fa`

### GENERATE the VCF:
  `bcftools mpileup --max-depth 100000 --no-BAQ -q 0 -Q 0 -f ~/genomic-data/HG00096_exome/hs37d5.fa HG00096.mapped.ILLUMINA.bwa.GBR.exome.20120522.bam.cram | bcftools call -mv -Oz -o HG00096.exome.vcf.gz`
  `tabix -p vcf HG00096.exome.vcf.gz`

### RUNNING varlap
  RUST COMMAND:
    `cargo build -r`
    `time ./target/release/rust-varlap -v /home/jch/genomic-data/HG00096_exome/HG00096.exome.vcf.gz -b /home/jch/genomic-data/HG00096_exome/HG00096.mapped.ILLUMINA.bwa.GBR.exome.20120522.bam.cram --varclass snv -o /home/jch/genomic-data/HG00096_exome/exome-rust-snv.csv --fasta-file /home/jch/genomic-data/HG00096_exome/hs37d5.fa`
  PYTHON COMMAND
    `source varlap_dev/bin/activate`
    `time varlap --varclass SNV --format VCF -- /home/jch/genomic-data/HG00096_exome/HG00096.mapped.ILLUMINA.bwa.GBR.exome.20120522.bam < /home/jch/genomic-data/HG00096_exome/HG00096.exome.vcf > exome-python-snv.varlap.csv`

### COMPARING CSV outputs
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
' /home/jch/genomic-data/HG00096_exome/exome-rust-snv.csv exome-python-snv.varlap.csv > exome-snv.diff.txt
```

NOTE ERROR FOR CRAM:
[W::cram_decode_slice] Slice ends beyond reference end at #23:58890239-59373566
[W::cram_decode_slice] Slice ends beyond reference end at #76:2-182896
[W::cram_decode_slice] Slice ends beyond reference end at #83:261949-547496

There was some difference when comparing the CRAM rust version to the BAM python version; the avg_base_qual scores were different; however, when using the BAM rust version the outputs are the same
- Most likely due to lossy CRAM files?: (The CRAM files the 1000 Genomes project distributes are lossy cram files which reduce the base quality scores using the Illumina 8-bin compression scheme as described in the lossy compression section on the cram usage page)[https://www.internationalgenome.org/category/cram/]

### Trying Lossless BAM -> CRAM conversion

`samtools index -c HG00096.mapped.ILLUMINA.bwa.GBR.exome.20120522.lossless.bam.cram`

Running the rust varlap version again results in the same output as the python version; therefore, most likely reason for discrepancy is that the CRAM files are lossy; also the errors when decoding slice though I am not as sure what the finer details are re: this error.


## NA18623

https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/data/NA18623/alignment/

https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/data/NA18623/alignment/NA18623.chrom20.ILLUMINA.bwa.CHB.low_coverage.20130415.bam

https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/data/NA18623/alignment/NA18623.chrom20.ILLUMINA.bwa.CHB.low_coverage.20130415.bam.bai

`bcftools mpileup --max-depth 100000 --no-BAQ -f ~/genomic-data/HG00096_exome/hs37d5.fa NA18623.chrom20.ILLUMINA.bwa.CHB.low_coverage.20130415.bam | bcftools call -mv -Ov -o NA18623.chrom20.vcf`
`tabix -p vcf NA18623.chrom20.vcf`

### RUNNING varlap
RUST COMMAND:
  `cargo build -r`
  `time ./target/release/rust-varlap -v ~/genomic-data/NA18623/NA18623.chrom20.sorted.vcf -r ~/genomic-data/NA18623/NA18623.chrom20.ILLUMINA.bwa.CHB.low_coverage.20130415.bam --varclass snv -o ~/genomic-data/NA18623/NA18623-chr20-rust-snv.csv`
PYTHON COMMAND
  `source varlap_dev/bin/activate`
  `time varlap --varclass SNV --format VCF -- ~/genomic-data/NA18623/NA18623.chrom20.ILLUMINA.bwa.CHB.low_coverage.20130415.bam < ~/genomic-data/NA18623/NA18623.chrom20.sorted.vcf > ~/genomic-data/NA18623/NA18623-chr20-python-snv.varlap.csv`
