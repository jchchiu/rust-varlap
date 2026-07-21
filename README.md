# rust-varlap

rust-varlap is a rewrite of [varlap](https://github.com/bjpop/varlap) from python into rust.

## Assumptions
- The input variant files (vcf, csv, tsv) are sorted by variant chromosome and position in ascending order
- Variants are not filtered or interpreted
- The input read files (bam, cram) are sorted by read chromosome and position in ascending order

## Changes
- The algorithm now iterates over reads instead of variants
- Multiple bam inputs are not currently supported (Can be implemented if needed; however, output will most likely be separate csvs instead of a combined csv)
- Changed header output field labels:
    - e.g. from [label + ' ' + 'ref avg nm'] to [label + ' ' + 'ref_avg_nm']
    - Easier to split header if necessary (split by ' ' will separate bam label and statistic field)
- CRAM files are now supported
- gzipped input variant files are supported

## Notes
- For csv files, multiple alt alleles are not supported (e.g. A, T). Please separate them into multiple rows (i.e. row for A, row for T)

## To-do
- ~~Fix csv/tsv parsing (figure out best idiomatic way in rust to reduce resources/minimize friction when attempting multithreading/duplicate code)~~
- ~~Add support for dealing with multiple chromosomes~~
- ~~Fix output csv header so it matches varlap exactly (?)~~
    - See changes for further info
- ~~Improve error handling and add logging~~
- Add unit tests/integration tests
- Add region scanning mode (split into subcommands, 'allele' and 'region')
- Implement multithreading
- Deal with paired reads (implement something similar to mosdepth?)
