# rust-varlap

rust-varlap is a rewrite of [varlap](https://github.com/bjpop/varlap) from python into rust.

## Assumptions
- The input variant files (vcf, csv, tsv) are sorted by variant chromosome and position in ascending order
- Variants are not filtered or interpreted
- The input read files (bam, cram) are sorted by read chromosome and position in ascending order

## Changes
- The algorithm now iterates over reads instead of variants
- Multiple bam inputs are not currently supported
- Changed header output field labels:
    - e.g. from [label + ' ' + 'ref avg nm'] to [label + ' ' + 'ref_avg_nm']
    - Easier to split header if necessary (split by ' ' will separate bam label and statistic field)

## To-do
- **Fix csv/tsv parsing (figure out best idiomatic way in rust to reduce resources/minimize friction when attempting multithreading/duplicate code)**
- **Add support for dealing with multiple chromosomes**
- Fix output csv header so it matches varlap exactly (?)
- Improve error handling and add logging
- Add unit tests/integration tests
- Add region scanning mode (split into subcommands, 'allele' and 'region')
- Implement multithreading
- Deal with paired reads (implement something similar to mosdepth?)
