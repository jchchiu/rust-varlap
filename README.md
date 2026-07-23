# rust-varlap

rust-varlap is a rewrite of [varlap](https://github.com/bjpop/varlap) from python into rust.

## Assumptions
- The input variant files (vcf, csv, tsv) are sorted by variant chromosome and position in ascending order
- Variants are not filtered or interpreted
- The input read files (bam, cram) are sorted by read chromosome and position in ascending order

## Changes
- The algorithm now iterates over reads instead of variants
    - Instead of calling pileup over variants we are now fetching reads for variants within a given range
    - Variants are binned by chromosome, and further binned if variants are some `n` distance apart (gap)
    - From preliminary testing, we have set the default gap size as 100 kb (100,000 bp) as it balances performance between files with sparsely and densely populated variants
    - An optional hyperparameter has can be used to control this gap size between bins; generally speaking, sparsely populated variants throughout chromosomes may perform better if gap size is small (e.g. `5 x mean read length`); conversely, for densely populated variants all throughout the chromosome performance may be better if gap size is large (e.g. 1 Mb or even more)
- Multiple bam inputs are not currently supported
    - Can be implemented if needed; however, output will most likely be separate csvs instead of a combined csv
- Changed header output field labels:
    - e.g. from [label + ' ' + 'ref avg nm'] to [label + ' ' + 'ref_avg_nm']
    - Easier to split header if necessary (split by ' ' will separate bam label and statistic field)
- CRAM files are now supported
- gzipped input variant files are now supported
- Region mode with bed files/Outliers mode is not currently supported

## Notes
- For csv files, multiple alt alleles are not supported (e.g. `A,T`). Please separate them into multiple rows (i.e. row for `A` alt, row for `T` alt)

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
- Add option to sort variants if unsorted
- Make ref/alt input optional (Still output stats that do not need this info)
