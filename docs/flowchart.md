```mermaid
flowchart TD
    subgraph Variant Handling

    var_file@{ shape: manual-file, label: "FILE INPUT \n Variant file \n (vcf, csv, tsv) (.gz)"}

    var_file --> var_parse@{ shape: rect, label: "VARIANT PARSER" }
    var_class@{ shape: manual-input, label: "USER INPUT \n Variant Class \n (Indel, Snv)"} --> var_parse
    var_parse -->|if variant passes QC| variants@{ shape: docs, label: "Variant info \n [Struct] \n {chrom, pos, ref, alt, vartype, features}" }
    variants -.->|check for unique chromosomes| chrom_unique@{ shape: bow-rect, label: "Chromosome info \n [Hashmap] \n {unique_chroms: first index in variants vector}" }
    variants -- add variant info to queue--> var_queue@{ shape: bow-rect, label: "Variants Vector \n [Vector Queue] \n {Variant info}" }

    end

    subgraph BAM Handling

    var_queue --> bam_parse@{ shape: rect, label: "BAM PARSER"}
    chrom_unique -.-> bam_parse
    var_class@{ shape: manual-input, label: "USER INPUT \n Variant Class \n (Indel, Snv)"}

    reads_file@{ shape: manual-file, label: "FILE INPUT \n Reads file \n (bam, cram)"}
    fasta_file@{ shape: manual-file, label: "FILE INPUT \n Fasta file \n (fasta/fa)"}

    reads_file --> bam_parse
    fasta_file -->|only required for cram| bam_parse
    bam_parse --> variant_conditional@{ shape: hex, label: "If read start > \n variant position" }
    variant_conditional --> remove_variant@{ shape: diamond, label: "Pop variant from queue" }
    remove_variant --> variant_conditional

    end


    scanning_mode@{ shape: manual-input, label: "USER INPUT \n Scanning Mode \n (Allele , Region)"} -- Region --> region_parse@{ shape: rect, label: "REGION PARSER" }
    bed_file@{ shape: manual-file, label: "FILE INPUT \n Region file \n (bed)"}
    bed_file --> region_parse
    scanning_mode -- Allele --> fetch_reads@{ shape: subproc, label: "Fetch reads for a chromosome"}
    region_parse --> fetch_reads

```

EXPANDED VERSION

```mermaid
flowchart TD

    bam_file@{ shape: manual-file, label: "BAM file"}
    var_file@{ shape: manual-file, label: "FILE INPUT \n Variant file \n (vcf, csv, tsv)"}
    var_file --> var_parse@{ shape: rect, label: "Variant parser" }
    var_parse --> file_checker@{ shape: subproc, label: "File format checker" }
    file_checker --> gzip_checker@{ shape: subproc, label: "Gzip checker" }
    gzip_checker --> var_queue@{ shape: bow-rect, label: "Variants Queue" }
    var_queue --> bam_parse@{ shape: rect, label: "Bam parser"}
    bam_file --> bam_parse
    mode@{ shape: manual-input, label: "USER INPUT \n (Indel, Snv)"} --> var_parse
```