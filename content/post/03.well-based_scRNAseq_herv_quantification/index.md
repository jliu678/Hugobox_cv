---
title: 🧬 A snakemake pipline-- quantify hERV, transgenes and genome from well-based single-cell RNAseq
summary: Well based scRNAseq has many more details than 10x scRNAseq to pay attention to in order to quantify features using Starsolo. Here, we discuss and present the snakemake pipline to quantify hERV, transgenes and genome from well-based single-cell RNAseq data  
date: 2025-03-26
authors:
  - admin
tags:
  - scRNAseq Workflow Algorithm
  - hERV
  - transgene
  - EM
  - scRNAseq
  - well-based
image:
  caption: 'Image credit: [**Heather Green on Unsplash**](https://unsplash.com)'
---

## Github source code
Please see [github repo](https://github.com/jliu678/snakemake-pipline_quantify-hERV-trangene_well-based-scRNAseq).

## What are hERVs and transgenes
Human Endogenous Retroviruses (hERVs) are ancient viral sequences embedded in the human genome. The roles of hERVs in gene regulation, immunity, development and cancer are under intense research. Transgenes including GFP, CRE, Luciferase, rtTA/Tet-On/Tet-Off and epitopically expressed genes are common targets to be quantified in sample from transgenic mouse models. 


## How to quantify them from Well-based scRNAseq 
The EM algorithm is well-suited for quantifying them as elucidated in [my previous blog](aaaaaaaaa). Briefly, EM algorithm is advantageous in dealing with  multimaping, which is commonly seen for hERV and transgene quantification, because:
- mathematically, EM is well suitable for estimating parameters in Gaussian Mixture Models, which is similar to the model describing multimapping
- biologically, EM can takes into consideration the mapping probability and mapping bias, though Starsolo dose not use bias correction

[Starsolo](https://github.com/alexdobin/STAR/blob/master/docs/STARsolo.md#all-parameters-that-control-starsolo-output-are-listed-again-below-with-defaults-and-short-descriptions) is particularly handy for this task because it:
- implements EM algorithm
- follows CellRanger logic for cell barcode whitelisting and UMI deduplication to produce nearly identical gene counts in the same format as the CellRanger does 
- is the foundation of alignment in CellRanger

## Why resembling CellRanger helps
The similarity of STARsolo with CellRanger is important. While CellRanger is not a gold standard, a substantial and growing body of biological interpretations is derived from data processed by it. Aligning with CellRanger enables more meaningful comparisons with such datasets. This is particularly crucial when data integration is required, as no current method can fully correct for technical artifacts without distorting true biological signals.

## Snakemake Workflow Management
Snakemake is used to manage a workflow as illustrated below. Briefly, I generated star index and wrangled all fastqs belonging to the same end reading and the same sample into one fastq file. The combined fastqs and star index were input to starsolo in rule `star_mapping`. After making sure both input fastq wrangling and rule `star_mapping` were done, I moved the input fastqs to a specified directory to safely carry out cleaning up.
{{< figure
  src="dag.svg"
  alt=""
  link=""
  caption=""
  class="ma0 w-75"
>}}

Compared with that for 10x scRNAseq, the Snakemake file to quantify hERV, transgenes and genome from well-based single-cell RNAseq data has more details to pay attention to. Below discussed the key unique lines of the Snakemake file for well-based single-cell RNAseq, and please see the common parts of the Snakemake files for both well-based and 10x scRNAseq data in [my previous blog](10xscrnaseq).  

### Align and count genes, hERV and transgenes
Rule `star_mapping` calls STARsolo with a detailed configuration for single-cell RNA-seq analysis, using a customized CB/UMI structure and enabling EM-based handling of multimappers.

The soloType for seqWell(well-based scRNAseq) is `CB_UMI_Complex` as excerpted below from [star-solo github](https://github.com/alexdobin/STAR/blob/master/docs/STARsolo.md#all-parameters-that-control-starsolo-output-are-listed-again-below-with-defaults-and-short-descriptions)

    - soloType input options:
        - CB_UMI_Simple   ... (a.k.a. Droplet) one UMI and one Cell Barcode of fixed length in read2, e.g. Drop-seq and 10X Chromium.
        - CB_UMI_Complex  ... one UMI of fixed length, but multiple Cell Barcodes of varying length, as well as adapters sequences are allowed in read2 only, e.g. inDrop.
       
Microwell-seq uses three separate barcodes with adaptors aside them rather than a single barcode. Generating their three white lists requires reading through the methods in [Mapping the Mouse Cell Atlas by Microwell-Seq](https://doi.org/10.1016/j.cell.2018.02.001). NOTE the results only showed 5' to 3' strands, while barcodes 2 and 3 were presented on the 3' to 5' strands and reversed when being transcribed onto the beads, we have to take the complement and reverse the resultant sequence. Below successfully accomplished the three white lists step by step:

- Download table S1 of paper [Mapping the Mouse Cell Atlas by Microwell-Seq](https://doi.org/10.1016/j.cell.2018.02.001) containing the sequences of well-specific barcodes. The 1st barcode is one of A1-96; the 2nd barcode is one of B1-96; the 3rd barcode is one of C1-96
- Take A1 `TTTAGGGATAACAGGGTAATAAGCAGTGGTATCAACGCAGAGTACGTTTTAGGCGACTCACTACAGGG` as example. We first remove adaptor whose sequence can be seen in the ["Step-by-step generation of Indexed-Beads-seqA:
" section of Microwell-seq github introduction](https://teichlab.github.io/scg_lib_structs/methods_html/Microwell-seq.html), which is the sequences beside [barcode1] in `5'-TTTAGGGATAACAGGGTAATAAGCAGTGGTATCAACGCAGAGTACGT[barcode1]CGACTCACTACAGGG -3'`
- After removing adaptor the barcode inside A1 is "TTTTAGG". That is all for A1.
- But for B1-96 and C1-96 barcodes, we need remove adaptor as above and further get the reverse complementary of them. Please note the direction of sequence listed in the table S1 downloaded above and that of the adaptor in  ["Step-by-step generation of Indexed-Beads-seqA:" section of Microwell-seq github introduction](https://teichlab.github.io/scg_lib_structs/methods_html/Microwell-seq.html). Also check directions of the raw reads in [sra of GSM3980130: Adult-Cervix1; Homo sapiens; RNA-Seq (SRR9843415)](https://trace.ncbi.nlm.nih.gov/Traces/?view=run_browser&page_size=10&acc=SRR9843415&display=reads)
- `soloCBposition` and `soloUMIposition` of Starsolo needs input in this pattern (0_0_0_5) to know the start and end of barcodes and UMI. See ["Library sequencing:" section of Microwell-seq github introduction](https://teichlab.github.io/scg_lib_structs/methods_html/Microwell-seq.html) and ["Complex barcodes:" section of Star-Solo github introduction](https://github.com/alexdobin/STAR/blob/master/docs/STARsolo.md#complex-barcodes). For example: 0_0_0_5 means the barcode starts at bp 0 and ends at bp 5.
- Because there are 3 seperate barcodes, the wl option in "config.yaml" needs to specify 3 different files corresponding to the 3 barcode segments.
- Importantly, in the `--readFilesIn` option, the 1st file has to be cDNA read, and the 2nd file has to be the barcode (cell+UMI) read.

The rest of the configurations are copied from the [STARsolo guide on github](https://github.com/alexdobin/STAR/blob/master/docs/STARsolo.md#all-parameters-that-control-starsolo-output-are-listed-again-below-with-defaults-and-short-descriptions). That gives
```python
rule star_mapping:
    input:
        "combined_reads/{input_reads}.r1.fastq.gz",
        "combined_reads/{input_reads}.r2.fastq.gz",
        idx=expand("{star_index}", star_index=config["star_index"]),
    output: 
        directory("mapped_reads/{input_reads}/")
    message:
        "star mapping..."
    threads: 24
    resources:
        mem = 64000,
        time = '4:00:00'
    shell:
        """STAR --genomeDir {input.idx} --runThreadN {threads} \
            --soloType CB_UMI_Complex \
            --soloCBposition 0_0_0_5 0_21_0_26 0_42_0_47 \
            --soloUMIposition 0_48_0_53 \
            --soloCBwhitelist {config[wl]} \
            --clipAdapterType CellRanger4 --outFilterScoreMin 10 \
            --soloCBmatchWLtype 1MM_multi \
            --soloUMIfiltering MultiGeneUMI_CR --soloUMIdedup 1MM_CR \
            --soloFeatures Gene GeneFull SJ --soloCellFilter EmptyDrops_CR \
            --soloMultiMappers EM \
            --readFilesIn {input[1]} {input[0]} --readFilesCommand 'gunzip -c'\
            --outFileNamePrefix {output}/"""
```

