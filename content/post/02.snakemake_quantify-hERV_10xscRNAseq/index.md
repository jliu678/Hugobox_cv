---
title: 🧬 A snakemake pipline-- quantify hERV, transgenes and genome from 10x scRNAseq
summary: A snakemake pipline to quantify hERV, transgenes and genome from 10x single-cell RNAseq data  
date: 2025-03-25
authors:
  - admin
tags:
  - Snakemake scRNAseq Workflow
  - EM algorithm
  - hERV
  - transgene
  - scRNAseq
  - 10X
image:
  caption: 'Image credit: [**The Creative Idea on Unsplash**](https://unsplash.com)'
---
## Github source code
Please see [github repo](https://github.com/jliu678/snakemake-pipline_quantify-hERV-trangene_10x-scRNAseq).

## What are hERVs and transgenes
Human Endogenous Retroviruses (hERVs) are ancient viral sequences embedded in the human genome. The roles of hERVs in gene regulation, immunity, development and cancer are under intense research. Transgenes including GFP, CRE, Luciferase, rtTA/Tet-On/Tet-Off and epitopically expressed genes are common targets to be quantified in sample from transgenic mouse models. 


## How to quantify them from 10x scRNAseq 
The EM algorithm is well-suited for quantifying them as elucidated in [my previous blog](aaaaaaaaa). Briefly, EM algorithm is advantageous in dealing with  multimaping, which is commonly seen for hERV and transgene quantification, because:
- mathematically, EM is well suitable for estimating parameters in Gaussian Mixture Models, which is similar to the model describing multimapping
- biologically, EM can takes into consideration the mapping probability and mapping bias, though Starsolo dose not use bias correction

[Starsolo](https://github.com/alexdobin/STAR/blob/master/docs/STARsolo.md#all-parameters-that-control-starsolo-output-are-listed-again-below-with-defaults-and-short-descriptions) is particularly handy for this task because it:
- implements EM algorithm
- produces nearly identical gene counts in the same format as the 10x CellRanger does by following CellRanger logic for cell barcode whitelisting and UMI deduplication
- is the foundation of alignment in CellRanger

## Why compatibility with CellRanger matters
The compatibility of STARsolo with CellRanger is important. While CellRanger is not a gold standard, a substantial and growing body of biological interpretations is derived from data processed by it. Aligning with CellRanger enables more meaningful comparisons with such datasets. This is particularly crucial when data integration is required, as no current method can fully correct for technical artifacts without distorting true biological signals.

## Snakemake Workflow Management
Snakemake is used to manage a workflow as illustrated below. Briefly, I generated star index and wrangled all fastqs belonging to a same end reading and a same sample into one fastq file. The combined fastqs and star index were input to starsolo in rule `star_mapping`. After making sure both input fastq wrangling and rule `star_mapping` were done, I moved the input fastqs to a specified directory to safely carry out cleaning up.
{{< figure
  src="dag.svg"
  alt=""
  link=""
  caption=""
  class="ma0 w-75"
>}}

Key lines of the Snakemake file are discussed below
### Import config.yaml
This line imports a configuration file in YAML format. The config file specify the genome and annotation files, the whitelist file(s), the star index directory and the CellBarcode/UMI lengths, etc. The cell barcode will be provided to `--soloCBwhitelist` option in STARsolo to:
- Match observed barcodes in your data,
- Correct sequencing errors (if close matches are allowed),
- Filter out barcodes that are not from real, expected cells.
```python
configfile: "config.yaml"
```
Part of values defined in the yaml are as below, other will be discussed otherwhere in the blog.
```yaml
genome: "combined.fa"
annotations: "combined.gtf"
wl: ["wl1.txt", "wl2.txt", "wl3.txt"]
star_index: "star/"
```

### Generates index for STARsolo.
 The name of the genome and annotation files, and the desired output directory are defined in the config file. Here we use `--genomeSAsparseD 3` to make the agreement between STARsolo and CellRanger even more perfect ([STARsolo link](https://github.com/alexdobin/STAR/blob/master/docs/STARsolo.md#how-to-make-starsolo-raw-gene-counts-almost-identical-to-cellrangers)). Output of this rule is required by `star_mapping`

```python
rule gen_star_idx:
    input: 
        expand("{genome}", genome=config["genome"]),
        expand("{annotations}", annotations=config["annotations"]),
    output:
        directory(config["star_index"])
    threads: 24
    resources:
        mem = 32000,
        time = "20:00:00",
        partition = "normal"
    message: 
        "generating star index..."
    shell:
        """STAR --runMode genomeGenerate --runThreadN {threads} \
            --genomeFastaFiles {input[0]} --sjdbGTFfile {input[1]} \
            --genomeSAsparseD 3 --genomeDir {output[0]}/"""
```

### Wrangle input fastqs
This rule combines reads from the `./raw_reads`. Because reads can be performed on different lanes, this rule combines them into one fa file, as STARsolo does not support multiple lanes per read. The outputs are produced in the `./combined_reads` dir. Output of this rule is required by `star_mapping`

```python
# use absolute dir to avoid instability
import os
cwd = os.getcwd() + "/raw_reads"
input_reads = [d for d in os.listdir(cwd) if os.path.isdir(os.path.join(cwd, d))]

rule combine_reads:
    input: 
        "raw_reads/{input_reads}/"
    output: 
        "combined_reads/{input_reads}.r1.fastq.gz",
        "combined_reads/{input_reads}.r2.fastq.gz",
    message:
        "combining reads..."
    run:
        shell("find -regextype sed -regex .*/{input}/.*[Rr]1.*\.fastq\.gz -exec cat {{}} + > {output[0]}")
        shell("find -regextype sed -regex .*/{input}/.*[Rr]2.*\.fastq\.gz -exec cat {{}} + > {output[1]}")
```

### Align and count genes, hERV and transgenes
Rule `star_mapping` calls STARsolo with a detailed configuration for single-cell RNA-seq analysis, using a customized CB/UMI structure and enabling EM-based handling of multimappers.

The soloType for 10x dropSeq is `CB_UMI_Simple` as excerpted below from [star-solo github](https://github.com/alexdobin/STAR/blob/master/docs/STARsolo.md#all-parameters-that-control-starsolo-output-are-listed-again-below-with-defaults-and-short-descriptions)

    - soloType input options:
        - CB_UMI_Simple   ... (a.k.a. Droplet) one UMI and one Cell Barcode of fixed length in read2, e.g. Drop-seq and 10X Chromium.
        - CB_UMI_Complex  ... one UMI of fixed length, but multiple Cell Barcodes of varying length, as well as adapters sequences are allowed in read2 only, e.g. inDrop.
       
10x barcode inclusion list (formerly barcode whitelist) can be downloaded [here](https://kb.10xgenomics.com/hc/en-us/articles/115004506263-What-is-a-barcode-inclusion-list-formerly-barcode-whitelist#:~:text=For%20example%2C%20there%20are%20roughly,AAACCTGAGAAAGTGG%20AAACCTGAGAACAACT%20AAACCTGAGAACAATC%20AAACCTGAGAACTCGG%20AAACCTGAGAACTGTA).

Importantly, in the --readFilesIn option, the 1st file has to be cDNA read, and the 2nd file has to be the barcode (cell+UMI) read, i.e.

10X data have a variety of CB and UMI lengths, I decided for them to be a configurable option. 
Please Ref to [barcode-geomery section of Star-Solo github](https://github.com/alexdobin/STAR/blob/master/docs/STARsolo.md#barcode-geometry) to quickly understand below arguments. 

1. 10x three prime: common values are below for 10x three prime, and accurate values can be found in 10x library protocol used to generate the data.
   ```yaml
   - soloCBstart 1
   - soloCBlen 16
   - soloUMIstart 17
   - soloUMIlen 10
   ```
   
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
                --soloType CB_UMI_Simple \
                --soloCBstart 1 --soloCBlen {config[cblen]} \
                --soloUMIstart {config[umistart]} --soloUMIlen {config[umilen]} \
                --soloBarcodeReadLength 1 --soloCBwhitelist {config[wl]} \
                --clipAdapterType CellRanger4 --outFilterScoreMin 10 \
                --soloCBmatchWLtype 1MM_multi_Nbase_pseudocounts \
                --soloUMIfiltering MultiGeneUMI_CR --soloUMIdedup 1MM_CR \
                --soloFeatures {soloFeaturesvalues} \
                --soloCellFilter EmptyDrops_CR \
                --soloMultiMappers EM \
                --readFilesIn {input[1]} {input[0]} --readFilesCommand 'gunzip -c'\
                --outFileNamePrefix {output}/"""
    ```
2. 10x five prime: common values are below for 10x five prime, and accurate values can be found in 10x library protocol used to generate the data.
   ```yaml
   - soloCBstart 1
   - soloCBlen 16
   - soloUMIstart 17
   - soloUMIlen 10
   ```

   Please read this [thread](https://github.com/alexdobin/STAR/issues/768) to achieve reasonable results for PE and SE alignment/quantification of the same 5' 2x150 bp samples. Altogether, the rule is
  
    ```python
    rule star_mapping:
      input:
          config["combined_reads_dir"] + "/{input_reads}.r1.fastq.gz",
          config["combined_reads_dir"] + "/{input_reads}.r2.fastq.gz",
          idx=expand("{star_index}", star_index=config["star_index"]),
      output: 
          directory(config["mapping_dir"] + "/{input_reads}/")
      message:
          "star mapping..."
      threads: 24
      resources:
          mem = 32000,
          time = '6:00:00'
      shell:
          """STAR --genomeDir {input.idx} --runThreadN {threads} \
              --soloType CB_UMI_Simple \
              --soloCBstart 1 --soloCBlen {config[cblen]} \
              --soloUMIstart {config[umistart]} --soloUMIlen {config[umilen]} \
              --soloBarcodeReadLength 1 --soloCBwhitelist {config[wl]} \
              --soloStrand Forward --outFilterScoreMin 30 \
              --soloCBmatchWLtype 1MM_multi_Nbase_pseudocounts \
              --soloUMIfiltering MultiGeneUMI_CR --soloUMIdedup 1MM_CR \
              --soloFeatures Gene GeneFull SJ Velocyto --soloCellFilter EmptyDrops_CR \
              --soloMultiMappers EM \
              --soloBarcodeMate 1   --clip5pNbases 39 0 \
              --readFilesIn {input[0]} {input[1]} --readFilesCommand 'gunzip -c'\
              --outFileNamePrefix {output}/"""
    ```