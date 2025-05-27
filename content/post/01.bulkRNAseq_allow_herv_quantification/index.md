---
title: 🧬 Quantify hERV and transgenes from Bulk-RNAseq
summary: Human Endogenous Retroviruses (hERVs) are ancient viral sequences embedded in the human genome. Transgenes are common in transgenic mouse models. To quantify them from sequencing reads, we need-- a) modify fasta and gtf files to include the their sequencing and annotation; b) a feature quantification algorithm to handle multimapping commonly seen for hERV and transgenes. Here I discussed the algorithms for feature quantification, and successfully quantified hERV and transgenes by implementing an EM algorithm.  
date: 2025-03-24
authors:
  - admin
tags:
  - hERV
  - transgene
  - feature quantification
  - algorithm
  - Expectation-Maximization (EM)
  - EM with MAP
image:
  caption: 'Image credit: [**Illustration by Round Icons on unsplash**](https://unsplash.com)'
---
## Github source code
Please see source codes in [github](https://github.com/jliu678/herv_project_siyi/tree/main) that can be deploied in Cloud Cluster Computational platform and local desktops, and allows for: 
- batch computaion with tunable batch size
- automatic environment setting up
- file preprocessing
- QC
- index buidling, alignment and feature counting by either `subread` or `salmon`  

The github repo also contains [usage examples using `SLURM` or `LSF` job scheduler](https://github.com/jliu678/herv_project_siyi/tree/main/example_usage) and [handy utility tools](https://github.com/jliu678/herv_project_siyi/tree/main/utils).
## What are hERVs and transgenes
Human Endogenous Retroviruses (hERVs) are ancient viral sequences embedded in the human genome. The roles of hERVs in gene regulation, immunity, development and cancer are under intense research. Transgenes including GFP, CRE, Luciferase, rtTA/Tet-On/Tet-Off and epitopically expressed genes are common targets to be quantified in sample from transgenic mouse models. 

## How to quantify them
To quantify them, we of course need modify fasta and gtf files to include the their sequencing and annotation. Also, because their sequencing reads are frequently repetitive or overlaping with host genome, we need a feature quantification algorithm to handle multimapping commonly seen for HERV and transgenes.


### Their fasta and gtf 
The transgenes can be included as chromosomes in fasta as below
```bash
> EGFP
ATGGTGAGCAAGGGCGAGGAGCTGTTCACCGGGGTGGTGCCCATCCTGGTCGAGCTGGACGGCG
ACGTAAACGGCCACAAGTTCAGCGTGTCCGGCGAGGGCGAGGGCGATGCCACCTACGGCAAGCT
GACCCTGAAGTTCATCTGCACCACCGGCAAGCTGCCCGTGCCCTGGCCCACCCTCGTGACCACC
CTGACCTACGGCGTGCAGTGCTTCAGCCGCTACCCCGACCACATGAAGCAGCACGACTTCTTCA
AGTCCGCCATGCCCGAAGGCTACGTCCAGGAGCGCACCATCTTCTTCAAGGACGACGGCAACTA
CAAGACCCGCGCCGAGGTGAAGTTCGAGGGCGACACCCTGGTGAACCGCATCGAGCTGAAGGGC
ATCGACTTCAAGGAGGACGGCAACATCCTGGGGCACAAGCTGGAGTACAACTACAACAGCCACA
ACGTCTATATCATGGCCGACAAGCAGAAGAACGGCATCAAGGTGAACTTCAAGATCCGCCACAA
CATCGAGGACGGCAGCGTGCAGCTCGCCGACCACTACCAGCAGAACACCCCCATCGGCGACGGC
CCCGTGCTGCTGCCCGACAACCACTACCTGAGCACCCAGTCCGCCCTGAGCAAAGACCCCAACG
AGAAGCGCGATCACATGGTCCTGCTGGAGTTCGTGACCGCCGCCGGGATCACTCTCGGCATGGA
CGAGCTGTACAAGTAA
```

The transgenes are formated as below in gtf
```bash
EGFP    self    exon    1       720     500     +       .       gene_id "custom_3"; gene_name "EGFP"
Cre     self    exon    1       1023    500     +       .       gene_id "custom_4"; gene_name "Cre"
ERT2    self    exon    1       931     500     +       .       gene_id "custom_5"; gene_name "ERT2"
```

The hERV sequences are already in reference genome, but their gtf need be added as below. We will talk about **[how to generate the hERV gtf matching the latest reference genome in another blog](willbecomposedsoon)**.
```bash
chr10   RepeatMasker    exon    3083835 3085549 4.5     -       .       gene_id "RepMasker_dc1fbd39-e43b-4ac4-90ef-1da20b1d6249"; gene_name "chr10:LTR/ERVL-MaLR:ORR1A0:3"
chr10   RepeatMasker    exon    3114480 3115385 6       -       .       gene_id "RepMasker_c49d8024-0b45-46a9-bcaa-c44f7a305a03"; gene_name "chr10:LTR/ERVK:RLTR13F:8"
chr10   RepeatMasker    exon    3122665 3123530 13.8    +       .       gene_id "RepMasker_5154950a-1e53-4528-a5cf-37d121e3e1a4"; gene_name "chr10:LTR/ERVK:RLTR13D6:11"
chr11   RepeatMasker    exon    3055489 3061982 25.9    +       .       gene_id "RepMasker_f6a34047-357d-4ed5-9744-eaf4c862be26"; gene_name "chr11:LTR/ERVK:IAPLTR3-int:1"
chr11   RepeatMasker    exon    3139464 3140302 11.7    +       .       gene_id "RepMasker_ad8b23ab-de54-4f73-be4d-f8405ad13e66"; gene_name "chr11:LTR/ERVK:RLTR33:3"
chr11   RepeatMasker    exon    3410867 3413285 16.9    +       .       gene_id "RepMasker_ec06076e-c8d0-4a94-8cfa-358f352bb1e8"; gene_name "chr11:LTR/ERV1:RLTR14-int:7"
chr12   RepeatMasker    exon    3052672 3059141 1.6     +       .       gene_id "RepMasker_8ea0db28-2cd0-4d7e-b151-b5b6480ee320"; gene_name "chr12:LTR/ERVK:IAPEz-int:1"
chr12   RepeatMasker    exon    3065049 3072337 27.3    -       .       gene_id "RepMasker_f3379b4a-cdfa-4b24-a4c9-5dacbdaf004c"; gene_name "chr12:LTR/ERV1:MMERGLN-int:2"
chr12   RepeatMasker    exon    3084084 3084938 5.8     -       .       gene_id "RepMasker_062a5737-d991-4220-9648-33c43ad467a2"; gene_name "chr12:LTR/ERVK:IAPEz-int:4"
```

### How tools deal with multimapping
Below list common feature count/qunatification tools emphasizing on the core algorithms they wrapped to deal with multimapping. [These common tools have similar quantification accuracy](https://bmcgenomics.biomedcentral.com/articles/10.1186/s12864-017-4002-1).

{{< table path="a.csv" header="true" caption="**Table 1. Multimapping options and algorithms of tools**" >}}

They by default may simply discard multimapping reads, which ignores most of hERVs. With carefully specified arguments, they either distribute read counts equally among certian features (fraction) or estimate using EM algorithm.

---
### Why EM algorithm for multimapping
"Fraction" algorithm overlooked how likely the reads can match with reference sequences. The biological prior information, if reasonably exploited, can significantly improve the multimapping quantification. Beatifully, multimapping by nature is similar to the problem of Gaussian Mixtures and EM has been proved to resolve Gaussian Mixture Model (GMM) well.

#### EM is proved to resolve GMM well
GMM describes a set of data points of which each comes from one of several Gaussian distributions, but we don't know which one. This missing information makes direct estimation (like via maximum likelihood) of the Gaussian parameters hard. EM solve this by iterating below E and M steps until converge:

- E-step (Expectation): Use initiallized parameter to compute the **posterior probability (sounds like bayesian, but it is not in standard EM)** that each data point was generated by each Gaussian (i.e., the responsibility γ for each component).

- M-step (Maximization): Update the parameters (π, μ, σ²) by performing **weighted MLE** using the expected values (responsibilities) computed in the E-step.

---
#### Multimapping is modeled as mixtures of distributions that account for mapping bias
[Salmon model ](https://www.nature.com/articles/nmeth.4197) describes a set of fragments (seqeuncing reads) of which each fragment *f* originates from one of the several transcripts *t* at a probability distribution {{< math >}} $P(t \mid f)$ {{< /math >}}, but we don't know which one. Thus EM goes:

##### *E-Step*

For each fragment {{< math >}} $f$ {{< /math >}} and transcript {{< math >}} $t$ {{< /math >}}:

{{< math >}}
$$
P(t \mid f) = \frac{\alpha_t \cdot P(f \mid t)}{\sum_{t'} \alpha_{t'} \cdot P(f \mid t')}
$$
{{< /math >}}

Where:

- {{< math >}} $\alpha_t$ {{< /math >}}: Current estimate of the relative abundance of transcript {{< math >}} $t$ {{< /math >}}.
- {{< math >}} $P(f \mid t)$ {{< /math >}}: Probability of observing fragment {{< math >}} $f$ {{< /math >}} given transcript {{< math >}} $t$ {{< /math >}}, incorporating bias corrections.

---

##### *Bias Corrections*

The fragment-transcript probability is modeled as:

{{< math >}}
$$
P(f \mid t) = P_{\text{len}}(l_f \mid t) \cdot P_{\text{pos}}(p_f \mid t) \cdot P_{\text{seq}}(s_f \mid t)
$$
{{< /math >}}

Where:

- {{< math >}} $l_f$ {{< /math >}}: Length of fragment {{< math >}} $f$ {{< /math >}}.
- {{< math >}} $p_f$ {{< /math >}}: Position of fragment {{< math >}} $f$ {{< /math >}} on transcript {{< math >}} $t$ {{< /math >}}.
- {{< math >}} $s_f$ {{< /math >}}: Sequence content of fragment {{< math >}} $f$ {{< /math >}}.
- {{< math >}} $P_{\text{len}}$ {{< /math >}}: Probability of fragment length, modeling the empirical fragment length distribution.
- {{< math >}} $P_{\text{pos}}$ {{< /math >}}: Positional bias model, capturing biases in fragment start locations.
- {{< math >}} $P_{\text{seq}}$ {{< /math >}}: Sequence-specific bias model, correcting for nucleotide composition effects.

Importantly, the bias models: {{< math >}} $P_{\text{len}}$ {{< /math >}}, {{< math >}} $P_{\text{pos}}$ {{< /math >}} and {{< math >}} $P_{\text{seq}}$ {{< /math >}} are not directly updated by the EM algorithm during transcript abundance estimation. Instead, these bias models are learned once, typically in a preprocessing or auxiliary phase before the EM iterations begin. The EM algorithm then uses the fixed values of these bias models to compute the probability {{< math >}} $P(f \mid t)$ {{< /math >}} in the E-step.

**The only parameters updated by EM are the relative abundances {{< math >}} $\alpha_t$ {{< /math >}} of the transcripts.**

---

##### *M-Step*

Update the abundance estimates:

{{< math >}}
$$
\alpha_t = \sum_f P(t \mid f)
$$
{{< /math >}}

Then normalize:

{{< math >}}
$$
\sum_t \alpha_t = 1
$$
{{< /math >}}

---

## Quantificaiton using Salmon

Now we see EM algorithm superior than fraction algorithms in quantifing multimappings since EM 
- is well suitable for estimating parameters in Gaussian Mixture Models, which is similar to the multimapping mathematical model
- takes multimapping mathematical models of mapping probability that accounts for the mapping bias

We implemented Salmon since it allows both EM and full bayesian EM. We will talk about [EM (with MLE), EM with MAP and bayesian EM](donesoon) as well as [Frequentist vs Bayesian](donesoon) in other blogs.

Please see source codes in [github](https://github.com/jliu678/herv_project_siyi/tree/main) that can be deploied in Cloud Cluster Computational platform and local desktops, and allows for: 
- batch computaion with tunable batch size
- automatic environment setting up
- file preprocessing
- QC
- index buidling, alignment and feature counting by either `subread` or `salmon`  

The github repo also contains [usage examples using `SLURM` or `LSF` job scheduler](https://github.com/jliu678/herv_project_siyi/tree/main/example_usage) and [handy utility tools](https://github.com/jliu678/herv_project_siyi/tree/main/utils).

---

### Excerpt scripts
Below excerpt scripts call each needed analysis step specified as `ANALYSIS_STEP` array, including file preprocessing, index buidling by `build_index()`, QC by `qc_all()`, alignment by `align_all()` if using `subread`, and feature counting by `count_all()`.
```bash
	for i in ${ANALYSIS_STEP[@]}; do
	  timed_print "$i-ing..."	
		case "$i" in 
			# for the first 2, do not run if they are in a child/forked process
			index) if [ $CHILD = false ]; then build_index; fi ;; 
			convert) if [ $CHILD = false ]; then get_pairs_all; fi ;;
			qc) qc_all ;;
			align) align_all ;;
			count) count_all ;;
		esac 
		timed_print "finished $i"
	done
```
Some of the core functions in the anylysis are below:
#### build index
```bash
salmon_build_index(){ # supports multiple transcripts
	# get dir name for the index, which is the transcript names seperated by a hyphen(-)
	# eg. for erv.fa and line.fa, you get: erv-line
	local index_name="salmon/$(IFS=-; echo "${TRANSCRIPTS[*]%.*}")" 
	if [ "$OVER_WRITE" = "true" ] || [ ! -d "${index_name}_index" ]; then 
		timed_print "building salmon index @: ${index_name}_index"

		# salmon requires a decoy, which in this case is the entire reference genome
		if [ "$OVER_WRITE" = "true" ] || [ ! -f "${index_name}_decoys.txt" ]; then
			timed_print "building decoys..."
			grep "^>" "${REF_GENOME}" | cut -d " " -f 1 > "${index_name}_decoys.txt"
			sed -i.bak -e 's/>//g' "${index_name}_decoys.txt"
		fi 

		# salmon requires all the transcripts being quantified to be concat into one file
		if [ "$OVER_WRITE" = "true" ] || [ ! -f "${index_name}_gentrome.fa" ]; then
			timed_print "building gentrome..."
			cat "${TRANSCRIPTS[@]}" "${REF_TRANSCRIPT}" "${REF_GENOME}" > "${index_name}.fa"
		fi

		# adding --gencode flag due to salmon doing some preprocessing for it
		# this flag does nothing for none gencode files (I think)
		salmon index -t ${index_name}.fa -d "${index_name}_decoys.txt" -i "${index_name}_index" --gencode
	fi
}

build_index() {
	case $ALIGN_METHOD in
		subread) subread_build_index ;;
		salmon) salmon_build_index ;;
	esac
}
```

#### QC
```bash
fastp_qc(){ 
	# results stored in /tmp/{name of the source}/qc/* and all end in *.qc.fq.gz
	if [ "$OVER_WRITE" = "true" ] || [ ! -f "tmp/${SOURCE}/qc/$1.qc.fq.gz" ]; then
		if [[ $# -eq 2 ]]; then # if pair ended
			./fastp -i "tmp/${SOURCE}/$1.fq.gz" -o "tmp/${SOURCE}/qc/$1.qc.fq.gz" \
				-I "tmp/${SOURCE}/$2.fq.gz" -O "tmp/${SOURCE}/qc/$2.qc.fq.gz" \
				-j "results/fastp/$1.json" -h "results/fastp/$1.html"
		elif [[ $# -eq 1 ]]; then # if single ended
			./fastp -i "tmp/${SOURCE}/$1.fq.gz" -o "tmp/${SOURCE}/qc/$1.qc.fq.gz" \
				-j "results/fastp/$1.json" -h "results/fastp/$1.html"
		fi	
	fi
}

qc_all(){ # only supports fastp as of now
	if [[ ! -d "tmp/${SOURCE}/qc" ]]; then 
		mkdir "tmp/${SOURCE}/qc"
	fi

	if [[ ! -d "results/${QC_METHOD}" ]]; then 
		mkdir "results/${QC_METHOD}"
	fi

	timed_print "qc-ing with $QC_METHOD"
	while IFS=, read -r r1 r2; do # loop over pair file, csv
		if [ "$OVER_WRITE" = "true" ] || [ ! -d "results/${ALIGN_METHOD}/$r1" ]; then 
			#if salmon quant already exists

			timed_print "qc-ing ${r1} and ${r2}"
			case $QC_METHOD in 
				fastp) fastp_qc $r1 $r2 ;;
			esac 
		fi
	done < $PAIR_FILE
}
```

#### quantification
```bash
salmon_quant() {
	# get the same index_name as specified in "salmon_build_index()"
	local index_name="salmon/$(IFS=-; echo "${TRANSCRIPTS[*]%.*}")"
	#salmon outputs into a dir, so we check for that instead
	#results stored in results/salmon/*
	if [ "$OVER_WRITE" = "true" ] || [ ! -d "results/salmon/$1" ]; then
		mkdir results/salmon/$1
		if [[ $# -eq 2 ]]; then # if pair ended
			salmon quant -i "${index_name}_index" -l A -1 "tmp/${SOURCE}/qc/$1.qc.fq.gz" -2 "tmp/${SOURCE}/qc/$2.qc.fq.gz" --validateMappings -o "results/salmon/$1"
			#sleep 1
		elif [[ $# -eq 1 ]]; then # if single ended
			salmon quant -i "${index_name}_index" -l A -r "tmp/${SOURCE}/qc/$1.qc.fq.gz" --validateMappings -o "results/salmon/$1"
			#sleep 1
		fi 
	fi
}

salmon_count() {
	while IFS=, read -r r1 r2; do
		salmon_quant $r1 $r2
	done < $PAIR_FILE
}

count_all(){
	if [ ! -d "results/${ALIGN_METHOD}" ]; then 
		mkdir "results/${ALIGN_METHOD}"
	fi
	
	case $ALIGN_METHOD in 
		subread) subread_count ;;
		salmon) salmon_count
	esac
}
```
