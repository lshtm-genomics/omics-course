# Transcriptomics Task

!!! Warning

    Be sure to use the `tb_transcriptomics_task` Codespace when initiating VM for this practical.

## Introduction

As we move from raw sequencing analysis into transcriptomic interpretation, it is important to apply your understanding of RNA-seq workflows to real Mycobacterium tuberculosis datasets. Gene expression analysis is a powerful approach for uncovering biological differences between bacterial lineages, particularly in relation to metabolism, stress response, and virulence.

In this practical, you will work with pre-processed gene count data derived from RNA sequencing experiments. These data represent transcriptional activity across multiple MTB isolates and lineages. Your objective is to identify patterns of differential gene expression that may explain phenotypic differences between strains, including potential adaptations linked to lineage evolution and survival strategies.

One sample has not yet been processed into counts, and you will generate this yourself using alignment and gene quantification tools. The remaining samples are provided as ready-to-use count files to ensure the analysis remains feasible within computational and storage constraints.

You will integrate all samples using a provided metadata file and perform differential expression analysis in R using DESeq2.

---

## Goal

By the end of this task you will:

1. Generate a gene count table for one sample (using HTSeq or featureCounts)
2. Combine it with provided count files
3. Perform differential expression analysis in R
4. Interpret lineage-specific expression patterns
5. Visualise results using heatmaps
---

## Input Data

You are provided with:

RNA-seq reads (one sample only)

- `*_1.fastq.gz`
- `*_2.fastq.gz`


This sample must be aligned to the reference genome and converted into a count table using HTSeq or featureCounts

You are also provided with the reference genome

- `tb.fasta`

and the gff file

- `tb.gff`

All other samples have been counted already due to space constraints

- `*counts.txt`

Finally, instead of manually inputting the lineage, you are provided with lineage information in:

- `tb_profiler.txt`

This needs to be read into R and used to construct the DESeq2 metadata table.


---

## Workflow Overview

You will complete the following stages:

1. Map one RNA-seq sample
2. Generate gene counts for that sample
3. Combine with provided count files
4. Build metadata-driven DESeq2 dataset
5. Perform differential expression analysis
6. Visualise results (heatmap + summary tables)

---

## Tips

???+ note "Tips for Mtb Transcriptomic Analysis"

    - Ensure sample names are consistent across FASTQ files, count files, and metadata before running DESeq2.

    - Always verify that the last column of featureCounts files is used as the expression count matrix.

    - Lineage is your primary experimental variable; treat it as a factor in DESeq2.

    - If samples do not align properly in R, check ordering before constructing the DESeq dataset.

    - Low-count genes add noise—filter them before running DESeq2.

    - Heatmaps work best with normalised counts (`counts(dds, normalized=TRUE)`).

    - PCA is a useful first step to confirm whether lineages separate biologically.

    - One mislabelled sample can break the entire analysis, so validate metadata carefully before proceeding.