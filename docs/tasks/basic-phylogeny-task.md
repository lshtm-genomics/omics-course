# Phylo Task

## Introduction

Welcome to the **Phylo analysis practical**, where you will reconstruct evolutionary relationships between Mycobacterium tuberculosis isolates starting from raw sequencing data.

This workflow mirrors real outbreak genomics pipelines used to:
- track transmission
- identify clusters of infection
- understand lineage structure
- compare genomic relationships between isolates

You will move from raw reads to a final phylogenetic tree and visualise it using iTOL.
---

## Goal

By the end of this task you will produce:
- consensus genomes for each isolate
- a core-genome alignment
- a phylogenetic tree
- an annotated iTOL visualisation

---

## Input Data

You are given paired-end sequencing files :

- `*_1.fastq.gz`
- `*_2.fastq.gz`

Each pair represents a TB isolate.

You are also given a reference genome:

- *Mycobacterium tuberculosis reference (tb.fasta)*

And given an annotation file to use in Itol at the end 

- `itol.txt`

These will need to be downloaded using 

- `wget -O - http://genomics.lshtm.ac.uk/phylo_task.tar.gz | tar -xvz`
---

## Workflow Overview

You will complete the following stages:

1. Read alignment
2. Variant detection
3. Consensus genome generation
4. Multi-sample genome collection
5. Phylogenetic reconstruction
6. Tree visualisation and annotation

---

## Stage 1 — Read Alignment

Each sample must be aligned to the reference genome.

### What you should think about:
- ensuring correct pairing of reads
- generating a sorted alignment file per sample
- indexing the alignment for downstream analysis

### Output expected:
- one BAM file per sample

---

## Stage 2 — Variant Detection

From each alignment, identify differences relative to the reference genome.

### What you should think about:
- detecting SNPs and small variants
- filtering low-confidence calls
- producing a variant file per sample

### Output expected:
- VCF file per sample

---

## Stage 3 — Consensus Genome Construction

Use variant information to reconstruct each isolate’s genome.

### What you should think about:
- Use bcftools on the vcfs to create a consensus with `bcftools consensus`
- generating a full-length FASTA per sample
- ensuring consistent naming across samples

### Output expected:
- one consensus FASTA per isolate

---


## Stage 4 — Combine Genomes

Merge all consensus genomes into a single dataset. Here we will use Parsnp, you can read the manual [here](https://harvest.readthedocs.io/en/latest/content/parsnp.html)


### What you should think about:
- move the generated fastas into one dir
- running parsnp on that one directory

### Output expected:
- multi-FASTA file containing all isolates

---


## Stage 5 — iTOL Visualisation

Upload your tree to **iTOL (Interactive Tree of Life)**. We have provided an annotation file and a script to generate the annotation file if interested.

### What to look for:
- clustering patterns
- lineage structure
- potential transmission clusters
- outliers or unusual branches

---

## Key Hints

- Keep sample names consistent across all files
- Low-quality samples may break downstream clustering
- Might need to remove some samples
- The tree reflects **genetic distance, not just lineage**
- Lineage labels are metadata overlays, not the result itself
- Parsnp builds a core-genome alignment, not just SNPs
- Interpretation is more important than execution

---

