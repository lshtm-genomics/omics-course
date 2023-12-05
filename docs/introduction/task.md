# Task

## Introduction

As we transition from our foundational lectures into practical application, it is crucial to evaluate your understanding and proficiency in the core bioinformatics processes of mapping, variant calling, and assembly. These tools are the cornerstone of genomic analysis and have wide-reaching implications in the field of infectious diseases, including tuberculosis. This test will simulate a real-world scenario where you will analyze raw sequencing data to extract meaningful insights that can inform treatment decisions and deepen our understanding of pathogen genomics.

So your task is to utilize your bioinformatics expertise to uncover the genetic factors that may explain the variability in clinical outcomes among these patients. By analyzing the sequencing data, you will identify crucial genetic variations—specifically, single nucleotide polymorphisms (SNPs) in drug-resistance genes and lineage-defining deletions. Your findings will not only classify the TB strains infecting each patient but also predict their resistance to commonly used medications, ultimately guiding more effective treatment strategies. This exercise will test your ability to apply genome mapping and variant analysis to real-world infectious disease challenges.

### Task

Six People have come into the local hospital and presented with symptoms of TB. Samples were taken, and were confirmed to be TB infection so antibiotics were given, however some patients are still not responding toward the antibiotics and are struggling to fight off infection. 

TB can also have certain deletions on the genome that can identify the lineage of the strain (check out TB lineages and locations across the world) such as an entire deletion on the PPE50 gene in lineage 1 strains and Rv0072 for lineage 2 strains. This is not the only identifies but it is what we will use for this test.

You will find the patients data within the following directory:

```
cd ~/data/tb/task

```

Your tasks are below:

1. Perform suitable research, what types of drug resistance are there with TB and what causes them. 
2. Decide whether to map with or without a reference. If you choose mapping you will need to find a suitable reference genome.
3. Call the variants and find the SNPs that are within the genes you found earlier that are known to cause drug resistance.
4. Identify the individuals that struggled to fight infection and what antibiotics were they using at the hospital.
5. Find the strain that comes with each patient by finding these deletions.