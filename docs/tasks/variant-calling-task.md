# Variant Discovery Task

!!! Warning

    Be sure to use the `task1` Codespace when initiating VM for this practical.

## Introduction

As we transition from our foundational lectures into practical application, it is crucial to evaluate your understanding and proficiency in the core bioinformatics processes of mapping and variant calling. These tools are the cornerstone of genomic analysis and have wide-reaching implications in the field of infectious diseases, including tuberculosis. This test will simulate a real-world scenario where you will analyse raw sequencing data to extract meaningful insights that can inform treatment decisions and deepen our understanding of pathogen genomics.

So your task is to utilise your bioinformatics expertise to uncover the genetic factors that may explain the variability in clinical outcomes among these patients. By analysing the sequencing data, you will identify crucial genetic variations—specifically, single nucleotide polymorphisms (SNPs) in drug-resistance genes and lineage-defining deletions. Your findings will not only classify the TB strains infecting each patient but also predict their resistance to commonly used medications, ultimately guiding more effective treatment strategies. This exercise will test your ability to apply genome mapping and variant analysis to real-world infectious disease challenges.

These genomes are larger than what you have been working on as we have only given you a snippet of the genome. Now you will work on the whole genome so it will take a bit of time!

### Task

Six People have come into the local hospital and presented with symptoms of TB. Samples were taken, and were confirmed to be TB infection so antibiotics were given, however some patients are still not responding toward the antibiotics and are struggling to fight off infection. 

TB can also have certain deletions on the genome that can identify the lineage of the strain (check out TB lineages and locations across the world) such as an entire deletion on the PPE50 gene in lineage 1 strains and Rv0072 for lineage 2 strains. This is not the only identifications but it is what we will use for this test.

We have provided some tools to run the task, but you may have to install other tools using conda `conda install **your packages**`

You will find your information in the `task codespace` and the tb directory.


```
conda activate task
cd ~/tb

```

Your tasks are below:

1. Read the sample list, and download the data from NCBI or ENA.
2. Map to the reference.
3. Call the variants to find the SNPs.
4. Use the drug resistance information in the csv file to identify which drug has been used.
5. Find the strain (lineage type) that comes with each patient by finding deletions in the txt file provided.

### Tips

???+ note "Tips for Mtb Genomic Analysis"

    - You can search the sample code on ENA or NCBI and use tools like wget to download your sequences.

    - You can google other variant callers but the ones we have used are acceptable

    - Check the drug resistance list to find locations however you wish.

    - You can use structural variant caller to identify SV's or you can look at the positions in IGV to find the lineage information.


