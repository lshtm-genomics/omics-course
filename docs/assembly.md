# Assembly Practical

## Introduction

The main challenge of sequencing a genome is determining how to arrange reads into chromosomes. We have already explored mapping, where the sequence data is aligned to a known reference genome. A complementary technique, where no reference is used or available, is called de novo assembly.

De novo assembly is the process of reconstruction of the sample genome sequence without comparison to other genomes. It follows a bottom-up strategy by which reads are overlapped and grouped into contigs. Contigs are joined into scaffolds covering, ideally, the whole of each chromosome in the organism. However de novo assembly from next generation sequence (NGS) data faces several challenges. Read lengths are short and therefore detectable overlap between reads is lower and repeat regions harder to resolve. Longer read lengths will overcome these limitations but this is technology limited. These issues can also be overcome by increasing the coverage depth, i.e. the number of reads over each part of the genome. The higher the coverage then the greater the chance of observing overlaps among reads to create larger contigs and being able to span short repeat regions. In addition, the modern hardware commonly output paired-end reads. These are two reads that are separated by a gap of known size. They enable the resolution of repeats greater than the read length by employing their expected separation and orientation as location constraints. Sequencing errors add difficulty since algorithms must allow certain mismatches when overlapping reads and joining regions, possibly leading to discarding true overlaps and false positives (Miller, Koren, & Sutton, 2010)

Despite all mentioned limitations, the high coverage currently achieved, growing read lengths (e.g. 100 read depth, >150bp – Illumina HiSeq2500) and paired-end information make it feasible to obtain assemblies from small genomes (e.g. bacterial) with relatively low fragmentation. Current assemblers employ graph theory to represent sequences and their overlaps as a set of nodes and edges, being classified into three main groups: greedy, Overlap/Layout/Consensus (OLC) and de Bruijn graph assemblers


!!! question
    This is a question?

!!! info
    this is osme info