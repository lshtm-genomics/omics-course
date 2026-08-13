# Methylation Task: Lineage Assignment

!!! Warning
    Use the **methylation** Codespace for this task. It's the same one from the methylation practicala as it already contains the data you need.

## Introduction

In this task you'll combine two datasets, the Nanopore barcodes and the PacBio strains, into a single phylogeny, and use it to work out the **lineages of the PacBio samples**.

The Nanopore samples already have lineage assignments (from `tb-profiler`). The PacBio samples don't. By placing both sets of sequences on one tree, you can infer each PacBio strain's lineage from the Nanopore samples it clusters with.

## Objective

**Assign a lineage to each PacBio sample** by building a combined phylogeny and reading lineages across from the Nanopore samples.

## Starting point

```bash
cd data/methylation/tree/
```

Have a look at what's already there before you start (`ls`).

## Steps

You'll need to work out the exact commands yourself, that's the task, but here's the route:

1. **Combine the sequences.** You have separate FASTA files for the Nanopore and PacBio samples. Merge them into one multi-FASTA.
2. **Align.** The combined sequences need to be aligned together with **mafft** before you can build a tree.
3. **Build the tree.** Infer a phylogeny from the alignment (the same tree-building tool you used in the phylogenetics practical works here).
4. **Annotate in iTOL.** Upload the tree, then bring in the annotations:
    - the **lineage** of each sample, and
    - which platform (**Nanopore** or **PacBio**) each sample came from.
   You have an annotation file for each dataset, you'll need to **combine the two annotation files into one** so every tip on the tree is labelled.
5. **Read off the answer.** Find where each PacBio sample sits and assign it the lineage of the Nanopore clade it falls within.

## Hints

- The tree-building and mafft steps mirror what you did in the phylogenetics practical, reuse those commands.
- For the annotations: both files are the same format (tip label + value). Combining them is just a matter of getting all tips into one file with consistent columns, watch that the tip labels match the names in your tree exactly, or iTOL won't attach them.
- A PacBio sample sitting *inside* a Nanopore lineage clade takes that lineage. One sitting on its own long branch, between clades, is ambiguous, note it rather than forcing an assignment.