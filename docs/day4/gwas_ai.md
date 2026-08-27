# TB Genomics: GWAS & Machine Learning


## What you're doing today

You have whole-genome sequencing data from 1,000 simulated *M. tuberculosis*
samples, each with a rifampicin resistance status and lineage assignment.
Over this session you'll:

1. Look at the population structure hiding in the genomes (which samples are
   genetically related to each other, and why that matters).
2. Run a genome-wide association study (GWAS) two different ways, and watch
   what happens to the results when you do — and don't — correct for that
   population structure.
3. Train machine learning models to predict resistance directly from
   genotypes, and check whether they "discover" the same variant the GWAS
   found.

By the end you'll have handled the two central problems in genomic
prediction: **confounding** (does an association reflect biology, or just
shared ancestry?) and **overfitting** (does a model's apparent skill
survive contact with new data?).

## Setup

```bash
conda activate gwas_ai
cd ~/GWAS_AI/tb
```

You should see `matched_geno_aligned.csv` (the genotype matrix — 1,000
samples x several thousand variants) and `matched_pheno_aligned.csv`
(sample IDs, lineage, sublineage, and drug resistance status) in this
folder, plus the Python scripts you'll run below.

> **Why "aligned"?** Genotype and phenotype files can easily drift out of
> sync — different sample order, samples present in one file but not the
> other. That's a silent failure mode: it won't crash your code, it'll
> just quietly pair the wrong genotype with the wrong phenotype and give
> you meaningless results. Normally the first step of any pipeline like
> this is to explicitly check and fix that (a script called
> `align_geno_pheno.py` does exactly this) — that step has already been
> done for you here, which is why your files are named `..._aligned.csv`.
> Worth remembering for any dataset you work with outside this course:
> never assume two files are in matching sample order just because they
> came from the same place.

Let's take a look at what these files look like. 

```bash
head -n 5 "matched_geno_aligned.csv"
```

> You should see a complex read out. If you look closely you should see
> many 0 and 1s that are separated by a ','. This is where the term 'csv'
> file comes from - a comma separated value. These are the genotypes. 
> They are encoded in a binary format. You may wonder how we got from a 
> VCF to these genotypes. Well we designated 1 if 80% of the reads supported
> the alternative allele (variant allele frequency). We designated 0 if 80%
> of the reads supported a reference allele. If in between we assign as missing.
> Sometimes people may assign heterozygous calls as 0.5, we have kept them missing
> as M.tb is a haploid organism and they are rare. Sometimes we also have low depth
> or missing variants - we can impute these using a variety of statistical methods
> (not covered here), but for now we have just used the reference genotype '0'.

> VCFs can be converted to a table using:

```
bcftools query -f '%CHROM/t%POS/t%REF/t%ALT[/t%GT]/n' your.vcf.gz > genotype_table.tsv
bcftools query -l your.vcf.gz  > sample_names.txt

```

> Simply it prints out the chromosome, position, reference allele, alternative allele
> and the genotypes (0 or 1) per sample. We also outputted the sample names using 'query'.
> We have already run this for you and made a clean csv file for you to work with. 
> This important before doing any modelling.

Now let's look at the phenotypes. 

```bash
head -n 5 "matched_pheno_aligned.csv"

```

> Here we have the same csv format, except for the phenotypic data. Here we have
> each sample ID, whether the sample is resistant (1) or susceptible (0) to
> isoniazid or rifampicin and then also the lineage of each TB sample.

---

## Part 1 — Population structure 

*M. tuberculosis* reproduces clonally — essentially no recombination — so
closely related strains stay genetically similar across their *entire*
genome, not just at a few loci. This is very different from something
like human GWAS, where recombination each generation breaks correlations
down over distance.

Let's have a look at population structure using our data.

```bash
python3 population_structure.py \
    --geno_matrix matched_geno_aligned.csv \
    --prefix tb \
    --n_pcs 5 \
    --colour_by matched_pheno_aligned.csv \
    --colour_col sublineage --top_n_categories 20
```


This does two things: builds a pairwise genetic distance matrix between all
1,000 samples, and runs PCA (Principal Component Analysis) on the genotype
matrix to summarise that structure in a handful of numbers per sample (the
"PCs"). Open **`tb_pca_scatter.png`**:

![PCA of genotype matrix](tb_pca_scatter.png)

Each point is one sample, coloured by sublineage. Samples cluster tightly
by lineage — this is the population structure. Hold onto this plot; it's
the reason the next two sections give different-looking answers to the
same question.

> **What's a PC, in one sentence?** A principal component is a new,
> made-up axis that captures as much of the variation between samples as
> possible. PC1 alone here explains ~31% of all the genetic variation
> across the genome — that's an enormous amount for a single number to
> capture, and it's telling you the dominant signal in your genotypes is
> "which lineage is this sample from," not anything to do with individual
> variants.

> You will notice throughout this module we will be using lots of python 
> scripts. Python is a popular programming language for machine learning. 
> Let's open up the population_structure.py scripts to understand how these 
> files are structured. 

!!! question "Exercise"
=== "Question 1"

    Open up population_structure.py - what do you think 'import', 'def', 'def main(args)' and 'parser' mean?
    Feel free to run an internet search to help you.

=== "Answer 1"

    In Python, these terms import external libraries (import), define reusable functions (def), set up the script's
    main execution function (def main(args)), and handle command-line inputs passed by the user (parser).
    

This is how many scripts and bioinformatic tools work, but don't worry too much for now- you can just apply the scripts.

---

## Part 2 — A naive GWAS

Even though we have observed strong population structure. We should 
first think about what happens if we directly test every variant in the genome for 
association with rifampicin resistance. Here we use ordinary logistic regression.
This performs one test per variant, with no correction for the fact that samples 
aren't independent of each other.

```bash
python3 logistic_gwas.py --geno_matrix matched_geno_aligned.csv --metadata matched_pheno_aligned.csv \
    --drug rifampicin --prefix logistic --covariates tb_pcs.csv
```

This runs the test twice — once with no correction at all, and once with
the 5 PCs from Part 1 added as covariates — so you can compare directly.

```bash
python3 compare_covariate_correction.py \
    --assoc_with_covars logistic_with_covars.assoc.txt \
    --assoc_without_covars logistic_no_covars.assoc.txt \
    --prefix logistic_logistic --p_col p_lrt
```

Open **`logistic_logistic_qq_comparison.png`**:

![QQ plot comparing with and without PC correction](logistic_logistic_qq_comparison.png)

This is a **QQ plot**. If none of your variants were truly associated with
resistance, the p-values across thousands of tests should be roughly
uniformly distributed — plotting observed vs. expected significance would
trace the diagonal dashed line. Points sitting *above* the diagonal mean
you're seeing more significant results than chance alone would produce.

<details>
<summary><b>❓ Discussion: what does the reduced inflation (λ) actually mean?</b></summary>

The number printed above each line — **λ (lambda), the genomic inflation
factor** — summarises how far off the diagonal the bulk of your points
sit. λ = 1.0 means your test statistics behave exactly as expected under
the null; λ > 1 means they're systematically inflated.

**Without PCs, λ = 2.86.** That's a lot of inflation — it means that
across the genome, far more variants look "significant" than should, by
chance, given how many tests you ran. This is a direct consequence of the
population structure you saw in Part 1: if resistance rates differ even
slightly between lineages (for reasons that have nothing to do with any
specific variant — different treatment history, different geography,
whatever), then *any* variant that happens to be more common in one
lineage than another will show a spurious statistical association with
resistance, purely by riding along with lineage. With thousands of
lineage-tagging variants in your genome, that produces genome-wide,
systematic inflation — not just one or two false positives, but a
pervasive bias.

**With PCs, λ drops to 1.21.** Adding the 5 PCs as covariates lets the
model "explain away" the resistance-rate differences that track lineage,
so a variant now has to explain something *beyond* lineage to look
significant. λ isn't all the way down to 1.0, which makes sense — 5 PCs
capture broad structure but not all fine-scale relatedness, and this
naive logistic regression has no other correction mechanism at all (more
on that in Part 4).

**Why this matters practically:** if you'd trusted the uncorrected p-values
and picked out your "top hits," most of them would likely be
lineage-tagging noise, not real drug-resistance biology — you'd be chasing
false leads.
</details>

---

## Part 3 — Which variant actually stands out?

```bash
python3 manhattan_plot.py \
    --assoc logistic_with_covars.assoc.txt \
    --p_col p_lrt \
    --prefix lr_gwas_tb
```

Open **`lr_gwas_tb_manhattan.png`**:

![Manhattan plot, logistic regression with PC correction](lr_gwas_tb_manhattan.png)

This is a **Manhattan plot** — every tested variant plotted by genomic
position (x-axis) against -log10(p-value) (y-axis), so more significant
hits sit higher up. The dashed lines are significance thresholds corrected
for the fact that you ran thousands of simultaneous tests (see the box
below).

One variant towers over everything else, dozens of orders of magnitude
past the significance threshold.

<details>
<summary><b>❓ Discussion: which marker is this, and what's known about it?</b></summary>

The labelled point is **`rpoB_761155_C_T`** — a C→T mutation at position
761,155 in the standard H37Rv reference genome, inside the *rpoB* gene.

**Your task:** search for "rpoB 761155" or "rpoB S450L rifampicin
resistance" and see what comes up.

**What you should find:** this is one of the best-characterised
drug-resistance mutations in all of clinical microbiology. *rpoB* encodes
the beta subunit of RNA polymerase, the enzyme rifampicin targets. Codon
450 sits inside the **Rifampicin Resistance-Determining Region (RRDR)**, an
81-base-pair stretch (codons 426–452) where the overwhelming majority of
resistance-conferring mutations cluster. The specific substitution here —
**S450L** (serine to leucine at codon 450, sometimes written S531L using
older *E. coli*-based numbering) — is consistently reported as the single
most common rifampicin-resistance mutation worldwide, frequently cited as
accounting for something in the range of 40-50% of resistant isolates on
its own. Structurally, it disrupts the rifampicin-binding pocket on RNA
polymerase, meaning the drug can no longer block transcription the way it
does in a susceptible strain.

The fact that a purely statistical GWAS — with no prior biological
knowledge fed in — lands squarely on the exact variant decades of
molecular microbiology has independently identified as the primary cause
of rifampicin resistance is itself worth sitting with. That's the whole
premise of GWAS: let the data point you to the biology, rather than the
other way around.
</details>

---

## Part 4 — A "proper" GWAS with GEMMA (~15 min)

`logistic_gwas.py` treated every sample as statistically independent —
which you now know is false. GEMMA fits a **linear mixed model (LMM)**
instead, which explicitly models relatedness between samples via a
**kinship matrix**. Let's run it two ways — without and with the same PCs
— using a relaxed minor allele frequency threshold (`0.001`) appropriate
for this dataset's rarer variants.

```bash
# WITHOUT PCs
python3 convert_to_gemma.py \
    --geno_matrix matched_geno_aligned.csv --metadata matched_pheno_aligned.csv \
    --drug rifampicin --prefix gwas_no_pcs
bash ./run_gemma.sh gwas_no_pcs 0.001 0.05

# WITH PCs
python3 convert_to_gemma.py \
    --geno_matrix matched_geno_aligned.csv --metadata matched_pheno_aligned.csv \
    --drug rifampicin --prefix gwas_with_pcs --covariates tb_pcs.csv
bash ./run_gemma.sh gwas_with_pcs 0.001 0.05

python3 manhattan_plot.py \
    --assoc gemma_output/gwas_no_pcs_assoc.assoc.txt \
    --p_col p_wald \
    --prefix gwas_no_pcs_maf001

python3 manhattan_plot.py \
    --assoc gemma_output/gwas_with_pcs_assoc.assoc.txt \
    --p_col p_wald \
    --prefix gwas_with_pcs_maf001
```

**Without PCs:**
![GEMMA Manhattan plot, no PCs](gwas_no_pcs_maf001_manhattan.png)

**With PCs:**
![GEMMA Manhattan plot, with PCs](gwas_with_pcs_maf001_manhattan.png)

<details>
<summary><b>❓ Discussion: why do these two plots look so similar, when Part 2/3's did not?</b></summary>

In Part 2, adding PCs visibly reduced genome-wide inflation (λ 2.86 → 1.21).
Here, running GEMMA with vs. without the *same* PCs barely changes
anything.

**Your task:** look up "linear mixed model kinship matrix population
structure GWAS" and think about what GEMMA is doing differently from plain
logistic regression.

**What you should find:** GEMMA doesn't just optionally accept covariates
— every GEMMA run first builds a **kinship (relatedness) matrix** from the
genotypes themselves (that's what `run_gemma.sh` does before the
association step), and every association test is fit as a mixed model
that includes this kinship matrix as a **random effect**, whether or not
you pass it any PCs at all. In effect, GEMMA is already asking "does this
variant's association with resistance hold up *after* accounting for how
related every pair of samples is to each other?" — which captures
population structure directly and comprehensively, rather than through a
handful of summary axes.

This is the real reason your naive-logistic-vs-GEMMA comparison matters
more than your GEMMA-with-vs-without-PCs comparison: the kinship
correction is doing most of the structure-correction work by default,
before PCs ever enter the picture. PCs on top of kinship provide some
additional refinement, but the big jump in rigor happened when you moved
from "independent samples" logistic regression to a kinship-aware mixed
model in the first place.
</details>

---

## Part 5 — Key concepts before the machine learning section

A few terms you'll see in the next part, in plain language:

- **Train/test split**: you don't evaluate a model on the same data it
  learned from — that's like grading a student's exam using the answer
  key they wrote it with. You hold out a chunk of samples (the *test set*)
  the model never sees during training, and check performance only on
  those.
- **Cross-validation (CV)**: instead of relying on a single train/test
  split (which could be lucky or unlucky by chance), you split the
  *training* data into several folds, train on some, validate on the rest,
  and repeat — rotating which fold is held out each time. This gives you a
  mean performance estimate *and* a sense of how much it varies, which a
  single split can't.
- **Overfitting**: when a model learns the training data too well —
  including its noise and quirks — rather than the general pattern. The
  signature is a big gap between training performance (often near-perfect)
  and test/CV performance (much lower).
- **Sensitivity (recall)**: of all the samples that are *actually*
  resistant, what fraction did the model correctly flag? Low sensitivity
  means you're missing real resistant cases.
- **Specificity**: of all the samples that are *actually* susceptible,
  what fraction did the model correctly flag? Low specificity means false
  alarms.
- **AUC (Area Under the ROC Curve)**: a single number (0.5 = random
  guessing, 1.0 = perfect) summarising how well the model ranks resistant
  samples above susceptible ones, across every possible decision
  threshold at once.
- **Feature importance**: for tree-based models, a score per variant
  reflecting how much that variant contributed to the model's decisions.
  High importance is not proof of biological causation on its own — but
  when it lines up with an independent GWAS hit (as you're about to see),
  that's a strong corroborating signal.

---

## Part 6 — Splitting the data and your first model (~10 min)

```bash
python split_data.py --geno_matrix matched_geno_aligned.csv --metadata matched_pheno_aligned.csv \
    --test_size 0.2 --random_state 123 --stratify_col rifampicin --prefix ml
```

`--stratify_col rifampicin` keeps the ratio of resistant/susceptible
samples the same in both the train and test sets — important because
resistance is relatively rare, and a plain random split could by chance
dump most resistant samples into one side.

```bash
python3 make_predictions.py --geno_train ml_train_geno.csv --geno_test ml_test_geno.csv \
    --pheno_train ml_train_pheno.csv --pheno_test ml_test_pheno.csv \
    --model Quick --drug rifampicin --prefix ml --cv_folds 5
```

> `--model Quick` is a single, **unconstrained decision tree** — no depth
> limit, no pruning. It's called "Quick" because it trains almost
> instantly (unlike the ensemble/grid-searched models later), which also
> makes it the easiest model to overfit and the clearest one to diagnose.

Open `ml_rifampicin_Quick_results.csv`:

| Train AUC | Test AUC | Sensitivity | Specificity | F1 | CV Val AUC (mean ± std) |
|---|---|---|---|---|---|
| 1.000 | 0.929 | 0.857 | 1.000 | 0.923 | 0.957 ± 0.029 |

And the top features driving the model, from
`ml_rifampicin_Quick_feature_importance.csv`:

| Rank | Feature | Importance |
|---|---|---|
| 1 | **rpoB_761155_C_T** | **0.896** |
| 2 | Chromosome_4273354_A_G | 0.062 |
| 3 | Chromosome_4355823_A_G | 0.034 |
| 4 | Chromosome_4294111_A_G | 0.005 |
| 5 | Chromosome_4382121_A_G | 0.003 |
| 6-10 | (5 more variants) | 0.000 |

<details>
<summary><b>❓ Discussion: are these results good? Does the model agree with the GWAS?</b></summary>

**Is Sens/Spec good?** Specificity is perfect (1.000) — every susceptible
test sample was correctly identified. Sensitivity (0.857) is lower, but
look at the test set size: at this resistance prevalence, the test set
likely contains well under 20 actually-resistant samples, so a single
misclassified case can shift sensitivity by several percentage points.
Compare the single-split Test AUC (0.929) against the 5-fold CV AUC
(0.957 ± 0.029, averaged over many more samples across folds) — CV is the
more stable, trustworthy estimate of how this model actually performs.

**Train AUC = 1.000** is a red flag on its own — a perfect score on
training data usually means the model has memorised it rather than learned
a generalisable rule. That the CV and test AUC are still both high (~0.93-
0.96) tells you the model *also* generalises well here — but don't assume
a perfect training score is automatically fine; check the next section for
a case where it isn't as reassuring.

**Does feature importance match the GWAS?** Yes, directly:
`rpoB_761155_C_T` — the exact variant that stood out in Part 3's Manhattan
plot — is also by far the most important feature here (0.896, more than
10x the next-highest variant). Two completely different statistical
approaches (a genome-wide association test, and a decision tree learning
to classify samples) independently converged on the same answer. That
convergence is much stronger evidence than either method alone.
</details>

Now look at the overfitting diagnostics:

```bash
python3 overfitting_diagnostics.py --geno_train ml_train_geno.csv \
    --pheno_train ml_train_pheno.csv \
    --drug rifampicin --prefix ml --cv_folds 5
```

Open **`ml_validation_curve.png`**:

![Validation curve](ml_validation_curve.png)

This is a **validation curve** — it re-trains a decision tree at every
possible depth limit (`max_depth`) and plots training vs. cross-validation
AUC at each one, all using the training data only. Notice:

- **Train AUC (blue)** shoots to ~1.0 almost immediately and stays there
  — the tree can always fit the training data perfectly given enough
  depth. This line alone tells you nothing about generalisation.
- **CV AUC (red)** peaks around `max_depth=3`, then *drops* and plateaus
  from depth 4 onward — more complexity past that point doesn't help, and
  slightly hurts. The unconstrained tree we ran above (far right,
  "None") is sitting past that peak.
- The pink shaded band is the spread across CV folds — it's wide, another
  reminder that any single performance number here comes with real
  uncertainty at this sample size.

The takeaway: train AUC is *never* a reliable stopping signal on its own;
you need the CV (or held-out test) curve to see where the model actually
stops improving.

---

## Part 7 — Comparing model types (~10 min)

`Quick` was a single decision tree. Now try a **regularised logistic
regression** — a fundamentally different kind of model that fits a
weighted linear combination of all variants, with a penalty that
discourages the model from relying too heavily on any one feature.

```bash
python3 make_predictions.py --geno_train ml_train_geno.csv --geno_test ml_test_geno.csv \
    --pheno_train ml_train_pheno.csv --pheno_test ml_test_pheno.csv \
    --model LRL2 --drug rifampicin --prefix ml --cv_folds 5
```

> **LRL1 vs LRL2** — both are logistic regression with a regularisation
> penalty on the coefficients, but the penalty is different:
> - **LRL1 (Lasso, L1 penalty)** can shrink irrelevant coefficients all
>   the way to exactly zero — effectively performing automatic feature
>   selection. With thousands of variants and only one or two truly
>   causal, this is often the more appropriate choice.
> - **LRL2 (Ridge, L2 penalty)** shrinks coefficients smoothly toward zero
>   but rarely to *exactly* zero — it spreads weight thinly across many
>   correlated variants rather than picking one. With thousands of raw SNP
>   features and relatively few training samples, this can make LRL2
>   noticeably more prone to overfitting on noise than a tree-based model
>   or LRL1 here — worth checking directly by comparing its results table
>   against `Quick`'s.

Open the new `ml_rifampicin_LRL2_results.csv` and compare it side by side
with the `Quick` table above. Does LRL2 do better, worse, or about the
same? Does its overfitting gap (Train AUC − Test AUC) look larger or
smaller than the decision tree's? That comparison is the point of running
more than one model type — no single algorithm is automatically "best,"
and the right diagnostic (train/test gap, CV spread, feature importance)
tells you *why* one wins over another for this particular dataset.

---

## Wrap-up

By now you've seen the same underlying signal — the `rpoB` S450L
mutation — surface independently through a naive GWAS, a proper
mixed-model GWAS, and a machine learning model's feature importance
ranking, while also seeing exactly how population structure can generate
convincing-looking false positives if you don't correct for it, and how a
model can look perfect on training data while telling you very little
about how it'll perform on a new sample. Those two lessons — check for
confounding, check for overfitting — generalise well beyond this specific
dataset.
