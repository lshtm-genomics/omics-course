# Machine Learning for Antimicrobial Resistance Prediction
### A Practical Workshop - Please use the ml-amr-workshop codespace
*Mycobacterium tuberculosis · Genotype-to-Phenotype Prediction*

---

## Contents

1. [Overview and Learning Objectives](#1-overview-and-learning-objectives)
2. [Scientific Background](#2-scientific-background)
3. [Software Used](#3-software-used)
4. [Running the Pipeline](#4-running-the-pipeline)
5. [Hyperparameter Experiments](#5-hyperparameter-experiments)
6. [Interpreting Your Results](#6-interpreting-your-results)
7. [Phenotype File Format](#7-phenotype-file-format)
8. [Discussion Questions](#8-discussion-questions)
9. [Quick Reference](#9-quick-reference)

---

## 1. Overview and Learning Objectives

This workshop teaches you to build, evaluate, and interpret machine learning classifiers that predict rifampicin resistance in *Mycobacterium tuberculosis* directly from whole-genome sequencing (WGS) data. You will work with a real genotype-phenotype dataset, process it end-to-end, and critically assess how modelling choices affect predictive performance.

**By the end of this workshop you will be able to:**

- Encode raw VCF genotype data into a binary feature matrix suitable for ML
- Split data correctly to avoid data leakage
- Train logistic regression and decision tree classifiers with different hyperparameters
- Interpret sensitivity, specificity, F1, and AUC in a clinical AMR context
- Explain the effect of regularisation strength and tree depth on model behaviour
- Extract and interpret feature importances as biological hypotheses
- Compare model runs visually and draw evidence-based conclusions

---

## 2. Scientific Background

### 2.1 Why predict AMR from genotype?

Culture-based drug susceptibility testing (DST) for *M. tuberculosis* takes 2–6 weeks using solid media, or 1–2 weeks using liquid culture systems such as MGIT. WGS can be completed in 24–48 hours and, when paired with a trained classifier, can predict resistance before culture results are available. This has direct clinical impact: faster appropriate therapy, reduced transmission, and better stewardship.

Genotypic prediction works because resistance in *M. tuberculosis* is almost entirely driven by chromosomal mutations — there is no horizontal gene transfer of resistance determinants. ML approaches are valuable because they can learn from patterns across thousands of variants simultaneously, including variants of uncertain significance not captured by curated catalogues such as WHO-TBProfiler.

### 2.2 Rifampicin resistance in *M. tuberculosis*

Rifampicin (RIF) is one of the most important first-line anti-TB drugs and a cornerstone of standard short-course therapy. It works by binding to the β subunit of RNA polymerase, blocking transcription. Resistance arises almost exclusively through mutations in an 81 bp region of *rpoB* (Rv0667) known as the **rifampicin resistance-determining region (RRDR)**, spanning codons 426–452 (H37Rv numbering).

Key *rpoB* codons and their clinical significance:

| Codon | Common mutation | Resistance level |
|---|---|---|
| S450 | S450L | High-level — most prevalent globally |
| H445 | H445Y / H445D / H445R | High-level |
| D435 | D435V | High-level |
| L452 | L452P | Variable |
| I491F | I491F | Low-level — may be missed by some assays |

Rifampicin resistance is a strong proxy for multidrug-resistant TB (MDR-TB), since most MDR strains are also isoniazid-resistant. Detecting it rapidly is therefore of major public health importance.

> **💡 Why is *rpoB* such a good target for ML?**
> Because resistance is so strongly concentrated in the RRDR, a well-trained model should assign high importance to a small number of *rpoB* positions. If your top features are dominated by RRDR variants, this is a good sign the model has learned genuine biology rather than noise.

### 2.3 The feature matrix

Each sample (row) is a sequenced *M. tuberculosis* isolate. Each feature (column) is a SNP, identified by `CHROM_POS_ALT`. Genotypes are encoded numerically:

```
0/0  →  0.0   (homozygous reference — no mutation)
0/1  →  0.5   (heterozygous — mixed infection or within-host diversity)
1/1  →  1.0   (homozygous alternate — mutation present)
./.  →  NaN   (missing data — imputed during training)
```

*M. tuberculosis* is haploid, so heterozygous calls (0.5) typically indicate mixed infection rather than true heterozygosity. A 0.5 call at an *rpoB* RRDR position may indicate a mixed susceptible/resistant infection, which has implications for both phenotyping and clinical interpretation.

### 2.4 Logistic Regression

Logistic regression models the log-odds of the outcome (resistant = 1) as a linear combination of input features:

```
log( p / (1-p) ) = β₀ + β₁x₁ + β₂x₂ + ... + βₙxₙ
```

**Regularisation** penalises large coefficients to prevent overfitting. Two types are available:

- **L2 (Ridge):** adds the sum of squared coefficients to the loss. Shrinks all coefficients towards zero but rarely to exactly zero. Good default.
- **L1 (Lasso):** adds the sum of absolute coefficients. Drives many coefficients to exactly zero, performing implicit feature selection. Particularly relevant here — since resistance signal is concentrated in *rpoB*, L1 should produce a sparse model with high weight on a handful of RRDR variants.

**`C`** is the inverse of regularisation strength. A small C (e.g. 0.01) applies strong regularisation — a simpler model with higher bias. A large C (e.g. 100) applies little regularisation — a more complex model that may overfit. Experiment with C across several orders of magnitude.

> **💡 Clinical relevance of LR coefficients**
> Because each coefficient β corresponds directly to one SNP, a positive β means that SNP is associated with resistance and a negative β with susceptibility. For rifampicin, you would expect large positive coefficients at known *rpoB* RRDR positions. Any large coefficients outside *rpoB* are worth scrutinising — they may reflect population structure rather than direct resistance mechanisms.

### 2.5 Decision Trees

A decision tree recursively partitions the feature space using binary splits on individual SNPs. At each node the algorithm selects the feature and threshold that best separates resistant from susceptible samples according to an impurity criterion (Gini impurity or information gain / entropy).

**Key hyperparameters:**

- **`max_depth`:** controls how many levels the tree can grow. A depth of 1 is a single split (stump); unlimited depth means the tree can perfectly memorise the training data (overfitting). For rifampicin, a very shallow tree (depth 2–3) may already achieve high performance, since a single *rpoB* RRDR mutation is often sufficient to confer resistance.
- **`min_samples_split`:** the minimum number of samples required to split a node. Higher values prevent splits on small, potentially noisy groups.
- **`min_samples_leaf`:** the minimum number of samples that must remain in a leaf. Increasing this smooths the decision boundaries.
- **`criterion`:** Gini impurity and entropy (information gain) almost always give equivalent results in practice.

> **💡 Overfitting vs underfitting**
> Run the same model on the training and test sets and compare performance. A model that scores 0.98 on training but 0.62 on test is overfitting. A model that scores 0.65 on both is underfitting. For rifampicin, even simple models tend to perform well — so pay close attention to sensitivity specifically. Are all resistant samples being caught?

### 2.6 Performance metrics in an AMR context

Accuracy alone is a poor metric for AMR prediction because resistance is often the minority class. A classifier that always predicts susceptible can achieve high accuracy while being clinically useless. The following metrics give a more complete picture:

| Metric | Formula | Interpretation |
|---|---|---|
| Sensitivity (Recall) | TP / (TP + FN) | Proportion of true resistants correctly identified |
| Specificity | TN / (TN + FP) | Proportion of true susceptibles correctly identified |
| F1 Score | 2 × (Prec × Rec) / (Prec + Rec) | Harmonic mean of precision and recall — useful for imbalanced classes |
| AUC | Area under ROC curve | Probability that the model ranks a random resistant above a susceptible sample |
| Accuracy | Correct / Total | Overall correctness — misleading with class imbalance |

For rifampicin resistance prediction, **sensitivity is the priority** — a missed resistant case means a patient receives an ineffective regimen and may deteriorate or transmit MDR-TB. However, specificity still matters: falsely labelling a susceptible patient as resistant leads to unnecessary use of toxic, expensive second-line drugs.

---

## 3. Software Used

The environment is pre-installed in your Codespace. The following tools are used in this workshop:

| Software | Purpose |
|---|---|
| **bcftools** | Indexing and querying the multi-sample VCF |
| **Python 3** | All scripting |
| **pandas** | Data loading and manipulation |
| **numpy** | Numerical operations and missing value handling |
| **scikit-learn** | Model training, imputation, and evaluation metrics |
| **matplotlib** | Visualisation plots |

---

## 4. Running the Pipeline

### File structure

```
machine_learning/
├── data/
│   ├── samples.vcf.gz          # Multi-sample VCF (provided)
│   └── phenotypes.csv          # Drug susceptibility labels (provided)
├── 1_encode_features.py
├── 2_split_data.py
├── 3_train_evaluate.py
├── 4_visualise.py
└── plots/                      # Created automatically by script 4
```

---

### Step 1 — Feature Encoding

Convert the multi-sample VCF into a binary SNP matrix. This uses bcftools to query raw genotype calls and encodes them numerically in Python.

```bash
conda activate ml-amr
python 1_encode_features.py \
    --vcf data/samples.vcf.gz \
    --out_csv features.csv
```

**Expected output:** `features.csv` — rows = samples, columns = SNPs named `CHROM_POS_ALT`.

**Check:** open the file and verify the index column contains your sample IDs and that values are 0, 0.5, 1, or NaN.

---

### Step 2 — Data Splitting

Split the feature matrix and phenotype labels into 80% training and 20% test sets. The random seed ensures reproducibility.

```bash
python 2_split_data.py \
    --geno_matrix features.csv \
    --metadata data/phenotypes.csv \
    --seed 42
```

**Output files:** `train_geno.csv`, `train_pheno.csv`, `test_geno.csv`, `test_pheno.csv`

> **💡 Why fix the random seed?**
> With a fixed seed the split is identical every time you run the script, so results are reproducible and comparable across students. Try changing `--seed` to a different value and re-running — does model performance change? This is a useful exercise in understanding variance due to data splitting.

---

### Step 3 — Train, Evaluate, and Extract Feature Importances

This is the main script. It trains a model, evaluates it on the held-out test set, prints metrics to the terminal, and saves results to CSV — all without saving the model to disk.

#### Logistic Regression — L2 penalty (default)

```bash
python 3_train_evaluate.py \
    --geno_train  train_geno.csv \
    --pheno_train train_pheno.csv \
    --geno_test   test_geno.csv \
    --pheno_test  test_pheno.csv \
    --model   LR \
    --drug    Rifampicin \
    --C       1.0 \
    --penalty l2 \
    --prefix  run_LR_C1
```

#### Logistic Regression — L1 penalty with strong regularisation

```bash
python 3_train_evaluate.py \
    --geno_train  train_geno.csv \
    --pheno_train train_pheno.csv \
    --geno_test   test_geno.csv \
    --pheno_test  test_pheno.csv \
    --model   LR \
    --drug    Rifampicin \
    --C       0.01 \
    --penalty l1 \
    --prefix  run_LR_L1_C001
```

#### Decision Tree — shallow (depth 3)

```bash
python 3_train_evaluate.py \
    --geno_train  train_geno.csv \
    --pheno_train train_pheno.csv \
    --geno_test   test_geno.csv \
    --pheno_test  test_pheno.csv \
    --model     DT \
    --drug      Rifampicin \
    --max_depth 3 \
    --prefix    run_DT_depth3
```

#### Decision Tree — no depth limit

```bash
python 3_train_evaluate.py \
    --geno_train  train_geno.csv \
    --pheno_train train_pheno.csv \
    --geno_test   test_geno.csv \
    --pheno_test  test_pheno.csv \
    --model     DT \
    --drug      Rifampicin \
    --max_depth 0 \
    --prefix    run_DT_nolimit
```

> **Note:** `--max_depth 0` means no limit (fully grown tree). You may observe overfitting — high training performance, lower test performance.

---

### Step 4 — Visualise and Compare

Pass the metrics and feature importance CSVs from Step 3 to produce comparison plots. Include as many runs as you like.

```bash
python 4_visualise.py \
    --metrics \
        run_LR_C1_Rifampicin_LR_metrics.csv \
        run_LR_L1_C001_Rifampicin_LR_metrics.csv \
        run_DT_depth3_Rifampicin_DT_metrics.csv \
        run_DT_nolimit_Rifampicin_DT_metrics.csv \
    --fi \
        run_LR_C1_Rifampicin_LR_feature_importance.csv \
        run_LR_L1_C001_Rifampicin_LR_feature_importance.csv \
        run_DT_depth3_Rifampicin_DT_feature_importance.csv \
        run_DT_nolimit_Rifampicin_DT_feature_importance.csv \
    --out_dir plots/
```

**Output:** `plots/performance_comparison.png`, `plots/performance_radar.png`, `plots/feature_importance_comparison.png`

---

## 5. Hyperparameter Experiments

The table below summarises the hyperparameters available, values to try, and the expected effect. Run multiple combinations and compare with script 4.

| Parameter | Values to try | Effect |
|---|---|---|
| `--penalty` | `l1`, `l2` | L1 = sparse coefficients (feature selection); L2 = shrinks all coefficients |
| `--C` | `0.01, 0.1, 1, 10, 100` | Low C = strong regularisation; High C = fits training data more closely |
| `--max_iter` | `100, 1000, 10000` | Increase if you see a ConvergenceWarning |
| `--max_depth` (DT) | `2, 3, 5, 10, 0` | Small = underfitting; 0 = unlimited = likely overfitting |
| `--min_samples_split` | `2, 10, 20, 50` | Higher values prevent splits on small noisy groups |
| `--min_samples_leaf` | `1, 5, 10, 20` | Higher values create smoother decision boundaries |
| `--criterion` | `gini`, `entropy` | Usually equivalent; entropy slightly more expensive |

**Suggested experiment sequence:**

1. Baseline LR: `--C 1.0 --penalty l2`
2. Strong regularisation: `--C 0.01` — does F1 drop? Does overfitting reduce?
3. Switch to L1: same C, `--penalty l1` — how many features have non-zero coefficients? Do they map to *rpoB*?
4. Baseline DT: `--max_depth 5`
5. Shallow DT: `--max_depth 2` — can a 2-level tree predict rifampicin resistance well? Why?
6. Unlimited DT: `--max_depth 0` — compare training vs test performance explicitly
7. Visualise all runs together with script 4 and discuss

---

## 6. Interpreting Your Results

### 6.1 Reading the terminal output

Script 3 prints a results block for each run:

```
============================================================
  Drug   : Rifampicin
  Model  : Logistic Regression | penalty=l2 | C=1.0
============================================================
  Accuracy   : 0.961
  Sensitivity: 0.943
  Specificity: 0.971
  F1 Score   : 0.952
  AUC        : 0.957

  Top 10 features:
  Feature                        Importance
  NC_000962.3_761155_T            1.842
  NC_000962.3_761139_C            1.204
  NC_000962.3_761161_A            0.987
  ...
============================================================
```

- **For LR:** Importance = the model coefficient. Positive = associated with resistance; negative = associated with susceptibility. Magnitude indicates strength of association.
- **For DT:** Importance = mean decrease in Gini impurity contributed by that feature across all splits. Values sum to 1.0 across all features.

### 6.2 Validating top features against known *rpoB* resistance loci

Cross-reference your top features against the *rpoB* RRDR. In *M. tuberculosis* H37Rv (NC_000962.3), the RRDR spans approximately positions **761,110–761,190**. Key positions to recognise:

| H37Rv position (approx.) | Codon | Common mutation |
|---|---|---|
| ~761,155 | S450 | S450L — most common RIF resistance mutation globally |
| ~761,139 | H445 | H445Y / H445D |
| ~761,101 | D435 | D435V |
| ~761,161 | L452 | L452P |

If your top features fall within this region, the model has learned genuine biology. If high-importance features map outside *rpoB*, consider two explanations:

- **Population structure:** variants that differ between lineages may correlate with resistance without causing it, if resistance prevalence differs between lineages in your dataset. This is why lineage filtering before training is important.
- **Overfitting to noise:** increase regularisation and re-run.

> **💡 A useful sanity check**
> Run with L1 regularisation and a moderately strong C (e.g. 0.1). Count how many features have non-zero coefficients. For rifampicin, a well-regularised model may achieve high performance with only a handful of features — most of which should be *rpoB* RRDR variants. If you see dozens of non-zero coefficients at non-*rpoB* positions, try reducing C further.

### 6.3 Common pitfalls

- **Perfect training performance, poor test performance:** classic overfitting. Increase regularisation (lower C) or constrain tree depth.
- **Low sensitivity with high specificity:** the model is predicting mostly susceptible. Check class balance — if resistant samples are rare, the model may not have enough examples to learn from.
- **Top features change completely between runs:** model instability, often due to correlated features (linkage disequilibrium between *rpoB* variants). L1 regularisation or increasing `--min_samples_leaf` can help.
- **ConvergenceWarning for LR:** increase `--max_iter` (try 10000).
- **Very high performance on all metrics:** this is expected for rifampicin — the RRDR signal is strong. The interesting scientific question is whether a very simple model achieves comparable performance to a complex one, and which features it uses.

---

## 7. Phenotype File Format

The phenotype file must be a CSV with the following structure:

```
sample_id,Rifampicin
ERR1234567,1
ERR1234568,0
ERR1234569,1
ERR1234570,0
ERR1234571,
```

- **First column:** sample identifiers — must match exactly the sample names in the VCF
- **`Rifampicin` column:** binary — `1` = resistant, `0` = susceptible
- **Empty cells:** allowed — samples with missing labels are automatically excluded
- **Column name:** must match exactly what you pass to `--drug` (i.e. `--drug Rifampicin`)

> **💡 Sample ID mismatch is the most common error**
> If you get a KeyError or empty results, run `bcftools query -l samples.vcf.gz` and compare the output with the first column of your phenotype file. Even a single trailing space or character difference will cause a mismatch.

---

## 8. Discussion Questions

Work through these questions with your group after completing the experiments:

1. How does changing `C` in logistic regression affect model performance and the number of features with non-zero coefficients when using L1? At what value of C does performance start to drop noticeably, and why?

2. Compare the top-10 features from logistic regression and the decision tree. Are they the same features? Do they all map to the *rpoB* RRDR? If not, what might explain the differences?

3. A shallow decision tree (depth 2–3) may achieve very high performance for rifampicin. What does this tell you about the genetic architecture of rifampicin resistance compared to a polygenic trait? How would you expect results to differ if you ran the same pipeline on a more complex phenotype?

4. For clinical deployment of a WGS-based resistance prediction tool, what sensitivity threshold would you set as a minimum? How does this trade off against specificity, and who should make that decision — the microbiologist, the clinician, or the public health team?

5. The dataset has been filtered to a single *M. tuberculosis* lineage. Why might this matter? What problems could arise if you trained on one lineage and predicted on isolates from another, even if overall accuracy remained high?

6. Heterozygous calls (encoded as 0.5) at *rpoB* RRDR positions likely indicate mixed infection — a susceptible majority population with a resistant minority. How should these samples be labelled phenotypically, and how would you handle them differently in the pipeline?

7. Rifampicin resistance prediction from WGS is already implemented in tools like TBProfiler and Mykrobe using curated mutation catalogues. What is the advantage of an ML approach? In what scenarios might a catalogue-based approach outperform it?

---

## 9. Quick Reference

### All command-line flags — script 3

```
Required:
  --geno_train          Training genotype CSV
  --pheno_train         Training phenotype CSV
  --geno_test           Test genotype CSV
  --pheno_test          Test phenotype CSV
  --model               LR or DT
  --drug                Drug name (must match phenotype column header exactly)
  --prefix              Output file prefix

Shared optional:
  --top_n               Number of top features to report (default: 10)

Logistic Regression:
  --penalty             l1 or l2 (default: l2)
  --C                   Inverse regularisation strength (default: 1.0)
  --max_iter            Solver iterations (default: 1000)

Decision Tree:
  --max_depth           Tree depth, 0 = unlimited (default: 0)
  --min_samples_split   Min samples to split a node (default: 2)
  --min_samples_leaf    Min samples at a leaf (default: 1)
  --criterion           gini or entropy (default: gini)
```

### Output files from script 3

| File | Contents |
|---|---|
| `{prefix}_{drug}_{model}_metrics.csv` | Accuracy, Sensitivity, Specificity, F1, AUC |
| `{prefix}_{drug}_{model}_feature_importance.csv` | Top N features and their importance scores |

### *rpoB* RRDR quick reference (NC_000962.3)

| Position range | Region |
|---|---|
| 761,110 – 761,190 | Rifampicin resistance-determining region (RRDR) |
| ~761,155 | Codon S450 — S450L most common resistance mutation globally |
| ~761,139 | Codon H445 — H445Y / H445D common |
| ~761,101 | Codon D435 — D435V common |
| ~761,161 | Codon L452 — L452P variable resistance level |
