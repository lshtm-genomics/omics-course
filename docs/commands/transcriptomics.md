```
conda activate rnaseq
cd ~/data/transcriptomics


bwa index H37Rv.fa

bwa mem -t 2 H37Rv.fa Mtb_L1_1.fastq.gz Mtb_L1_2.fastq.gz | samtools sort -@ 2 - -o Mapping_Mtb/Mtb_L1.bam
samtools index Mapping_Mtb/Mtb_L1.bam
bwa mem -t 2 H37Rv.fa Mtb_L4_1.fastq.gz Mtb_L4_2.fastq.gz | samtools sort -@ 2 - -o Mapping_Mtb/Mtb_L4.bam
samtools index Mapping_Mtb/Mtb_L4.bam

python -m HTSeq.scripts.count -f bam -r pos -s reverse -t gene ./Mapping_Mtb/Mtb_L1.bam Mtb.gtf > ./Mapping_Mtb/Mtb_L1_htseq_count.txt
python -m HTSeq.scripts.count -f bam -r pos -s reverse -t gene ./Mapping_Mtb/Mtb_L4.bam Mtb.gtf > ./Mapping_Mtb/Mtb_L4_htseq_count.txt

cd ~/data/transcriptomics/Mapping_Mtb/HTSeqCounts

R
library(DESeq2)
library(gplots)

directory <- "~/data/transcriptomics/Mapping_Mtb/HTSeqCounts/" 

sampleFiles <- grep("Mtb", list.files(directory), value = TRUE)

sampleCondition <- c("l4","l1","l1","l4","l1","l4")

sampleTable <- data.frame(
    sampleName = sampleFiles,
    fileName = sampleFiles,
    condition = sampleCondition
)

dds <- DESeqDataSetFromHTSeqCount(
    sampleTable = sampleTable,
    directory = directory,
    design = ~ condition
)

keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]

dds$condition <- factor(dds$condition, levels = c("l1","l4"))
dds <- DESeq(dds)
res <- results(dds)

summary(res)
res <- res[!is.na(res$padj),]
resOrdered <- res[order(res$pvalue),]
resOrdered

counts_heatmap <- counts(dds, normalized = TRUE)
idx <- rownames(resOrdered)[1:17]
counts_heatmap <- counts_heatmap[rownames(counts_heatmap)%in%idx,]

counts_heatmap

colnames(counts_heatmap) <- c("L4_1","L1_2","L1_3","L4_4","L1_5","L4_6")
heatmap.2(as.matrix(counts_heatmap), scale="row", col=greenred(75), Rowv=NA, dendrogram = "col", trace="none", density.info = "none")


```