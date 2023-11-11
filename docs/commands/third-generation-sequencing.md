```
conda activate nanopore
cd ~/data/nanopore_activity/basecalling/raw_fast5_reads/ 
guppy_basecaller --config dna_r9.4.1_450bps_fast.cfg --trim_adapters --compress_fastq  --input_path ~/data/nanopore_activity/basecalling/raw_fast5_reads  --save_path ~/data/nanopore_activity/basecalling/fastq 

cd ~/data/nanopore_activity/basecalling/fastq/pass
pycoQC -f ~/data/nanopore_activity/basecalling/fastq/sequencing_summary.txt -o pycoqc_results.html
firefox ~/data/nanopore_activity/basecalling/fastq/pass/pycoqc_results.html &

porechop -i ~/data/nanopore_activity/basecalling/fastq/pass/*.fastq.gz -o basecalled_reads.porechop.fastq

cd ~/data/nanopore_activity/kraken
kraken --db ~/data/nanopore_activity/kraken/KDB/ --output temp.krak ~/data/nanopore_activity/kraken/basecalled_reads.porechop.fastq
rcf -k temp.krak -o rcf.html
firefox rcf.html &

cd ~/data/nanopore_activity/mapping
minimap2 -ax map-ont ./reference.fasta basecalled_reads.fastq | samtools view -q 15 -b | samtools sort -o alignment.bam

cd ~/data/nanopore_activity/variant_calling 
samtools depth ~/data/nanopore_activity/mapping/alignment.bam > depth_statistics

R
data<-read.table("depth_statistics")
plot(data$V3,type="l",xlab="Reference Position", ylab="read Depth")
quit()

bcftools mpileup -q 8 -B -I -Ou -f reference.fasta ~/data/nanopore_activity/mapping/alignment.bam | bcftools call -mv -Oz -o calls.vcf.gz
bcftools index calls.vcf.gz

bcftools consensus -f reference.fasta calls.vcf.gz -o consensus_sequence.fasta

cd ~/data/nanopore_activity/phylogenetics
aliview zika_dataset.fasta
cat zika_dataset.fasta consensus_sequence.fasta > zika_all.fasta
mafft zika_all.fasta > zika_all_aligned.fasta
iqtree -m GTR+G -s ./zika_all_aligned.fasta -bb 1000 -nt AUTO

figtree ./zika_all_aligned.fasta.treefile


```