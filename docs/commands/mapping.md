``` bash
conda activate mapping
cd ~/data/tb/
bwa index ~/data/tb/tb.fasta

# sample1
trimmomatic PE sample1_1.fastq.gz sample1_2.fastq.gz -baseout sample1.fastq LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:36
bwa mem -R "@RG\tID:sample1\tSM:sample1\tPL:Illumina" ~/data/tb/tb.fasta sample1_1P.fastq sample1_2P.fastq | samtools view -b - | samtools sort -o sample1.bam -
samtools index sample1.bam
rm sample1_1P.fastq sample1_1U.fastq sample1_2P.fastq sample1_2U.fastq

# sample2 
trimmomatic PE sample2_1.fastq.gz sample2_2.fastq.gz -baseout sample2.fastq LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:36
bwa mem -R "@RG\tID:sample2\tSM:sample2\tPL:Illumina" ~/data/tb/tb.fasta sample2_1P.fastq sample2_2P.fastq | samtools view -b - | samtools sort -o sample2.bam -
samtools index sample2.bam
rm sample2_1P.fastq sample2_1U.fastq sample2_2P.fastq sample2_2U.fastq

cd ~/data/malaria/
bwa index ~/data/malaria/Pf3D7_05.fasta
bwa mem ~/data/malaria/Pf3D7_05.fasta ~/data/malaria/IT.Chr5_1.fastq.gz ~/data/malaria/IT.Chr5_2.fastq.gz | samtools view -b - | samtools sort -o IT.Chr5.bam -
samtools index IT.Chr5.bam
samtools view IT.Chr5.bam | grep "IL39_6014:8:61:7451:18170"

```