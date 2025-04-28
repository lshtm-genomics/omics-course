```
conda activate assembly
cd ~/assembly/assembly
ls tb_ILL

jellyfish count -C -m 21 -s 100M -t 10 -o reads.jf <(zcat tb_ILL/*.fastq.gz)

jellyfish histo reads.jf > reads.histo

spades.py -1  tb_ILL/sample1_1.fastq.gz -2  tb_ILL/sample1_2.fastq.gz -o short/ -k 55 --isolate -t 4

quast -r tbdb.fasta -o quast/short short/scaffolds.fasta 

conda deactivate

conda activate busco

busco -i short/scaffolds.fasta -l bacteria_odb10 -o busco_output -m genome -c 4

conda deactivate

conda activate assembly

nano gapclosing.txt

GapCloser -a short/scaffolds.fasta -b gapclosing.txt -o short/scaffolds_gapClosed.fasta


flye --nano-hq tb_ONT/sample1_ONT.fastq.gz --genome-size 4.1m --threads 8 --read-error 0.06 -o long/

ntLink scaffold target=long/assembly.fasta reads=tb_ONT/sample1_ONT.fastq.gz G=500 rounds=3 t=4

mkdir long/tgs_gapcloser

tgsgapcloser --scaff long/assembly.fasta.k32.w100.z1000.stitch.abyss-scaffold.fa --reads tb_ONT/sample1_ONT.fastq.gz --output long/tgs_gapcloser --ne

cp long/tgs_gapcloser.scaff_seqs long/tgs_gapcloser.scaff_seqs.fa

minimap2 -t 4 -x map-ont long/tgs_gapcloser.scaff_seqs.fa tb_ONT/sample1_ONT.fastq.gz > racon.paf

racon -u --no-trimming -t 4 tb_ONT/sample1_ONT.fastq.gz racon.paf long/tgs_gapcloser.scaff_seqs.fa > long/final_assembly.fa

spades.py -1 tb_ILL/sample1_1.fastq.gz -2 tb_ILL/sample1_2.fastq.gz --nanopore tb_ONT/sample1_ONT.fastq.gz -t 4 -o hybrid/spades/

masurca masurca_config
bash assemble.sh

quast -o quast/hybrid/masurca CA.mr.55.17.15.0.02/primary.genome.scf.fasta
```