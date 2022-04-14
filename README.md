# P.chondrillina-project-scripts

## Genome assembly

1, Base calling

guppy_basecaller -i {path to directory of fast5 files} -s {output directory} -c dna_r9.4.1_450bps_sup.cfg --device auto —-min_qscore 7

NanoPlot --summary sequencing_summary.txt -o {output folder}

Version: guppy_basecaller 5.0.11, NanoPlot 1.32.1

2, Flye genome assembly

flye --nano-raw longreads.fastq -o {output directory} -t {num}

Version: Flye 2.8.1

3, long reads polishing with racon and medaka

minimap2 -ax map-ont -t {num} assembly.fasta longreads.fastq > longread_map.sam

racon -t {num} -u -m 8 -x -6 -g -8 -w 500 longreads.fastq longread_map.sam assembly.fasta > racon_assembly.fasta

medaka consensus -i longreads.fastq -d racon_assembly.fasta -o {output directory} -m r941_min_sup_g507

Version: minimap2 2.2.2, racon 1.4.20, medaka 1.4.4

4, trim short read

trimmomatic PE -threads {num} -phred33 shortread1.fq.gz shortread2.fq.gz output_shortread1_paired.fq.gz output_shortread2_unpaired.fq.gz
output_shortread2_paired.fq.gz  output_shortread2_unpaired.fq.gz ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 HEADCROP:3 MINLEN:36

5, short reads polishing with nextPolish

set the path to short reads

realpath output_shortread1_paired.fq.gz output_shortread2_paired.fq.gz >sgs.fofn

create a file to set up the polishing steps:
nano run.cfg

Here is the content of run.cfg
***********************************************
[General]

job_type = local

job_prefix = nextPolish

task = 1212

rewrite = yes

rerun = 3

parallel_jobs = 4

multithread_jobs = 12

genome ={path to assembly file}

genome_size = auto

workdir ={working directory}

polish_options = -p {multithread_jobs}

[sgs_option] #optional

sgs_fofn = sgs.fofn

sgs_options = -max_depth 100 -minimap2

***********************************************
run the polihing steps:
nextPolish run2.cfg

Version: nextPolish 1.4.0

5, short reads polishing with Hapo-g

python3 hapog.py -g assembly.fasta --pe1 output_shortread1_paired.fq.gz --pe2 output_shortread2_paired.fq.gz  -o 
{outputpath} -t {num}

Version: Hapo-g 1.2 

6, statistic summary of assembly

ststs.sh in=assembly.fasta

Version: bbmap 38.93

7, use BUSCO to test completeness of assembly

busco --offline -f -m genome -l basidiomycota_odb10 -c {num} -i assembly.fasta -o {output folder}

8, blast on contigs

download nt databese:

mkdir -p nt \
wget "ftp://ftp.ncbi.nlm.nih.gov/blast/db/nt.??.tar.gz" -P nt/ && \
        for file in nt/*.tar.gz; \
            do tar xf $file -C nt && rm $file; \
        done

blastn on assembly:

DB = nt

blastn -db $DB/nt \
       -query assembly.fasta \
       -outfmt "6 qseqid staxids bitscore std" \
       -max_target_seqs 10 \
       -max_hsps 1 \
       -evalue 1e-25 \
       -num_threads {num} \
       -out {output name}
       
       
download mitochondrial database:

wget https://ftp.ncbi.nlm.nih.gov/blast/db/mito.tar.gz
tar –xvzf mito.tar.gz

blastn on mitochondrial contigs:

DB = mito
blastn -db $DB/nt \
       -query assembly.fasta \
       -outfmt "6 qseqid staxids bitscore std" \
       -max_target_seqs 10 \
       -max_hsps 1 \
       -evalue 1e-25 \
       -num_threads {num} \
       -out {output name}
       
Version: Blast 2.12.0

9, extract subset sequences from assembly:

seqtk subseq asssembly.fast contig_names.lst > subset.fasta

Version: suqtk 1.3

10, mapping reads against assembly without secondary alignments:

For nanopore longread data:

minimap2 -ax map-ont --secondary=no -t {num} assembly.fasta longreads.fastq | samtools sort -@ {num} -O BAM -o {name of mapping file}.bam

For paired-end shortread data:

minimap2 -ax sr --secondary=no -t {num} assembly.fasta shortread1.fq shortread2.fq | samtools sort -@ {num} -O BAM -o {name of mapping file}.bam

Version: minimap2 2.17 samtools 1.3.1

11, Use blobtools to visualize the taxonomic clssification of contigs, GC content, and coverage of assembly

