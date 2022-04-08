# P.chondrillina-project-scripts

## Genome assembly

1, Base calling

guppy_basecaller -i {path to directory of fast5 files} -s {output directory} -c dna_r9.4.1_450bps_sup.cfg --device auto —-min_qscore 7

Version: guppy_ba

basecaller 5.0.11

2, Flye genome assembly

flye --nano-raw longreads.fastq -o {output directory} -t {num}

Version: Flye 2.8.1

3, long reads polishing 

minimap2 -ax map-ont -t {num} assembly.fasta longreads.fastq > longread_map.sam

racon -t {num} -u -m 8 -x -6 -g -8 -w 500 longreads.fastq longread_map.sam assembly.fasta > racon_assembly.fasta

medaka consensus -i longreads.fastq -d racon_assembly.fasta -o {output directory} -m r941_min_sup_g507

Version: minimap2 2.2.2, racon 1.4.20, medaka 1.4.4

4, short reads polishing

set the path to short reads

realpath forward_paired-end_shortreads.fq,gz reverse_paired-end_shortreads.fq.gz >sgs.fofn

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
