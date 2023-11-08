#!/bin/bash


NPROC=$(nproc)
#outdir=minimap2.proG.contigs
input=fastq
index=/mnt/d/data/databases/EMU/species_taxid.fasta

mkdir $outdir

## Run minimap2
for fastq in $(ls $input); do
	base=$(basename $input/$fastq .fastq)
	echo 'started mapping ' $base ' to ' $index ' on:'
	date
	emu abundance --type map-ont --db /mnt/i/data/databases/EMU -N 5 --keep-files --keep-counts --keep-read-assignments --output-unclassified --threads 4 /mnt/i/data/Gemini/BactMix/fastq/$fastq
	## --MD -ax
	#samtools view ${outdir}/${base}.sam -b -@ ${NPROC} -o ${outdir}/${base}.bam
	#samtools sort  -@ ${NPROC} -o ${outdir}/${base}.bam ${outdir}/${base}'.bam' #-T reads.tmp 
	#samtools index -@ ${NPROC} ${outdir}/${base}.bam 
	echo 'mapping done, on:'
	date
	#
done


