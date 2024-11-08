#!/bin/bash

star_path=$1
path_to_gtf=$2
genome_dir=$3
path_to_genome_assembly=$4

threads=16
#the sjdbOverhang should be readlength-1 (eg 49 in case of 50 bp reads).

#Generating the genome
$1/STAR --runMode genomeGenerate --sjdbGTFfile $path_to_gtf --genomeDir $genome_dir= --genomeFastaFiles $path_to_genome_assembly= --sjdbOverhang 100 --runThreadN 36

out=$6
now=`date`

#aligning to genome and generating counts
while read inp
do
	suffix="_R1.fastq.gz"
	bn=$(basename $inp $suffix)
	path=$(dirname $inp)
	fastq1=$path/$bn"_R1.fastq.gz"
	fastq2=$path/$bn"_R2.fastq.gz"
	echo $bn
		(
		echo "starts $now"
		#STAR two pass alignment mode
		$star_path/STAR --genomeDir $genome_dir --runThreadN $threads --readFilesIn $fastq1 $fastq2 --readFilesCommand "gunzip -c"  --outSAMtype BAM SortedByCoordinate --twopassMode Basic --quantMode GeneCounts --outFileNamePrefix $out/$bn"."
		echo "ends $now"
		) 1>>$out$bn".pipeline.log" 2>>$out$bn".pipeline.log" 
done < $5

#collating counts data
cd $out
paste *.ReadsPerGene.out.tab > counts.txt