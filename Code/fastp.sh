#!/bin/bash

now=`date`

while read inp
do
	suffix="_R1.fastq.gz"
	bn=$(basename $inp $suffix)
	path=$(dirname $inp)
	fastq1=$path/$bn"_R1.fastq.gz"
	fastq2=$path/$bn"_R2.fastq.gz"
	
	out1=$path/Trimmed/$bn"_trimmed_R1.fastq.gz"
	out2=$path/Trimmed/$bn"_trimmed_R2.fastq.gz"
	
	echo $bn
		(
		echo "starts $now"

		$1/fastp -i $fastq1 -I $fastq2 -o $out1 -O $out2

		echo "ends $now"
		) 1>>$out$bn".pipeline.log" 2>>$out$bn".pipeline.log" 
done < $2
