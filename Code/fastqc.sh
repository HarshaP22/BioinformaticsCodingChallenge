#!/bin/bash

mkdir $3/FastQC
$1/fastqc $2/*.fastq.gz -t 24 -o $3/FastQC
