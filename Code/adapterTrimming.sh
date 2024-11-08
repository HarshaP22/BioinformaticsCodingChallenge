#!/bin/bash

$1/fastp -i /home/bio3/current_data/bioinformatics/harsha/bioinformatics/Heart_ZT0_1_R1.fastq.gz -I /home/bio3/current_data/bioinformatics/harsha/bioinformatics/Heart_ZT0_1_R2.fastq.gz -o /home/bio3/current_data/bioinformatics/harsha/bioinformatics/Trimmed/Heart_ZT0_1_trimmed_R1.fastq.gz -O /home/bio3/current_data/bioinformatics/harsha/bioinformatics/Trimmed/Heart_ZT0_1_trimmed_R2.fastq.gz 

bash /home/bio3/current_data/bioinformatics/harsha/fastp.sh /home/bio3/current_data/bioinformatics/harsha/bioinformatics/list_fastq.txt
