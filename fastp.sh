#!/bin/bash

#Trimming away low quality reads and adapters sequences with fastp
for i in $(more list.txt)
do

PATH=/home/alejandra/Documentos/DR_MURA/Data

fastp\
 -i $PATH/${i}_1_sequence.txt.gz -I $PATH/${i}_2_sequence.txt.gz \
 -o ${i}_1_sequence.txt.out.gz -O ${i}_2_sequence.txt.out.gz

done
