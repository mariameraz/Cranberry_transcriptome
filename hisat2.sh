#!/bin/bash

#Aligning reads to the genome reference with hisat2

for i in $(more list.txt)
do
    INDEX_PATH=/home/alejandra/Documentos/DR_MURA/Hisat2/Index

hisat2 -x $INDEX_PATH/V_macrocarpo_Stevens_v1\
 -1 ${i}_1_sequence.txt.out.gz\
 -2 ${i}_2_sequence.txt.out.gz -S {i}.sam 

done
