#PBS -N hisat_alingment
#PBS -q default
#PBS -l nodes=1:ppn=8,walltime=999:999:00,mem=30gb,vmem=30gb
#PBS -V

cd $PBS_O_WORKDIR

module load hisat2/2.1.0

INDEX_PATH=/PATH/index
FASTA_PATH=/PATH/Fastp

for i in $(more list.txt)
do

hisat2 -p 8 --dta -q $INDEX_PATH/macrocarpon.idx -1 $FASTA_PATH/${i}_1_sequence.txt.out.gz -2 $FASTA_PATH/${i}_2_sequence.txt.out.gz -S ${i}.sam --new-summary

done
