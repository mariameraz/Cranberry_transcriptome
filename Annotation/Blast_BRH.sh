#PBS -N blast_brh
#PBS -q default
#PBS -l nodes=1:ppn=1,walltime=99999:999:00,mem=30gb,vmem=30gb
#PBS -V

cd $PBS_O_WORKDIR

module load ncbi-blast+/2.6.0

IDX_PATH=/ABS_PATH_TO_IDX_FOLDER

blastp -query ../macrocarpon.proteins.fasta\
 -db $IDX_PATH/uniprot_fasta_ubh.fasta -outfmt 7\
 -max_target_seqs 1 -evalue 1e-5 -out BRH.tab #BRH: "Best reciprocal (bidirectional) hit
 
 #Get ids
 grep 'sp|' BRH.tab | cut -f 1,2 > ids.txt
