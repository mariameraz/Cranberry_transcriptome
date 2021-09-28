#PBS -N blast_ubh
#PBS -q default
#PBS -l nodes=1:ppn=1,walltime=99999:999:00,mem=30gb,vmem=30gb
#PBS -V

cd $PBS_O_WORKDIR

module load ncbi-blast+/2.6.0

IDX_PATH=/ABS_PATH_TO_IDX_FOLDER

blastp -query ../uniprot.fasta\
 -db $IDX_PATH/macrocarpon.proteins.fasta -outfmt 7\
 -max_target_seqs 10 -evalue 1e-5 -out uniprot_vs_macro_UBH.tab #UBH:"Undirectional best hits".
 
#Filter best uniprot hits IDs
grep 'sp|' macrocarpon_to_uniprot.tab | cut -f 2 > uniprot_vs_macro_UBH.txt
