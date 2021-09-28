#PBS -N blast_idx
#PBS -q default
#PBS -l nodes=1:ppn=1,walltime=999:999:00,mem=30gb,vmem=30gb
#PBS -V

cd $PBS_O_WORKDIR

module load ncbi-blast+/2.6.0

makeblastdb -in macrocarpon.proteins.fasta -dbtype prot

