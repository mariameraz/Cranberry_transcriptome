#PBS -N hisat.idx
#PBS -q default
#PBS -l nodes=1:ppn=1,walltime=999:999:00,mem=30gb,vmem=30gb
#PBS -V

cd $PBS_O_WORKDIR

module load hisat2/2.1.0

hisat2-build -p 16 -f V_macrocarpon_Stevens_v1.fasta macrocarpon.idx
