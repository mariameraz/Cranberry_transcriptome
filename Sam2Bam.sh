#PBS -N sam2bam
#PBS -q default
#PBS -l nodes=1:ppn=1,walltime=999:999:00,mem=30gb,vmem=30gb
#PBS -V

cd $PBS_O_WORKDIR

module load samtools/1.9

for i in $(more list.txt)
do

#sam to bam file
samtools view -bS ${i}.sam > ./bamfiles/${i}.bam

#sorting bam files
samtools sort ./bamfiles/${i}.bam -o ./bamfiles/${i}.sorted.bam

done
