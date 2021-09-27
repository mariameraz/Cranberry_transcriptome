#PBS -N fastp
#PBS -l nodes=1:ppn=1,mem=32gb,vmem=32gb,walltime=250:00:00
#PBS -q default
#PBS -V

cd $PBS_O_WORKDIR

module load fastp/0.20

for i in $(more list.txt)
do

/data/software/fastp/fastp -q 30 -i ${i}_1_sequence.txt.gz -I ${i}_2_sequence.txt.gz\
 -o ${i}_1_sequence.txt.out.gz -O ${i}_2_sequence.txt.out.gz\
 -h {i}.html -j ${i}.json

done
