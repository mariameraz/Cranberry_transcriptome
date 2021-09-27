#PBS -N featurecounts
#PBS -q default
#PBS -l nodes=1:ppn=1,walltime=999:999:00,mem=30gb,vmem=30gb
#PBS -V

cd $PBS_O_WORKDIR

module load subread/1.5.1

featureCounts -T 4 -t exon -g gene_id -a V_macrocarpon_Stevens_v1-geneModels.gtf -o counts.cranberry.txt *.sorted.bam
