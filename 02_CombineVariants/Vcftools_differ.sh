#PBS -N catVariants.compare
#PBS -q ZJU
#PBS -l nodes=1:ppn=1
#PBS -l walltime=00:60:00
cd /data/szj/software/gatk3.8

vcftools --vcf /data/szj/1000C/5GenotypeGVCFs.vcf --diff /data/szj/software/gatk3.8/interval_test/300GATK/011vcfs_intervals_result/5all_interVals_GenotypeGVCFs.vcf --diff-site --out /data/szj/software/gatk3.8/interval_test/300GATK/02_CombineVariants/Diff.site