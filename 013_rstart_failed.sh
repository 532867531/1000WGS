cd /data/szj/software/gatk3.8/interval_test/300GATK/01_GenotypeGVCFs
for onefile in `ls *.PBS`;do
	#echo ${onefile};
	file=${onefile/GenotypeGVCFs/}
	file=${file/.intervals.PBS/}
	
	outvcffile="/data/szj/software/gatk3.8/interval_test/300GATK/011vcfs_intervals_result/5${file}.intervalsGenotypeGVCFs.vcf.idx"
	if [ ! -f  $outvcffile ];then
		dos2unix $onefile $onefile;
		echo $outvcffile;
		qsub $onefile && sleep 1s;
	fi
done