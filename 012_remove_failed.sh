cd /data/szj/software/gatk3.8/interval_test/300GATK/011vcfs_intervals_result
for onefile in `ls *.vcf`;do
	if [ ! -f "${onefile}.idx" ];then
		echo ${onefile};
		rm -rf ${onefile}
	fi
done