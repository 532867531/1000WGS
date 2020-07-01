#parameters
vcf_output_dir="/data/szj/software/gatk3.8/interval_test/300GATK/011vcfs_intervals_result/"
interval_dir="/data/szj/software/gatk3.8/interval_test/300GATK/00intervals/"

setwd("D:/华大/1000WGS/300GATK/01_GenotypeGVCFs")
model=read.delim(file = "D:/华大/1000WGS/300GATK/01_GenotypeGVCFs/GenotypeGVCFs-538.model",stringsAsFactors = F,header = F,col.names = F)
model=read.delim(file = "D:/华大/1000WGS/300GATK/01_GenotypeGVCFs/GenotypeGVCFs-1000.model",stringsAsFactors = F,header = F,col.names = F)

(files=list.files(path = "D:/华大/1000WGS/300GATK/00intervals",pattern = "file\\d+.intervals"))
for(onefile in files){
  the_model=model
  for(rowindex in c(1:nrow(model))){
    the_model[rowindex,"FALSE."]=sub(pattern = "_RANGE_",replacement = onefile,x = the_model[rowindex,"FALSE."])
    the_model[rowindex,"FALSE."]=sub(pattern = "_VCFOUTPUTDIR_",replacement = vcf_output_dir,x = the_model[rowindex,"FALSE."])
    the_model[rowindex,"FALSE."]=sub(pattern = "_INTERVALDIR_",replacement = interval_dir,x = the_model[rowindex,"FALSE."])
  }
  write.table(x = the_model,file = paste("GenotypeGVCFs",onefile,".PBS",sep = ""),row.names = F,col.names = F,quote = F)
}

