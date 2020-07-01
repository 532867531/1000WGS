setwd("D:/华大/1000WGS/300GATK/01_GenotypeGVCFs/log")
(files=list.files(pattern = "//d+$"))

data=read.delim(file = "D:/华大/1000WGS/300GATK/01_GenotypeGVCFs/log/log.log",stringsAsFactors = F,header = F)
rows=c()
for(rowindex in c(1:nrow(data))){
  if(grepl(pattern = "Total runtime",x = data[rowindex,"V1"])){
    rows=c(rows,data[rowindex,"V1"])
  }
}

write.table(x = rows,file = "log1.log",row.names = F,col.names = F,quote = F)
