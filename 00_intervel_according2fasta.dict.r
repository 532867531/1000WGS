setwd("D:/����/1000WGS/300GATK")
list.files()

##��ȡ�ο��������ֵ��ļ�dict
dict=read.delim(file = "C:/Users/Administrator/Documents/Homo_sapiens_assembly38.dict",stringsAsFactors = F,header = F)
dict$V3=as.numeric(sub(pattern = "LN:",replacement = "",x = dict$V3))
dict$V2=sub(pattern = "SN:",replacement = "",x = dict$V2)
##ѡȡ��Ҫ��contig
dict=dict[c(2:3367),]
dict=dict[c(2:26),]

##���մ��interval���б��Լ�interval ��С
intervals=list()
length_intervals=60000000

##�����ܹ���Ҫ�ָ��bp����
sums=sum(dict$V3)

dataframe1=data.frame(CUM=cumsum(dict$V3),CHR=dict$V2,LEN=dict$V3,stringsAsFactors = F)##ʵ���ۻ�milestone
dataframe2=data.frame(CUM=c(c(1:(sums/length_intervals))*length_intervals,sums),LEN=NA,CHR=NA,stringsAsFactors = F)##�����ۻ�milestone

dataframe=rbind(dataframe1,dataframe2)
dataframe=dataframe[order(dataframe$CUM),]


minus=0
for(rowindex in c(1:nrow(dataframe))){
  if(grepl(pattern = "chr",x = dataframe[rowindex,"CHR"])){
    minus=dataframe[rowindex,"CUM"]                           ##���ڸ���contigʱ���������0
    dataframe[rowindex,"CUTOFF"]=dataframe[rowindex,"LEN"]
  }else{
    dataframe[rowindex,"CUTOFF"]=dataframe[rowindex,"CUM"]-minus  ##���ڸ���contigʱ���������0
  }
}
##�����Ⱦɫ��
label=NA
for(rowindex in c(nrow(dataframe):1)){
  if(grepl(pattern = "chr",x = dataframe[rowindex,"CHR"])){
    label=dataframe[rowindex,"CHR"]
  }else{
    dataframe[rowindex,"CHR"]=label
  }
}
dataframe=dataframe[which(!is.na(dataframe$CHR)),]

##��contig�ۺ�milestone�ڵ�
data=aggregate(x = dataframe,by = list(CHR=dataframe$CHR),FUN = function(x){paste(x,collapse = "|")})
rownames(data)=data$CHR
data=data[unique(dataframe$CHR),]

intervals=c()
intervals_dataframe=na.omit(data.frame(CHR=NA,BEGIN=NA,END=NA,stringsAsFactors = F))


##����ÿ��Ⱦɫ��interval
for(rowindex in c(1:nrow(data))){
  CHR=data[rowindex,"CHR"]
  CUTOFF=as.numeric(c(0,unlist(stringi::stri_split(str = data[rowindex,"CUTOFF"],regex = "\\|"))))
  for(index in c(2:length(CUTOFF))){
      intervals=c(intervals,paste(CHR,CUTOFF[index-1]+1,CUTOFF[index],sep = "-"))
      intervals_dataframe[nrow(intervals_dataframe)+1,c("CHR","BEGIN","END")]=c(CHR,CUTOFF[index-1]+1,CUTOFF[index])
  }
}
clipr::write_clip(content = paste(intervals,sep = "\n"))

##д���ļ�
intervals_dataframe$TOwrite=paste(paste(intervals_dataframe$CHR,intervals_dataframe$BEGIN,sep = ":"),intervals_dataframe$END,sep = "-")
intervals_dataframe$LENGTH=as.numeric(intervals_dataframe$END)-as.numeric(intervals_dataframe$BEGIN)+1
bps=0
fileIndex=1
for(rowindex in c(1:nrow(intervals_dataframe))){
  write.table(x = intervals_dataframe[rowindex,"TOwrite"],file = paste("file",fileIndex,".intervals",sep = ""),append = T,quote = F,row.names = F,col.names = F)
  bps=bps+as.numeric(intervals_dataframe[rowindex,"LENGTH"])
  if(bps>=length_intervals){
    fileIndex=fileIndex+1
    bps=0
  }  
}
shell.exec(getwd())
