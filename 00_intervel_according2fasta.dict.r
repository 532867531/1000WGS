setwd("D:/华大/1000WGS/300GATK")
list.files()

##读取参考基因组字典文件dict
dict=read.delim(file = "C:/Users/Administrator/Documents/Homo_sapiens_assembly38.dict",stringsAsFactors = F,header = F)
dict$V3=as.numeric(sub(pattern = "LN:",replacement = "",x = dict$V3))
dict$V2=sub(pattern = "SN:",replacement = "",x = dict$V2)
##选取需要的contig
dict=dict[c(2:3367),]
dict=dict[c(2:26),]

##最终存放interval的列表以及interval 大小
intervals=list()
length_intervals=60000000

##汇总总共需要分割的bp长度
sums=sum(dict$V3)

dataframe1=data.frame(CUM=cumsum(dict$V3),CHR=dict$V2,LEN=dict$V3,stringsAsFactors = F)##实际累积milestone
dataframe2=data.frame(CUM=c(c(1:(sums/length_intervals))*length_intervals,sums),LEN=NA,CHR=NA,stringsAsFactors = F)##参照累积milestone

dataframe=rbind(dataframe1,dataframe2)
dataframe=dataframe[order(dataframe$CUM),]


minus=0
for(rowindex in c(1:nrow(dataframe))){
  if(grepl(pattern = "chr",x = dataframe[rowindex,"CHR"])){
    minus=dataframe[rowindex,"CUM"]                           ##用于更换contig时坐标相对置0
    dataframe[rowindex,"CUTOFF"]=dataframe[rowindex,"LEN"]
  }else{
    dataframe[rowindex,"CUTOFF"]=dataframe[rowindex,"CUM"]-minus  ##用于更换contig时坐标相对置0
  }
}
##倒序标染色体
label=NA
for(rowindex in c(nrow(dataframe):1)){
  if(grepl(pattern = "chr",x = dataframe[rowindex,"CHR"])){
    label=dataframe[rowindex,"CHR"]
  }else{
    dataframe[rowindex,"CHR"]=label
  }
}
dataframe=dataframe[which(!is.na(dataframe$CHR)),]

##分contig聚合milestone节点
data=aggregate(x = dataframe,by = list(CHR=dataframe$CHR),FUN = function(x){paste(x,collapse = "|")})
rownames(data)=data$CHR
data=data[unique(dataframe$CHR),]

intervals=c()
intervals_dataframe=na.omit(data.frame(CHR=NA,BEGIN=NA,END=NA,stringsAsFactors = F))


##按照每个染色体interval
for(rowindex in c(1:nrow(data))){
  CHR=data[rowindex,"CHR"]
  CUTOFF=as.numeric(c(0,unlist(stringi::stri_split(str = data[rowindex,"CUTOFF"],regex = "\\|"))))
  for(index in c(2:length(CUTOFF))){
      intervals=c(intervals,paste(CHR,CUTOFF[index-1]+1,CUTOFF[index],sep = "-"))
      intervals_dataframe[nrow(intervals_dataframe)+1,c("CHR","BEGIN","END")]=c(CHR,CUTOFF[index-1]+1,CUTOFF[index])
  }
}
clipr::write_clip(content = paste(intervals,sep = "\n"))

##写入文件
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
