source("E:/华大/表型/xiaomin/曼哈顿图绘制/qqman.r")
setwd("E:/华大/表型/xiaomin/pheno_new_normal")
(files=list.files(path = getwd(),pattern = ".*logistic$",full.names = F))

dataframe_files=data.frame(files=files,stringsAsFactors = FALSE)
dataframe_files$keyword=stringi::stri_extract(str = dataframe_files$files,regex = "[a-zA-Z0-9]+")
data_gwas_catalouge=as.data.frame(readxl::read_excel(path = "D:/gwas/all_pha2/3966.disease.snp.xlsx",sheet = "Sheet1"))
data_gwas_catalouge[,"CHR_BP"]=paste(data_gwas_catalouge$Chr_id,data_gwas_catalouge$Chr_pos,sep = "_")

(diseases_gwas_catalouge=unlist(stringi::stri_split(str = unique(unlist(data_gwas_catalouge$NAME_FILE_LOCAL)),regex = "\\|")))
##读取gwas_catalogue的文档
for(rowindex in c(1:nrow(dataframe_files))){
  keyword=dataframe_files[rowindex,"keyword"]
  if(!keyword %in% diseases_gwas_catalouge){
    message(paste(keyword,"不存在已知位点！"))
    next
  }
    file_data=dataframe_files[rowindex,"files"]
    data=as.data.frame(readr::read_delim(file = file_data,delim = "\t"))
    
    ##提取SNP位基因名称
    data1=apply(X = data,MARGIN = 1,FUN = function(x){
      x[["LABEL"]]=stringi::stri_split(str = x[["SNP"]],regex = ":")[[1]][6]
      x[["LABEL"]]=stringi::stri_replace_all(str = x[["LABEL"]],replacement ="" ,regex ="\\(.*?\\)" )
      x[["LABEL"]]=stringi::stri_replace_all(str = x[["LABEL"]],replacement ="" ,regex ="\\(.*" )
      x
    })
    data=as.data.frame(t(data1))
    
  ##
    data[,"CHR_BP"]=paste(data$CHR,data$BP,sep = "_")
  ##筛选GWAS CATALOGUE存在的
    data_gwas_catalouge_selected=data_gwas_catalouge[which(grepl(pattern = keyword,x = data_gwas_catalouge$NAME_FILE_LOCAL)),c("CHR_BP")]  
  
  ####<<<<<<<<<<<<<<<<<<<定义function检测是否在50kb范围内
  findtheNear = function(CHR_BP1,CHR_BP2){##CHR_BP2为已知的位点list,CHR_BP1为待验证位点的list
    listreturn=c()
    ##CHR_BP1转成DATAFARME
    CHR_BP1=data.frame(CHR_BP1=CHR_BP1,CHR1=NA,BP1=NA,isok=FALSE,stringsAsFactors = FALSE)
    CHR_BP1$CHR1=as.numeric(stringi::stri_extract(str = CHR_BP1$CHR_BP1,regex = "^\\d+"))
    CHR_BP1$BP1=as.numeric(stringi::stri_extract(str = CHR_BP1$CHR_BP1,regex = "\\d+$"))
    for(one in CHR_BP2){
      CHR2=as.numeric(stringi::stri_extract(str = one,regex = "^\\d+")) 
      BP2=as.numeric(stringi::stri_extract(str = one,regex = "\\d+$"))
      ##验证的化为TRUE
      CHR_BP1[intersect(which(CHR_BP1$CHR1==CHR2),which(abs(CHR_BP1$BP1-BP2)<=50000)),"isok"]=TRUE
    }
    ##返回值处理
    listreturn=unique(CHR_BP1$CHR_BP1[which(CHR_BP1$isok==TRUE)])
    return(listreturn)
  } 
  ####<<<<<<<<<<<<<<<<<<<定义function检测是否在50kb范围内
  for(colindex in c(1:ncol(data))){
    data[,colindex]=as.character(data[,colindex])
  }
    data[,"COLOR"]="black"
    data$CHR=as.numeric(as.character(data$CHR))
    data$BP=as.numeric(data$BP)
    data$P=as.numeric(as.character(data$P))
    data$CHR_BP=as.character(data$CHR_BP)
    data$LABEL=as.character(data$LABEL)
    data[which(data$CHR_BP %in% findtheNear(data$CHR_BP,data_gwas_catalouge_selected)),"COLOR"]="green3"
    ##绘图
    str(data)
    library(calibrate)
    setwd("E:/华大/表型/xiaomin/pheno_transform")
    tiff(filename = paste(file_data,".tiff",sep = ""),width = 2800,height = 1200,res = 200)
    data=data[which(!is.na(data$P)),]
    qqman_(x = data,chr = "CHR",bp = "BP",p = "P",snp = "CHR_BP",col = "COLOR",highlight = data[which(data$COLOR=="green3"),c("CHR_BP")],annotateTop = TRUE,annotatePval = 8)
    dev.off()
  
}
