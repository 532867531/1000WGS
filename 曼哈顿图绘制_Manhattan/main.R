source("E:/华大/表型/xiaomin/曼哈顿图绘制/qqman.r")
#qqman_(x = qqman::gwasResults,suggestiveline = 2,genomewideline = 10,highlight = T,annotatePval = 0.05)

##绘制BMI
rm(list = ls())
setwd("E:/华大/表型/xiaomin/pheno_transform")
source("E:/华大/表型/xiaomin/曼哈顿图绘制/qqman.r")
#data_bmi=as.data.frame(readxl::read_excel(path = "bmi.assoc.linear.sig.xlsx"))
data_bmi=as.data.frame(readxl::read_excel(path = "bmi.adjust.assoc.linear.sig.xlsx"))
data_bmi[,"CHR_BP"]=paste(data_bmi$CHR,data_bmi$BP,sep = "_")
##读取Meta-analysis_Locke_et_al+UKBiobank_2018_UPDATED.txt
#setwd("E:/华大/表型/summary-statistics-ukbiobank/GIANT_GWAS/BMI and Height GIANT and UK BioBank Meta-analysis Summary Statistics")
#data_uk_bmi=readr::read_delim(delim="\t",file = "Meta-analysis_Locke_et_al+UKBiobank_2018_UPDATED.txt")
data_uk_bmi=readr::read_delim(delim="\t",file = "E:/华大/表型/summary-statistics-ukbiobank/asian/BMI/29273807/BMI_Eastern_Asian.txt")
colnames(data_uk_bmi)[10]="P1"
data_uk_bmi2=readr::read_delim(delim="\t",file = "E:/华大/表型/summary-statistics-ukbiobank/asian/BMI/29273807/BMI_South_Asian.txt")
colnames(data_uk_bmi2)
colnames(data_uk_bmi2)[10]="P2"

head(data_uk_bmi)

colnames(data_uk_bmi)=paste(colnames(data_uk_bmi),"-")
colnames(data_uk_bmi2)=paste(colnames(data_uk_bmi2),"-")

data_uk_bmi[,"CHR_BP"]=paste(data_uk_bmi$`CHR -`,data_uk_bmi$`POS -`,sep = "_")
data_uk_bmi2[,"CHR_BP"]=paste(data_uk_bmi2$`CHR -`,data_uk_bmi2$`POS -`,sep = "_")
##merge
data_bmi_=merge(x = data_bmi,y=data_uk_bmi,by.x = "CHR_BP",by.y = "CHR_BP",all.x = T)
data_bmi_=merge(x = data_bmi_,y=data_uk_bmi2,by.x = "CHR_BP",by.y = "CHR_BP",all.x = T)
colnames(data_bmi_)

data_bmi_[,"COLOR"]="black"
data_bmi_[c(which(data_bmi_$`P1 -`<=0.01),which(data_bmi_$`P2 -`<=0.01)),"COLOR"]="red"
data_bmi_[c(which(data_bmi_$`P1 -`<=1e-05),which(data_bmi_$`P2 -`<=1e-05)),"COLOR"]="green"
##绘图
data_bmi_$SNP=stringi::stri_extract(str = data_bmi_$SNP,regex = "rs\\d+")
colnames(data_bmi_)
library(calibrate)
setwd("E:/华大/表型/xiaomin/pheno_transform")
#tiff(filename = "bmi.assoc.linear.sig.xlsx-BMI-South_North_Asian.tiff",width = 2800,height = 1200,res = 200)
tiff(filename = "bmi.adjust.assoc.linear.sig.xlsx-BMI-South_North_Asian.tiff",width = 2800,height = 1200,res = 200)
qqman_(x = data_bmi_,chr = "CHR",bp = "BP",p = "P",snp = "SNP",col = "COLOR",highlight = data_bmi_[which(data_bmi_$COLOR!="black"),c("SNP")],annotateTop = TRUE,annotatePval = 8)
dev.off()





##绘制Arthritis
rm(list = ls())
source("E:/华大/表型/xiaomin/曼哈顿图绘制/qqman.r")
setwd("E:/华大/表型/xiaomin/pheno_transform")
data_bmi=as.data.frame(readxl::read_excel(path = "arthritis.adjust.assoc.logistic.sig.xlsx"))
#data_bmi=as.data.frame(readxl::read_excel(path = ""))
data_bmi[,"CHR_BP"]=paste(data_bmi$CHR,data_bmi$BP,sep = "_")
##读取Meta-analysis_Locke_et_al+UKBiobank_2018_UPDATED.txt
data_uk_bmi=readr::read_delim(delim="\t",file = "E:/华大/表型/summary-statistics-ukbiobank/asian/arthritis/24390342/RA_GWASmeta_Asian_v2.txt")
head(data_uk_bmi)
data_uk_bmi=as.data.frame(data_uk_bmi)
colnames(data_uk_bmi)=paste(colnames(data_uk_bmi),"-")
data_uk_bmi[,"CHR_BP"]=paste(data_uk_bmi$`Chr -`,data_uk_bmi[,3],sep = "_")
##merge
data_bmi_=merge(x = data_bmi,y=data_uk_bmi,by.x = "CHR_BP",by.y = "CHR_BP",all.x = T)
data_bmi_[,"COLOR"]="black"
data_bmi_[which(data_bmi_$`P-val -`<=0.01),"COLOR"]="red"
data_bmi_[which(data_bmi_$`P-val -`<=1e-05),"COLOR"]="green"
##绘图
data_bmi_$SNP=stringi::stri_extract(str = data_bmi_$SNP,regex = "rs\\d+")
colnames(data_bmi_)
library(calibrate)
setwd("E:/华大/表型/xiaomin/pheno_transform")
tiff(filename = "arthritis.adjust.assoc.logistic.sig.xlsx-RA_GWASmeta_Asian_v2.txt.tiff",width = 2800,height = 1200,res = 200)
qqman_(x = data_bmi_,chr = "CHR",bp = "BP",p = "P",snp = "SNP",col = "COLOR",highlight = data_bmi_[which(data_bmi_$COLOR!="black"),c("SNP")],annotateTop = TRUE,annotatePval = 8)
dev.off()




##绘制glucoma
rm(list = ls())
setwd("E:/华大/表型/xiaomin/pheno_transform")
source("E:/华大/表型/xiaomin/曼哈顿图绘制/qqman.r")
data_bmi=as.data.frame(readxl::read_excel(path = "glaucoma.adjust.assoc.logistic.sig.xlsx"))
data_bmi[,"CHR_BP"]=paste(data_bmi$CHR,data_bmi$BP,sep = ":")
data_bmi[,"CHR_BP_rs"]=stringi::stri_extract(str = data_bmi$SNP,regex = "rs\\d+")
setwd("E:/华大/表型/summary-statistics-ukbiobank/asian/glaucoma")
list.files()
data_uk_bmi1=readr::read_delim(delim="\t",file = "Meta_CA_age_sex_asians_MAF0.01_20160807_newGC.1tbl")
data_uk_bmi1=as.data.frame(data_uk_bmi1)
data_uk_bmi1=data_uk_bmi1[,c("SNP","P-value")]
head(data_uk_bmi1)
colnames(data_uk_bmi1)=paste(colnames(data_uk_bmi1),"-")
##
data_uk_bmi2=readr::read_delim(delim="\t",file = "Meta_IOP_age_sex_asians_MAF0.01_20160806_newGC.1tbl")
data_uk_bmi2=as.data.frame(data_uk_bmi2)
head(data_uk_bmi2)
data_uk_bmi2=data_uk_bmi2[,c("SNP","P-value")]
colnames(data_uk_bmi2)=paste(colnames(data_uk_bmi2),"-")
##
data_uk_bmi3=readr::read_delim(delim="\t",file = "Meta_VCDR_age_sex_asians_MAF0.01_20160808_newGC.1tbl")
data_uk_bmi3=as.data.frame(data_uk_bmi3)
head(data_uk_bmi3)
data_uk_bmi3=data_uk_bmi3[,c("SNP","P-value")]
colnames(data_uk_bmi3)=paste(colnames(data_uk_bmi3),"-")
##
data_uk_bmi4=readr::read_delim(delim="\t",file = "Meta_CA_age_sex_asians_MAF0.01_20160807_newGC.1tbl")
data_uk_bmi4=as.data.frame(data_uk_bmi4)
head(data_uk_bmi4)
data_uk_bmi4=data_uk_bmi4[,c("SNP","P-value")]
colnames(data_uk_bmi4)=paste(colnames(data_uk_bmi4),"-")
##merge
data_bmi_=merge(x = data_bmi,y=data_uk_bmi1,by.x = "CHR_BP_rs",by.y = "SNP -",all.x = T)
data_bmi_=merge(x = data_bmi_,y=data_uk_bmi1,by.x = "CHR_BP",by.y = "SNP -",all.x = T)

data_bmi_=merge(x = data_bmi_,y=data_uk_bmi2,by.x = "CHR_BP_rs",by.y = "SNP -",all.x = T)
data_bmi_=merge(x = data_bmi_,y=data_uk_bmi2,by.x = "CHR_BP",by.y = "SNP -",all.x = T)

data_bmi_=merge(x = data_bmi_,y=data_uk_bmi3,by.x = "CHR_BP_rs",by.y = "SNP -",all.x = T)
data_bmi_=merge(x = data_bmi_,y=data_uk_bmi3,by.x = "CHR_BP",by.y = "SNP -",all.x = T)

data_bmi_=merge(x = data_bmi_,y=data_uk_bmi4,by.x = "CHR_BP_rs",by.y = "SNP -",all.x = T)
data_bmi_=merge(x = data_bmi_,y=data_uk_bmi4,by.x = "CHR_BP",by.y = "SNP -",all.x = T)

data_bmi_[,"COLOR"]="black"
data_bmi_[c(which(data_bmi_$`P-value -.x`<=0.01),which(data_bmi_$`P-value -.y`<=0.01),which(data_bmi_$`P-value -.x.1`<=0.01),which(data_bmi_$`P-value -.y.1`<=0.01),which(data_bmi_$`P-value -.x.2`<=0.01),which(data_bmi_$`P-value -.y.2`<=0.01),which(data_bmi_$`P-value -.x.3`<=0.01),which(data_bmi_$`P-value -.y.3`<=0.01)),"COLOR"]="red"
data_bmi_[c(which(data_bmi_$`P-value -.x`<=1e-05),which(data_bmi_$`P-value -.y`<=1e-05),which(data_bmi_$`P-value -.x.1`<=1e-05),which(data_bmi_$`P-value -.y.1`<=1e-05),which(data_bmi_$`P-value -.x.2`<=1e-05),which(data_bmi_$`P-value -.y.2`<=1e-05),which(data_bmi_$`P-value -.x.3`<=1e-05),which(data_bmi_$`P-value -.y.3`<=1e-05)),"COLOR"]="green"
##绘图
colnames(data_bmi_)
library(calibrate)
setwd("E:/华大/表型/xiaomin/pheno_transform")
data_bmi_$SNP=stringi::stri_extract(str = data_bmi_$SNP,regex = "rs\\d+")
data_bmi_=data_bmi_[which(!is.na(data_bmi_$P)),]
tiff(filename = "glaucoma.adjust.assoc.logistic.sig.xlsx.tiff",width = 2800,height = 1200,res = 200)
qqman_(x = data_bmi_,chr = "CHR",bp = "BP",p = "P",snp = "SNP",col = "COLOR",highlight = data_bmi_[which(data_bmi_$COLOR!="black"),c("SNP")],annotateTop = TRUE,annotatePval = 8)
dev.off()
shell.exec(getwd())








##绘制tubercolusis???
rm(list = ls())
source("E:/华大/表型/xiaomin/曼哈顿图绘制/qqman.r")
setwd("E:/华大/表型/xiaomin/pheno_transform")
data_bmi=as.data.frame(readxl::read_excel(path = "tuber.adjust.assoc.logistic.sig.xlsx"))
data_bmi[,"CHR_BP"]=paste(data_bmi$CHR,data_bmi$BP,sep = "_")
##读取Meta-analysis_Locke_et_al+UKBiobank_2018_UPDATED.txt
data_uk_bmi=readr::read_delim(delim="\t",file = "E:/华大/表型/summary-statistics-ukbiobank/asian/tuberculosis/ZhengR_30287856_GCST006619/tb_summary.txt")
data_uk_bmi=as.data.frame(data_uk_bmi)
head(data_uk_bmi)
colnames(data_uk_bmi)=paste(colnames(data_uk_bmi),"-")
data_uk_bmi[,"CHR_BP"]=paste(data_uk_bmi$`chr -`,data_uk_bmi$`pos -`,sep = "_")
##merge
data_bmi_=merge(x = data_bmi,y=data_uk_bmi,by.x = "CHR_BP",by.y = "CHR_BP",all.x = T)
data_bmi_[,"COLOR"]="black"
data_bmi_[which(data_bmi_$`pvalue -`<=0.01),"COLOR"]="red"
data_bmi_[which(data_bmi_$`pvalue -`<=1e-05),"COLOR"]="green"
##绘图
colnames(data_bmi_)
library(calibrate)
setwd("E:/华大/表型/xiaomin/pheno_transform")
data_bmi_$SNP=stringi::stri_extract(str = data_bmi_$SNP,regex = "rs\\d+")
tiff(filename = "tuber.adjust.assoc.logistic.sig.xlsx.tiff",width = 2800,height = 1200,res = 200)
data_bmi_=data_bmi_[which(!is.na(data_bmi_$P)),]
qqman_(x = data_bmi_,chr = "CHR",bp = "BP",p = "P",snp = "CHR_BP",col = "COLOR",highlight = data_bmi_[which(data_bmi_$COLOR!="black"),c("CHR_BP")],annotateTop = TRUE,annotatePval = 8)
#qqman_(x = data_bmi_,chr = "CHR",bp = "BP",p = "P",snp = "SNP",col = "COLOR",annotateTop = TRUE,annotatePval = 8)
dev.off()
shell.exec(getwd())






##绘制PD
rm(list = ls())
setwd("E:/华大/表型/xiaomin/pheno_transform")
source("E:/华大/表型/xiaomin/曼哈顿图绘制/qqman.r")
data_bmi=as.data.frame(readxl::read_excel(path = "parkinson.assoc.logistic.sig.xlsx"))
data_bmi=as.data.frame(readxl::read_excel(path = "parkinson.adjust.assoc.logistic.sig.xlsx"))
data_bmi=as.data.frame(readxl::read_excel(path = "parkinson.adjust.assoc.linear.sig.xlsx"))

data_bmi[,"CHR_BP"]=paste(data_bmi$CHR,data_bmi$BP,sep = "_")
##读取Meta-analysis_Locke_et_al+UKBiobank_2018_UPDATED.txt
data_uk_bmi=readr::read_delim(delim="\t",file = "E:/华大/表型/summary-statistics-ukbiobank/PD/PD_meta_analysis_PDGENE_PDWBS/PD_meta_analysis_PDGENE_PDWBS.txt")
head(data_uk_bmi)
data_uk_bmi=as.data.frame(data_uk_bmi)
colnames(data_uk_bmi)=paste(colnames(data_uk_bmi),"-")
data_uk_bmi[,"CHR_BP"]=paste(data_uk_bmi$`CHR -`,data_uk_bmi$`BP -`,sep = "_")
##merge
data_bmi_=merge(x = data_bmi,y=data_uk_bmi,by.x = "CHR_BP",by.y = "CHR_BP",all.x = T)
data_bmi_[,"COLOR"]="black"
data_bmi_[which(data_bmi_$`P.META -`<=0.01),"COLOR"]="red"
data_bmi_[which(data_bmi_$`P.META -`<=1e-05),"COLOR"]="pink"
##绘图
colnames(data_bmi_)
library(calibrate)
setwd("E:/华大/表型/xiaomin/pheno_transform")
tiff(filename = "parkinson.assoc.logistic.sig.tiff",width = 2800,height = 1200,res = 200)
tiff(filename = "parkinson.adjust.assoc.logistic.sig.tiff",width = 2800,height = 1200,res = 200)
tiff(filename = "parkinson.adjust.assoc.linear.sig.tiff",width = 2800,height = 1200,res = 200)
  head(data_bmi_)
data_bmi_$SNP=stringi::stri_extract(str = data_bmi_$SNP,regex = "rs\\d+")
data_bmi_=data_bmi_[which(!is.na(data_bmi_$P)),]
qqman_(x = data_bmi_,chr = "CHR",bp = "BP",p = "P",snp = "SNP",col = "COLOR",highlight = data_bmi_[which(data_bmi_$COLOR!="black"),c("SNP")],annotateTop = TRUE,annotatePval = 8)
dev.off()
shell.exec(getwd())





##绘制Stroke
rm(list = ls())
source("E:/华大/表型/xiaomin/曼哈顿图绘制/qqman.r")
setwd("E:/华大/表型/xiaomin/pheno_transform")
data_bmi=as.data.frame(readxl::read_excel(path = "stroke.adjust.assoc.logistic.sig.xlsx"))

data_bmi[,"CHR_BP_rs"]=stringi::stri_extract(str = data_bmi$SNP,regex = "rs\\d+")
data_bmi[,"SNP"]=stringi::stri_extract(str = data_bmi$SNP,regex = "rs\\d+")

data_bmi[,"CHR_BP"]=paste(data_bmi$CHR,data_bmi$BP,sep = "_")
##读取Meta-analysis_Locke_et_al+UKBiobank_2018_UPDATED.txt
data_uk_bmi1=readr::read_delim(delim="\t",file = "E:/华大/表型/summary-statistics-ukbiobank/4818561.Malik.2016/METAANALYSIS1_CE.TBL")
data_uk_bmi1=as.data.frame(data_uk_bmi1)
data_uk_bmi1=data_uk_bmi1[,c("MarkerName","P-value")]
head(data_uk_bmi1)
colnames(data_uk_bmi1)=paste(colnames(data_uk_bmi1),"1")
##
data_uk_bmi2=readr::read_delim(delim="\t",file = "E:/华大/表型/summary-statistics-ukbiobank/4818561.Malik.2016/METAANALYSIS1_IS.TBL")
head(data_uk_bmi2)
data_uk_bmi2=as.data.frame(data_uk_bmi2)
data_uk_bmi2=data_uk_bmi2[,c("MarkerName","P-value")]
colnames(data_uk_bmi2)=paste(colnames(data_uk_bmi2),"2")
##
data_uk_bmi3=readr::read_delim(delim="\t",file = "E:/华大/表型/summary-statistics-ukbiobank/4818561.Malik.2016/METAANALYSIS1_LVD.TBL")
head(data_uk_bmi3)
data_uk_bmi3=as.data.frame(data_uk_bmi3)
data_uk_bmi3=data_uk_bmi3[,c("MarkerName","P-value")]
colnames(data_uk_bmi3)=paste(colnames(data_uk_bmi3),"3")
##
data_uk_bmi4=readr::read_delim(delim="\t",file = "E:/华大/表型/summary-statistics-ukbiobank/4818561.Malik.2016/METAANALYSIS1_SVD.TBL")
head(data_uk_bmi4)
data_uk_bmi4=as.data.frame(data_uk_bmi4)
data_uk_bmi4=data_uk_bmi4[,c("MarkerName","P-value")]
colnames(data_uk_bmi4)=paste(colnames(data_uk_bmi4),"4")
 ##merge
data_bmi_=merge(x = data_bmi,y=data_uk_bmi1,by.x = "CHR_BP_rs",by.y = "MarkerName 1",all.x = T)
#data_bmi_=merge(x = data_bmi,y=data_uk_bmi1,by.x = "CHR_BP",by.y = "MarkerName",all.x = T)

data_bmi_=merge(x = data_bmi_,y=data_uk_bmi2,by.x = "CHR_BP_rs",by.y = "MarkerName 2",all.x = T)
#data_bmi_=merge(x = data_bmi_,y=data_uk_bmi2,by.x = "CHR_BP",by.y = "MarkerName",all.x = T)

data_bmi_=merge(x = data_bmi_,y=data_uk_bmi3,by.x = "CHR_BP_rs",by.y = "MarkerName 3",all.x = T)
#data_bmi_=merge(x = data_bmi_,y=data_uk_bmi3,by.x = "CHR_BP",by.y = "MarkerName",all.x = T)

data_bmi_=merge(x = data_bmi_,y=data_uk_bmi4,by.x = "CHR_BP_rs",by.y = "MarkerName 4",all.x = T)
#data_bmi_=merge(x = data_bmi_,y=data_uk_bmi4,by.x = "CHR_BP",by.y = "MarkerName",all.x = T)##绘图
colnames(data_bmi_)

data_bmi_[,"COLOR"]="black"
data_bmi_[c(which(data_bmi_$`P-value 1`<=0.01),which(data_bmi_$`P-value 2`<=0.01),which(data_bmi_$`P-value 3`<=0.01),which(data_bmi_$`P-value 4`<=0.01),which(data_bmi_$`P-value -.x.2`<=0.01),which(data_bmi_$`P-value -.y.2`<=0.01),which(data_bmi_$`P-value -.x.3`<=0.01),which(data_bmi_$`P-value -.y.3`<=0.01)),"COLOR"]="red"
data_bmi_[c(which(data_bmi_$`P-value 1`<=1e-05),which(data_bmi_$`P-value 2`<=1e-05),which(data_bmi_$`P-value 3`<=1e-05),which(data_bmi_$`P-value 4`<=1e-05),which(data_bmi_$`P-value -.x.2`<=1e-05),which(data_bmi_$`P-value -.y.2`<=1e-05),which(data_bmi_$`P-value -.x.3`<=1e-05),which(data_bmi_$`P-value -.y.3`<=1e-05)),"COLOR"]="green"

library(calibrate)
setwd("E:/华大/表型/xiaomin/pheno_transform")
tiff(filename = "stroke.adjust.assoc.logistic.sig.xlsx-Malik.2016.tiff",width = 2800,height = 1200,res = 200)
head(data_bmi_)
data_bmi_$SNP=stringi::stri_extract(str = data_bmi_$SNP,regex = "rs\\d+")
data_bmi_=data_bmi_[which(!is.na(data_bmi_$P)),]
qqman_(x = data_bmi_,chr = "CHR",bp = "BP",p = "P",snp = "SNP",col = "COLOR",highlight = data_bmi_[which(data_bmi_$COLOR!="black"),c("SNP")],annotateTop = TRUE,annotatePval = 8)
dev.off()
shell.exec(getwd())









##绘制Heart
rm(list = ls())
source("E:/华大/表型/xiaomin/曼哈顿图绘制/qqman.r")
setwd("E:/华大/表型/xiaomin/pheno_transform")
data_bmi=as.data.frame(readxl::read_excel(path = "heart.adjust.assoc.logistic.sig.xlsx"))

data_bmi[,"CHR_BP_rs"]=stringi::stri_extract(str = data_bmi$SNP,regex = "rs\\d+")
data_bmi[,"SNP"]=stringi::stri_extract(str = data_bmi$SNP,regex = "rs\\d+")

data_bmi[,"CHR_BP"]=paste(data_bmi$CHR,data_bmi$BP,sep = "_")
##读取Meta-analysis_Locke_et_al+UKBiobank_2018_UPDATED.txt
data_uk_bmi1=readr::read_delim(delim="\t",file = "E:/华大/表型/summary-statistics-ukbiobank/Heart failure GWAS/2018.HRC.GWAS.UKBB/HF_HRC_GWAS_UKBB_EUR.txt")
data_uk_bmi1=as.data.frame(data_uk_bmi1)
head(data_uk_bmi1)
colnames(data_uk_bmi1)=paste(colnames(data_uk_bmi1),"1")
data_uk_bmi1[,"CHR_BP"]=paste(data_uk_bmi1$`chr 1`,data_uk_bmi1$`pos 1`,sep = "_")
##
data_uk_bmi2=readr::read_delim(delim="\t",file = "E:/华大/表型/summary-statistics-ukbiobank/Heart failure GWAS/2018.HRC.GWAS.UKBB/NICM_HRC_GWAS_UKBB_EUR.txt")
data_uk_bmi2=as.data.frame(data_uk_bmi2)
head(data_uk_bmi2)
colnames(data_uk_bmi2)=paste(colnames(data_uk_bmi2),"2")
data_uk_bmi2[,"CHR_BP"]=paste(data_uk_bmi2$`chr 2`,data_uk_bmi2$`pos 2`,sep = "_")
##
##merge
data_bmi_=merge(x = data_bmi,y=data_uk_bmi1,by.x = "CHR_BP",by.y = "CHR_BP",all.x = T)
#data_bmi_=merge(x = data_bmi,y=data_uk_bmi1,by.x = "CHR_BP",by.y = "MarkerName",all.x = T)

data_bmi_=merge(x = data_bmi_,y=data_uk_bmi2,by.x = "CHR_BP",by.y = "CHR_BP",all.x = T)
#data_bmi_=merge(x = data_bmi_,y=data_uk_bmi2,by.x = "CHR_BP",by.y = "MarkerName",all.x = T)


data_bmi_[,"COLOR"]="black"
data_bmi_[c(which(data_bmi_$`p.value 1`<=0.01),which(data_bmi_$`p.value 2`<=0.01)),"COLOR"]="red"
data_bmi_[c(which(data_bmi_$`p.value 1`<=1e-05),which(data_bmi_$`p.value 2`<=1e-05)),"COLOR"]="green"

library(calibrate)
setwd("E:/华大/表型/xiaomin/pheno_transform")
tiff(filename = "heart.adjust.assoc.logistic.sig.xlsx.2018.HRC.tiff",width = 2800,height = 1200,res = 200)
data_bmi_$SNP=stringi::stri_extract(str = data_bmi_$SNP,regex = "rs\\d+")
head(data_bmi_)
data_bmi_=data_bmi_[which(!is.na(data_bmi_$P)),]
qqman_(x = data_bmi_,chr = "CHR",bp = "BP",p = "P",snp = "SNP",col = "COLOR",highlight = data_bmi_[which(data_bmi_$COLOR!="black"),c("SNP")],annotateTop = TRUE,annotatePval = 8)
dev.off()
shell.exec(getwd())



##绘制cataract
rm(list = ls())
source("E:/华大/表型/xiaomin/曼哈顿图绘制/qqman.r")
setwd("E:/华大/表型/xiaomin/pheno_transform")
data_bmi=as.data.frame(readxl::read_excel(path = "cataract.adjust.assoc.logistic.sig.xlsx"))
data_bmi[,"CHR_BP"]=paste(data_bmi$CHR,data_bmi$BP,sep = "_")
##读取Meta-analysis_Locke_et_al+UKBiobank_2018_UPDATED.txt
data_uk_bmi=readr::read_delim(delim="\t",file = "E:/华大/表型/summary-statistics-ukbiobank/Cataract-6148_4_logistic.EUR.sumstats.MACfilt.txt/6148_4_logistic.EUR.sumstats.MACfilt.txt")
head(data_uk_bmi)
data_uk_bmi=as.data.frame(data_uk_bmi)
colnames(data_uk_bmi)=paste(colnames(data_uk_bmi),"-")
data_uk_bmi[,"CHR_BP"]=paste(data_uk_bmi$`CHR -`,data_uk_bmi$`BP -`,sep = "_")
##merge
data_bmi_=merge(x = data_bmi,y=data_uk_bmi,by.x = "CHR_BP",by.y = "CHR_BP",all.x = T)
data_bmi_[,"COLOR"]="black"
data_bmi_[which(data_bmi_$`P -`<=0.01),"COLOR"]="red"
data_bmi_[which(data_bmi_$`P -`<=1e-05),"COLOR"]="green"
##绘图
colnames(data_bmi_)
library(calibrate)
setwd("E:/华大/表型/xiaomin/pheno_transform")
data_bmi_$SNP=stringi::stri_extract(str = data_bmi_$SNP,regex = "rs\\d+")
tiff(filename = "cataract.adjust.assoc.logistic.sig.xlsx.tiff",width = 2800,height = 1200,res = 200)
data_bmi_=data_bmi_[which(!is.na(data_bmi_$P)),]
qqman_(x = data_bmi_,chr = "CHR",bp = "BP",p = "P",snp = "SNP",col = "COLOR",highlight = data_bmi_[which(data_bmi_$COLOR!="black"),c("SNP")],annotateTop = TRUE,annotatePval = 8)
dev.off()
shell.exec(getwd())


##绘制epilepsy
rm(list = ls())
source("E:/华大/表型/xiaomin/曼哈顿图绘制/qqman.r")
setwd("E:/华大/表型/xiaomin/pheno_transform")
data_bmi=as.data.frame(readxl::read_excel(path = "epilepsy.adjust.assoc.logistic.sig.xlsx"))
data_bmi[,"CHR_BP"]=paste(data_bmi$CHR,data_bmi$BP,sep = "_")
##读取Meta-analysis_Locke_et_al+UKBiobank_2018_UPDATED.txt
data_uk_bmi=readr::read_delim(delim="\t",file = "E:/华大/表型/summary-statistics-ukbiobank/epilepsy/all_epilepsy_METAL/all_epilepsy_METAL")
head(data_uk_bmi)
data_uk_bmi=as.data.frame(data_uk_bmi)
colnames(data_uk_bmi)=paste(colnames(data_uk_bmi),"-")
data_uk_bmi[,"CHR_BP"]=paste(data_uk_bmi$`CHR -`,data_uk_bmi$`BP -`,sep = "_")
##merge
data_bmi_=merge(x = data_bmi,y=data_uk_bmi,by.x = "CHR_BP",by.y = "CHR_BP",all.x = T)
data_bmi_[,"COLOR"]="black"
data_bmi_[which(data_bmi_$`P-value -`<=0.01),"COLOR"]="red"
data_bmi_[which(data_bmi_$`P-value -`<=1e-05),"COLOR"]="green"
##绘图
colnames(data_bmi_)
library(calibrate)
setwd("E:/华大/表型/xiaomin/pheno_transform")
data_bmi_$SNP=stringi::stri_extract(str = data_bmi_$SNP,regex = "rs\\d+")
tiff(filename = "epilepsy.adjust.assoc.logistic.sig.xlsx.tiff",width = 2800,height = 1200,res = 200)
data_bmi_=data_bmi_[which(!is.na(data_bmi_$P)),]
qqman_(x = data_bmi_,chr = "CHR",bp = "BP",p = "P",snp = "SNP",col = "COLOR",highlight = data_bmi_[which(data_bmi_$COLOR!="black"),c("SNP")],annotateTop = TRUE,annotatePval = 8)
dev.off()
shell.exec(getwd())



##绘制T2D
rm(list = ls())
source("E:/华大/表型/xiaomin/曼哈顿图绘制/qqman.r")
setwd("E:/华大/表型/xiaomin/pheno_transform")
data_bmi=as.data.frame(readxl::read_excel(path = "diabete.adjust.assoc.logistic.sig.xlsx"))
data_bmi[,"CHR_BP"]=paste(data_bmi$CHR,data_bmi$BP,sep = "_")
##读取Meta-analysis_Locke_et_al+UKBiobank_2018_UPDATED.txt
data_uk_bmi=readr::read_delim(delim="\t",file = "E:/华大/表型/summary-statistics-ukbiobank/UKBB/disease_T2D.sumstats/disease_T2D.sumstats")
head(data_uk_bmi)
data_uk_bmi=as.data.frame(data_uk_bmi)
colnames(data_uk_bmi)=paste(colnames(data_uk_bmi),"-")
data_uk_bmi[,"CHR_BP"]=paste(data_uk_bmi$`CHR -`,data_uk_bmi$`POS -`,sep = "_")
##merge
data_bmi_=merge(x = data_bmi,y=data_uk_bmi,by.x = "CHR_BP",by.y = "CHR_BP",all.x = T)
data_bmi_[,"COLOR"]="black"
data_bmi_[which(data_bmi_$`P -`<=0.01),"COLOR"]="red"
data_bmi_[which(data_bmi_$`P -`<=1e-05),"COLOR"]="green"
##绘图
colnames(data_bmi_)
library(calibrate)
setwd("E:/华大/表型/xiaomin/pheno_transform")
data_bmi_$SNP=stringi::stri_extract(str = data_bmi_$SNP,regex = "rs\\d+")
tiff(filename = "diabete.adjust.assoc.logistic.sig.xlsx-disease_T2D.sumstats.tiff",width = 2800,height = 1200,res = 200)
data_bmi_=data_bmi_[which(!is.na(data_bmi_$P)),]
qqman_(x = data_bmi_,chr = "CHR",bp = "BP",p = "P",snp = "SNP",col = "COLOR",highlight = data_bmi_[which(data_bmi_$COLOR!="black"),c("SNP")],annotateTop = TRUE,annotatePval = 8)
dev.off()
shell.exec(getwd())






##绘制Hypertension
rm(list = ls())
source("E:/华大/表型/xiaomin/曼哈顿图绘制/qqman.r")
setwd("E:/华大/表型/xiaomin/pheno_transform")
data_bmi=as.data.frame(readxl::read_excel(path = "hyper_r.adjust.assoc.logistic.sig.xlsx"))
#data_bmi=as.data.frame(readxl::read_excel(path = "hyper_r.adjust2.assoc.logistic.sig.xlsx"))
data_bmi[,"CHR_BP"]=paste(data_bmi$CHR,data_bmi$BP,sep = "_")
##读取Meta-analysis_Locke_et_al+UKBiobank_2018_UPDATED.txt
data_uk_bmi=readr::read_delim(delim="\t",file = "E:/华大/表型/summary-statistics-ukbiobank/UKBB/hypertension/41204_I10_logistic.EUR.sumstats.MACfilt.txt/41204_I10_logistic.EUR.sumstats.MACfilt.txt")
head(data_uk_bmi)
data_uk_bmi=as.data.frame(data_uk_bmi)

colnames(data_uk_bmi)=paste(colnames(data_uk_bmi),"-")
data_uk_bmi[,"CHR_BP"]=paste(data_uk_bmi$`CHR -`,data_uk_bmi$`BP -`,sep = "_")
##merge
data_bmi_=merge(x = data_bmi,y=data_uk_bmi,by.x = "CHR_BP",by.y = "CHR_BP",all.x = T)
data_bmi_[,"COLOR"]="black"
data_bmi_[which(data_bmi_$`P -`<=0.01),"COLOR"]="red"
data_bmi_[which(data_bmi_$`P -`<=1e-05),"COLOR"]="green"
##绘图
colnames(data_bmi_)
library(calibrate)
setwd("E:/华大/表型/xiaomin/pheno_transform")
data_bmi_$SNP=stringi::stri_extract(str = data_bmi_$SNP,regex = "rs\\d+")
tiff(filename = "hyper_r.adjust.assoc.logistic.sig.xlsx-41204_I10_logistic.tiff",width = 2800,height = 1200,res = 200)
#tiff(filename = "hyper_r.adjust2.assoc.logistic.sig.xlsx-41204_I10_logistic.tiff",width = 2800,height = 1200,res = 200)
data_bmi_=data_bmi_[which(!is.na(data_bmi_$P)),]
qqman_(x = data_bmi_,chr = "CHR",bp = "BP",p = "P",snp = "SNP",col = "COLOR",highlight = data_bmi_[which(data_bmi_$COLOR!="black"),c("SNP")],annotateTop = TRUE,annotatePval = 8)
dev.off()
shell.exec(getwd())


##绘制Cholelith
rm(list = ls())
source("E:/华大/表型/xiaomin/曼哈顿图绘制/qqman.r")
setwd("E:/华大/表型/xiaomin/pheno_transform")
data_bmi=as.data.frame(readxl::read_excel(path = "cholelith.adjust.assoc.logistic.sig.xlsx"))
data_bmi[,"CHR_BP"]=paste(data_bmi$CHR,data_bmi$BP,sep = "_")
##读取Meta-analysis_Locke_et_al+UKBiobank_2018_UPDATED.txt
data_uk_bmi=readr::read_delim(delim="\t",file = "E:/华大/表型/summary-statistics-ukbiobank/UKBB/Cholelithiasis_41202_K80_logistic.EUR.sumstats.MACfilt.txt/41202_K80_logistic.EUR.sumstats.MACfilt.txt")
data_uk_bmi=as.data.frame(data_uk_bmi)
head(data_uk_bmi)


colnames(data_uk_bmi)=paste(colnames(data_uk_bmi),"-")
data_uk_bmi[,"CHR_BP"]=paste(data_uk_bmi$`CHR -`,data_uk_bmi$`BP -`,sep = "_")
##merge
data_bmi_=merge(x = data_bmi,y=data_uk_bmi,by.x = "CHR_BP",by.y = "CHR_BP",all.x = T)
data_bmi_[,"COLOR"]="black"
data_bmi_[which(data_bmi_$`P -`<=0.01),"COLOR"]="red"
data_bmi_[which(data_bmi_$`P -`<=1e-05),"COLOR"]="green"
##绘图
colnames(data_bmi_)
library(calibrate)
setwd("E:/华大/表型/xiaomin/pheno_transform")
data_bmi_$SNP=stringi::stri_extract(str = data_bmi_$SNP,regex = "rs\\d+")
tiff(filename = "cholelith.adjust.assoc.logistic.sig.xlsx-Cholelithiasis_41202_K80.tiff",width = 2800,height = 1200,res = 200)
data_bmi_=data_bmi_[which(!is.na(data_bmi_$P)),]
qqman_(x = data_bmi_,chr = "CHR",bp = "BP",p = "P",snp = "SNP",col = "COLOR",highlight = data_bmi_[which(data_bmi_$COLOR!="black"),c("SNP")],annotateTop = TRUE,annotatePval = 8)
dev.off()
shell.exec(getwd())





##绘制respiratory-asthma
rm(list = ls())
source("E:/华大/表型/xiaomin/曼哈顿图绘制/qqman.r")
setwd("E:/华大/表型/xiaomin/pheno_transform")
data_bmi=as.data.frame(readxl::read_excel(path = "respiratory.adjust.assoc.logistic.sig.xlsx"))
data_bmi[,"CHR_BP"]=paste(data_bmi$CHR,data_bmi$BP,sep = "_")
##读取Meta-analysis_Locke_et_al+UKBiobank_2018_UPDATED.txt
data_uk_bmi=readr::read_delim(delim="\t",file = "E:/华大/表型/summary-statistics-ukbiobank/UKBB/UKBB_409K/disease_ASTHMA_DIAGNOSED.sumstats/disease_ASTHMA_DIAGNOSED.sumstats")
data_uk_bmi=as.data.frame(data_uk_bmi)
head(data_uk_bmi)

colnames(data_uk_bmi)=paste(colnames(data_uk_bmi),"-")
data_uk_bmi[,"CHR_BP"]=paste(data_uk_bmi$`CHR -`,data_uk_bmi$`POS -`,sep = "_")
##merge
data_bmi_=merge(x = data_bmi,y=data_uk_bmi,by.x = "CHR_BP",by.y = "CHR_BP",all.x = T)
data_bmi_[,"COLOR"]="black"
data_bmi_[which(data_bmi_$`P -`<=0.01),"COLOR"]="red"
data_bmi_[which(data_bmi_$`P -`<=1e-05),"COLOR"]="green"
##绘图
colnames(data_bmi_)
library(calibrate)
setwd("E:/华大/表型/xiaomin/pheno_transform")
data_bmi_$SNP=stringi::stri_extract(str = data_bmi_$SNP,regex = "rs\\d+")
tiff(filename = "respiratory.adjust.assoc.logistic.sig.xlsx-disease_ASTHMA_DIAGNOSED.sumstats.tiff",width = 2800,height = 1200,res = 200)
data_bmi_=data_bmi_[which(!is.na(data_bmi_$P)),]
qqman_(x = data_bmi_,chr = "CHR",bp = "BP",p = "P",snp = "SNP",col = "COLOR",highlight = data_bmi_[which(data_bmi_$COLOR!="black"),c("SNP")],annotateTop = TRUE,annotatePval = 8)
dev.off()
shell.exec(getwd())



##绘制respiratory-asthma2
rm(list = ls())
source("E:/华大/表型/xiaomin/曼哈顿图绘制/qqman.r")
setwd("E:/华大/表型/xiaomin/pheno_transform")
data_bmi=as.data.frame(readxl::read_excel(path = "respiratory.adjust.assoc.logistic.sig.xlsx"))
data_bmi[,"CHR_BP"]=paste(data_bmi$CHR,data_bmi$BP,sep = "_")
##读取Meta-analysis_Locke_et_al+UKBiobank_2018_UPDATED.txt
#data_uk_bmi=readr::read_delim(delim="\t",file = "E:/华大/表型/summary-statistics-ukbiobank/UKBB/respiratory/Asthma-41204_J45_logistic.EUR.sumstats.MACfilt.txt/41204_J45_logistic.EUR.sumstats.MACfilt.txt")
data_uk_bmi=readr::read_delim(delim="\t",file = "E:/华大/表型/summary-statistics-ukbiobank/UKBB/respiratory/Asthma-41204_J45_logistic.EUR.sumstats.MACfilt.txt/41204_J45_logistic.EUR.sumstats.MACfilt.txt")
data_uk_bmi=as.data.frame(data_uk_bmi)
head(data_uk_bmi)

colnames(data_uk_bmi)=paste(colnames(data_uk_bmi),"-")
data_uk_bmi[,"CHR_BP"]=paste(data_uk_bmi$`CHR -`,data_uk_bmi$`BP -`,sep = "_")
##merge
data_bmi_=merge(x = data_bmi,y=data_uk_bmi,by.x = "CHR_BP",by.y = "CHR_BP",all.x = T)
data_bmi_[,"COLOR"]="black"
data_bmi_[which(data_bmi_$`P -`<=0.01),"COLOR"]="red"
data_bmi_[which(data_bmi_$`P -`<=1e-05),"COLOR"]="green"
##绘图
colnames(data_bmi_)
library(calibrate)
setwd("E:/华大/表型/xiaomin/pheno_transform")
data_bmi_$SNP=stringi::stri_extract(str = data_bmi_$SNP,regex = "rs\\d+")
tiff(filename = "respiratory.adjust.assoc.logistic.sig.xlsx-Asthma-41204_J45.tiff",width = 2800,height = 1200,res = 200)
data_bmi_=data_bmi_[which(!is.na(data_bmi_$P)),]
qqman_(x = data_bmi_,chr = "CHR",bp = "BP",p = "P",snp = "SNP",col = "COLOR",highlight = data_bmi_[which(data_bmi_$COLOR!="black"),c("SNP")],annotateTop = TRUE,annotatePval = 8)
dev.off()
shell.exec(getwd())



##绘制respiratory-asthma-TCGA
rm(list = ls())
source("E:/华大/表型/xiaomin/曼哈顿图绘制/qqman.r")
setwd("E:/华大/表型/xiaomin/pheno_transform")
data_bmi=as.data.frame(readxl::read_excel(path = "respiratory.adjust.assoc.logistic.sig.xlsx"))
data_bmi[,"CHR_BP"]=paste(data_bmi$CHR,data_bmi$BP,sep = "_")
##读取Meta-analysis_Locke_et_al+UKBiobank_2018_UPDATED.txt
data_uk_bmi=readr::read_delim(delim="\t",file = "E:/华大/表型/summary-statistics-ukbiobank/UKBB/respiratory/TAGC_meta-analyses_results_for_asthma_risk/TAGC meta-analyses results for asthma risk/TAGC_Multiancestry_and_European-Ancestry_Meta-analyses_Results.tsv")
data_uk_bmi=as.data.frame(data_uk_bmi)
head(data_uk_bmi)

colnames(data_uk_bmi)=paste(colnames(data_uk_bmi),"-")
data_uk_bmi[,"CHR_BP"]=paste(data_uk_bmi$`chr -`,data_uk_bmi$`position -`,sep = "_")
##merge
data_bmi_=merge(x = data_bmi,y=data_uk_bmi,by.x = "CHR_BP",by.y = "CHR_BP",all.x = T)
data_bmi_[,"COLOR"]="black"
data_bmi_[which(data_bmi_$`Multiancestry_pval_fix -`<=0.01),"COLOR"]="red"
data_bmi_[which(data_bmi_$`Multiancestry_pval_fix -`<=1e-05),"COLOR"]="green"
##绘图
colnames(data_bmi_)
library(calibrate)
setwd("E:/华大/表型/xiaomin/pheno_transform")
data_bmi_$SNP=stringi::stri_extract(str = data_bmi_$SNP,regex = "rs\\d+")
tiff(filename = "respiratory.adjust.assoc.logistic.sig.xlsx-TAGC.tiff",width = 2800,height = 1200,res = 200)
data_bmi_=data_bmi_[which(!is.na(data_bmi_$P)),]
qqman_(x = data_bmi_,chr = "CHR",bp = "BP",p = "P",snp = "SNP",col = "COLOR",highlight = data_bmi_[which(data_bmi_$COLOR!="black"),c("SNP")],annotateTop = TRUE,annotatePval = 8)
dev.off()
shell.exec(getwd())
