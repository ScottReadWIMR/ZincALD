# E-MTAB-7751
## Get data, run once
# library(ArrayExpress)
# library(limma)
# options("timeout"=500)
# mexp<-getAE("E-MTAB-7751",type="raw",path='data') #all subsequent fail, had to get manifest from Illumina
# library(limma)
# targets<-readTargets("data/E-MTAB-7751.sdrf.txt")
# rownames(targets)<-gsub("data/|\\.idat","",targets$Array.Data.File)
# idat.files <- dir("data",pattern="idat",full.names=T)
# x<-read.idat(idat.files,bgxfile='data/HumanWG-6_V2_0_R4_11223189_A.bgx')
# boxplot(log2(x$E),range=0,ylab="log2 intensity")
# y <- neqc(x)
# colnames(y)<-gsub("data/|\\.idat","",colnames(y))
# table(colnames(y) %in% rownames(targets))
# boxplot(log2(y$E),range=0,ylab="log2 intensity",las=2)
# y$samples<-targets[colnames(y),]
# plotMDS(y,labels=targets$Characteristics.disease.)
# saveRDS(y,file="r_objects/211105_E-MTAB-7751_data.RDS")
# Nuked All raw data exept for manifest file

## Plot
y<-readRDS(file="r_objects/211105_E-MTAB-7751_data.RDS")
library(tidyverse)
library(limma)

exp<-y$E
rownames(exp)=make.unique(y$genes$Symbol)

goi<-c("PLK1",'NR1I2','ENTPD1','COL1A1','RORC','BRD4')

cand<-exp[unlist(lapply(goi, function(X) grep(X,rownames(exp),value=T))),]

library(tidyverse)
pdat<-data.frame(cbind(y$samples,t(cand))) %>% pivot_longer(cols=-(1:ncol(y$samples)),names_to='gene',values_to="log_exp")

ggplot(pdat,aes(color=Characteristics.infect.,x=Characteristics.infect.,y=log_exp)) + geom_boxplot() + guides(x=guide_axis(angle=45)) + facet_wrap(~gene,scale='free_y') + theme_bw()

ggplot(filter(pdat,gene=="PLK1"),aes(color=Characteristics.infect.,x=Characteristics.infect.,y=log_exp,shape=Characteristics.sex.,alpha=as.numeric(Characteristics.age.))) + geom_jitter(width=.2) + guides(x=guide_axis(angle=45)) + facet_wrap(~Characteristics.disease.) + theme_bw() + ggtitle("Recovered PLK1 expression from E-MTAB-7751",'Virus genotype-dependent gene expression in progressive liver disease prior to the development of advanced cirrhosis')

########### DATA ACCESS 12/11/21
library(GEOquery)
 library(affyPLM)
options(timeout="6000")
###HCV
HCVids=c(
	"GSE20119", #no PLK1
	"GSE15331", #seems Zscore/normalized, 2 colour!!!!
	"GSE14323",
	"GSE11190",
	"GSE168186",
	"GSE54648",
	"GSE54648",
	"GSE70779"
)
HCVdats<-lapply(HCVids,getGEO)
names(HCVdats)<-HCVids
saveRDS(HCVdats,paste0("r_objects/",format(Sys.time(), '%y%m%d_%H'),"_HCV_GEOdata.RDS"))

HCVids2=c(
	"GSE54101", #no PLK1
	"GSE54100", #seems Zscore/normalized, 2 colour!!!!
	"GSE15654",
	"GSE34798",
	"GSE7123"
)
HCVdats2<-lapply(HCVids2,getGEO)
names(HCVdats2)<-HCVids2

HCVdatsFull<-c(HCVdats,HCVdats2)
saveRDS(HCVdatsFull,paste0("r_objects/",format(Sys.time(), '%y%m%d_%H'),"_HCVall_GEOdata.RDS"))

HCVdats<-readRDS('r_objects/211116_09_HCVall_GEOdata.RDS')
library(tidyverse)

testID<-function(obj=HCVdats,id){
	require(affyPLM)
	message("number of subrecord(s): ", length(obj[id]))
	print(head(pData(obj[[id]][[1]]),2))
	print(head(fData(obj[[id]][[1]]),2))
	plotDensity(obj[[id]][[1]])
}
#prototype
GSE15331p<-data.frame(sid=paste0(pData(HCVdats[[2]][[1]])$`sample:ch1`,"_vs_",pData(HCVdats[[2]][[1]])$`sample:ch2`),PLK1=exprs(HCVdats[[2]][[1]])[ grep("PLK1",fData(HCVdats[[2]][[1]])$GENE_SYMBOL,ignore.case=T),]) %>% ggplot(aes(x=sid,color=sid,y=PLK1)) + geom_boxplot() + geom_jitter(width=.2) + guides(x=guide_axis(angle=45)) + ggtitle("PLK1 expression (GSE15331)",experimentData(HCVdats[[2]][[1]])@title)

#cycle
id="GSE14323"
 testID(id=id)
 voi<-c("characteristics_ch1","Gene Symbol")
 grep("PLK1",fData(HCVdats[[id]][[1]])[,voi[2]],value=T)
GSE14323p<-data.frame(sid=pData(HCVdats[[id]][[1]])[,voi[1]],PLK1=exprs(HCVdats[[id]][[1]])[grep("PLK1",fData(HCVdats[[id]][[1]])[,voi[2]]),]) %>% ggplot(aes(x=sid,color=sid,y=PLK1)) + geom_boxplot() + geom_jitter(width=.2) + guides(x=guide_axis(angle=45)) + ggtitle(paste0("PLK1 expression (",id,")"),experimentData(HCVdats[[id]][[1]])@title)


id=names(HCVdats)[4]
 testID(id=id)
 voi<-c("treatment status:ch1","Gene Symbol")
 grep("PLK1",fData(HCVdats[[id]][[1]])[,voi[2]],value=T)
 plotDensities(log2(exprs(HCVdats[[id]][[1]])))

GSE11190p<-data.frame(sid=pData(HCVdats[[id]][[1]])[,c("title",voi[1])],PLK1=log2(exprs(HCVdats[[id]][[1]]))[grep("PLK1",fData(HCVdats[[id]][[1]])[,voi[2]])[1],]) %>% mutate(id= gsub("-Pat.[0-9][0-9]|-Pat.[0-9]|-[0-9][0-9]$|[0-9]$","",unlist(lapply(strsplit(sid.title,"\\ \\["), function(X) X[1])))) %>% ggplot(aes(x=id,y=PLK1)) + geom_boxplot() + geom_jitter(width=.2,inherit.aes=F,aes(x=id, y=PLK1,color=sid.characteristics_ch1)) + guides(x=guide_axis(angle=45)) + ggtitle(paste0("PLK1 expression (",id,")"),experimentData(HCVdats[[id]][[1]])@title)




id=names(HCVdats)[5]
 testID(id=id) # no data

 [[1]] CHIP
 [[2]] PRb meth
 [[3]]# RNAseq GSM5133412-GSM5133439

fnam<-paste0("data/GSE168186/", dir('data/GSE168186/',pattern="quant"))
snam<-gsub("data/GSE168186/","",unlist(lapply(strsplit(fnam,"_"), function(X) X[1])))
rdat<-lapply(fnam,data.table::fread)
names(rdat)<-snam
tpms<-data.frame(row.names=rdat[[1]]$Name,do.call(cbind,lapply(rdat, function(X) X[,4])))
colnames(tpms)<-gsub('.TPM','',colnames(tpms))

#finding PLK
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
tx<-transcripts(TxDb.Hsapiens.UCSC.hg38.knownGene,columns=c("GENEID","tx_name"))
plk<-tx[unlist(lapply(tx$GENEID, function(X) '5347' %in% X))]$tx_name
plkexp<-tpms[plk,]
# ENST00000300093.9 is the proper one

pGSE168186<-data.frame(pData(HCVdats[[id]][[3]])[colnames(tpms),c('disease stage:ch1','cell line:ch1','source_name_ch1')],PLK1=as.numeric(tpms['ENST00000300093.9',])) %>% ggplot(aes(x=source_name_ch1,y=PLK1)) + geom_jitter(width=.1,inherit.aes=F,aes(x=source_name_ch1, y=PLK1,color=disease.stage.ch1)) + guides(x=guide_axis(angle=45)) + ggtitle(paste0("PLK1 expression (",id,") NOTE TPMS this time"),experimentData(HCVdats[[id]][[3]])@title) + theme_bw() + ylab("PLK1 main transcript TPMs")




id=names(HCVdats)[6]
 testID(id=id) #normalised
 voi<-c("cell population:ch1","ILMN_Gene")
 grep("PLK1",fData(HCVdats[[id]][[1]])[,voi[2]],value=T)
 plotDensities(log2(exprs(HCVdats[[id]][[1]]))) | plotDensities(exprs(HCVdats[[id]][[1]])) #looks normalized

pGSE54648<-data.frame(sid=pData(HCVdats[[id]][[1]])[,c("title",'days post infection:ch1','donor id:ch1',voi[1])],PLK1=exprs(HCVdats[[id]][[1]])[grep("PLK1",fData(HCVdats[[id]][[1]])[,voi[2]])[1],]) %>% ggplot(aes(x=sid.cell.population.ch1,y=PLK1)) + geom_boxplot() + geom_jitter(width=.2,inherit.aes=F,aes(x=sid.cell.population.ch1, y=PLK1,color=sid.donor.id.ch1)) + guides(x=guide_axis(angle=45)) + ggtitle(paste0("PLK1 expression (",id,")"),experimentData(HCVdats[[id]][[1]])@title) + facet_wrap(~sid.days.post.infection.ch1)


id=names(HCVdats)[8]
id
 testID(id=id) #logged but not normalized
 plotDensities(normalize(HCVdats[[id]][[1]])) # ok

 voi<-c("characteristics_ch1","GB_ACC")
 grep("NM_005030",fData(HCVdats[[id]][[1]])[,voi[2]],value=T)
pGSE70779<-data.frame(sid=pData(HCVdats[[id]][[1]])[,voi[1]],PLK1=exprs(normalize(HCVdats[[id]][[1]]))[grep("NM_005030",fData(HCVdats[[id]][[1]])[,voi[2]]),]) %>% ggplot(aes(x=sid,color=sid,y=PLK1)) + geom_boxplot() + geom_jitter(width=.2) + guides(x=guide_axis(angle=45)) + ggtitle(paste0("PLK1 expression (",id,")"),experimentData(HCVdats[[id]][[1]])@title)




id=names(HCVdats)[9] 
id
 testID(id=id) #Not logged but  normalised
plotDensities(log2(exprs(HCVdats[[id]][[1]]))) # ok

 voi<-c("disease state:ch1","Symbol") # add prognostic prediction:ch1
 grep("PLK1",fData(HCVdats[[id]][[1]])[,voi[2]],value=T)
pGSE54101<-data.frame(sid=pData(HCVdats[[id]][[1]])[,c("title",'prognostic prediction:ch1',voi[1])],PLK1=log2(exprs(HCVdats[[id]][[1]]))[grep("PLK1",fData(HCVdats[[id]][[1]])[,voi[2]]),]) %>% ggplot(aes(x=sid.prognostic.prediction.ch1,color=sid.prognostic.prediction.ch1,y=PLK1)) + geom_boxplot() + geom_jitter(width=.2) + guides(x=guide_axis(angle=45)) + ggtitle(paste0("PLK1 expression (",id,")"),experimentData(HCVdats[[id]][[1]])@title)


id=names(HCVdats)[10] 
id
 testID(id=id) #logged maybe not normalised
plotDensities(normalize(exprs(HCVdats[[id]][[1]]))) # ok

 voi<-c("prognostic prediction:ch11","gene_symbol") # add prognostic prediction:ch1
 grep("PLK1",fData(HCVdats[[id]][[1]])[,voi[2]],value=T)
# PLK1 not in 186 gene array


id=names(HCVdats)[11]
id
 testID(id=id) #not logged but normalized
 plotDensities(log2(exprs(HCVdats[[id]][[1]]))) # ok

 voi<-c("prediction (p<0.05):ch1","ILMN_Gene")
 grep("PLK1",fData(HCVdats[[id]][[1]])[,voi[2]],value=T)
pGSE15654<-
data.frame(sid=pData(HCVdats[[id]][[1]])[,voi[1]],PLK1=log2(exprs(HCVdats[[id]][[1]]))[grep("PLK1",fData(HCVdats[[id]][[1]])[,voi[2]]),]) %>% ggplot(aes(x=sid,color=sid,y=PLK1)) + geom_boxplot() + geom_jitter(width=.2) + guides(x=guide_axis(angle=45)) + ggtitle(paste0("PLK1 expression (",id,")"),experimentData(HCVdats[[id]][[1]])@title)



id=names(HCVdats)[12]
id
 testID(id=id) #normalized, TWO COLOR

 voi<-c("characteristics_ch2","GENE_SYMBOL")
 grep("PLK1",fData(HCVdats[[id]][[1]])[,voi[2]],value=T)
pGSE34798<-data.frame(sid=paste0(pData(HCVdats[[id]][[1]])$`disease_state:ch1`,"_vs_",pData(HCVdats[[id]][[1]])$`disease_state:ch2`),PLK1=exprs(HCVdats[[id]][[1]])[grep("PLK1",fData(HCVdats[[id]][[1]])[,voi[2]]),]) %>% ggplot(aes(x=sid,color=sid,y=PLK1)) + geom_boxplot() + geom_jitter(width=.2) + guides(x=guide_axis(angle=45)) + ggtitle(paste0("PLK1 expression (",id,")"),experimentData(HCVdats[[id]][[1]])@title)



id=names(HCVdats)[13]
id
 testID(id=id) # not normalized or logged

 plotDensities(log2(exprs(normalize(HCVdats[[id]][[1]])))) # ok
 voi<-c("characteristics_ch1","Gene Symbol")
 grep("PLK1",fData(HCVdats[[id]][[1]])[,voi[2]],value=T)
pGSE7123<-data.frame(row.names=rownames(pData(HCVdats[[id]][[1]])),do.call(rbind,lapply(strsplit(pData(HCVdats[[id]][[1]])[,voi[1]],"_"),function(X) return(c(id=X[which(X=="ID")+1],time=X[which(X=="TimePt")+1],gt=X[which(X=="Race")+1],respo=X[which(X=="Resp")+1])))),PLK1=log2(exprs(normalize((HCVdats[[id]][[1]]))))[grep("PLK1",fData(HCVdats[[id]][[1]])[,voi[2]]),]) %>% ggplot(aes(x=respo,fill=respo,y=PLK1)) + geom_boxplot() + geom_jitter(width=.2,aes(col=as.numeric(time))) + guides(x=guide_axis(angle=45)) + ggtitle(paste0("PLK1 expression (",id,")"),experimentData(HCVdats[[id]][[1]])@title)
pGSE7123




###HBV
HBVids=c(
	"GSE83148",
	"GSE84044",
	"GSE66698"
)
HBVdats<-lapply(HBVids,getGEO)
names(HBVdats)<-HBVids
saveRDS(HBVdats,paste0("r_objects/",format(Sys.time(), '%y%m%d_%H'),"_HBV_GEOdata.RDS"))


id=names(HBVdats)[1]
id
 testID(obj=HBVdats,id=id)
 voi<-c("title","Gene Symbol")
 grep("PLK1",fData(HBVdats[[id]][[1]])[,voi[2]],value=T)
pGSE83148<-data.frame(sid=pData(HBVdats[[id]][[1]])[,voi[1]],pData(HBVdats[[id]][[1]])[,c('alt:ch1','ast:ch1','hbv-dna:ch1')],PLK1=exprs(HBVdats[[id]][[1]])[grep("PLK1",fData(HBVdats[[id]][[1]])[,voi[2]])[1],]) %>% mutate(id=unlist(lapply(strsplit(sid,"\\-"), function(X) X[1]))) %>%  ggplot(aes(x=id,fill=id,y=PLK1)) + geom_boxplot() + geom_jitter(width=.1,inherit.aes=F,aes(x=id,color=alt.ch1,shape=hbv.dna.ch1,y=PLK1)) + guides(x=guide_axis(angle=45)) + ggtitle(paste0("PLK1 expression (",id,")"),experimentData(HBVdats[[id]][[1]])@title) + scale_fill_manual(values=c("white",'grey','grey30'))
pGSE83148




id=names(HBVdats)[2]
id
 testID(obj=HBVdats,id=id) #logged and norm :)
 voi<-c("title","Gene Symbol")
tmp<-data.frame(sid=pData(HBVdats[[id]][[1]])[,voi[1]],pData(HBVdats[[id]][[1]])[,c('age:ch1','batch:ch1','gender:ch1',"scheuer score g:ch1","scheuer score s:ch1","tissue:ch1")],PLK1=exprs(HBVdats[[id]][[1]])[grep("PLK1",fData(HBVdats[[id]][[1]])[,voi[2]])[1],]) %>% mutate(ageOver55=age.ch1>55)
vars<-make.names(c('gender:ch1',"scheuer score g:ch1","scheuer score s:ch1","ageOver55"))
p<-lapply(vars, function(X) ggplot(tmp,aes_string(x=X,color=X,y="PLK1")) + geom_boxplot() + geom_jitter(width=.1)  + guides(x=guide_axis(angle=45)) + ggtitle(paste0("PLK1 expression (",id," vs ",X,")"),experimentData(HBVdats[[id]][[1]])@title))
cowplot::plot_grid(plotlist=p)


id=names(HBVdats)[3]
id
 testID(obj=HBVdats,id=id) #logged and norm :)
 voi<-c("title","Gene Symbol")
 grep("PLK1",fData(HBVdats[[id]][[1]])[,voi[2]],value=T)
tmp<-data.frame(sid=pData(HBVdats[[id]][[1]])[,voi[1]],pData(HBVdats[[id]][[1]])[,c('age:ch1',"alt:ch1","gender:ch1","pre or post treatment:ch1","response:ch1","viral titer in serum:ch1")],PLK1=exprs(HBVdats[[id]][[1]])[grep("PLK1",fData(HBVdats[[id]][[1]])[,voi[2]])[1],])

bvars<-make.names(c("gender:ch1","pre or post treatment:ch1"))
cvars<-make.names(c('age:ch1',"alt:ch1","viral titer in serum:ch1","PLK1"))
p1<-lapply(bvars, function(X) ggplot(tmp,aes_string(x=X,color=X,y="PLK1")) + geom_boxplot() + geom_jitter(width=.1)  + guides(x=guide_axis(angle=45)) + ggtitle(paste0("PLK1 expression (",id," vs ",X,")"),experimentData(HBVdats[[id]][[1]])@title))
pheatmap::pheatmap(data.frame(row.names=rownames(tmp),apply(tmp[,cvars],2, as.numeric)),scale='column',main="Continuous Vars, Scale by Column")
p1[[1]] | p1[[2]]

###NAFLD
NAFLDids=c(
	"GSE66676",
	"GSE49541",
	"GSE834521"
)
NAFLDdats<-lapply(NAFLDids,function(X) tryCatch(getGEO(X), error=function(e) return(NA)))
names(NAFLDdats)<-NAFLDids

load("data/E-MEXP-3291.eSet.r") # from https://www.ebi.ac.uk/arrayexpress/files/E-MEXP-3291/E-MEXP-3291.eSet.r
plotDensities(study)
NAFLDdats$E_MEXP_3291<-rma(study)
#add gene annotations
library(affycoretools)
library('hugene10sttranscriptcluster.db')
NAFLDdats$E_MEXP_3291<-annotateEset(NAFLDdats$E_MEXP_3291,hugene10sttranscriptcluster.db,columns = c("PROBEID", "ENTREZID", "SYMBOL", "GENENAME"), multivals = "first")
plotDensities(NAFLDdats$E_MEXP_3291)

saveRDS(NAFLDdats,paste0("r_objects/",format(Sys.time(), '%y%m%d_%H'),"_NAFLD_GEOdata.RDS"))

id=names(NAFLDdats)[1]
id
 testID(obj=NAFLDdats,id=id) #rownorm 
 voi<-c("title","gene_assignment")
 grep("NM_005030",fData(NAFLDdats[[id]][[1]])[,voi[2]],value=T) ##1
tmp<-data.frame(sid=pData(NAFLDdats[[id]][[1]])[,voi[1]],pData(NAFLDdats[[id]][[1]])[,c("age (yrs):ch1",'bmi:ch1','cholesterol:ch1','gender:ch1','hdl:ch1','histology:ch1','ldl:ch1','sample type:ch1','triglycerides:ch1')],PLK1=exprs(NAFLDdats[[id]][[1]])[grep("NM_005030",fData(NAFLDdats[[id]][[1]])[,voi[2]])[1],])

bvars<-make.names(c('gender:ch1','histology:ch1'))
cvars<-make.names(c("age (yrs):ch1",'bmi:ch1','cholesterol:ch1','hdl:ch1','ldl:ch1','triglycerides:ch1',"PLK1"))
p1<-lapply(bvars, function(X) ggplot(tmp,aes_string(x=X,color=X,y="PLK1")) + geom_boxplot() + geom_jitter(width=.1)  + guides(x=guide_axis(angle=45)) + ggtitle(paste0("PLK1 expression (",id," vs ",X,")"),experimentData(NAFLDdats[[id]][[1]])@title))
pheatmap::pheatmap(data.frame(row.names=rownames(tmp),apply(tmp[,cvars],2, as.numeric)),scale='column',main="Continuous Vars, Scale by Column",show_rownames=F)
p1[[1]] / p1[[2]] 



id=names(NAFLDdats)[2]
id
 testID(obj=NAFLDdats,id=id) #logged and norm
 voi<-c("Stage:ch1","Gene Symbol")
 grep("PLK1",fData(NAFLDdats[[id]][[1]])[,voi[2]],value=T) ##1

pGSE49541<-data.frame(sid=pData(NAFLDdats[[id]][[1]])[,voi[1]],PLK1=exprs(NAFLDdats[[id]][[1]])[grep("PLK1",fData(NAFLDdats[[id]][[1]])[,voi[2]])[1],])  %>%  ggplot(aes(x=sid,fill=sid,y=PLK1)) + geom_boxplot() + geom_jitter(width=.1) + guides(x=guide_axis(angle=45)) + ggtitle(paste0("PLK1 expression (",id,")"),experimentData(NAFLDdats[[id]][[1]])@title) 
+ scale_fill_manual(values=c("white",'grey','grey30'))
pGSE49541



id=names(NAFLDdats)[4]
id
message("number of subrecord(s): ", length(obj[id]))
print(head(pData(NAFLDdats[[id]]),2))
print(head(fData(NAFLDdats[[id]]),2))
plotDensity(NAFLDdats[[id]])

 voi<-c("Factor.Value.DISEASE_STATE.","SYMBOL")
 grep("PLK1",fData(NAFLDdats[[id]])[,voi[2]],value=T) 

pE_MEXP_3291<-data.frame(sid=pData(NAFLDdats[[id]])[,voi[1]],PLK1=exprs(NAFLDdats[[id]])[grep("PLK1",fData(NAFLDdats[[id]])[,voi[2]])[1],]) %>% mutate(id=unlist(lapply(strsplit(sid,"\\-"), function(X) X[1]))) %>%  ggplot(aes(x=sid,fill=sid,y=PLK1)) + geom_boxplot() + geom_jitter(width=.1) + guides(x=guide_axis(angle=45)) + ggtitle(paste0("PLK1 expression (",id,")"),"Transcription profiling by array of H. sapiens liver cells to investigate global gene expression associated with the progression of NAFLD")
pE_MEXP_3291



###HCC
HCCids=c(
	"GSE44074",
	"GSE94660",
	"GSE121248",
	"GSE84402",
	"GSE25097",
	"GSE14520"
)

HCCdats<-lapply(HCCids,getGEO)
names(HCCdats)<-HCCids
saveRDS(HCCdats,paste0("r_objects/",format(Sys.time(), '%y%m%d_%H'),"_HCC_GEOdata.RDS"))


##1
id=names(HCCdats)[1]
 testID(HCCdats,id=id) # normalized
 voi<-c("specimen group:ch1","Symbol") #c("source_name_ch1" "characteristics_ch1","source_name_ch2","characteristics_ch2")
 grep("PLK1",fData(HCCdats[[id]][[1]])[,voi[2]],value=T)

pGSE44074<-data.frame(pData(HCCdats[[id]][[1]])[,c(voi[1],"title","source_name_ch1", "characteristics_ch1","source_name_ch2","characteristics_ch2")],PLK1=exprs(HCCdats[[id]][[1]])[grep("PLK1",fData(HCCdats[[id]][[1]])[,voi[2]])[1],]) %>% mutate(id= paste0(gsub("\\ |[0-9]","",title)," VS reference")) %>% ggplot(aes(x=specimen.group.ch1,y=PLK1,color=specimen.group.ch1)) + geom_boxplot() + geom_jitter(width=.1) + guides(x=guide_axis(angle=45)) + ggtitle(paste0("PLK1 expression (",id,")"),experimentData(HCCdats[[id]][[1]])@title) + geom_hline(yintercept=0)
pGSE44074



##2
id=names(HCCdats)[2]
id
 testID(HCCdats,id=id) # normalized
 voi<-c("disease state:ch1","Symbol") #to add'tissue:ch1'

dat<-data.frame(data.table::fread("data/GSE94660_RPKM_normalized.txt.gz"),row.names=1)
plotDensities(dat)#ok i guess
table(colnames(dat)==pData(HCCdats[[id]][[1]])$description) #ok

pGSE946604<-data.frame(pData(HCCdats[[id]][[1]])[,c(voi[1],"title",'tissue:ch1')],PLK1_normRPKM=as.numeric(dat['PLK1',])) %>% ggplot(aes(x=paste0(disease.state.ch1,"\n",tissue.ch1),y=PLK1_normRPKM,color=paste0(disease.state.ch1,"\n",tissue.ch1))) + geom_boxplot() + geom_jitter(width=.1) + guides(x=guide_axis(angle=45)) + ggtitle(paste0("PLK1 expression (",id,")"),experimentData(HCCdats[[id]][[1]])@title) + geom_hline(yintercept=0)
pGSE946604


sdats<-lapply(HCCdats, function(X) pData(X[[1]]))
openxlsx::write.xlsx(sdats,file="HCC_sampleAnnotations.xlsx")
#helpful!

#####template
id=names(HCCdats)[1]
 testID(HCCdats,id=id)
 voi<-c("treatment status:ch1","Gene Symbol")
 grep("PLK1",fData(HCCdats[[id]][[1]])[,voi[2]],value=T)
 plotDensities(log2(exprs(HCCdats[[id]][[1]])))

p<-data.frame(sid=pData(HCCdats[[id]][[1]])[,c("title",voi[1])],PLK1=log2(exprs(HCCdats[[id]][[1]]))[grep("PLK1",fData(HCCdats[[id]][[1]])[,voi[2]])[1],]) %>% mutate(id= gsub("-Pat.[0-9][0-9]|-Pat.[0-9]|-[0-9][0-9]$|[0-9]$","",unlist(lapply(strsplit(sid.title,"\\ \\["), function(X) X[1])))) %>% ggplot(aes(x=id,y=PLK1)) + geom_boxplot() + geom_jitter(width=.2,inherit.aes=F,aes(x=id, y=PLK1,color=sid.characteristics_ch1)) + guides(x=guide_axis(angle=45)) + ggtitle(paste0("PLK1 expression (",id,")"),experimentData(HCCdats[[id]][[1]])@title)
#########
##############