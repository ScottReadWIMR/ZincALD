library(GEOquery)
library(limma)
library(tidyverse)
GSEs<-c(
"GSE103580",
"GSE94397",
"GSE94399"
)
rdat<-lapply(unique(GSEs), getGEO)
names(rdat)<-unique(GSEs)
saveRDS(rdat,file="210610_drunkenLiver/210610_AlcLiv_geoData.RDS")

keepcol=c("characteristics_ch1","disease status:ch1","characteristics_ch1.1","outcome at 6 months:ch1")

expdat<-lapply(rdat, function(X){
	ExpressionSet(assayData=log2(exprs(affyPLM::normalize.ExpressionSet.quantiles(X[[1]]))),phenoData=AnnotatedDataFrame(pData(X[[1]])[,colnames(pData(X[[1]])) %in% keepcol]),featureData=AnnotatedDataFrame(fData(X[[1]])[,c(1,14:16,20,21)]),experimentData=experimentData(X[[1]]))
})
saveRDS(exp,file="210610_drunkenLiver/210610_AlcLiv_CleanData.RDS")


X=expdat[[1]]
keepgenes<-!is.na(fData(X)[,1]) & rowSums(exprs(X)>2)>3
condvec<-make.names(pData(X)[,3])
cond<-factor(condvec,levels=unique(condvec)[c(2,3,1)])
des<-model.matrix(~cond)
f<-lmFit(X[keepgenes,],des)
f2<-eBayes(f)
outs=list(DEG=topTable(f2,n=1000000),gene_anno=data.frame(included=keepgenes,fData(X)),samdat=pData(X),expresssion=data.frame(exprs(X[keepgenes,])))
openxlsx::write.xlsx(outs,file="210610_drunkenLiver/210610_GSE103580_DEG_noCutoff_vsMild.xlsx",asTable=T,rowNames=T)



#mdsplot
pdf("210610_drunkenLiver/210610_MDSplots.pdf",width=16,height=10)
for (i in 1:length(expdat)) plotMDS(expdat[[i]],col=as.numeric(as.factor(pData(expdat[[i]])[,3])),main=paste(names(expdat)[i],experimentData(expdat[[i]])@title))
dev.off()




X=expdat[[2]]
keepgenes<-!is.na(fData(X)[,1]) & rowSums(exprs(X)>2)>3
condvec<-make.names(pData(X)[,3])
cond<-factor(condvec,levels=unique(condvec)[c(2,1)])
des<-model.matrix(~cond)
f<-lmFit(X[keepgenes,],des)
f2<-eBayes(f)
outs=list(DEG=topTable(f2,n=1000000),gene_anno=data.frame(included=keepgenes,fData(X)),samdat=pData(X),expresssion=data.frame(exprs(X[keepgenes,])))
openxlsx::write.xlsx(outs,file="210610_drunkenLiver/210610_GSE94397_DEG_noCutoff_vsAlive.xlsx",asTable=T,rowNames=T)

subset(outs$DEG,grepl("metallothi",Gene.Title)) %>% ggplot(aes(x=Gene.Symbol,y=logFC,color=-log10(P.Value))) + geom_point() + ggtitle("GSE94397 log foldchange of metallothionein genes") + ylab("logFC (dead/trans vs. alive")


X=expdat[[3]]
keepgenes<-!is.na(fData(X)[,1]) & rowSums(exprs(X)>2)>3
condvec<-make.names(pData(X)[,3])
cond<-factor(condvec,levels=unique(condvec)[c(2,1)])
des<-model.matrix(~cond)
f<-lmFit(X[keepgenes,],des)
f2<-eBayes(f)
outs=list(DEG=topTable(f2,n=1000000),gene_anno=data.frame(included=keepgenes,fData(X)),samdat=pData(X),expresssion=data.frame(exprs(X[keepgenes,])))
openxlsx::write.xlsx(outs,file="210610_drunkenLiver/210610_GSE94399_DEG_noCutoff_vsAlive.xlsx",asTable=T,rowNames=T)

subset(outs$DEG,grepl("metallothi",Gene.Title)) %>% ggplot(aes(x=Gene.Symbol,y=logFC,color=-log10(P.Value))) + geom_point() + ggtitle("GSE94399 log foldchange of metallothionein genes") + ylab("logFC (dead/trans vs. alive")

#where GSE9439[79] is respective DEG from last 2
pp<-full_join(GSE94397,GSE94399[,-c(2:6)],by="ID",suffix=c("_GSE94397","_GSE94399"))
ggplot(pp,aes(logFC_GSE94397,logFC_GSE94399)) + geom_point(color='lightgrey') + geom_point(data=subset(pp,grepl("metallothi",Gene.Title))) + geom_text_repel(data=subset(pp,grepl("metallothi",Gene.Title)),aes(label=Gene.Symbol),max.overlaps=Inf) + geom_hline(yintercept=0) + geom_vline(xintercept=0)