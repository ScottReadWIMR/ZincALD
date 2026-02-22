library(GEOquery)
library(limma)
library(tidyverse)
library(edgeR)
# GSEs<-c(
# "GSE28619"
# )
# rdat<-lapply(unique(GSEs), getGEO)
# names(rdat)<-unique(GSEs)

keepcol=c("characteristics_ch1","disease status:ch1","characteristics_ch1.1","outcome at 6 months:ch1","liver sample group:ch1")

#"GSE155907"
pd<-pData(getGEO("GSE155907")[[1]])
pd<-pd[,colnames(pd) %in% keepcol] %>% mutate(characteristics_ch1.1=factor(characteristics_ch1.1,levels=unique(characteristics_ch1.1)))

edat<-data.table::fread("GSE155907_alcHep/GSE155907_GeneCounts.tsv.gz")
expdat<-data.frame(edat[,-c(1,3,4,5,6)],row.names=1)
genedat<-data.frame(edat[,c(2,3,1,4,5,6)],row.names=1)

GSE155907<-DGEList(counts=expdat,genes=genedat,samples=pd[,colnames(pd) %in% keepcol])

X=getGEO("GSE28619")
# data is already logged!
#GSE28619<-ExpressionSet(assayData=log2(exprs(affyPLM::normalize.ExpressionSet.quantiles(X[[1]]))),phenoData=AnnotatedDataFrame(pData(X[[1]])[,colnames(pData(X[[1]])) %in% keepcol]),featureData=AnnotatedDataFrame(fData(X[[1]])[,c(1:13)]),experimentData=experimentData(X[[1]]))

#run once
#saveRDS(list(GSE155907=GSE155907,GSE28619=GSE28619),file="210610_drunkenLiver/210707_AlcLiv_CleanData.RDS")


GSE155907 = GSE155907[rowSums(cpm(GSE155907) >1) >=1,]

design=model.matrix(~characteristics_ch1.1,data=GSE155907$samples)
y<-estimateDisp(GSE155907,design)
fit <- glmQLFit(y, design)
qlf_comp<-glmQLFTest(fit)
att_comp<-data.frame(topTags(qlf_comp,n=1000000))

outs=list(DEG=data.frame(att_comp),gene_anno=data.frame(GSE155907$genes),samdat=data.frame(GSE155907$samples),expression=data.frame(voom(GSE155907)$E))
#openxlsx::write.xlsx(outs,file="210610_drunkenLiver/210707_GSE155907_DEG_noCutoff_vsControl.xlsx",asTable=T,rowNames=T)
#subset(outs$DEG,grepl("^MT[0-9]",GeneName)) %>% ggplot(aes(x=GeneName,y=logFC,color=-log10(FDR))) + geom_point() + ggtitle("GSE155907 log foldchange of maybe metallothionein genes") + ylab("logFC (alc hep vs. alive")


# data is already logged AND IF YOU FUCK UP GETGEO, that dataset is ruined for the whole session
#X=GSE28619
X=getGEO("GSE28619")
Z=X[[1]]
Y=ExpressionSet(assayData=exprs(X[[1]]),phenoData=AnnotatedDataFrame(pData(X[[1]])[,colnames(pData(X[[1]])) %in% keepcol]),featureData=AnnotatedDataFrame(fData(X[[1]])[,c(1:13)]),experimentData=experimentData(X[[1]]))
keepgene<-!is.na(fData(Y)[,1]) & rowSums(exprs(Y)>2)>3
condvec<-gsub("X..","",make.names(pData(Y)[,2]))
cond<-factor(condvec,levels=unique(condvec))
des<-model.matrix(~cond)
f<-lmFit(Y[keepgene,],des)
f2<-eBayes(f)
outs=list(DEG=topTable(f2,n=1000000),gene_anno=data.frame(included=keepgene,fData(Y)),samdat=pData(Y),expresssion=data.frame(exprs(Y[keepgene,])))
#openxlsx::write.xlsx(outs,file="210610_drunkenLiver/210708_GSE28619_NOLOG_DEG_noCutoff_vsControl.xlsx",asTable=T,rowNames=T)

#
# subset(outs$DEG,grepl("metallothi",Gene.Title)) %>% ggplot(aes(x=Gene.Symbol,y=logFC,color=-log10(P.Value))) + geom_point() + ggtitle("GSE28619 log foldchange of metallothionein genes") + ylab("logFC (dead/trans vs. alive")

# subset(fData(X[[1]])[,1:13],grepl("CD24",`Gene Symbol`,ignore.case=T))

# subset(outs$DEG,grepl("^CD24",Gene.Symbol)) %>% ggplot(aes(x=Gene.Symbol,y=logFC,color=-log10(P.Value))) + geom_point() + ggtitle("GSE155907 log foldchange of maybe metallothionein genes") + ylab("logFC (alc hep vs. alive")


#mdsplot
pdf("210610_drunkenLiver/210707_MDSplots.pdf",width=16,height=10)
plotMDS(GSE155907,col=as.numeric(GSE155907$samples$characteristics_ch1.1),main="GSE155907: Super Enhancer Regulation of Cytokine-Induced Chemokine Production in Alcoholic Hepatitis")
plotMDS(GSE28619[keepgenes,],col=as.numeric(as.factor(pData(GSE28619)[,2])),main=paste("GSE28619",experimentData(GSE28619)@title))
dev.off()


