library(GEOquery)
library(limma)
library(tidyverse)
library(edgeR)
library(data.table)
library(Glimma)
library(biomaRt)
library(sva)
GSEs<-c( 
"GSE103580", #various stages of alcoholic liver disease
"GSE94397",	# sever alcoholic hepatitis patients (derivation cohort)
"GSE94399"	# sever alcoholic hepatitis
 )
options(timeout=6000)
 rdat<-lapply(unique(GSEs), getGEO)
names(rdat)<-unique(GSEs)

par(mfrow=c(1,3)) ;  lapply(rdat, function(X) plotDensities(log2(exprs(X[[1]]))))

 sdatlist<-lapply(rdat, function(X) pData(X[[1]])[c('title','geo_accession','characteristics_ch1.1',"characteristics_ch1")])
for(i in 1:length(sdatlist)) sdatlist[[i]]$list=names(sdatlist)[i]
colnames(sdatlist[[1]])<-colnames(sdatlist[[1]])[c(1:2,4,3,5)]

lapply(sdatlist,head)

comPat<-rbindlist(lapply(sdatlist, function(X) X[,colnames(sdatlist[[2]])] %>% rownames_to_column("sample_id")))

comDat<-do.call(cbind,lapply(rdat, function(X) exprs(X[[1]])))

comQDat<-do.call(cbind,lapply(rdat, function(X) exprs(affyPLM::normalize.ExpressionSet.quantiles(X[[1]]))))

comBat<-sva::ComBat(log2(comDat),batch=comPat$list)

comQBat<-sva::ComBat(log2(comQDat),batch=comPat$list)
#based on https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0156594

par(mfrow=c(2,4))
plotMDS(log2(comDat),col=as.numeric(as.factor(comPat$list)),main="input data, by batch")
plotMDS(log2(comDat),col=as.numeric(as.factor(comPat$characteristics_ch1.1)),main="input data, by sample type")

plotMDS(log2(comQDat),col=as.numeric(as.factor(comPat$list)),main="within batch norm input data, by batch")
plotMDS(log2(comQDat),col=as.numeric(as.factor(comPat$characteristics_ch1.1)),main="within batch norm input data, by sample type")

plotMDS(comBat,col=as.numeric(as.factor(comPat$list)),main="combat data, by batch")
plotMDS(comBat,col=as.numeric(as.factor(comPat$characteristics_ch1.1)),main="combat data, by sample type")

plotMDS(comQBat,col=as.numeric(as.factor(comPat$list)),main="combat within batch norm data, by batch")
plotMDS(comQBat,col=as.numeric(as.factor(comPat$characteristics_ch1.1)),main="combat within batch norm data, by sample type")
dev.copy2pdf(file="reports/220204_BatchNormArrays_ALD.pdf")

par(mfrow=c(2,2))
 plotDensities(log2(comDat),main="raw")
 plotDensities(comBat,main="combat")
 plotDensities(log2(comQDat),main="within batch norm")
 plotDensities(comQBat,main="combat within batch norm")
 dev.copy2pdf(file="reports/220204_BatchNormArrays_ALD_densities.pdf")

saveRDS(rdat,file="r_objects/220204_AlcHep_Unify_rdat.RDS")

forScott=list(
	expression=comQBat,
	sampleAnnotation=comPat,
	geneAnnotation=fData(rdat[[1]][[1]])[,1:27]
)
saveRDS(forScott,file="r_objects/220204_AlcHep_Unify_ScottDat.RDS")

openxlsx::write.xlsx(forScott,file='reports/220204_Normalized_BatchCorrected_ALD_expression.xlsx',rowNames=c(T,F,F))