library(tidyverse)
library(limma)
library(edgeR)
 library(affyPLM)
 library(data.table)
 library(patchwork)
 gdat<-readRDS("../ruby/data/mirSmash/IDtab.RDS")

  ##functions
  load(file="r_objects/220114_helperFunctions.rdata")

alcLivdat=list()
indat<-readRDS("/Users/brian/Dropbox (Sydney Uni)/projects/scott/GEOstuff/210610_drunkenLiver/210707_AlcLiv_CleanData.RDS")



#workspace
i=2

names(indat)[i]
names(indat)[i] %in% ls()
names(indat)[i] %in% names(alcLivdat)

dim(indat[[i]])
plotDensities(exprs(indat[[i]]))
plotDensities(log2(exprs(indat[[i]])))
plotDensities(log2(exprs(normalize(indat[[i]][[1]]))))

ex=exprs(indat[[i]])
plotDensities(ex)

ga<-data.frame(fData(indat[[i]]))
sa<-data.frame(pData(indat[[i]]))
ti=experimentData(indat[[i]])@title

head(ga,2)
gaid="GeneName"
head(unique(ga[,gaid]))

head(sa,2)
sa$scat<-gsub("disease state: ","",sa$characteristics_ch1.1)
said=c("scat")
unique(sa[,said])
saindex<-c(1,2)
unique(sa[,said])[saindex]
table(sa[,said])


#done
 #alcLivdat
 #
 #
 #
 #

ID=names(indat)[i]
#ID="PRJNA512027"
ID
GSE155907 <-mkobj()
testplot(GSE155907)
rm(sa,ga,ex,ti,i,said,gaid,saindex)



alcLivdat=list()
alcLivdat$GSE28619=GSE28619
alcLivdat$GSE155907=GSE155907

alcLivdat<-alcLivdat[order(names(alcLivdat))]

pl<-readRDS("r_objects/220113_LiverPlotTool.RDS")
p<-pl(alcLivdat,"PITX1")
cowplot::plot_grid(plotlist=p$expPlot) + p$dePlot


saveRDS(alcLivdat,file="r_objects/220119_AlcLiver_data.RDS")


# doing seq
i=1
plotDensities(cpm(indat[[i]],log=T))
ex<-cpm(indat[[1]],log=T)
ga<-indat[[i]]$genes
sa=indat[[i]]$samples
ti='Super Enhancer Regulation of Cytokine-Induced Chemokine Production in Alcoholic Hepatitis'


#new one GSE142530

id='	PRJNA597328, AKA 		SRP238597'#not there
#get data
 library("recount3")
 human_projects <- available_projects()
 
 proj_info <- subset(
    human_projects,
    project == "SRP238597" & project_type == "data_sources"
)
proj_info


#get other data
#https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE142530&format=file&file=GSE142530%5FAnnoted%2DRNAseq%2Dwith%2DSampleIDs%2Ecsv%2Egz

library(data.table)
library(tidyverse)
library(GEOquery)
library(edgeR)
dat<-fread("data/GSE142530_Annoted-RNAseq-with-SampleIDs.csv.gz",data.table=F)
sdat<-data.frame(row.names=colnames(dat)[-c(1,2)],title=do.call(rbind,strsplit(t(dat[1,-c(1,2)]),"_")))
load("r_objects/220114_helperFunctions.rdata")


dat$gene_name[dat$gene_name==""]=NA
dat$gene_name[is.na(dat$gene_name)]=dat$gene_id[is.na(dat$gene_name)]

counts<-data.frame(row.names=make.unique(dat[-1,]$gene_name),sapply(dat[-1,-c(1:2)],as.integer))

keep<-rowSums(cpm(counts,log=T)>1)>1 #pretty good but 2 was not inclusive so 1
#keep<-rowSums(counts>1)>1 #kinda shit
plotDensities(cpm(counts[keep,],log=T))

geodat<-getGEO('GSE142530')
pData(geodat[[1]])$description.1 %in% rownames(sdat)
pData(geodat[[1]])$description.1 == rownames(sdat)

samples=sdat %>% mutate(description.1=rownames(.)) %>% left_join(mutate(pData(geodat[[1]]),sid=make.names(`disease state:ch1`))) %>% column_to_rownames('description.1')
#filter out crappy sample
skeep<-!samples$title.1=="not.used"

#check
table(rownames(samples)==colnames(counts))
d<-DGEList(counts=counts[keep,skeep],samples=samples[skeep,])
d = calcNormFactors(d)
plotMDS(d,labels=d$samples$`disease.state.ch1`)

 gdat<-readRDS("../ruby/data/mirSmash/IDtab.RDS")
ID="GSE142530"
ex=cpm(counts[keep,skeep],log=T)
ga=data.frame(ex) %>% select(NULL) %>% mutate(gene_name=rownames(.))
gaid="gene_name"
sa<-samples[skeep,]
said="sid"
saindex=c(3,1,2)
ti=experimentData(geodat[[1]])@title

GSE142530=mkobj()

alcLivdat<-readRDS(file="r_objects/220119_AlcLiver_data.RDS")

alcLivdat$GSE142530=GSE142530

alcLivdat<-alcLivdat[order(names(alcLivdat))]

pl<-readRDS("r_objects/220113_LiverPlotTool.RDS")
p<-pl(alcLivdat,"AKR1B10")
cowplot::plot_grid(plotlist=p$expPlot) + p$dePlot

#done once for scott, manual mkobj

tmp<- mutate(data.frame(ex),gene_symbol=make.unique(ga[,gaid])) %>%   filter(gene_symbol %in% gdat$hgnc_symbol & !gene_symbol=="" & !is.na(gene_symbol)) %>%  pivot_longer(cols=-gene_symbol,values_to="Log_Expr",names_to="sample_id") %>%  left_join(data.frame(sample_id=make.names(rownames(sa)),vars=sa[,said])) %>%  mutate(vars=factor(vars,levels=unique(vars)[saindex]))

 #de
expdat= mutate(data.frame(ex),gene_symbol=make.unique(ga[,gaid])) %>%  filter(gene_symbol %in% gdat$hgnc_symbol & !gene_symbol=="" & !is.na(gene_symbol))  %>% data.frame(row.names=NULL) %>% column_to_rownames("gene_symbol")

cond= paste0(sa[,said])
cond<-factor(cond,levels=unique(cond)[saindex])

des<-model.matrix(~cond)
fit<-lmFit(expdat,des)
fit2<-eBayes(fit)
out<-lapply(2:ncol(des), function(i) topTable(fit2,coef=i,n=100000))
names(out)<-make.names(colnames(des)[-1])
for(x in 1:length(out)) {
out[[x]]$list=names(out)[x]
out[[x]]$gene_symbol=rownames(out[[x]])
}

out$exp=ex
out$samples=cbind(sa,des[,-1])

openxlsx::write.xlsx(lapply(out,data.frame),file="../scott/reports/220714_GSE142530_data_fixed.xlsx",rowNames=T,asTable=T)

saveRDS(alcLivdat,file="r_objects/220714_AlcLiver_data.RDS")
