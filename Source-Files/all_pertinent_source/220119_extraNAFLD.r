  
  
library(tidyverse)
library(limma)
library(edgeR)
 library(affyPLM)
 library(data.table)
 library(patchwork)
 gdat<-readRDS("../ruby/data/mirSmash/IDtab.RDS")

  ##functions
  load(file="r_objects/220114_helperFunctions.rdata")

nafldat<-readRDS("r_objects/220114_NAFLD_data.RDS")


indat<-readRDS("../iram/r_objects/211202_14_FIB_GEOdata.RDS")


GSE167523 is seq
GSE163211 is 800 features
GSE151158 is 618 features
GSE49542 is meth
GSE49541 is already in

#workspace
i=7

ID=names(indat)[i]
names(indat)[i] %in% ls()
names(indat)[i] %in% names(nafldat)

dim(indat[[i]][[1]])
plotDensities(exprs(indat[[i]][[1]]))
plotDensities(log2(exprs(indat[[i]][[1]])))
plotDensities(log2(exprs(normalize(indat[[i]][[1]]))))

ex=exprs(indat[[i]][[1]])
plotDensities(ex)

ga<-data.frame(fData(indat[[i]][[1]]))
sa<-data.frame(pData(indat[[i]][[1]]))
ti=experimentData(indat[[i]][[1]])@title

head(ga,2)
gaid="gene_name"
head(unique(ga[,gaid]))

head(sa,2)
said=c("disease.subtype.ch1")
unique(sa[,said])
saindex<-c(1,2)
unique(sa[,said])[saindex]
table(sa[,said])


#done
 #GSE46300
 #PRJNA512027
 #GSE135251
 #
 #

ID=names(indat)[i]
ID
#or ID="PRJNA512027"
GSE167523 <-mkobj()
testplot(GSE167523)
rm(sa,ga,ex,ti,i,said,gaid,saindex)

nafldat$GSE46300=GSE46300
nafldat$PRJNA512027=PRJNA512027
nafldat$GSE135251=GSE135251_
nafldat$GSE167523=GSE167523


nafldat<-nafldat[order(names(nafldat))]

pl<-readRDS("r_objects/220113_LiverPlotTool.RDS")
p<-pl(nafldat,"PITX1")
cowplot::plot_grid(plotlist=p$expPlot) + p$dePlot

#nafldat$GSE46300$de<-filter(nafldat$GSE46300$de,!list=='none.vsCont') #coz only 2 

saveRDS(nafldat,file="r_objects/220119_NAFLD_data.RDS")


#extra from IRAM added above
id='PRJNA512027, AKA 	SRP174668'#there
#get data
 library("recount3")
 human_projects <- available_projects()
 
 proj_info <- subset(
    human_projects,
    project == "SRP174668" & project_type == "data_sources"
)
proj_info
rse_gene <- create_rse(proj_info)

#get data obj
counts1 <- compute_read_counts(rse_gene)
plotDensities(edgeR::cpm(counts1,log=T)) # good
rse_gene_expanded <-    expand_sra_attributes(rse_gene)

#make data for mine
sa<-colData(rse_gene_expanded)[, ncol(colData(rse_gene)):ncol(colData(rse_gene_expanded))] 
 ex=edgeR::cpm(counts1,log=T)
 ga=rowData(rse_gene)
 table(rownames(ga)==rownames(ex))
 ti="Nonalcoholic steatohepatitis (NASH) is strongly associated with obesity and type 2 diabetes. Transcriptomic Profiling of Obesity-Related Nonalcoholic Steatohepatitis Reveals a Core Set of Fibrosis-Specific Genes"


GSE167523 added too

indat=readRDS("../mark/r_objects/211119_13_NAFLD_GEOdata.RDS")
 i=1

rdat<-data.table::fread("../iram/data/GSE167523_Raw_gene_counts_matrix.txt.gz")
counts<-data.frame(row.names=make.unique(rdat$Gene),select(rdat,-Gene))
plotDensities(edgeR::cpm(counts,log=T)) # fine

 ex=edgeR::cpm(counts,log=T)
 ga=data.frame(ex) %>% select(NULL) %>% mutate(gene_name=rownames(.))
 table(rownames(ga)==rownames(ex))

sa<-data.frame(pData(indat[[i]][[1]]))
table(sa$description == colnames(ex))
colnames(ex)=rownames(sa)

ti=experimentData(indat[[i]][[1]])@title



## FROM SCOTT

load('~/Dropbox (Sydney Uni)/projects/scott/GEOstuff/forTheKiddies/210809_processedData.rdata')
Y=GSE48452 NAFLD,array#in already
"GSE48452" %in% names(nafldat) #in already

d==GSE135251 #NALFD,seq not in , PRJNA558102, AKA SRP217231' #NOT in RECOUNT 
"GSE135251" %in% names(nafldat) #so do from obj

plotDensities(cpm(d,log=T)) # great

#make obj
ex=cpm(d,log=T)
ga=d$genes
sa=d$samples
ti="Transcriptomic profiling across the nonalcoholic fatty liver disease spectrum reveals gene signatures for steatohepatitis and fibrosis"

#insert to top
