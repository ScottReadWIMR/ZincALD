
library(GEOquery)
library(limma)
library(tidyverse)
library(edgeR)
library(data.table)
library(Glimma)
library(biomaRt)
library(patchwork)

options(timeout=6000)

#nashid='GSE135251'
load("GEOstuff/forTheKiddies/210809_processedData.rdata")
GSE135251_exprs=cpm(d,log=T)
rownames(GSE135251_exprs)<-make.unique(d$genes$SYMBOL)
pGSE135251<-cbind(d$samples,PLK1=GSE135251_exprs["PLK1",rownames(d$samples)]) %>% ggplot(aes(x=group.ch1,y=PLK1)) + geom_boxplot() + geom_jitter(width=.1) + ggtitle('GSE135251',"just liver stuff")

#ithernash GSE48452
#finding PLK1 
fData(Y)[grep("PLK1",fData(Y)$gene_assignment)[1],]
pGSE48452<-data.frame(pData(Y),PLK1=exprs(Y)["7994109",rownames(pData(Y))]) %>% ggplot(aes(x=group.ch1,y=PLK1)) + geom_boxplot() + geom_jitter(width=.1) + ggtitle('GSE48452',Y@experimentData@title)



#hoda
GSE41804<-getGEO('GSE41804')
exp_GSE41804<-exprs(GSE41804[[1]])
rownames(exp_GSE41804)<-make.names(make.unique(fData(GSE41804[[1]])$`Gene Symbol`))

pGSE41804<-data.frame(pData(GSE41804[[1]])[,c("title","genotype:ch1","tissue:ch1")],PLK1=exp_GSE41804["PLK1",rownames(pData(GSE41804[[1]]))]) %>% mutate(gt=unlist(lapply(strsplit(genotype.ch1," "),function(X) X[1])),tissue.ch1=factor(tissue.ch1,levels=rev(unique(tissue.ch1))),id=gsub("[0-9]","",title)) %>% ggplot(aes(x=id,y=PLK1)) + geom_boxplot() + geom_point(position=position_jitter(width=.1)) + ggtitle('GSE41804',GSE41804[[1]]@experimentData@title)


#ueda
GSE44074<-getGEO('GSE44074')
exp_GSE44074<-exprs(GSE44074[[1]])
rownames(exp_GSE44074)<-make.names(make.unique(fData(GSE44074[[1]])$`Symbol`))

table(colnames(exp_GSE44074)==rownames(pData(GSE44074[[1]])))

p44074<-data.frame(pData(GSE44074[[1]])[,c("title","specimen group:ch1","description")],PLK1=exp_GSE44074["PLK1",rownames(pData(GSE44074[[1]]))])  %>% ggplot(aes(x=specimen.group.ch1,y=PLK1)) + geom_boxplot() + geom_point(position=position_jitter(width=.1)) + guides(x=guide_axis(angle=45)) + ggtitle('GSE44074',GSE44074[[1]]@experimentData@title)


(pGSE41804 + p44074) / (pGSE135251 + pGSE48452)

#nault
#wget https://ftp.ncbi.nlm.nih.gov/geo/series/GSE129nnn/GSE129322/suppl/GSE129322_non-normalized.txt.gz
GSE129322<-fread("~/projects/scott/GSE129322_non-normalized.txt.gz")

#PLK not there sadface

########################
save(GSE44074,GSE41804,d,Y,pGSE41804,p44074,pGSE135251,pGSE48452,file="r_objects/210928_someLiverViral.rdata")
#######################

Tcga

library(TCGAretriever)
all_studies <- get_cancer_studies()
interest <- grepl("liver", all_studies[,2],ignore.case=T)

all_studies[interest, ]

keep<-

my_csid <- "lihc_tcga_pan_can_atlas_2018"

 lihc_tcga<-get_genetic_profiles(csid = my_csid)

 lihc_cas<-get_case_lists(csid = my_csid)
 
 lihc_tcga_pan_can_atlas_2018_rna_seq_v2_mrna
 lihc_tcga_pan_can_atlas_2018_rna_seq_v2_mrna

  q_genes<-c("PLK1","GAPDH")
  q_cases <- "lihc_tcga_pan_can_atlas_2018_rna_seq_v2_mrna"
  rna_prf="lihc_tcga_pan_can_atlas_2018_rna_seq_v2_mrna"

  lihc_RNA<-TCGAretriever::get_profile_data(case_id = q_cases, gprofile_id = rna_prf, glist = q_genes)

rownames(lihc_RNA) <- lihc_RNA$COMMON
lihc_RNA <- lihc_RNA[2, -c(1,2)]

keepID<-intersect(rownames(clinDat), colnames(lihc_RNA))

clinDat<-get_clinical_data(case_id = q_cases)
rownames(clinDat)<-clinDat$CASE_ID


keepID<-intersect(rownames(clinDat), colnames(lihc_RNA))

tcga<-data.frame(clinDat[keepID,],PLK1=as.numeric(lihc_RNA["PLK1",keepID]))