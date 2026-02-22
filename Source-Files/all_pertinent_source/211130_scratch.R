library(GEOquery)
 library(affyPLM)
 library(tidyverse)
 library(limma)
 library(ggpubr)
 library(patchwork)
options(timeout="6000")
###HCV
FIBids=c(
	"GSE167523",
"GSE163211",
"GSE151158",
"GSE46300",
"GSE49542"
)# (both RNA and methylation)

PRJNA512027 is a SRA 


FIBdats<-lapply(FIBids,getGEO)
names(FIBdats)<-FIBids


saveRDS(FIBdats,paste0("r_objects/",format(Sys.time(), '%y%m%d_%H'),"_FIB_GEOdata.RDS"))



testID<-function(obj=FIBdats,id){
	require(affyPLM)
	message(id, " : number of subrecord(s): ", length(obj[id]))
	print(experimentData(obj[[id]][[1]])@title)
	print(head(pData(obj[[id]][[1]]),2))
	if(nrow(exprs(obj[[id]][[1]]))==0) message("no expData: could be seq, google: ",id,"\n\n") else {
	print(head(fData(obj[[id]][[1]]),2))
	plotDensity(obj[[id]][[1]])
	}
}

for(i in 1:length(FIBdats)) testID(obj=FIBdats,id=names(FIBdats)[i])


#1
id=names(FIBdats)[1]
id
 testID(id=id)
 voi<-c("disease subtype:ch1",NA)

rdat<-data.table::fread("data/GSE167523_Raw_gene_counts_matrix.txt.gz")
counts<-data.frame(row.names=make.unique(rdat$Gene),select(rdat,-Gene))
cpms<-edgeR::cpm(counts,log=T)
 grep("RARRES1",rownames(cpms),value=T)
table(colnames(cpms) == pData(FIBdats[[id]][[1]])$description) #fine
plotDensities(cpms)#ok i guess

GSE167523p<-data.frame(row.names=pData(FIBdats[[id]][[1]])$description,sid=pData(FIBdats[[id]][[1]])[,voi[1]],RARRES1=as.numeric(cpms["RARRES1",])) %>% ggplot(aes(x=sid,color=sid,y=RARRES1)) + geom_boxplot() + geom_jitter(width=.2) + guides(x=guide_axis(angle=45)) + ggtitle(paste0("RARRES1 expression (",id,")"),experimentData(FIBdats[[id]][[1]])@title) + stat_compare_means(method='wilcox.test',paired=F)
GSE167523p

#2
id=names(FIBdats)[2]
id
 testID(id=id)
 voi<-c("characteristics_ch1.8","ORF")
 grep("RARRES1",fData(FIBdats[[id]][[1]])[,voi[2]],value=T) #doesnt exist
plotDensities(log2(exprs(normalize(FIBdats[[id]][[1]],)))) #possibly overdone?


#3
id=names(FIBdats)[3]
id
 testID(id=id)
 voi<-c("source_name_ch1","ID")
 grep("RARRES1",fData(FIBdats[[id]][[1]])[,voi[2]],value=T)
  grep("RAR",rownames(fData(FIBdats[[id]][[1]])),value=T)

#4 
id=names(FIBdats)[4]
id
 testID(id=id)
 voi<-c("source_name_ch1","ILMN_Gene")
 grep("RARRES1",fData(FIBdats[[id]][[1]])[,voi[2]],value=T) #three probes,
 filter(fData(FIBdats[[id]][[1]]),ILMN_Gene=="RARRES1") # went with ALL isoforms (2, Probe_Type=A(see https://support.illumina.com/content/dam/illumina-support/documents/documentation/software_documentation/genomestudio/genomestudio-2011-1/genomestudio-gx-module-v1-0-user-guide-11319121-a.pdf))

plotDensities(log2(exprs(FIBdats[[id]][[1]]))) #dist weird, illumina amirite?: logging breaks it... will see what raw data does

pGSE46300<-data.frame(sid=pData(FIBdats[[id]][[1]])[,voi[1]],RARRES1=exprs(FIBdats[[id]][[1]])[grep("RARRES1",fData(FIBdats[[id]][[1]])[,voi[2]])[2],]) %>% ggplot(aes(x=sid,color=sid,y=RARRES1)) + geom_boxplot() + geom_jitter(width=.2) + guides(x=guide_axis(angle=45)) + ggtitle(paste0("RARRES1 expression (",id,")"),experimentData(FIBdats[[id]][[1]])@title) + stat_compare_means(method='wilcox.test',paired=F)



#5
id=names(FIBdats)[5]
id
 testID(id=id)
 voi<-c("characteristics_ch1.1","UCSC_RefGene_Name")# or characteristics_ch1.2
 grep("RARRES1",fData(FIBdats[[id]][[1]])[,voi[2]],value=T)
 filter(fData(FIBdats[[id]][[1]]),grepl("RARRES1",UCSC_RefGene_Name))

mdat<-data.frame(exprs(FIBdats[[id]][[1]]))[grepl("RARRES1",fData(FIBdats[[id]][[1]])$UCSC_RefGene_Name),]
pdat<- filter(fData(FIBdats[[id]][[1]]),grepl("RARRES1",UCSC_RefGene_Name))
sdat<-pData(FIBdats[[id]][[1]])

GSE49542m<-cbind(sdat,t(mdat)) %>% pivot_longer(cols=-c(1:ncol(sdat)),values_to="M.score",names_to="ProbeID") %>% left_join(y=rownames_to_column(pdat,"ProbeID")) %>% ggplot(aes(x=as.factor(as.numeric(MAPINFO)),color=characteristics_ch1.2,y=M.score)) + geom_boxplot()  + stat_compare_means(method='wilcox.test',paired=F,label = "p.signif",show.legend=T) + geom_text(data=pdat,inherit.aes=F,aes(y=-10,x=as.factor(as.numeric(MAPINFO)),label=Relation_to_UCSC_CpG_Island),color='black',angle=90) + ggtitle(paste0("RARRES1 Methylation (",id,")"),experimentData(FIBdats[[id]][[1]])@title) + xlab("Location (chr3)")
GSE49542m

#RNASeq is previous one
GSE49541<-getGEO('GSE49541')

FIBdats$GSE49541=GSE49541
saveRDS(FIBdats,paste0("r_objects/",format(Sys.time(), '%y%m%d_%H'),"_FIB_GEOdata.RDS"))

id="GSE49541"
id
 testID(id=id)
 voi<-c("Stage:ch1","Gene Symbol")
 grep("RARRES1",fData(FIBdats[[id]][[1]])[,voi[2]],value=T) #more than 1
 filter(fData(FIBdats[[id]][[1]]),grepl("RARRES1",`Gene Symbol`)) # pick consensus for sequence type (221872_at),3)


GSE49541p<-data.frame(sid=pData(FIBdats[[id]][[1]])[,voi[1]],RARRES1=exprs(FIBdats[[id]][[1]])[grep("RARRES1",fData(FIBdats[[id]][[1]])[,voi[2]])[3],]) %>% ggplot(aes(x=sid,color=sid,y=RARRES1)) + geom_boxplot() + geom_jitter(width=.2) + guides(x=guide_axis(angle=45)) + ggtitle(paste0("RARRES1 expression (",id,")"),experimentData(FIBdats[[id]][[1]])@title) + stat_compare_means(method='wilcox.test',paired=F)
GSE49541p + GSE49542m + plot_layout(widths=c(1,6))


### summary plots
GSE49541p + GSE49542m + plot_layout(widths=c(1,6))
pGSE46300
GSE167523p


### generic
id=names(FIBdats)[1]
id
 testID(id=id)
 voi<-c("characteristics_ch1","Gene Symbol")
 grep("RARRES1",fData(FIBdats[[id]][[1]])[,voi[2]],value=T)

p<-data.frame(sid=pData(FIBdats[[id]][[1]])[,voi[1]],RARRES1=exprs(FIBdats[[id]][[1]])[grep("RARRES1",fData(FIBdats[[id]][[1]])[,voi[2]]),]) %>% ggplot(aes(x=sid,color=sid,y=RARRES1)) + geom_boxplot() + geom_jitter(width=.2) + guides(x=guide_axis(angle=45)) + ggtitle(paste0("RARRES1 expression (",id,")"),experimentData(FIBdats[[id]][[1]])@title) + stat_compare_means(method='wilcox.test',paired=F)


# RECOUNT3 FTW
id='PRJNA512027, AKA SRP174668'

 library("recount3")
 human_projects <- available_projects()
 
 proj_info <- subset(
    human_projects,
    project == "SRP174668" & project_type == "data_sources"
)

rse_gene_SRP174668 <- create_rse(proj_info)
rse_gene_SRP174668

FIBdats$PRJNA512027<-rse_gene_SRP174668
#saveRDS(FIBdats,paste0("r_objects/",format(Sys.time(), '%y%m%d_%H'),"_FIB_GEOdata.RDS")) # r_objects/211202_14_FIB_GEOdata.RDS

assay(rse_gene_SRP174668, "counts") <- transform_counts(rse_gene_SRP174668)
plotDensities(edgeR::cpm(assay(rse_gene_SRP174668),log=T)) # good

rse_gene_SRP174668_expanded <-    expand_sra_attributes(rse_gene_SRP174668)
colData(rse_gene_SRP174668_expanded)[, ncol(colData(rse_gene_SRP174668)):ncol(colData(rse_gene_SRP174668_expanded))]

voi=c("sra_attribute.isolate","gene_name")

PRJNA512027p<-data.frame(sid=colData(rse_gene_SRP174668_expanded)[,voi[1]],RARRES1=edgeR::cpm(assay(rse_gene_SRP174668_expanded),log=T)[grep("RARRES1",as(rowRanges(rse_gene_SRP174668_expanded),'data.frame')[,voi[2]]),]) %>% mutate(sid=factor(sid,levels=unique(sid)[c(1,5,6,3,2,8,4,7)])) %>% ggplot(aes(x=sid,color=sid,y=RARRES1)) + geom_boxplot() + geom_jitter(width=.2) + guides(x=guide_axis(angle=45)) + ggtitle(paste0("RARRES1 expression (",id,")"),'Nonalcoholic steatohepatitis (NASH) is strongly associated with obesity and type 2 diabetes') + stat_compare_means(method='anova',paired=F)

#summary data for delivery
sumdat<-list(
	GSE167523plot=GSE167523p$data,
	GSE46300plot=pGSE46300$data,
	GSE49542m_plot=GSE49542m$data,
	GSE49541e_plot=GSE49541p$data,
	PRJNA512027_plot=PRJNA512027p$data
)
openxlsx::write.xlsx(sumdat,'reports/211202_Fibrosis_plotData.xlsx',rowNames=T,asTable=T)

sdat<-lapply(FIBdats[-c(2,3)], function(X) pData(X[[1]]))
sdat$PRJNA512027<-colData(rse_gene_SRP174668_expanded)[, ncol(colData(rse_gene_SRP174668)):ncol(colData(rse_gene_SRP174668_expanded))]
openxlsx::write.xlsx(lapply(sdat,data.frame),'reports/211202_Fibrosis_availableAnnotations.xlsx',rowNames=T,asTable=T)