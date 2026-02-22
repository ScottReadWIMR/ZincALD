
HCVdats=readRDS("../mark/r_objects/211116_09_HCVall_GEOdata.RDS")
naf=readRDS("../mark/r_objects/211119_13_NAFLD_GEOdata.RDS")
hcc=readRDS("../mark/r_objects/211112_12_HCC_GEOdata.RDS")
hbv=readRDS("../mark/r_objects/211112_11_HBV_GEOdata.RDS")
alcLiv=readRDS("/Users/brian/Dropbox (Sydney Uni)/projects/scott/GEOstuff/210610_drunkenLiver/210707_AlcLiv_CleanData.RDS")

HCVdats[[1]]

hcv_MD<-lapply(HCVdats, tryCatch(function(X) pData(X[[1]]),error=function(e) NA))
naf_MD<-lapply(naf, function(X) tryCatch(pData(X[[1]]),error=function(e) NA))
hcc_MD<-lapply(hcc, function(X) tryCatch(pData(X[[1]]),error=function(e) NA))
hbv_MD<-lapply(hbv, function(X) tryCatch(pData(X[[1]]),error=function(e) NA))
alcLiv_MD<-lapply(alcLiv, function(X) tryCatch(pData(X[[1]]),error=function(e) NA))

indat_naf<-readRDS("../iram/r_objects/211202_14_FIB_GEOdata.RDS")
load('~/Dropbox (Sydney Uni)/projects/scott/GEOstuff/forTheKiddies/210809_processedData.rdata')

#PRJNA512027
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

naf_MD[[4]]<-pData(naf[[4]])
naf_MD$GSE167523<-pData(indat_naf$GSE167523[[1]])
naf_MD$GSE46300<-pData(indat_naf$GSE46300[[1]])
naf_MD$GSE48452<-pData(Y)
naf_MD$GSE135251=d$samples
naf_MD$PRJNA512027<-sa

alcLiv_MD[[1]]<-alcLiv[[1]]$samples
alcLiv_MD[[2]]<-pData(alcLiv[[2]])



lapply(hcv_MD, function(X) X[1:2,1:2])
lapply(naf_MD[-3], function(X) X[1:2,1:2])
lapply(hcc_MD, function(X) X[1:2,1:2])
lapply(hbv_MD, function(X) X[1:2,1:2])
lapply(alcLiv_MD, function(X) X[1:2,1:2])

openxlsx::write.xlsx(hcv_MD,"HCV_studies_metadata.xlsx")
openxlsx::write.xlsx(naf_MD[-3],"NAFLD_studies_metadata.xlsx")
openxlsx::write.xlsx(hcc_MD,"HCC_studies_metadata.xlsx")
openxlsx::write.xlsx(hbv_MD,"HBV_studies_metadata.xlsx")
openxlsx::write.xlsx(alcLiv_MD,"ALD_studies_metadata.xlsx")

hist(as.numeric(naf_MD$GSE66676$`bmi:ch1`),breaks=100,xlim=c(0,100)) ; abline(v=25,col='red')
 hist(as.numeric(naf_MD$GSE48452$`bmi:ch1`),breaks=100,xlim=c(0,100)) ; abline(v=25,col='red')

