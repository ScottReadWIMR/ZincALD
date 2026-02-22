1. Concerning patients with chronic hepatitis C, PLK1 expression is higher in HCV-infected liver than uninfected liver (GSE15331, Peng 2009, BMC Genomics).
2. Further, PLK1 expression in the liver is reduced following HCV cure with DAA therapy (GSE70779, Meissner 2016, J Viral Hep).
3. Finally, in patients with early (Child-Pugh A), HCV-induced cirrhosis, PLK1 is part of a 218 gene expression signature that predicts predicts death, progression to advanced cirrhosis and development of HCC over 10 years (GSE15654, Hoshida 2013, Gastroenterology). This 218 gene signature (which includes PLK1) was originally identified in patients with HCV, cirrhosis and HCC, and predicted overall death (Hoshida 2008, NEJM).
4. Looking at GSE15654, can we identify individual patients in GSE15654 who developed HCC?? It would be nice if they had higher PLK1 at baseline than those who did not develop HCC.
5. I tried doing this on the GEO2R website, grouping by “1” or “0” in HCC column. From a bit of looking around I think the PLK1 probe ID is ILMN_1736176? If so, looking at the full table from this simple GEO2R analysis, it looks like p value is 0.0835519, which would not be significant. Is that correct? Is there a better way to look at this??
Similar to hepatitis C, patients with chronic hepatitis B have higher PLK1 expression in the liver than uninfected patients (GSE83148, Zhou 2017, Liver Int))
 

library(tidyverse)
 library(affyPLM)
 library(ggpubr)
HCVdats<-readRDS('r_objects/211116_09_HCVall_GEOdata.RDS')
HCVdats

p1<-data.frame(sid=paste0(pData(HCVdats[[2]][[1]])$`sample:ch1`,"_vs_",pData(HCVdats[[2]][[1]])$`sample:ch2`),PLK1=exprs(HCVdats[[2]][[1]])[ grep("PLK1",fData(HCVdats[[2]][[1]])$GENE_SYMBOL,ignore.case=T),]) %>% filter(sid %in% c("HCV positive_vs_Reference","HCV negative_vs_Reference")) %>% ggplot(aes(x=sid,y=PLK1)) + geom_boxplot() + geom_jitter(width=.2) + guides(x=guide_axis(angle=45))  + stat_compare_means(method='anova')+ ggtitle("PLK1 expression (GSE15331)",experimentData(HCVdats[[2]][[1]])@title) + theme_bw(base_size=12) + ylab("PLK Expression (Ratio)") + xlab("")

id=names(HCVdats)[8]
voi<-c("characteristics_ch1","GB_ACC")
p2<-data.frame(sid=pData(HCVdats[[id]][[1]])[,voi[1]],PLK1=exprs(normalize(HCVdats[[id]][[1]]))[grep("NM_005030",fData(HCVdats[[id]][[1]])[,voi[2]]),]) %>% ggplot(aes(x=sid,y=PLK1)) + geom_boxplot() + stat_compare_means(method='anova') + geom_jitter(width=.2) + guides(x=guide_axis(angle=45)) + stat_compare_means(method='anova') + ggtitle(paste0("PLK1 expression (",id,")"),experimentData(HCVdats[[id]][[1]])@title)+ theme_bw(base_size=12) + ylab("PLK Expression (Log intensity)") + xlab("")


id=names(HCVdats)[11]
voi<-c("prediction:ch1","ILMN_Gene")
data.frame(sid=pData(HCVdats[[id]][[1]])[,voi[1]],PLK1=log2(exprs(HCVdats[[id]][[1]]))[grep("PLK1",fData(HCVdats[[id]][[1]])[,voi[2]]),],hcc=pData(HCVdats[[id]][[1]])[,"hcc:ch1"])  %>% ggplot(aes(x=sid,color=hcc,y=PLK1)) + geom_boxplot()  + guides(x=guide_axis(angle=45)) + stat_compare_means(method='anova') + ggtitle(paste0("PLK1 expression (",id,")"),experimentData(HCVdats[[id]][[1]])@title) | data.frame(sid=pData(HCVdats[[id]][[1]])[,voi[1]],PLK1=log2(exprs(HCVdats[[id]][[1]]))[grep("PLK1",fData(HCVdats[[id]][[1]])[,voi[2]]),],hcc=pData(HCVdats[[id]][[1]])[,"hcc:ch1"])  %>% ggplot(aes(x=hcc,color=hcc,y=PLK1)) + geom_boxplot()  + guides(x=guide_axis(angle=45)) + stat_compare_means(method='anova') + ggtitle(paste0("PLK1 expression (",id,")"),experimentData(HCVdats[[id]][[1]])@title) | data.frame(sid=pData(HCVdats[[id]][[1]])[,voi[1]],PLK1=log2(exprs(HCVdats[[id]][[1]]))[grep("PLK1",fData(HCVdats[[id]][[1]])[,voi[2]]),],platelet_lessthan_100K=pData(HCVdats[[id]][[1]])[,"platelet < 100,000/mm3:ch1"])  %>% ggplot(aes(x=sid,color=platelet_lessthan_100K,y=PLK1)) + geom_boxplot()  + guides(x=guide_axis(angle=45)) + stat_compare_means(method='anova') + ggtitle(paste0("PLK1 expression (",id,")"),experimentData(HCVdats[[id]][[1]])@title) 

p3<-data.frame(sid=pData(HCVdats[[id]][[1]])[,voi[1]],PLK1=log2(exprs(HCVdats[[id]][[1]]))[grep("PLK1",fData(HCVdats[[id]][[1]])[,voi[2]]),],hcc=pData(HCVdats[[id]][[1]])[,"hcc:ch1"])  %>% ggplot(aes(x=sid,y=PLK1)) + geom_boxplot()  + geom_jitter(width=.1) + guides(x=guide_axis(angle=45)) + stat_compare_means(method='anova') + ggtitle(paste0("PLK1 expression (",id,")"),experimentData(HCVdats[[id]][[1]])@title) + theme_bw(base_size=12) + ylab("PLK Expression (Log intensity)") + xlab("")

keepcol=colnames(pData(HCVdats[[id]][[1]]))[c(45:57)]
	resdat<-data.frame(pData(HCVdats[[id]][[1]])[,keepcol],PLK=log2(exprs(HCVdats[[id]][[1]]))[grep("PLK1",fData(HCVdats[[id]][[1]])[,voi[2]]),])

 fit<-lm(PLK~bilirubin...1.0mg.dl.ch1+child.ch1+death.ch1+decomp.ch1+hcc.ch1+platelet...100.000.mm3.ch1+prediction.ch1,data=resdat)
 summary(fit)

 fit2<-lm(PLK~hcc.ch1+platelet...100.000.mm3.ch1+prediction.ch1,data=resdat)
 summary(fit2)

ggplot(resdat,aes(x=PLK,y=as.integer(days.to.hcc.ch1),color=prediction.ch1)) + geom_point() + geom_smooth(method="lm")


 HBVdats<-readRDS("r_objects/211112_11_HBV_GEOdata.RDS")

id=names(HBVdats)[1]
voi<-c("source_name_ch1","Gene Symbol")
p4<-data.frame(sid=pData(HBVdats[[id]][[1]])[,voi[1]],pData(HBVdats[[id]][[1]])[,c('alt:ch1','ast:ch1','hbv-dna:ch1')],PLK1=exprs(HBVdats[[id]][[1]])[grep("PLK1",fData(HBVdats[[id]][[1]])[,voi[2]])[1],]) %>% mutate(id=unlist(lapply(strsplit(sid,"\\-"), function(X) X[1]))) %>%  ggplot(aes(x=id,y=PLK1)) + geom_boxplot() + geom_jitter(width=.1,inherit.aes=F,aes(x=id,y=PLK1)) + stat_compare_means(method='anova')+ guides(x=guide_axis(angle=45)) + ggtitle(paste0("PLK1 expression (",id,")"),experimentData(HBVdats[[id]][[1]])@title) + scale_fill_manual(values=c("white",'grey','grey30')) + theme(legend.position=0) + theme_bw(base_size=12) + ylab("PLK Expression (Log intensity)") + xlab("")

p1 + p2 + p3 + p4 +plot_layout(nrow=2)
ggsave('reports/220303_PaperImage.pdf')