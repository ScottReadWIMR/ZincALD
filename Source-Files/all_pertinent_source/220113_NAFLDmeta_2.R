library(tidyverse)
library(limma)
library(edgeR)
 library(affyPLM)
 library(data.table)
 library(patchwork)
 gdat<-readRDS("../ruby/data/mirSmash/IDtab.RDS")

 ##functions

mkobj<-function(){
	require(tidyverse)
	require(limma)
#var names for future ex=ex,ga=ga,gaid=gaid,sa=sa,said=said,saindex=saindex
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
deres<-rbindlist(out) %>% mutate(list=factor(paste0(gsub("cond|cond_|x_cond","",list),".vsCont"),levels=unique(paste0(gsub("cond|cond_|x_cond","",list),".vsCont"))))
deres$study=ID

return(
	list(
	pl= list(
		expdat= tmp,
		title= ti
		),
	de=deres
))
}

testplot<-function(obj){
	require(tidyverse)
	require(patchwork)
print(
(obj$pl$expdat %>% filter(gene_symbol==obj$de$gene_symbol[1]) %>% ggplot(aes(x=vars,y=Log_Expr,color=vars)) + geom_boxplot() + geom_jitter()) | (filter(obj$de, gene_symbol==obj$de$gene_symbol[1]) %>% ggplot(aes(x=list,y=logFC)) + geom_bar(stat="identity"))
)
}

#mark
 indat=readRDS("../mark/r_objects/211119_13_NAFLD_GEOdata.RDS")
#scott 
 load('~/Dropbox (Sydney Uni)/projects/scott/GEOstuff/forTheKiddies/210809_processedData.rdata')
#d=GSE135251 NALFD,seq
#Y=GSE48452 NAFLD,array

#add array to inlist
indat$GSE48452=list(ESet=Y)
indat$E_MEXP_3291=list(ESet=indat$E_MEXP_3291)

#GSE834521 doesnt exist
# workspace
i=5
names(indat)[i]
names(indat)[i] %in% ls()
dim(indat[[i]][[1]])
plotDensities(indat[[i]][[1]])
plotDensities(log2(exprs(indat[[i]][[1]])))
plotDensities(log2(exprs(normalize(indat[[i]][[1]]))))

ex=exprs(indat[[i]][[1]])
plotDensities(ex)

ga<-data.frame(fData(indat[[i]][[1]]))
sa<-data.frame(pData(indat[[i]][[1]]))
ti=experimentData(indat[[i]][[1]])@title

head(ga,2)
gaid="lu"
head(unique(ga[,gaid]))

head(sa,2)
said=c("group.ch1")
unique(sa[,said])
saindex<-c(1,3,4,2)


#done
 #GSE66676
 
 #GSE49541
 #E_MEXP_3291
 #GSE48452

ID=names(indat)[i]
ID
GSE48452<-mkobj()
testplot( GSE48452)
rm(sa,ga,ex,ti,i,said,gaid,saindex)


NAFLD<-list(
	E_MEXP_3291=E_MEXP_3291,
	GSE48452=GSE48452,
	GSE49541=GSE49541,
	GSE66676=GSE66676
)
NAFLD<-NAFLD[order(names(NAFLD))]

pl<-readRDS("r_objects/220113_LiverPlotTool.RDS")
p<-pl(NAFLD,"TP53")
cowplot::plot_grid(plotlist=p$expPlot) / p$dePlot

#saving funcs
saveRDS(NAFLD,file="r_objects/220114_NAFLD_data.RDS")
save(mkobj,testplot,file="r_objects/220114_helperFunctions.rdata")

# manual
tmp<- mutate(data.frame(ex),gene_symbol=make.unique(ga[,gaid])) %>%  
 filter(gene_symbol %in% gdat$hgnc_symbol & !gene_symbol=="" & !is.na(gene_symbol)) %>% 
 pivot_longer(cols=-gene_symbol,values_to="Log_Expr",names_to="sample_id") %>% 
 left_join(data.frame(sample_id=make.names(rownames(sa)),vars=sa[,said])) %>% 
 mutate(vars=factor(vars,levels=unique(vars)[saindex]))

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
deres<-rbindlist(out)
deres$study=names(indat)[i]
return(
	list(
	pl= list(
		expdat= tmp,
		title= ti
		),
	de=deres
))
}

 =list(
	pl= list(
		expdat= tmp,
		title= ti
		),
	de=deres
)
testplot( )
rm(sa,ga,ex,tmp,ti,i,said,gaid,saindex)
 