library(tidyverse)
library(limma)
library(edgeR)
 library(affyPLM)
 library(data.table)
 gdat<-readRDS("../ruby/data/mirSmash/IDtab.RDS")


# Mark Read 
#/Users/brian/Dropbox (Sydney Uni)/projects/mark/source/211101_DataProcessing.r

ddir<-"../mark/r_objects/"

HCVdats=readRDS("../mark/r_objects/211116_09_HCVall_GEOdata.RDS")
naf=readRDS("../mark/r_objects/211119_13_NAFLD_GEOdata.RDS")
hcc=readRDS("../mark/r_objects/211112_12_HCC_GEOdata.RDS")
hbv=readRDS("../mark/r_objects/211112_11_HBV_GEOdata.RDS")


# scott 
#scott/GEOstuff/210610_drunkenLiver/210707_datawrangling.r

alcLiv=readRDS("/Users/brian/Dropbox (Sydney Uni)/projects/scott/GEOstuff/210610_drunkenLiver/210707_AlcLiv_CleanData.RDS")

#'~/Dropbox (Sydney Uni)/projects/scott/GEOstuff/forTheKiddies/210809_dataWrangle.r'

load('~/Dropbox (Sydney Uni)/projects/scott/GEOstuff/forTheKiddies/210809_processedData.rdata')

d=GSE135251 NALFD,seq
Y=GSE48452 NAFLD,array


#hcv
 names(HCVdats)

#GSE20119 is crap
#GSE168186 is Seq (GSE168178)
#GSE54100 is crap

# workspace
i=14
ID=names(HCVdats)[i]
names(HCVdats)[i] %in% ls()
dim(HCVdats[[i]][[1]])
plotDensities(HCVdats[[i]][[1]])
plotDensities(log2(exprs(HCVdats[[i]][[1]])))
plotDensities(log2(exprs(normalize(HCVdats[[i]][[1]])))

ex=log2(exprs(normalize(HCVdats[[i]][[1]])))

ga<-data.frame(fData(HCVdats[[i]][[1]]))
sa<-data.frame(pData(HCVdats[[i]][[1]]))
ti=experimentData(HCVdats[[i]][[1]])@title

head(ga,2)
gaid="Gene.Symbol"
head(unique(ga[,gaid]))

head(sa,2)
said=c("characteristics_ch1")
unique(sa[,said])
saindex<-c(2:1)



tmp<- mutate(data.frame(ex),gene_symbol=make.unique(ga[,gaid])) %>%  
 filter(gene_symbol %in% gdat$hgnc_symbol & !gene_symbol=="" & !is.na(gene_symbol)) %>% 
 pivot_longer(cols=-gene_symbol,values_to="Log_Expr",names_to="sample_id") %>% 
 left_join(data.frame(sample_id=rownames(sa),vars=sa[,said])) %>% 
 mutate(vars=factor(vars,levels=unique(vars)[saindex]))

 #de
expdat= mutate(data.frame(ex),gene_symbol=make.unique(ga[,gaid])) %>%  filter(gene_symbol %in% gdat$hgnc_symbol & !gene_symbol=="" & !is.na(gene_symbol))  %>% data.frame(row.names=NULL) %>% column_to_rownames("gene_symbol")

cond= paste0(sa[,said])
cond
unique(cond)
cond<-factor(cond,levels=unique(cond)[saindex])
cond

i
des<-model.matrix(~cond)
fit<-lmFit(expdat,des)
fit2<-eBayes(fit)
out<-lapply(2:ncol(des), function(i) topTable(fit2,coef=i,n=100000))
names(out)<-make.names(colnames(des)[-1])
	for(x in 1:length(out)) {
		out[[x]]$list=names(out)[x]
		out[[x]]$gene_symbol=rownames(out[[x]])
	}
deres<-rbindlist(out[-3])
i
deres$study=names(HCVdats)[i]

names(HCVdats)[i]


 =list(
	pl= list(
		expdat= tmp,
		title= ti
		),
	de=deres
)
testplot( )
rm(sa,ga,ex,tmp,ti,i,said,gaid,saindex)


GSE15331<-list(
	pl=list(
		expdat= mutate(data.frame(exprs(HCVdats$GSE15331[[1]])),gene_symbol=make.unique(fData(HCVdats$GSE15331[[1]])$GENE_SYMBOL)) %>%  filter(gene_symbol %in% gdat$hgnc_symbol & !gene_symbol=="" & !is.na(gene_symbol)) %>% pivot_longer(cols=-gene_symbol,values_to="Log_Expr",names_to="sample_id") %>% left_join(data.frame(sample_id=rownames(pData(HCVdats$GSE15331[[1]])),vars=paste0(pData(HCVdats$GSE15331[[1]])$`sample:ch1`,".v.",pData(HCVdats$GSE15331[[1]])$`sample:ch2`))),
		title=experimentData(HCVdats$GSE15331[[1]])@title
	),
	de=deres
)
testplot(GSE15331)


GSE14323 =list(
	pl= list(
		expdat= tmp,
		title= ti	
		),
	de=deres
)
testplot(GSE14323)
rm(sa,ga,ex,tmp,ti,i)


GSE20119 =list(
	pl= list(
		expdat= tmp,
		title= ti
		),
	de=deres
)
testplot(GSE20119 )
rm(sa,ga,ex,tmp,ti,i)

GSE54648 =list(
	pl= list(
		expdat= tmp,
		title= ti
		),
	de=deres
)
testplot( GSE54648)
rm(sa,ga,ex,tmp,ti,i,said,gaid,saindex)


 GSE70779=list(
	pl= list(
		expdat= tmp,
		title= ti
		),
	de=deres
)
testplot(GSE70779 )
rm(sa,ga,ex,tmp,ti,i,said,gaid,saindex)

GSE54101 =list(
	pl= list(
		expdat= tmp,
		title= ti
		),
	de=deres
)
testplot(GSE54101 )
rm(sa,ga,ex,tmp,ti,i,said,gaid,saindex)


 GSE15654=list(
	pl= list(
		expdat= tmp,
		title= ti
		),
	de=deres
)
testplot( GSE15654)
rm(sa,ga,ex,tmp,ti,i,said,gaid,saindex)

GSE34798 =list(
	pl= list(
		expdat= tmp,
		title= ti
		),
	de=deres
)
testplot(GSE34798 )
rm(sa,ga,ex,tmp,ti,i,said,gaid,saindex)



GSE7123 =list(
	pl= list(
		expdat= tmp,
		title= ti
		),
	de=deres
)
testplot( GSE7123)
rm(sa,ga,ex,tmp,ti,i,said,gaid,saindex)

id<-grep("GSE",ls(),value=T)
d[order(as.numeric(gsub("GSE","",id)))]

hcv=list(
	GSE7123=GSE7123,
	GSE14323=GSE14323,
	GSE15331=GSE15331,
	GSE15654=GSE15654,
	GSE20119=GSE20119,
	GSE34798=GSE34798,
	GSE54101=GSE54101,
	GSE54648=GSE54648,
	GSE70779=GSE70779
)
saveRDS(hcv,file="r_objects/220113_HCV_data.RDS")

extGene=function(bigObj,gene="TP53",titleTrim=50,textSize=12){
	require(tidyverse)
	require(ggpubr)
	require(janitor)
	ints<-lapply(bigObj, function(X){
		return(
			list(
				ex=filter(X$pl$expdat,gene_symbol==gene),
				ti=X$pl$title,
				de=filter(X$de,gene_symbol==gene) %>% mutate(title=X$pl$title)
			)
		)
	}
	)
	return(
		list(
			expPlot=lapply(ints, function(Y){
				ggplot(Y$ex,aes(x=vars,color=vars,y=Log_Expr)) + geom_boxplot(outlier.shape=NA)+ geom_jitter(width=.1) + ggtitle(paste(gene,":", Y$de$study[1]),substring(Y$ti,1,titleTrim)) + theme_minimal(base_size=textSize) + ylab("Expression (log2 intensity)") + stat_compare_means(method='anova') + guides(x=guide_axis(angle=20))
			}),
			dePlot=rbindlist(lapply(ints,function(Z) Z$de)) %>% ggplot(aes(x=list,y=logFC,fill=-log10(adj.P.Val))) + geom_bar(stat='identity')+ facet_wrap(~study,nrow=1,scale="free_x") + scale_fill_viridis_c(option="H") + ggtitle(gene, "foldchange in basic comparisons") + theme_bw(base_size=textSize+4) + geom_hline(yintercept=0) + guides(x=guide_axis(angle=20))
		)
	)
}

p<-extGene(gene="PLK1",titleTrim=100)

 cowplot::plot_grid(plotlist=p$expPlot)
 p$de

 saveRDS(extGene,file="r_objects/220113_LiverPlotTool.RDS")