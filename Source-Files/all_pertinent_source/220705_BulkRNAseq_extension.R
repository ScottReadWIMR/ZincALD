

#setup

expdat<-readRDS(file="r_objects/220630_Complete_proc.RDS")

doMerPlot<-function(obj,varname=NA,subt=NULL){
	require(tidyverse)
	require(ggpubr)
	if(is.na(varname)) {
	print(head(colnames(obj$md)))
	stop()
	}
	p<-data.frame(obj$md,MERTK=obj$expr["MERTK",]) %>% ggplot(aes(x=eval(parse(text=varname)),y=MERTK,color=eval(parse(text=varname)))) + geom_boxplot() + geom_jitter(width=.1,color='black',size=.2) + xlab(NULL) + ggtitle(obj$title,subt) + stat_compare_means(method="wilcox",paired=F) + theme_classic(base_size=14) + scale_color_viridis_d(option="H") + theme(legend.position=0,plot.title = element_text(size = 12,face="bold"))
	print(p)
	return(invisible(p$data))
}

# confirm mann whitney for GSE66494
doMerPlot(expdat$GSE66494,make.names('disease status:ch1'),"GSE66494")


#GSE110147, we need only control and IPF (first and last group).
expdat$GSE110147$md

keep<-expdat$GSE110147$md$`disease state:ch1` %in% c("Idiopathic pulmonary fibrosis","Normal control")

t_GSE110147<-expdat$GSE110147
t_GSE110147$md<-t_GSE110147$md[keep,]
t_GSE110147$expr<-t_GSE110147$expr[,keep]
table(rownames(t_GSE110147$md)==colnames(t_GSE110147$expr))

doMerPlot(t_GSE110147,make.names('disease state:ch1'),"GSE110147")

# from liver meta GSE135251 (NAFLD), GSE49541 (NAFLD)

 NAFLD=readRDS(file="~/projects/liver_meta/r_objects/220119_NAFLD_data.RDS")

 pl=readRDS(file="~/projects/liver_meta/r_objects/220113_LiverPlotTool.RDS")

 p_GSE135251<-filter(NAFLD$GSE135251$pl$expdat,gene_symbol=="MERTK") %>% ggplot(aes(x=vars,y=Log_Expr,color=vars)) + geom_boxplot() + geom_jitter(width=.1,color='black',size=.2) + xlab(NULL) + ggtitle(NAFLD$GSE135251$pl$title,"GSE135251") + stat_compare_means(method="anova") + theme_classic(base_size=14) + scale_color_viridis_d(option="H") + theme(legend.position=0,plot.title = element_text(size = 12,face="bold")) + ylab('MERTK')
p_GSE135251


p_GSE49541<-filter(NAFLD$GSE49541$pl$expdat,gene_symbol=="MERTK") %>% ggplot(aes(x=vars,y=Log_Expr,color=vars)) + geom_boxplot() + geom_jitter(width=.1,color='black',size=.2) + xlab(NULL) + ggtitle(NAFLD$GSE49541$pl$title,"GSE49541") + stat_compare_means(method="anova") + theme_classic(base_size=14) + scale_color_viridis_d(option="H") + theme(legend.position=0,plot.title = element_text(size = 12,face="bold")) + ylab('MERTK')
p_GSE49541


 nafldat<-list(
	GSE135251=p_GSE135251$data,
	GSE49541=p_GSE49541$data
 )

 openxlsx::write.xlsx(nafldat,file="reports/220706_two_NAFL_datasets.xlsx")