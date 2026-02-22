 suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(patchwork))
suppressPackageStartupMessages(library(DT))
suppressPackageStartupMessages(library(ggpubr))
 dat<-list(
    HCV=readRDS(file="~/projects/liver_meta/r_objects/220113_HCV_data.RDS"),
    NAFLD=readRDS(file="~/projects/liver_meta/r_objects/220119_NAFLD_data.RDS"),
    HCC=readRDS(file="~/projects/liver_meta/r_objects/220114_HCC_data.RDS"),
    HBV=readRDS(file="~/projects/liver_meta/r_objects/220114_HBV_data.RDS"),
    ALD=readRDS(file="~/projects/liver_meta/r_objects/220119_AlcLiver_data.RDS")
)

PPARa (PPAR-alpha)
FASN (Fatty acid synthase)
SREBP-1c
SREBP2
HNGCR (HMG CoA reductase)
LDLR (Low density lipoprotein receptor
VLCAD
LDHA (Lactate Dehydrogenase A)
PDK1, PDK2, PDK3 (pyruvate dehydrogenase kinase isoforms)
MTTP (microsomal triglyceride transfer protein)
ACAT1, ACAT2

cand<-c("PPARA","FASN","SREBF1","SREBF2","HMGCR","LDLR","ACADVL","LDHA","PDK1","PDK2","PDK3","MTTP","ACAT1","ACAT2")

cand2<-c("LDHA","PDK1","PDK2","PDK3")

pl<-function(bigObj,gene="TP53",titleTrim=50,textSize=12){
    require(tidyverse)
    require(ggpubr)
    require(janitor)
    ints<-lapply(bigObj, function(X){
        return(
            list(
                ex=dplyr::filter(X$pl$expdat,gene_symbol==gene),
                ti=X$pl$title,
                de=dplyr::filter(X$de,gene_symbol==gene) %>% mutate(title=X$pl$title)
            )
        )
    }
    )
    return(
        list(
            expPlot=lapply(ints, function(Y){
                ggplot(Y$ex,aes(x=vars,color=vars,y=Log_Expr)) + geom_boxplot(outlier.shape=NA)+ geom_jitter(width=.1) + ggtitle(paste(gene,":", Y$de$study[1]),substring(Y$ti,1,titleTrim)) + theme_minimal(base_size=textSize) + ylab("Expression (log2 intensity)") + stat_compare_means(method='anova') + guides(x=guide_axis(angle=20)) + theme(legend.position = 0)
            }),
            dePlot=rbindlist(lapply(ints,function(Z) Z$de)) %>% ggplot(aes(x=list,y=logFC,fill=-log10(adj.P.Val))) + geom_bar(stat='identity')+ facet_wrap(~study,nrow=1,scale="free_x") + scale_fill_viridis_c(option="H") + ggtitle(gene, "foldchange in basic comparisons") + theme_bw(base_size=textSize+4) + geom_hline(yintercept=0) + guides(x=guide_axis(angle=20))
        )
    )
}

pdf(file="plots/231222_mark_HCV_candExp.pdf",width=16,height=16)
lapply(cand, function(X) wrap_plots(pl(dat$HCV,X)[[1]],textSize=10))
dev.off()

pdf(file="plots/231222_mark_HBV_candExp.pdf",width=16,height=6)
lapply(cand, function(X) wrap_plots(pl(dat$HBV,X)[[1]],textSize=10))
dev.off()

lapply(cand2, function(X) wrap_plots(pl(dat$HCV,X)[[1]],textSize=10))
lapply(cand2, function(X) wrap_plots(pl(dat$HBV,X)[[1]],textSize=10))



# for clarifications for paper
pdatC<-lapply(cand2, function(X) pl(dat$HCV,X)[[1]])
pdatB<-lapply(cand2, function(X) pl(dat$HBV,X)[[1]])
names(pdatB)<-cand2


ppdat<-data.table::rbindlist(lapply(c(cand2), function(X) pdatB[[X]]$GSE83148$data))

ggplot(ppdat,aes(vars,Log_Expr)) + geom_boxplot(outlier.shape=NA) + geom_jitter(width=.1,size=.1) + facet_wrap(~gene_symbol,scale="free_y") + ggpubr::stat_compare_means(method="anova",size=2) + ggprism::theme_prism(base_size=10) + guides(x=guide_axis(angle=22))
ggsave("plots/240226_GSE83148_4cand.pdf")