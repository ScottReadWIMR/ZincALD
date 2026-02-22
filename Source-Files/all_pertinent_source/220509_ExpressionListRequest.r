suppressPackageStartupMessages(library(shiny))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(patchwork))
suppressPackageStartupMessages(library(DT))

#best look
HBV=readRDS(file="~/projects/liver_meta/r_objects/220114_HBV_data.RDS")$GSE83148
goi=c("HSD17B13","MBOAT7","TMARC1","FGF21","GCKR","MERTK","SERPINA1","TLL1","LEPR","PYGO1","PPP1R3B","LIPA","PCSK7","SOD2","UCP2")
goi<-goi[goi %in% HBV$de$gene_symbol]

ggplot(filter(HBV$de, gene_symbol %in% goi) %>% arrange(gene_symbol) %>% mutate(gene_symbol=factor(gene_symbol, levels=rev(unique(gene_symbol)))), aes(y=logFC,x=gene_symbol,fill=-log10(adj.P.Val))) + geom_bar(stat='identity') + scale_fill_viridis_c(option="H") + coord_flip() + facet_wrap(~study) + theme_classic(base_size=16) + ggplot(filter(HBV$pl$expdat, gene_symbol %in% goi),aes(x=vars,color=vars,y=Log_Expr)) + geom_boxplot() + geom_jitter(width=.1) + facet_wrap(~gene_symbol,scale='free_y',nrow=3) + stat_compare_means(method="wilcox",paired=F) + theme_classic(base_size=16) + guides(x=guide_axis(angle=45)) + theme(legend.position=0) + plot_layout(widths=c(2,6))