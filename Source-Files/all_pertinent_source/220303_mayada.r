220303_mayada.r

suppressPackageStartupMessages(library(shiny))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(patchwork))
suppressPackageStartupMessages(library(DT))
NAFLD=readRDS(file="~/projects/liver_meta/r_objects/220119_NAFLD_data.RDS")
pl=readRDS(file="~/projects/liver_meta/r_objects/220113_LiverPlotTool.RDS")

poi<-c("E_MEXP_3291",'GSE46300','GSE66676','PRJNA512027')
gene="RPLP0"
inp<-pl(NAFLD[poi],gene,1000)

p1<-inp$expPlot$E_MEXP_3291$data %>% filter(as.character(vars) %in% c("normal","Steatosis")) %>%  ggplot(aes(x=vars,color=vars,y=Log_Expr)) + geom_boxplot(outlier.shape=NA)+ geom_jitter(width=.1) + theme_minimal(base_size=12) + ylab("Expression (log2 intensity)") + stat_compare_means(method='anova') + guides(x=guide_axis(angle=20)) + theme(legend.position = 0)  + ggtitle(inp$expPlot$E_MEXP_3291$labels$title)


p2<-inp$expPlot$GSE46300$data %>% filter(as.character(vars) %in% c("low","high")) %>%  ggplot(aes(x=vars,color=vars,y=Log_Expr)) + geom_boxplot(outlier.shape=NA)+ geom_jitter(width=.1) + theme_minimal(base_size=12) + ylab("Expression (log2 intensity)") + stat_compare_means(method='anova') + guides(x=guide_axis(angle=20)) + theme(legend.position = 0)  + ggtitle(inp$expPlot$GSE46300$labels$title)



p3<-inp$expPlot$GSE66676$data %>% filter(as.character(vars) %in% c("No NAFLD","NAFLD not NASH")) %>%  ggplot(aes(x=vars,color=vars,y=Log_Expr)) + geom_boxplot(outlier.shape=NA)+ geom_jitter(width=.1) + theme_minimal(base_size=12) + ylab("Expression (log2 intensity)") + stat_compare_means(method='anova') + guides(x=guide_axis(angle=20)) + theme(legend.position = 0)  + ggtitle(inp$expPlot$GSE66676$labels$title)


p4<-inp$expPlot$PRJNA512027$data %>% filter(as.character(vars) %in% c("NORMAL","STEATOSIS 2","STEATOSIS 3")) %>% mutate(vars=ifelse(vars=="NORMAL","Normal","Steatosis 2 & 3 combined")) %>%  ggplot(aes(x=vars,color=vars,y=Log_Expr)) + geom_boxplot(outlier.shape=NA)+ geom_jitter(width=.1) + theme_minimal(base_size=12) + ylab("Expression (log2 intensity)") + stat_compare_means(method='anova') + guides(x=guide_axis(angle=20)) + theme(legend.position = 0)  + ggtitle(inp$expPlot$PRJNA512027$labels$title)

p1+p2+p3+p4+plot_layout(nrow=1)
ggsave("reports/220303_mayada_v1.pdf")
dat<-list(
	E_MEXP_3291=p1$data,
	GSE46300=p2$data,
	GSE66676=p3$data,
	PRJNA512027=p4$data
)
openxlsx::write.xlsx(dat,"reports/220303_mayada_v1.xlsx")


