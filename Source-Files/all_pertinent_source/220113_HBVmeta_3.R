library(tidyverse)
library(limma)
library(edgeR)
 library(affyPLM)
 library(data.table)
 library(patchwork)

  ##functions
  load(file="r_objects/220114_helperFunctions.rdata")

  #data
  indat=readRDS("../mark/r_objects/211112_11_HBV_GEOdata.RDS")
length(indat)

# workspace
i=3
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
gaid="Gene.Symbol"
head(unique(ga[,gaid]))

head(sa,2)
said=c("varoi")
unique(sa[,said])
saindex<-c(1:2)
unique(sa[,said])[saindex]
table(sa[,said])


#done
 #GSE83148
 #GSE84044
 #GSE66698
 #
 #

names(indat)[i]
<-mkobj()
testplot( )
rm(sa,ga,ex,ti,i,said,gaid,saindex)

HBV=list(
	GSE66698=GSE66698,
	GSE83148=GSE83148,
	GSE84044=GSE84044
)

pl<-readRDS("r_objects/220113_LiverPlotTool.RDS")
p<-pl(HBV,"PLK1")
cowplot::plot_grid(plotlist=p$expPlot) | p$dePlot
# manual
saveRDS(HBV,file="r_objects/220114_HBV_data.RDS")