
library(tidyverse)
library(limma)
library(edgeR)
 library(affyPLM)
 library(data.table)
 library(patchwork)
 gdat<-readRDS("../ruby/data/mirSmash/IDtab.RDS")

  ##functions
  load(file="r_objects/220114_helperFunctions.rdata")

 #data
indat=readRDS("../mark/r_objects/211112_12_HCC_GEOdata.RDS")
length(indat)

#GSE94660 doesnt exist, probably seq

#workspace
i=6
names(indat)[i]
names(indat)[i] %in% ls()
dim(indat[[i]][[1]])
plotDensities(exprs(indat[[i]][[1]]))
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
said=c("scat")
unique(sa[,said])
saindex<-c(2,1)
unique(sa[,said])[saindex]
table(sa[,said])


#done
 #GSE44074
 #GSE121248
 #GSE84402
 #GSE25097
 #GSE14520

names(indat)[i]
GSE14520 <-mkobj()
testplot(GSE14520)
rm(sa,ga,ex,ti,i,said,gaid,saindex)


HCC=list(
	GSE14520=GSE14520,
	GSE25097=GSE25097,
	GSE44074=GSE44074,
	GSE84402=GSE84402,
	GSE121248=GSE121248
	)

pl<-readRDS("r_objects/220113_LiverPlotTool.RDS")
p<-pl(HCC,"PLK1")
cowplot::plot_grid(plotlist=p$expPlot,nrow=1) / p$dePlot
# manual
saveRDS(HCC,file="r_objects/220114_HCC_data.RDS")


### noticed that they didnt line up and it annoyed me so i made the following to make sure that it works better, cleared backlog and changed mkobj in NAFLDto reflect this
fixList<-function(bigObj){
	bigObj<-bigObj[order(names(bigObj))]
	for(i in 1:length(bigObj)) bigObj[[i]]$de = mutate(bigObj[[i]]$de,list=factor(paste0(gsub("cond|cond_|x_cond","",list),".vsCont"),levels=unique(paste0(gsub("cond|cond_|x_cond","",list),".vsCont"))))
	return(bigObj)
	}

HCC2<-fixList(HCC)
p2<-extGene(HCC2,"PLK1")
cowplot::plot_grid(plotlist=p2$expPlot,nrow=1) / p2$dePlot
saveRDS(HCC2,file="r_objects/220114_HCC_data.RDS")

 dir("r_objects")
220113_HCV_data.RDS: done
220114_HBV_data: done


dat1<-readRDS("r_objects/220113_HCV_data.RDS")
dev.set(2)
p1<-extGene(dat1,"PLK1")
cowplot::plot_grid(plotlist=p1$expPlot,nrow=3) | p1$dePlot
dat2<-fixList(dat1)
dev.set(3)
p2<-extGene(dat2,"PLK1")
cowplot::plot_grid(plotlist=p2$expPlot,nrow=3) | p2$dePlot
saveRDS(dat2,file="r_objects/220114_HCC_data.RDS")
rm(dat2,dat1)
dir("r_objects")

