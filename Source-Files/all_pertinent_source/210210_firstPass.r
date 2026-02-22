library(RGSEA)
library(GEOquery)
library(limma)
library(tidyverse)

# #get studies of interest only run once
GSEs<-c(
	"GSE76510",
	"GSE2964",
	"GSE2111",
	"GSE6960",
	"GSE60159",
	"GSE25167",
	"GSE152703",
	"GSE99435"
)

 rdat<-lapply(unique(GSEs), getGEO)
names(rdat)<-unique(GSEs)
saveRDS(rdat,file="210210_zinc/210210_ZN_rawData.RDS")


#make esets
if (!exists("rdat")) rdat<-readRDS("210210_zinc/210210_ZN_rawData.RDS")
esets<-lapply(rdat[-7], function(X) {
	if (max(exprs(X[[1]]))>100)  {
		exprs(X[[1]])<-log2(exprs(X[[1]]))
		return(affyPLM::normalize.ExpressionSet.quantiles(X[[1]])) } else {
	affyPLM::normalize.ExpressionSet.quantiles(X[[1]])
		}
})



library(readxl)    
read_excel_allsheets <- function(filename, tibble = FALSE) {
  # I prefer straight data.frames
  # but if you like tidyverse tibbles (the default with read_excel)
  # then just pass tibble = TRUE
  sheets <- readxl::excel_sheets(filename)
  x <- lapply(sheets, function(X) readxl::read_excel(filename, sheet = X,col_names=))
  if(!tibble) x <- lapply(x, as.data.frame)
  names(x) <- sheets
  x
}
an<-read_excel_allsheets("210210_zinc/210210_esetAno.xlsx")
an<-lapply(an,function(x) {
	print(x)
	data.frame(x, stringsAsFactors=F,row.names=1)
})


tts<-lapply(names(esets[-6]), function(X){
	message(X)
tmp<-esets[[X]][,rownames(an[[X]])]
fData(tmp)<-fData(tmp)#[,c(1:2,grep("symbol|title|name,ILMN_Gene",colnames(fData(tmp)),ignore.case=T))]
des<-model.matrix(~0+cond,data=an[[X]])
cont<-makeContrasts(condZinc-condControl,levels=des)
f<-lmFit(tmp,des)
f2<-contrasts.fit(f,cont)
f3<-eBayes(f2)
return(list(tt=topTable(f3,number=100000,p.value=1),dat=exprs(tmp),an=an[[X]],title= experimentData(tmp)@title))
})
names(tts)<-names(esets[-6])

X="GSE25167"
tmp<-esets[[X]][,rownames(an[[X]])]
fData(tmp)<-fData(tmp)#[,c(1:2,grep("symbol|title|name",colnames(fData(tmp)),ignore.case=T))]
des<-model.matrix(~line+cond,data=an[[X]])
f<-lmFit(tmp,des)
f3<-eBayes(f)

tts$GSE25167=list(tt=topTable(f3,number=100000,p.value=1,coef=3),dat=exprs(tmp),an=an[[X]],title=experimentData(tmp)@title)

# GSE152703 failed, will need to get manually
#11th: its a custom lawl DL matrix from GEO page
GSE152703<-read.table("GSE152703_MTF1_Data_Matrix.txt",header=T,sep="\t",stringsAsFactors=F)
gan<-data.frame(do.call(rbind,strsplit(GSE152703[,1],"_")),row.names=GSE152703[,1])
counts<-data.frame(GSE152703,row.names=1)
colnames(gan)<- c("SYMBOL","NOTSURE")

anno<-data.frame(pData(rdat$GSE152703[[1]])[,c(1,2,37,38)],row.names=gsub(" ",".",gsub("rep|_P$","",unlist(lapply(strsplit(as.character(pData(rdat$GSE152703[[1]])$title),"\\-"), function(X) paste0(X[1],".",X[2]))))))

table(rownames(anno) == colnames(counts))
library(edgeR)
d<-DGEList(counts=counts,genes=gan,samples=anno)
skeep<-c("GSM4624109",
"GSM4624110",	
"GSM4624111",	
"GSM4624115",
"GSM4624116",
"GSM4624117")

y<-d[,rownames(d$samples)[d$samples$geo_accession %in% skeep]]
cond<-factor(gsub("ZnSO4","Zinc",y$samples$chemical.ch1))
rep<-rep(c(1:3),2)
des<-model.matrix(~rep+cond,data=an[[X]])
v<-voom(y,des)
f<-lmFit(v,des)
f2<-eBayes(f)
topTable(f2,coef=3)

tts$GSE152703=list(tt=arrange(topTable(f2,number=100000,p.value=1,coef=3),P.Value),dat=v$E,data.frame(y$samples,cond,rep),title=experimentData(rdat$GSE152703[[1]])@title)

#post talk, do paired for GSE2111

an$GSE2111$rep=c(1:4,1:4)

X="GSE2111"
tmp<-esets[[X]][,rownames(an[[X]])]
fData(tmp)<-fData(tmp)#[,c(1:2,grep("symbol|title|name",colnames(fData(tmp)),ignore.case=T))]
des<-model.matrix(~rep+cond,data=an[[X]])
f<-lmFit(tmp,des)
f3<-eBayes(f)

tts$GSE2111=list(tt=topTable(f3,number=100000,p.value=1,coef=3),dat=exprs(tmp),an=an[[X]],title=experimentData(tmp)@title)

saveRDS(tts,file="210210_zinc/210215_MostDE.RDS")

plotHeat<-function(thingy,num,tit='heatmap'){
print(pheatmap::pheatmap(thingy$dat[rownames(thingy$tt[1:num,]),],annotation_col=thingy$an,main=paste0(tit,thingy$title),annotation_row=thingy$tt[1:num,c('logFC','adj.P.Val','P.Value')],fontsize=6))
}

pdf("210210_zinc/210210_top50CandHeatmaps.pdf",height=8,width=8)
lapply(names(tts),function(X) plotHeat(thingy=tts[[X]],num=50,tit=X))
dev.off()

ttt<-lapply(tts[1:7],function(X) X$tt[X$tt$adj.P.Val<0.05,])
ttt$GSE152703<-tts$GSE152703$tt[tts$GSE152703$tt$P.Value<0.05,]
openxlsx::write.xlsx(ttt,file="210210_zinc/210210_DEG_padj0.05.xlsx")


### for 15/2/21
#dont subset
ttt<-lapply(tts,function(X) X$tt)

#countobs in other studies
symvec<-unlist(lapply(ttt,function(X) grep("symbol",colnames(X), ignore.case=T,value=T)))
symlist<-lapply(1:length(symvec), function(x) ttt[[x]][,symvec[x]])
attt<-lapply(1:length(symvec), function(X) {
	print(X)
tmp<-data.frame(do.call(cbind,lapply(c(1:length(symvec))[-X], function(y) {
	tv<-symlist[[X]]
	tv[tv==""]="NoSymbol"
	tv %in% symlist[[y]]
})))
colnames(tmp)<-paste0("SymbolIn_",names(symvec)[-X])
tmp$total_studies=rowSums(tmp)
return(data.frame(cbind(ttt[[X]],tmp)))
})
names(attt)<-names(symvec)

openxlsx::write.xlsx(attt,file="210210_zinc/210215_DEG_NoCutoff_CountObs.xlsx")

#### doing new dataset
oi<-getGEO("GSE108923")
sdat<-pData(oi[[1]])

# garb, download from geo (FPKM)
GSE108923_exp<-data.frame(read.table('210210_zinc/GSE108923_Expression_data.txt',header=T),row.names=1)

gdat<-GSE108923_exp[,1:5]
edat<-GSE108923_exp[,-c(1:5)]
pdat<-data.frame(row.names=colnames(edat),do.call(rbind,strsplit(colnames(edat),"_")))[,1:2]
colnames(pdat)<-c("condition","replicate")

plotDensities(edat) # holy hell nothing is expressed

#keep<-rowSums(edat>0)>3 #trim by expressed in 3 or more
keep<-rownames(edat) #realised this would make it different to everything else

eset<-ExpressionSet(assayData=data.matrix(edat[keep,]),phenoData=AnnotatedDataFrame(pdat),featureData=AnnotatedDataFrame(gdat[keep,-3]))

des<-model.matrix(~0+condition,data=pData(eset))
cont<-makeContrasts(
	norVSdep=conditionFull-conditionDepleted,
	repVSdep=conditionRepleted-conditionDepleted,
	norVSrep=conditionFull-conditionRepleted,
	levels=des
	)
f<-lmFit(eset,des)
f2<-contrasts.fit(f,cont)
f3<-eBayes(f2)

 res<-decideTests(f3)
 vennDiagram(res,include=c("up","down"),counts.col=c("red","blue"))
 head(res[rowSums(abs(res))==1,])

tts<-lapply(colnames(cont), function(X) topTable(f3,number=100000,p.value=1,coef=X))
names(tts)<-colnames(cont)

GSE108923=list(tt=tts,dat=exprs(eset),an=pData(eset),title=experimentData(oi[[1]])@title )
saveRDS(GSE108923,file="210210_zinc/210222_GSE108923.RDS")

#aside query
pdf("reports/210309_cand_GSE108923_redone.pdf")
par(mar=c(10,3,2,2)) ; barplot(data.matrix(edat[rownames(gdat)[gdat$GeneSymbol=="MT2A"],]),las=2,main="MT2A")
par(mar=c(10,3,2,2)) ; barplot(data.matrix(edat[rownames(gdat)[gdat$GeneSymbol=="MT1X"],]),las=2,main="MT1X")
par(mar=c(10,3,2,2)) ; barplot(data.matrix(edat[rownames(gdat)[gdat$GeneSymbol=="MT1F"],]),las=2,main="MT1F")
dev.off()

#xref to other lists
symvec<-unlist(lapply(ttt,function(X) grep("symbol",colnames(X), ignore.case=T,value=T)))
symlist<-lapply(1:length(symvec), function(x) ttt[[x]][,symvec[x]])

sattt<-lapply(1:length(GSE108923$tt), function(X) {
	print(X)
tmp<-data.frame(do.call(cbind,lapply(c(1:length(symvec)), function(y) {
	tv<-as.character(GSE108923$tt[[X]]$GeneSymbol)
	tv[tv==""]="NoSymbol"
	tv %in% symlist[[y]]
})))
colnames(tmp)<-paste0("SymbolIn_",names(symvec))
tmp$total_studies=rowSums(tmp)
return(data.frame(cbind(GSE108923$tt[[X]],tmp)))
})
names(sattt)<-names(GSE108923$tt)

openxlsx::write.xlsx(sattt,file="210210_zinc/210309_DEG_REDONE_GSE108923_NoCutoff_CountObs.xlsx")


# make plotting method
symexp<-lapply(1:length(symvec), function(x) ttt[[x]][,c(symvec[x],"logFC","adj.P.Val")])
names(symexp)<-names(symvec)
symexp=c(symexp,lapply(tts,function(X) X[,c(1,5,9)]))

for(i in 1:length(symexp)) {
	colnames(symexp[[i]])[1]="SYMBOL"
	symexp[[i]]$list=names(symexp)[i]
}

pldat<-data.table::rbindlist(symexp) %>% mutate(list=factor(list,levels=unique(list)))

plG<-function(symbol='MT1E'){
	require(tidyverse)
	if (length(symbol)>1) {
		print(subset(pldat,SYMBOL %in% symbol) %>% ggplot(aes(x=list,y=logFC,color=adj.P.Val<0.05))  + geom_jitter(width=.1) + ggtitle("multigene") + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) + geom_hline(yintercept=0) + facet_wrap(~SYMBOL,scale='free_y'))
	} else {
	print(subset(pldat,SYMBOL==symbol) %>% ggplot(aes(x=list,y=logFC,color=adj.P.Val<0.05))  + geom_jitter(width=.1) + ggtitle(symbol) + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) + geom_hline(yintercept=0))
	}
}

 save(pldat,plG,file="210210_zinc/210222_plotMeth.rdata")

 plG(grep("^MT1",pldat$SYMBOL,value=T))