# 116353-1	
#InTEAM Consortium - Alcoholic Hepatitis (phs001807.v1.p1)
#Disease-Specific (Alcoholism, IRB, NPU, RD) (phs001807.v1.p1.c1), JAAMH
#PI: Thomas Tu , Project #31612:  Mutational changes leading to Hepatocellular Carcinoma

cd /scratch/RDS-FMH-Brian-RW/scott

mkdir downloads fastq logs source



# save 36112 key this way
sshfs -o IdentityFile=~/.ssh/id_rsa_artemis bglo7514@hpc.sydney.edu.au:/scratch/RDS-FMH-Brian-RW/scott /Users/briglo/mtpt


cd /scratch/RDS-FMH-Brian-RW/scott
module load sratoolkit
 vdb-config --import prj_31612_D32850.ngc


 /usr/local/asperaconnect/ibm-aspera-connect-3.11.0.5-linux-g2.12-64.sh
 ~/.aspera/connect/bin/ascp -QTr -l 300M -k 1 -i  ~/.aspera/connect/etc/asperaweb_id_dsa.openssh  -W ADE611BB4D17F416A96D6B7746C5A3A2978A030F46BA0428D8A8082D1769051D47E02A0B8B5C34DD6B863D82D52643639D dbtest@gap-upload.ncbi.nlm.nih.gov:data/instant/briangloss/89662 ./downloads

vdb-decrypt 89662

module load sratoolkit
prefetch downloads/cart_DAR110735_202205200000.krt #going

ls  /scratch/DPM_HCC/dbgap_116351/downloads/sra/sra

cd /scratch/RDS-FMH-Brian-RW/scott
module load sratoolkit/3.0.0
 vdb-config --import downloads/prj_31612_D32850.ngc

 vdb-config --ignore-protected-repositories
prefetch --ngc downloads/prj_31612_D32850.ngc --cart downloads/cart_DAR116353_202210122206.krt -O /scratch/RDS-FMH-Brian-RW/scott/downloads/sra


#PBS -P DPM_HCC
#PBS -N dumpnzip
#PBS -j oe
#PBS -l select=1:ncpus=2:mem=12GB
#PBS -l walltime=4:00:00
#PBS -q defaultQ

cd /scratch/RDS-FMH-Brian-RW/scott/downloads/sra/"$snam"

module load sratoolkit/3.0.0

if [ ! -f "$snam".sra ] #; then echo "yay" ; fi
then 
    fastq-dump --split-3 --ngc /scratch/RDS-FMH-Brian-RW/scott/downloads/prj_31612_D32850.ngc -O /scratch/RDS-FMH-Brian-RW/scott/fastq/ "$snam" 
else 
    fastq-dump --split-3 --ngc /scratch/RDS-FMH-Brian-RW/scott/downloads/prj_31612_D32850.ngc -O /scratch/RDS-FMH-Brian-RW/scott/fastq/ "$snam".sra
fi

gzip /scratch/RDS-FMH-Brian-RW/scott/fastq/"$snam"_1*fastq 
 gzip /scratch/RDS-FMH-Brian-RW/scott/fastq/"$snam"_2*fastq # sucks... job finishes regardless of background process... do serially....

#patch with screen -r 106481.pts-37.login3
#for i in `ls *fastq | head -n 19` ; do echo "$i" ; gzip -f "$i" ; done

cd /scratch/RDS-FMH-Brian-RW/scott
# prototype | tail -n +2 | head -n 1 and | tail -n +3 | head -n 1
for i in `cut -d "," -f 1 downloads/SraRunTable.txt | tail -n +4` ; do echo "$i"_unpack >> logs/jerbs.txt ; qsub -o logs/dumpNzip_"$i".log -v "snam="$i"" source/SRAdumpNzip.sh | tee -a logs/jerbs.txt ; done

cd /scratch/RDS-FMH-Brian-RW/scott
for i in `cut -d "," -f 1 downloads/SraRunTable.txt` ; do ls fastq/"$i"* ; done


#SRR8937563 failed
# 2022-10-13T21:02:46 fastq-dump.3.0.0 err: data unexpected while reading file within network system module - Cannot KStblHttpFileTimedReadChunked
# 2022-10-13T21:11:33 fastq-dump.3.0.0 err: item not found while constructing within virtual database module - the path 'SRR8937563' cannot be opened as database or table
# fastq-dump quit with error code 3

# try by itself
cd /scratch/RDS-FMH-Brian-RW/scott
module load sratoolkit/3.0.0
prefetch --ngc downloads/prj_31612_D32850.ngc -O /scratch/RDS-FMH-Brian-RW/scott/downloads SRR8937563

# 2022-10-16T22:00:10 prefetch.3.0.0 err: data unexpected while reading file within network system module - Cannot KStblHttpFileRead
# 2022-10-16T22:01:51 prefetch.3.0.0 err: data unexpected while reading file within network system module - Cannot KStblHttpFileRead
# 2022-10-16T22:01:52 prefetch.3.0.0 fatal: SIGNAL - Segmentation fault

#try local coz https://github.com/ncbi/sra-tools/issues/634
prefetch --ngc ~/mtpt/downloads/prj_31612_D32850.ngc -O ~/Downloads SRR8937563
# same.. emailed sra...
# on the off chance its vpn related do local
prefetch --ngc ~/Downloads/prj_31612_D32850.ngc -O ~/Downloads SRR8937563

#waiting a week
i=SRR8937563
echo "$i"_unpack_postStuff >> logs/jerbs.txt ; qsub -o logs/redo_dumpNzip_"$i".log -v "snam="$i"" source/SRAdumpNzip.sh | tee -a logs/jerbs.txt
#6536924 looks ok, gzip failed KEKW fixrd in post

#! /bin/bash
#PBS -P Brian
#PBS -N star_mapping  
#PBS -j oe                
#PBS -l select=1:ncpus=12:mem=64GB   
#PBS -l walltime=16:00:00
#PBS -q defaultQ
# usage qsub -v "snam=test" runstar.sh 

module load star/2.7.1a
module load trimgalore
module load bowtie2

fastadir=/scratch/RDS-FMH-Brian-RW/scott/fastq
outdir=/scratch/RDS-FMH-Brian-RW/scott
genomedir=/project/RDS-FMH-WIMR_RNA-RW/briglo_tmp/joey/STAR_hg38_ercc/STAR
FQS=/project/RDS-FMH-WIMR_RNA-RW/briglo_tmp/joey/fastq_screen

if [ ! -d "$fastadir"/trimgalore ]; then mkdir "$fastadir"/trimgalore ; fi

trim_galore --paired --fastqc -o "$fastadir"/trimgalore  "$fastadir"/"$snam"_*.gz

"$FQS"/FastQ-Screen-0.14.1/fastq_screen --threads 4 --conf "$FQS"/FastQ_Screen_Genomes/fastq_screen.conf --outdir "$outdir"/fastqscreen --subset 1000000 "$fastadir"/"$snam"_1.fastq.gz

STAR --runMode alignReads \
    --runThreadN 10 \
    --readFilesIn "$fastadir"/trimgalore/"$snam"_*fq.gz \
    --readFilesCommand zcat \
    --outFileNamePrefix "$outdir"/starout/"$snam"_ \
    --genomeDir "$genomedir" \
    --readFilesType Fastx \
    --chimSegmentMin 20 \
    --chimOutType Junctions \
    --chimScoreJunctionNonGTAG -1 \
    --quantMode GeneCounts \
    --outReadsUnmapped Fastx \
    --outSAMtype BAM SortedByCoordinate \
    --bamRemoveDuplicatesType UniqueIdentical 


cd /scratch/RDS-FMH-Brian-RW/scott
# prototype | tail -n +2 | head -n 1  ## Run tail -n +3
for i in `cut -d "," -f 1 downloads/SraRunTable.txt | tail -n +3` ; do echo "$i"_alignCount >> logs/jerbs.txt ; qsub -o logs/alignCount_"$i".log -v "snam="$i"" source/runstar.sh | tee -a logs/jerbs.txt ; done
# 6529530 all g except for sample that failed


#waiting a week
i=SRR8937563
echo "$i"_alignCount2_ >> logs/jerbs.txt ; qsub -o logs/redo2_alignCount_"$i".log -v "snam="$i"" source/runstar.sh | tee -a logs/jerbs.txt
#6537290- didnt nitice gzip failed KEKW redo after manual 6540104

# QC (local)
sshfs -o IdentityFile=~/.ssh/id_rsa_artemis bglo7514@hpc.sydney.edu.au:/scratch/RDS-FMH-Brian-RW/scott /Users/briglo/mtpt

 multiqc ~/mtpt/logs ~/mtpt/fastq/trimgalore/ ~/mtpt/starout/ ~/mtpt/fastqscreen -o reports -n 221026_phs001807_qc




#R
library(readxl)    
library(tidyverse)
library(edgeR)
library(org.Hs.eg.db)
library(clusterProfiler)
library(Glimma)

cd("~/mtpt/starout")

fnam<-dir(pattern="ReadsPerGene.out.tab")
rdat<-lapply(fnam, function(X) data.frame(read.table(X),row.names=1))
names(rdat)<-paste0(gsub('_ReadsPerGene.out.tab','',fnam))

cd("~/projects/scott/")
saveRDS(rdat,"r_objects/221026_dbGAP_rawData.RDS")

# looks like stranded reverse
counts1<-do.call(cbind,lapply(rdat,function(X) X[,3]))
rownames(counts1)<-rownames(rdat[[1]])
#trim for failed sample
table(colSums(counts1)>0)
counts1<-counts1[,colSums(counts1)>0]

noint = rownames (counts1) %in% c("N_unmapped","N_multimapping","N_ambiguous",'N_noFeature')
cpms = cpm(counts1)
keep = rowSums(cpms >1) >=1 & !noint
counts = counts1[keep,]

#gene annotation
 gen<-readRDS('../joey/data/hg38_ercc_GRanges.RDS')
gen
names(gen)<-gen$gene_id
geneInfo<-data.frame(gen[rownames(counts)])[,c(1:3,5,10,12,14)]

hasEntrez<-bitr(geneInfo$gene_id, fromType="ENSEMBL", toType=c("ENTREZID", "SYMBOL"), OrgDb="org.Hs.eg.db")
 geneInfo$ENTREZID<-hasEntrez[match(geneInfo$gene_id,hasEntrez$ENSEMBL),"ENTREZID"]
geneInfo$manualBiotype=ifelse(geneInfo$gene_biotype=="protein_coding","protein_coding",ifelse(geneInfo$gene_biotype=="lncRNA",'lncRNA',ifelse(grepl("TR|IG",geneInfo$gene_biotype),"immune",ifelse(grepl("pseudogene",geneInfo$gene_biotype),"pseudo",ifelse(grepl("rRNA",geneInfo$gene_biotype),"rRNA",ifelse(grepl("miRNA|snoRNA|snRNA",geneInfo$gene_biotype),"small_rna","other"))))))


#sample annotation
read_excel_allsheets <- function(filename, tibble = FALSE) {
    # I prefer straight data.frames
    # but if you like tidyverse tibbles (the default with read_excel)
    # then just pass tibble = TRUE
    sheets <- readxl::excel_sheets(filename)
    x <- lapply(sheets, function(X) readxl::read_excel(filename, sheet = X))
    if(!tibble) x <- lapply(x, as.data.frame)
    names(x) <- sheets
    x
}

phe<-read_excel_allsheets("data/dbgap/phenoData.xlsx")
sample_anno<-left_join(phe[[1]],phe[[2]]) %>% left_join(phe[[3]]) %>% column_to_rownames("Run") %>% mutate(study_disease=replace_na(study_disease,"control"),study_disease=gsub("\\\\, ","_",study_disease)) 
%>% filter(rownames(.) %in% colnames(counts))

#make DGElist for edgeR
table(colnames(counts) %in% rownames(sample_anno))
table(rownames(counts) == geneInfo$gene_id)
rownames(counts)<-rownames(geneInfo) <- make.unique(make.names(geneInfo$gene_name))

d = DGEList(counts = counts[,rownames(sample_anno)], group = sample_anno$study_disease,samples=sample_anno,genes=geneInfo)
d = calcNormFactors(d)

glMDSPlot(d,groups=d$samples,folder='reports',main='dbGAP',html="221026_dbGAP_rnaSeq")
saveRDS(d, file="r_objects/221018_dbGAP_rnaSeqData.RDS")

forscot<-cbind(d$samples,t(cpm(d,log=T)))

openxlsx::write.xlsx(list(expr=forscot,gene_ano=d$genes),file="reports/rnaseq/221026_dbGAP_rnaSeq_expr.xlsx")


cond<-factor(d$samples$AffectionStatus)

#cond<-factor(d$samples$study_disease)

design<-model.matrix(~cond)
y<-estimateDisp(d,design)
fit <- glmQLFit(y, design)


qlf<-lapply(2:7, function(X) glmQLFTest(fit,coef = X))
names(qlf)=paste0("A_S_",2:7)
att<-lapply(qlf, function(X) data.frame(topTags(X,n=1000000)))


tt<-lapply(att, function(X) subset(X, FDR<0.05,abs(logFC)>1))
am<-unique(do.call(c,lapply(tt, function(X) rownames(X))))
df<-data.frame(row.names=am,do.call(cbind,lapply(tt, function(X) ifelse(am %in% rownames(X),1,0))))
library(eulerr)
plot(euler(df,shape="ellipse"),quantities=T)

ggplot(forscot,aes(AFFECTION_STATUS,AL138963.4)) + geom_boxplot()


att$alcHep_vs_cont<-to # set study_disease as cond
openxlsx::write.xlsx(att,file="reports/rnaseq/221026_dbGAP_rnaSeq_DGE_AFFECTION.xlsx",asTable=T,rowNames=T)
