BiocManager::install("org.Rn.eg.db", version = "3.8")
BiocManager::install("org.Hs.eg.db", version = "3.8")
BiocManager::install("clusterProfiler", version = "3.8")
BiocManager::install("rat2302.db", version = "3.8")
BiocManager::install("AnnotationDbi", version = "3.8")
BiocManager::install("clusterProfiler", version = "3.8")


# next --------------------------------------------------------------------
library("clusterProfiler") 
library("AnnotationDbi")
library("rat2302.db") 
library("DescTools")
library(org.Hs.eg.db)
library(org.Rn.eg.db)

setwd("/Users/onurtidin/Desktop/desktopvirtual/SIBEnrichmentAnalysis/exercise/dataExercise2")

# initial -----------------------------------------------------------------
require(clusterProfiler) 
rat <- read.table("rat_KD.txt", sep = "\t", header = T,stringsAsFactors=FALSE) 
pathway<-read.csv(file="Pathway.csv", stringsAsFactor=FALSE, header=TRUE)
ttestRat <- function(df, grp1, grp2) {
x = df[grp1]
y = df[grp2]
x = as.numeric(x)
y = as.numeric(y)
results = t.test(x, y)
results$p.value} 


# set score ---------------------------------------------------------------
# set score (those you get from a t-test or any other statistical test)
rawp <- apply(rat, 1, ttestRat, grp1 = c(2:7), grp2 = c(8:12))
names(rawp)	<-rat$row.names
sortedrawp	<-sort(rawp)
p_holm 		  <-p.adjust(sortedrawp,method="BH") 
names(p_holm)	<-names(sortedrawp)

SCORE		<-p_holm 
SCORE		<-sort(SCORE,decreasing=TRUE) 
head(SCORE) 

# get phenotype (term)
term2gene	<-data.frame(term=pathway$pathways,name=pathway$row.names)
head(term2gene) 



# next --------------------------------------------------------------------

# run GSEA

gsea.out <- GSEA(SCORE, TERM2GENE=term2gene, nPerm=10000, pvalueCutoff=1, pAdjustMethod = "BH")
gseaplot(gsea.out,"yes")
summary(gsea.out)
head(gsea.out)


# gene names --------------------------------------------------------------------

library(data.table)
library("annotate")
library("rat2302.db")
library("DescTools")
PROBES<- rat$row.names
OUT<-select(rat2302.db,keys= PROBES, columns=c("SYMBOL", "ENTREZID", "GENENAME"))
#ribosomal<-OUT[which(OUT$GENENAME %like% "ribosomal protein"),]
ribosomal<-OUT[grep("ribosomal protein",OUT$GENENAME),]
#ubiquitin<-OUT[whi˚ch(OUT$GENENAME %like% "ubiquitin-conjugating"),]


select(rat2302.db, c("1368587_at","1385248_a_at"), c("SYMBOL","ENTREZID", "GENENAME"))


# set score (those you get from a t-test or any other statistical test)
rawp <- apply(rat, 1, ttestRat, grp1 = c(2:7), grp2 = c(8:12))
names(rawp)<-rat$row.names
sortedrawp<-sort(rawp)
p_holm <- p.adjust(sortedrawp,method="BH")
names(p_holm)<-names(sortedrawp)
SCORE<-p_holm
SCORE=sort(SCORE,decreasing=TRUE)
head(SCORE)

# get phenotype (term)
term2gene<-data.frame(term="no",name=rat$row.names,stringsAsFactors=FALSE)
term2gene[which(term2gene$name %in% ribosomal$PROBEID),1]<-"yes"
head(term2gene)

# run GSEA
gsea.out<-GSEA(SCORE, TERM2GENE=term2gene, nPerm=10000, pvalueCutoff=1, pAdjustMethod = "BH")
gseaplot(gsea.out, "yes")
head(gsea.out)
summary(gsea.out)

