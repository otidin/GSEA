setwd(".../dataExercise3")
#source(https://bioconductor.org/biocLite.R)
biocLite("clusterProfiler")
require(clusterProfiler) 
 
# Get data
table<-read.csv("GSEA_data_input.csv")
case2.mat=matrix(c(
  length(which(table$scores<0.17 & table$case2=="PATHWAY")),
  length(which(table$scores<0.17 & table$case2=="NO")),
  length(which(table$scores>=0.17 & table$case2=="PATHWAY")),
  length(which(table$scores>=0.17 & table$case2=="NO")))
  ,nrow=2)
fisher.test(case2.mat)
df_case2 <- data.frame(table$gene.ID, table$scores, table$case2)
colnames(df_case2) <- c("ID","score","S")
head(df_case2)

SCORE=df_case2$score
names(SCORE)=df_case2$ID
SCORE=sort(SCORE,decreasing=TRUE)
head(SCORE)

term2gene_case2=data.frame(term=df_case2$S,name=df_case2$ID)
head(term2gene_case2)

gsea.out_case2<-GSEA(SCORE,
		TERM2GENE=term2gene_case2,
		nPerm=10000,
		pvalueCutoff=1,
		pAdjustMethod = "BH")

gseaplot(gsea.out_case2,"PATHWAY")


