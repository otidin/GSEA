setwd("/Users/lindadib/Work/SIB-BCF/work/SIB_tutoring/Dib_EA/2019/dataExercice4")



require(clusterProfiler)
library(org.Hs.eg.db)
keytypes(org.Hs.eg.db)

table       <-read.csv("HS_pvalues.csv")
SCORE       <-table$score
names(SCORE)<-table$gene.ENTREZ.ID
SCORE       <-sort(SCORE ,decreasing=T)

ego <- gseGO(geneList     = SCORE,
              OrgDb        = org.Hs.eg.db,
              ont          = "ALL",
              nPerm        = 1000,
	    pvalueCutoff = 1,
              verbose      = FALSE)

head(ego)             
library("DOSE")         
gseaplot(ego, geneSetID="GO:0048518")
dotplot(ego, showCategory=30)

#########
egoBP <- gseGO(geneList     = SCORE,
               OrgDb        = org.Hs.eg.db,
               ont          = "BP",
               nPerm        = 1000,
	       pvalueCutoff = 1,
               verbose      = FALSE)

head(egoBP)                   
gseaplot(egoBP, geneSetID="GO:0048518")
dotplot(egoBP, showCategory=30)
plotGOgraph(egoBP)
emapplot(egoBP)
