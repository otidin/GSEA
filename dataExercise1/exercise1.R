# library("clusterProfiler") 
# library("AnnotationDbi")
# library("rat2302.db")
# library(“DescTools”)
# library(org.Hs.eg.db)
# library(org.Rn.eg.db)

setwd('/Users/onurtidin/Desktop/desktopvirtual/SIBEnrichmentAnalysis/exercise/dataExercise1')
rat<-read.table("/Users/onurtidin/Desktop/desktopvirtual/SIBEnrichmentAnalysis/exercise/dataExercise1/rat_KD.txt")
dim(rat)
dimnames(rat)[[1]]
