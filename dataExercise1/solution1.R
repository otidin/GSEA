#Get data
setwd("/Users/onurtidin/Desktop/desktopvirtual/SIBEnrichmentAnalysis/exercise/dataExercise2")

#Question1
rat<-read.table("rat_KD.txt", sep = "\t", header = TRUE,stringsAsFactors=FALSE)
dimnames(rat)[[1]]<-rat[,1]
rowNb  <-which(rat[,1] == "1398751_at")
v1<- t.test(rat[rowNb,2:7], rat[rowNb,8:12], alternative = "less")
v1
v1<- t.test(rat[rowNb,2:7], rat[rowNb,8:12], alternative = "greater")
v1
v1<- t.test(rat[rowNb,2:7], rat[rowNb,8:12])
v1
#Question2
ttestRat <- function(df, grp1, grp2){
        x = df[grp1]
        y = df[grp2]
        x = as.numeric(x)
        y = as.numeric(y)
        results = t.test(x, y) 
        results$p.value}

rawp <- apply(rat, 1, ttestRat, grp1 = c(2:7), grp2 = c(8:12))
p_holm <- p.adjust(sort(rawp),method="BH") 
p_holm <- p.adjust(rawp,method="BH") 
hist(p_holm) 
