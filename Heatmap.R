
"Creating heatmaps of the human and guinea pig vaginal microbiota"
========================================================
  Embed Heatmaps, Infected samples:
  
  
#set working directory
setwd("../../data")
#install Packages needed
library("gplots", lib.loc="/Library/Frameworks/R.framework/Versions/3.1/Resources/library")
source("heatmapLib2[1].R")
table <- read.csv("ps_hu_gp.csv",header=T)
table2 <- read.csv("ps_gp_hm_meta.csv",header=T)
attach(table2)
nrow(table2)==length(species)
speciesf=factor(species)
levels(speciesf)
spColorTbl=c()
spColorTbl[levels(speciesf)[1]] ="red"
spColorTbl[levels(speciesf)[2]] = "blue"
row.names(table) <- table$SampleID
table <- table[,2:166]
table[1:10,1:10]
table <- as.matrix(table)
nr <- nrow(table)
hc <- hclust(dist(table),method="ward.D") # hierarchical clustering
nClrs <- 6
memb <- cutree(hc,k=nClrs)

totalRelAb <- apply(table, 2, sum)
totalRelAb[1:10]
o <- order(totalRelAb, decreasing=T)
table <- table[,o]
table <- as.matrix(table)
sideBars=cbind(spColorTbl[speciesf], colorTbl[memb])
colnames(sideBars)=c("Species", "Cluster")
o <- order(totalRelAb, decreasing=T)
table <- table[,o]
table <- as.matrix(table)

heatmap.2(table[,1:20],col=rainbow(50,start=1/6,end=0),cexCol=1,show.key=T)

heatmap2(table[,1:20],col=rainbow(50,start=1/6,end=0),Colv=NA,Rowv=as.dendrogram(hc),
         RowSideColors=sideBars,RowSideTitleCex=0.9,RowSideTitleLine=0.5,margins=c(13,15),
         cexCol=1.4,xlas=2, labRow=NA)

legend(0.8,0.2,legend=levels(factor(memb)),fill=colorTbl,title="Cluster",cex=0.6)

legend(0.8,0.3,legend=c("Human","Guinea pig"), pch=19, col=spColorTbl, title="Species",
       cex=0.4)

--------------------------------------------
  Embed Heatmaps, Non-infected samples:
  --------------------------------------------
#set working directory
setwd("../../data")
table <- read.csv("neg_hu_gp.csv",header=T)
table2 <- read.csv("neg_hu_gp_meta.csv",header=T)

attach(table2)
nrow(table2)==length(species)
speciesf=factor(species)
levels(speciesf)

spColorTbl=c()
spColorTbl[levels(speciesf)[1]] ="red"
spColorTbl[levels(speciesf)[2]] = "blue"

row.names(table) <- table$SampleID
table <- table[,2:166]
table[1:10,1:10]

table <- as.matrix(table)
nr <- nrow(table)
hc <- hclust(dist(table),method="ward.D") # hierarchical clustering
nClrs <- 7
memb <- cutree(hc,k=nClrs)

totalRelAb <- apply(table, 2, sum)
totalRelAb[1:10]

o <- order(totalRelAb, decreasing=T)
table <- table[,o]
table <- as.matrix(table)

sideBars=cbind(spColorTbl[speciesf], colorTbl[memb])
colnames(sideBars)=c("Species", "Cluster")

o <- order(totalRelAb, decreasing=T)
table <- table[,o]
table <- as.matrix(table)


heatmap2(table[,1:20],col=rainbow(50,start=1/6,end=0),Colv=NA,Rowv=as.dendrogram(hc),
         RowSideColors=sideBars,RowSideTitleCex=0.9,RowSideTitleLine=0.5,
         margins=c(13,15),cexCol=1.4,xlas=2, labRow=NA)
legend(0.8,0.2,legend=levels(factor(memb)),fill=colorTbl,title="Cluster",cex=0.6)
legend(0.8,0.3,legend=c("Human","Guinea pig"), pch=19, col=spColorTbl, title="Species", cex=0.4)