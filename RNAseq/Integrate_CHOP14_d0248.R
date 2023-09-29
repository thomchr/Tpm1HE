#get CHOP14 D0, 4, 8 with CHOP14 D2 in same file for PCA/heat map analysis

#as of 230118 cannot get colors to label anymore (new R version) but this might be the way to specify now: https://stackoverflow.com/questions/51161172/how-to-color-labels-of-dendogram-with-dendextend-and-heatmap-2-using-pre-defined

library(data.table)
library(dplyr)
library(gplots)
library("reshape2")
library(RColorBrewer)
library("ape")
library("pheatmap")

setwd("/Users/thomc/Dropbox/BulkRNASeq/Tpm1KO RNAseq STAR/chop14/")

df <- fread("CHOP14.Tpm1KO.Bulk.txt", header=T)
df2 <- fread("/Users/thomc/Dropbox/BulkRNASeq/Tpm1KO RNAseq STAR/chop10/CHOP14vs14KO.D2.Illumina.Tpm1KO.Bulk.txt", header=T)

x <- merge(df, df2, by.x=c("hgnc_symbol", "ENSG"), by.y=c("hgnc_symbol", "ENSG"))
#x[x$hgnc_symbol == "TPM1"]

write.table(x, "CHOP14.Tpm1KO.Bulk.D0248.txt", sep="\t", row.names=F, col.names=T, quote=F)
#put this into DEBrowser to filter and normalize -> saved output as CHOP14.Tpm1KO.D0248.afterbatchtable.tsv
#used debrowser to make a PCA separately


#Make a heat map from this



#or find overall similarity btw samples
#https://master.bioconductor.org/packages/release/workflows/vignettes/rnaseqGene/inst/doc/rnaseqGene.html
# that has good tutorial for PCA and heat map as well

library("vsn") #installed, had a lot of dependencies that failed but still loads ok
library("DESeq2") #vst is a DESeq2 function, force installed =  BiocManager::install("DESeq2", force = TRUE)

dds <- fread("CHOP14.Tpm1KO.Bulk.D0248.txt", header=T) #we recommend raw counts
dds <- dds[(dds$hgnc_symbol != ""),]
rownames(dds) <- dds$hgnc_symbol
dds <- as.matrix(dds[,-c(1,2)]) #drop hgnc_symbol and ENSG
vsd <- vst(dds, blind = FALSE) #transforms to a value

sampleDists <- dist(t(vsd)) #default calc for R dist. needs samples in rows and genes in columns. Default is method="euclidean" ??
sampleDists

#Heatmap of sample-to-sample distances using the variance stabilizing transformed values.
sampleDistMatrix <- data.matrix( sampleDists )
dim(sampleDistMatrix)
m <- sampleDistMatrix
#write.table(m, "CHOP14.Tpm1KO.Bulk.D0248.SampleDistMatrix.txt", sep="\t", row.names=T, col.names=T, quote=F)
#m <- as.matrix(fread("CHOP14.Tpm1KO.Bulk.D0248.SampleDistMatrix.txt", header=T))

#Labels - same ones used for all plots/data sets
WTd0 <- c('WT14-d0-2Aligned.sortedByCoord.out.bam',
          'WT14-d0-3Aligned.sortedByCoord.out.bam',
          'WT14-d0-4Aligned.sortedByCoord.out.bam',
          'WT14-d0-5Aligned.sortedByCoord.out.bam')
KOd0 <- c('KO9-d0-2Aligned.sortedByCoord.out.bam',
          'KO9-d0-3Aligned.sortedByCoord.out.bam',
          'KO9-d0-4Aligned.sortedByCoord.out.bam',
          'KO9-d0-5Aligned.sortedByCoord.out.bam')
WTd2 <- c('WT14-D2-1Aligned.sortedByCoord.out.bam',
          'WT14-D2-2Aligned.sortedByCoord.out.bam',
          'WT14-D2-3Aligned.sortedByCoord.out.bam',
          'WT14-D2-4Aligned.sortedByCoord.out.bam')
KOd2 <- c('T14KO9-D2-1Aligned.sortedByCoord.out.bam',
          'T14KO9-D2-2Aligned.sortedByCoord.out.bam',
          'T14KO9-D2-3Aligned.sortedByCoord.out.bam',
          'T14KO9-D2-4Aligned.sortedByCoord.out.bam')
WTd4 <- c('WT-d4-1Aligned.sortedByCoord.out.bam',
          'WT-d4-2Aligned.sortedByCoord.out.bam',
          'WT-d4-3Aligned.sortedByCoord.out.bam',
          'WT-d4-4Aligned.sortedByCoord.out.bam')
KOd4 <- c('KO9-d4-1Aligned.sortedByCoord.out.bam',
          'KO9-d4-2Aligned.sortedByCoord.out.bam',
          'KO9-d4-3Aligned.sortedByCoord.out.bam',
          'KO9-d4-5Aligned.sortedByCoord.out.bam')
WTd8 <- c('WT4-HPC-1Aligned.sortedByCoord.out.bam',
          'WT4-HPC-2Aligned.sortedByCoord.out.bam',
          'WT4-HPC-3Aligned.sortedByCoord.out.bam',
          'WT4-HPC-4Aligned.sortedByCoord.out.bam')
KO1d8 <- c('TPM1-KO1-1Aligned.sortedByCoord.out.bam',
           'TPM1-KO1-2Aligned.sortedByCoord.out.bam',
           'TPM1-KO1-3Aligned.sortedByCoord.out.bam',
           'TPM1-KO1-4Aligned.sortedByCoord.out.bam')
KO9d8 <- c('TPM1-KO9-1Aligned.sortedByCoord.out.bam',
           'TPM1-KO9-2Aligned.sortedByCoord.out.bam',
           'TPM1-KO9-3Aligned.sortedByCoord.out.bam',
           'TPM1-KO9-4Aligned.sortedByCoord.out.bam')

WTd0_col=sample(c("snow2"), length(WTd0), replace = TRUE, prob = NULL)
KOd0_col=sample(c("snow3"), length(KOd0), replace = TRUE, prob = NULL)
WTd2_col=sample(c("beige"), length(WTd2), replace = TRUE, prob = NULL)
KOd2_col=sample(c("bisque3"), length(KOd2), replace = TRUE, prob = NULL)
WTd4_col=sample(c("tomato"), length(WTd4), replace = TRUE, prob = NULL)
KOd4_col=sample(c("tomato4"), length(KOd4), replace = TRUE, prob = NULL)
WTd8_col=sample(c("lightblue1"), length(WTd8), replace = TRUE, prob = NULL)
KO1d8_col=sample(c("blue"), length(KO1d8), replace = TRUE, prob = NULL)
KO9d8_col=sample(c("mediumblue"), length(KO9d8), replace = TRUE, prob = NULL)

# separate cluster and create dendrogram
clust <- dist(as.matrix(m)) 
hc <- hclust(clust, method = "complete" )
dd <- as.dendrogram(hc)
pdf("CHOP14.clusteringtrial.pdf", height=10, width=5)
par(cex=0.7, mar=c(10,10,10,10))
wts = matrix(nrow = 1, ncol = nrow(m))
wts[1,] = 1 #set all default as 1
wts[1,21] = 2000 #WTd0
wts[1,36] = 1000 #WTd2
wts[1,17] = 500 # WTd4
wts[1,25] = 100 # WTd8
wts[1,13] = 2 #KO9d8
plot(reorder(dd, wts, agglo.FUN = mean), horiz=T, lwd = 2, main = "reordered")
dev.off()

####################################################################
####################################################################

#unrooted
pdf("CHOP14.unrooted.pdf", width=10, height=10)
phylo <- plot(as.phylo(hc), type = "fan", cex =0.5, tip.color = c(WTd0_col,KOd0_col,WTd2_col,KOd2_col,WTd4_col,KOd4_col,WTd8_col,KO1d8_col,KO9d8_col)) #no.margin=TRUE
dev.off()

#plots
pdf("CHOP14.wKO1HPCs.RNAseqSamples.Clustered.pdf")
par(oma=c(5,1,1,5))
dd <- as.dendrogram(reorder(dd, wts, agglo.FUN = mean), horiz=T)
heatmap <- heatmap.2(m, dendrogram="row", Rowv=dd, Colv=rev(dd), col=rev(brewer.pal(9, "Greys")), scale="none", margins=c(5,5), ColSideColors=matrix(c(WTd0_col,KOd0_col,WTd2_col,KOd2_col,WTd4_col,KOd4_col,WTd8_col,KO9d8_col,KO1d8_col), nrow=1, ncol=ncol(m)), cexRow=1, offsetRow=0, offsetCol=0, trace='none', density.info=c("none"), key=TRUE)
dev.off()



pdf("CHOP14.wKO1HPCs.legend.pdf")
heatmap <- heatmap.2(m, dendrogram="both", Rowv=dd, Colv=rev(dd), col=brewer.pal(9, "Greys"), scale="none", margins=c(1,1), ColSideColors=matrix(c(WTd0_col,KOd0_col,WTd2_col,KOd2_col,WTd4_col,KOd4_col,WTd8_col,KO9d8_col,KO1d8_col), nrow=1, ncol=ncol(m)), RowSideColors=matrix(c(WTd0_col,KOd0_col,WTd2_col,KOd2_col,WTd4_col,KOd4_col,WTd8_col,KO9d8_col,KO1d8_col), nrow=nrow(m), ncol=1), cexRow=0.5, cexCol=1, srtCol=90, offsetRow=-45, offsetCol=0, trace='none', density.info=c("none"), key=TRUE, symkey=FALSE, keysize=2, key.title="Dissimilarity", key.xlab="Dissimilarity") #, RowSideColorsSize=2, ColSideColorsSize=2) #, lmat = lmat, lwid = lwid, lhei = lhei) 

legend('topright', legend=c("WT iPSC", "KO iPSC", "WT d2", "KO d2", "WT d4", "KO d4", "WT HPC", "KO1 HPC", "KO9 HPC"), fill=c("snow2", "snow3", "beige", "bisque3", "tomato", "tomato4", "lightblue1", "blue", "mediumblue"), border=FALSE, bty="n", y.intersp = 0.8, cex=1)
dev.off()

####################################################################
####################################################################


#test 

sampleDistMatrix <- as.matrix( sampleDists )
dim(sampleDistMatrix)

rownames(sampleDistMatrix) # <- paste(vsd$condition, vsd$type, sep="-")
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)



####################################################################
####################################################################
#leave out KO1

dds <- fread("CHOP14.Tpm1KO.Bulk.D0248.txt", header=T) #we recommend raw counts
dds <- dds[(dds$hgnc_symbol != ""),]
rownames(dds) <- dds$hgnc_symbol
dds <- as.matrix(dds[,-c(1,2,11,12,13,14)]) #drop hgnc_symbol and ENSG, as well as KO1 HPCs
vsd <- vst(dds, blind = FALSE) #transforms to a value

sampleDists <- dist(t(vsd)) #default calc for R dist. needs samples in rows and genes in columns. Default is method="euclidean" ??
sampleDists

#Heatmap of sample-to-sample distances using the variance stabilizing transformed values.
sampleDistMatrix <- as.matrix( sampleDists )
dim(sampleDistMatrix)
m <- sampleDistMatrix

#Labels - same ones used for all plots/data sets
WTd0 <- c('WT14-d0-2Aligned.sortedByCoord.out.bam',
          'WT14-d0-3Aligned.sortedByCoord.out.bam',
          'WT14-d0-4Aligned.sortedByCoord.out.bam',
          'WT14-d0-5Aligned.sortedByCoord.out.bam')
KOd0 <- c('KO9-d0-2Aligned.sortedByCoord.out.bam',
          'KO9-d0-3Aligned.sortedByCoord.out.bam',
          'KO9-d0-4Aligned.sortedByCoord.out.bam',
          'KO9-d0-5Aligned.sortedByCoord.out.bam')
WTd2 <- c('WT14-D2-1Aligned.sortedByCoord.out.bam',
          'WT14-D2-2Aligned.sortedByCoord.out.bam',
          'WT14-D2-3Aligned.sortedByCoord.out.bam',
          'WT14-D2-4Aligned.sortedByCoord.out.bam')
KOd2 <- c('T14KO9-D2-1Aligned.sortedByCoord.out.bam',
          'T14KO9-D2-2Aligned.sortedByCoord.out.bam',
          'T14KO9-D2-3Aligned.sortedByCoord.out.bam',
          'T14KO9-D2-4Aligned.sortedByCoord.out.bam')
WTd4 <- c('WT4-HPC-1Aligned.sortedByCoord.out.bam',
          'WT4-HPC-2Aligned.sortedByCoord.out.bam',
          'WT4-HPC-3Aligned.sortedByCoord.out.bam',
          'WT4-HPC-4Aligned.sortedByCoord.out.bam')
KOd4 <- c('KO9-d4-1Aligned.sortedByCoord.out.bam',
          'KO9-d4-2Aligned.sortedByCoord.out.bam',
          'KO9-d4-3Aligned.sortedByCoord.out.bam',
          'KO9-d4-5Aligned.sortedByCoord.out.bam')
WTd8 <- c('WT4-HPC-1Aligned.sortedByCoord.out.bam',
          'WT4-HPC-2Aligned.sortedByCoord.out.bam',
          'WT4-HPC-3Aligned.sortedByCoord.out.bam',
          'WT4-HPC-4Aligned.sortedByCoord.out.bam')
#KO1d8 <- c('TPM1-KO1-1Aligned.sortedByCoord.out.bam',
#           'TPM1-KO1-2Aligned.sortedByCoord.out.bam',
#           'TPM1-KO1-3Aligned.sortedByCoord.out.bam',
#           'TPM1-KO1-4Aligned.sortedByCoord.out.bam')
KO9d8 <- c('TPM1-KO9-1Aligned.sortedByCoord.out.bam',
           'TPM1-KO9-2Aligned.sortedByCoord.out.bam',
           'TPM1-KO9-3Aligned.sortedByCoord.out.bam',
           'TPM1-KO9-4Aligned.sortedByCoord.out.bam')

WTd0_col=sample(c("snow2"), length(WTd0), replace = TRUE, prob = NULL)
KOd0_col=sample(c("snow3"), length(KOd0), replace = TRUE, prob = NULL)
WTd2_col=sample(c("beige"), length(WTd2), replace = TRUE, prob = NULL)
KOd2_col=sample(c("bisque3"), length(KOd2), replace = TRUE, prob = NULL)
WTd4_col=sample(c("tomato"), length(WTd4), replace = TRUE, prob = NULL)
KOd4_col=sample(c("tomato4"), length(KOd4), replace = TRUE, prob = NULL)
WTd8_col=sample(c("lightblue1"), length(WTd8), replace = TRUE, prob = NULL)
#KO1d8_col=sample(c("blue"), length(KO1d8), replace = TRUE, prob = NULL)
KO9d8_col=sample(c("mediumblue"), length(KO9d8), replace = TRUE, prob = NULL)



# separate cluster and create dendrogram
clust <- dist(as.matrix(m)) 
hc <- hclust(clust, method = "complete" )
dd <- as.dendrogram(hc)
pdf("CHOP14.clusteringtrial.pdf", width=5, height=10)
par(cex=0.7, mar=c(10,10,10,10))
wts = matrix(nrow = 1, ncol = nrow(m))
wts[1,] = 1 #set all default as 1
wts[1,17] = 2000 #WTd0
wts[1,32] = 1000 #WTd2
wts[1,13] = 500 # WTd4
wts[1,21] = 100 # WTd8
wts[1,19] = 2 #KO9d8
plot(reorder(dd, wts, agglo.FUN = mean), horiz=T, lwd = 2, main = "reordered")
dev.off()

####################################################################
####################################################################

#unrooted
jpeg("CHOP14.unrooted.jpg", width=400, height=2000, res=300, quality=100)
phylo <- plot(as.phylo(hc), type = "fan", cex =0.5, tip.color = c(WTd0_col,KOd0_col,WTd2_col,KOd2_col,WTd4_col,KOd4_col,WTd8_col,KO1d8_col,KO9d8_col)) #no.margin=TRUE
dev.off()

#plots
#pdf("CHOP14.RNAseqSamples.Clustered.pdf")
jpeg("CHOP14.RNAseqSamples.Clustered.jpg", width=2000, height=2000, res=300, quality=100)
par(oma=c(2,2,2,2))
dd <- as.dendrogram(reorder(dd, wts, agglo.FUN = mean), horiz=T)
heatmap <- heatmap.2(m, dendrogram="row", Rowv=dd, Colv=rev(dd), col=rev(brewer.pal(9, "Greys")), scale="none", margins=c(5,5), cexCol=0.25, cexRow=0.5, offsetRow=0, offsetCol=0, trace='none', density.info=c("none"), key=FALSE) #ColSideColors=matrix(c(WTd0_col,KOd0_col,WTd2_col,KOd2_col,WTd4_col,KOd4_col,WTd8_col,KO9d8_col,KO1d8_col), nrow=1, ncol=ncol(m)), 
dev.off()



jpeg("CHOP14.wKO1HPCs.legend.jpg", width=2000, height=2000, res=300, quality=100)
heatmap <- heatmap.2(m, dendrogram="row", Rowv=dd, Colv=rev(dd), col=rev(brewer.pal(9, "Greys")), scale="none", margins=c(9,9), ColSideColors=matrix(c(WTd0_col,KOd0_col,WTd2_col,KOd2_col,WTd4_col,KOd4_col,WTd8_col,KO9d8_col), nrow=1, ncol=ncol(m)), RowSideColors=matrix(c(WTd0_col,KOd0_col,WTd2_col,KOd2_col,WTd4_col,KOd4_col,WTd8_col,KO9d8_col), nrow=nrow(m), ncol=1), cexRow=0.5, cexCol=0.25, srtCol=90, offsetRow=-45, offsetCol=0, trace='none', density.info=c("none"), key=TRUE, symkey=FALSE, keysize=2, key.title="Dissimilarity", key.xlab="Dissimilarity") #, RowSideColorsSize=2, ColSideColorsSize=2) #, lmat = lmat, lwid = lwid, lhei = lhei) 

legend('topright', legend=c("WT iPSC", "KO iPSC", "WT d2", "KO d2", "WT d4", "KO d4", "WT HPC", "KO9 HPC"), fill=c("snow2", "snow3", "beige", "bisque3", "tomato", "tomato4", "lightblue1", "mediumblue"), border=FALSE, bty="n", y.intersp = 0.8, cex=0.5)
dev.off()
