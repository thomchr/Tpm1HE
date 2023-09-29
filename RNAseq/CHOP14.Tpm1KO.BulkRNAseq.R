#Install packages/ load packages

library(data.table)
#library(dplyr)
library(tidyverse)
library('biomaRt')
library(debrowser)

setwd("/Users/thomc/Dropbox/BulkRNASeq/Tpm1KO RNAseq STAR/chop14/")

###make a gene list file for analysis, takes STAR-aligned input, finds ENSG gene names, outputs tables
df <- readRDS("fc_chop14_TakarabulkRNA_CT.rds")
df <- as.data.frame(df$counts)
setDT(df, keep.rownames = "ENSG")

#convert ENSG to hgnc_symbol
mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))

genes <- df$ENSG
G_list <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id","hgnc_symbol"),values=genes,mart=mart)
df <- merge(df,G_list,by.x="ENSG",by.y="ensembl_gene_id")

#move hgnc_symbol to first column
out <- df %>% relocate(hgnc_symbol) 
out <- out[!duplicated(out$hgnc_symbol),] #remove duplicated gene names  #61K to 40K

write.table(out, "CHOP14.Tpm1KO.Bulk.txt", row.names=F, col.names=T, quote=F, sep="\t")


####analyze using DEBrowser
startDEBrowser() #brings up a separate window interface
#did all the DE analyses using this platform - can get significant DE genes, heatmap for individual comparisons, exported afterbatch.tsv which is normalized data for all genes, also have up/down DE genes for each two-sample comparison

####want to get a heatmap of all DE genes and cluster by sample to get an overall impression of how conserved DE genes are among each condition
library(ComplexHeatmap)
library("circlize")

####compile a list of all DE genes (p<0.05, logFC>1.5) from subcomparisons
#Eu_wt vs Eu_s
a <- fread("Euploid_GATA1wt vs Euploid_GATA1s/221116.Euploidwt_Euploids.GeneList.Up836.csv", header=T)
b <- fread("Euploid_GATA1wt vs Euploid_GATA1s/221116.Euploidwt_Euploids.GeneList.Down791.csv", header=T)
#test <- a[a$padj <= 0.05] #confirmed all rows/genes are significant

#Eu_wt vs T21_wt
c <- fread("Euploid_GATA1wt vs T21_GATA1wt/221116.Euploidwt_T21wt.GeneList.Up335.csv", header=T)
d <- fread("Euploid_GATA1wt vs T21_GATA1wt/221116.Euploidwt_T21wt.GeneList.Down418.csv", header=T)

#Eu_s vs T21_s
e <- fread("EuploidGATA1s_T21GATA1s/221116.Euploids_T21s.GeneList.Up1715.csv", header=T)
f <- fread("EuploidGATA1s_T21GATA1s/221116.Euploids_T21s.GeneList.Down2083.csv", header=T)

#T21_wt vs T21_s
g <- fread("T21_GATA1wt vs T21_GATA1s/221116.T21wt_T21s.GeneList.Up266.csv", header=T)
h <- fread("T21_GATA1wt vs T21_GATA1s/221116.T21wt_T21s.GeneList.Down770.csv", header=T)

#Eu_wt vs T21_s
i <- fread("Euploid_GATA1wt vs T21_GATA1s/Euwt_T21s.GeneList.Up1113.csv")
j <- fread("Euploid_GATA1wt vs T21_GATA1s/Euwt_T21s.GeneList.Down1471.csv")

delist <- bind_rows(a, b, c, d, e, f, g, h, i, j) #9798 DE genes
delist <- as.data.frame(delist %>% distinct(V1)) #4854 unique DE genes

#identify those genes in the afterbatch.tsv normalized data
norm <- fread("afterbatchtable.tsv", header=T)
DE.list <- inner_join(delist, norm, by=c("V1" = "ID"))

#make a clustered heatmap of those genes


colnames(DE.list) <- c("V1", "T21_GATA1wt",   #this is based on Kaoru.Metadata.txt file
                    "T21_GATA1wt",
                    "T21_GATA1wt",
                    "T21_GATA1wt",
                    "T21_GATA1wt",
                    "T21_GATA1wt",
                    "Euploid_GATA1s",
                    "Euploid_GATA1wt",
                    "Euploid_GATA1wt",
                    "Euploid_GATA1wt",
                    "Euploid_GATA1s",
                    "Euploid_GATA1wt",
                    "Euploid_GATA1wt",
                    "Euploid_GATA1wt",
                    "Euploid_GATA1s",
                    "Euploid_GATA1s",
                    "Euploid_GATA1s",
                    "Euploid_GATA1s",
                    "T21_GATA1s",
                    "T21_GATA1s",
                    "T21_GATA1s",
                    "T21_GATA1s",
                    "T21_GATA1s",
                    "T21_GATA1s")
row.names(DE.list) <- DE.list$V1 #set gene as rowname
DE.list <- DE.list[,-c(1)] #remove gene column
DE.list <- as.matrix(DE.list) # Final matrix should how rownames gene, and colnames samples

#Scale matrix by gene
scaled.DE <- t(apply(DE.list,1,scale))
colnames(scaled.DE) <- colnames(DE.list)

#Clustering
#Hierarchal clustering, two closest correlations grouped, iteratively using ward.d2 method with 1 - pearson correlation as the distance metric
clust<-hclust(as.dist(1-cor(t(scaled.DE), method="pearson")), method="ward.D2")

#Six clusters expected based on 
grp <- cutree(clust,k=12) #This object contains the genes in the group

col_fun1 = colorRamp2(c(-2,-1,0,1,2),c("blue","skyblue","white","lightcoral","red"),space="RGB")

#Prepare the first heatmap, scores normalized by maximum value
h_list = Heatmap(scaled.DE,
                 show_row_names=FALSE,
                 col=col_fun1,
                 cluster_rows=clust,
                 cluster_columns = TRUE,
                 #cluster_columns = FALSE,
                 #column_order=c(8:10,12:14,1:6,7,11,15:18,19:24),
                 heatmap_legend_param = list(title = "Scaled Expression", col_fun=col_fun1,color_bar="continuous"),
                 split = 12)+
  Heatmap(grp, name = "clusters", show_row_names = FALSE, width = unit (5, "mm"),
          col = rainbow(12)) #structure(names = c("1", "2", "3", "4","5","6","7","8","9","10","11","12"), 
                         # c(#"#999999","#56B4E9","#009E73" ,"#CC79A7","#F0E442" ,"#E69F00" , "#999999","#56B4E9","#009E73" ,"#CC79A7","#F0E442" ,"#E69F00")))+
  
  pdf("4condition.DEgenes.expression_heatmap.pdf", useDingbats=FALSE)
draw(h_list)
dev.off()

#pull identities for groups
write.table(grp, 'grp.txt', quote=F) #now can grep a given letter to search Ontology within specific clusters with something like awk -F " " '$2 ~ /5|6/' grp.txt






#Make a heat map for selected genes
library(ComplexHeatmap)
library("circlize")

df <- fread("CHOP14.Tpm1KO.Bulk.txt", header=T)
names(df) = gsub(pattern = "Aligned.sortedByCoord.out.bam", replacement = "", x = names(df))

genelist <- c("KIT", "KITLG", "CD34", "MECOM", "CDH5", "HOXA9", "PECAM", "RUNX1", "MLLT3", "EPOR", "GATA1", "GATA2", "MPL")
dfsub <- df %>% filter(hgnc_symbol %in% genelist)
dfsub <- dfsub %>% arrange(match(hgnc_symbol, genelist))
d4 <- dfsub[,c(1,19:22,7:10, 27:30, 15:18)]
row.names(d4) <- d4$hgnc_symbol
d4 <- d4[,-1]
mat <- as.matrix(d4)
row.names(mat) <- dfsub$hgnc_symbol

#Scale matrix by gene
scaled.mat <- t(apply(mat,1,scale))
colnames(scaled.mat) <- colnames(mat)
rownames(scaled.mat) <- rownames(mat)

col_fun = colorRamp2(c(-2,0,2),c("white","lightcoral","red3"),space="RGB")
Heatmap(scaled.mat, name = "mat", 
        col = col_fun,
        cluster_rows = FALSE,
        cluster_columns = FALSE, 
        row_names_side = "left")


png("HE.heatmap.png", width=1200, height=1000, res=300)
Heatmap(scaled.mat, name = "mat", 
        col = col_fun,
        cluster_rows = FALSE,
        cluster_columns = FALSE, 
        row_names_side = "left")
dev.off()
#################


#Scale matrix by gene
scaled.mat <- t(apply(mat,1,scale))
colnames(scaled.mat) <- colnames(scaled)

#Clustering
#Hierarchal clustering, two closest correlations grouped, iteratively using ward.d2 method with 1 - pearson correlation as the distance metric
clust<-hclust(as.dist(1-cor(t(scaled.mat), method="pearson")), method="ward.D2")

#Six clusters expected based on 
grp <- cutree(clust,k=2) #This object contains the genes in the group

col_fun1 = colorRamp2(c(-2,0,2),c("white","lightcoral","red"),space="RGB")

#Prepare the first heatmap, scores normalized by maximum value
h_list = Heatmap(scaled.mat,
                 show_row_names=FALSE,
                 col=col_fun1,
                 cluster_rows=clust,
                 #cluster_columns = TRUE,
                 cluster_columns = FALSE,
                 column_order=c(1:8),
                 heatmap_legend_param = list(title = "Scaled Expression", col_fun=col_fun1,color_bar="continuous"),
                 split = 12)+
  Heatmap(grp, name = "clusters", show_row_names = FALSE, width = unit (5, "mm"),
          col = rainbow(12)) #structure(names = c("1", "2", "3", "4","5","6","7","8","9","10","11","12"), 
# c(#"#999999","#56B4E9","#009E73" ,"#CC79A7","#F0E442" ,"#E69F00" , "#999999","#56B4E9","#009E73" ,"#CC79A7","#F0E442" ,"#E69F00")))+

pdf("4condition.DEgenes.expression_heatmap.pdf", useDingbats=FALSE)
draw(h_list)
dev.off()













# Attach the `pheatmap` library
library(pheatmap)

# install DESeq if necessary
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

# install and load package
BiocManager::install("DESeq")
library("DESeq")

# if you can't install DESeq, I have hosted the file at https://davetang.org/file/TagSeqExample.tab
# example_file <- "https://davetang.org/file/TagSeqExample.tab"

# load data and subset
example_file <- system.file ("extra/TagSeqExample.tab", package="DESeq")
data <- read.delim(example_file, header=T, row.names="gene")
data_subset <- as.matrix(data[rowSums(data)>50000,])

# create heatmap using pheatmap
pheatmap(data_subset)

# Create and store the annotated heatmap object
annot_colors=list(condition=c(Ad="#F0978D",P7="#63D7DE"))
heatmap_annotated <-
  pheatmap(
    selected_df,
    cluster_rows = TRUE,
    cluster_cols = TRUE,
    show_rownames = FALSE,
    annotation_names_col= FALSE,
    annotation_col = annotation_df,
    annotation_colors = annot_colors[1],
    main = "",
    colorRampPalette(c(
      "deepskyblue",
      "black",
      "yellow"
    ))(25),
    scale = "row" # Scale values with respect to genes (rows)
  )
