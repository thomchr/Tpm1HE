#subset RNA seq sample columns for GSEA

library(data.table)

setwd("/Users/thomc/Dropbox/BulkRNASeq/Tpm1KO RNAseq STAR/chop14/GSEA/")

df<-fread("CHOP14.Tpm1KO.Bulk.txt", header=T)
colnames(df)

#D0
x<-df[,c(1,2,23,24,25,26,3,4,5,6)]
names(x)[names(x) == 'ENSG'] <- 'Description'
x$Description<-"NA"
x<-x[(x$hgnc_symbol != ""),]
write.table(x, "CHOP14.Tpm1KO.D0.txt", quote=F, sep="\t", col.names=T, row.names=F)

#D2
df<-fread("/Users/thomc/Dropbox/BulkRNASeq/Tpm1KO RNAseq STAR/chop10/CHOP14vs14KO.D2.Illumina.Tpm1KO.Bulk.txt", header=T)
colnames(df)
x<-df[,c(1,2,7,8,9,10,3,4,5,6)]
names(x)[names(x) == 'ENSG'] <- 'Description'
x$Description<-"NA"
x<-x[(x$hgnc_symbol != ""),]
write.table(x, "CHOP14.Tpm1KO.D2.txt", quote=F, sep="\t", col.names=T, row.names=F)


#D4
x<-df[,c(1,2,19,20,21,22,7,8,9)]
names(x)[names(x) == 'ENSG'] <- 'Description'
x$Description<-"NA"
x<-x[(x$hgnc_symbol != ""),]
write.table(x, "CHOP14.Tpm1KO.D4.txt", quote=F, sep="\t", col.names=T, row.names=F)

#D8
x<-df[,c(1,2,27,28,29,30,15,16,17,18)]
names(x)[names(x) == 'ENSG'] <- 'Description'
x$Description<-"NA"
x<-x[(x$hgnc_symbol != ""),]
write.table(x, "CHOP14.Tpm1KO.D8.txt", quote=F, sep="\t", col.names=T, row.names=F)

