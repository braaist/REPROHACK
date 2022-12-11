#groupe 1 projet hackaton

##chargement de la librairie DESeq2
library(DESeq2)

##nettoyage de l'environnement
rm(list = objects())
graphics.off()

##chargement et restructuration des données 
args=commandArgs(trailingOnly = TRUE)
count_data<- read.csv(args, sep = "\t", header=T, comment.char = '#')
##met data
#count_data<- read.csv("out_count_tab.txt", sep = "\t", header=T, comment.char = '#')
count_data<-count_data[-c(2:6)]
rownames(count_data) <- count_data[,1]
count_data <- count_data[,-1]

##Construction des metadata
col2<-c("MUT","MUT","MUT","WT","WT","WT","WT","wt")
metadata<-data.frame(row.names = colnames(count_data) ,col2)
condition<-factor((metadata$col2))

## analyse sur les données filtrées 
data_seq<-DESeqDataSetFromMatrix(count_data, data.frame(condition) , ~condition)
dds<-DESeq(data_seq)
res<-results(dds)
summary(res)

##correction de la pvalue, on a 219 genes differentiellement exprimés.  83 surexprimé et 136 réprimés
sign_gene<-na.omit(res)
sign_gene<-sign_gene[sign_gene$padj<0.05,]
sign_gene_df<-as.data.frame(sign_gene)
summary(sign_gene)

##les best hits
bestreg<-sign_gene_df[order(sign_gene_df$padj),]
bestreg<-tail(bestreg, n=10)
write.table(bestreg,"bestreg_unfiltered.txt" )


##Analyse sur données filtrés
##Suprression des gènes ou il y a moins d'un read par echantillon, il reste 15353 gènes
kept_gene<-rowSums( counts(dds) >= 1 ) == 8
dds_filtered <- dds[kept_gene,]
res_filtered<-results(dds_filtered)
summary(res_filtered)

##correction pvalue alpha=0.05, on a 120 genes differentiellement exprimés, 66 surexprimé et 54 réprimés.
sign_gene_filtered<-na.omit(res_filtered)
sign_gene_filtered<-sign_gene_filtered[sign_gene_filtered$padj < 0.05,]
sign_gene_filtered_df<-as.data.frame(sign_gene_filtered)
summary(sign_gene_filtered)

##les best hits
bestreg_filtered<-sign_gene_filtered_df[order(sign_gene_filtered_df$padj),]
write.table(bestreg_filtered,"bestreg_filtered.txt" )

#Controle qualité du modèle

##On regarde l'histogramme des pvalues pour voir leur distribution et si ils ont une distribution uniforme avec un pic à 0
pval<-res$pvalue
png("hist_pval_unfiltered.png")
hist(res$pvalue)
dev.off()

pval<-res$pvalue
png("hist_pval_filtered.png")
hist(res_filtered$pvalue)
dev.off()

##visualisation of the dispersion

png("Maplot.png")
plotMA(dds)
dev.off()

png("Maplot_filtered.png")
plotMA(dds_filtered)
dev.off()


png("dispersionnplot.png")
plotDispEsts(dds)
dev.off()

png("dispersionnplot_filtered.png")
plotDispEsts(dds_filtered)
dev.off()

##le gène avec la plus petite pvalue et celle qui a la plus grande pour voir les differences d'expression dans les 2 conditionss
png("best_filtered.png")
plotCounts(dds,gene=which.min(res$padj), intgroup = "condition")
dev.off()

png("best_unfiltered.png")
plotCounts(dds_filtered,gene=which.min(res_filtered$padj), intgroup = "condition")
dev.off()

##Analyse en composante principal
##transormation des données de comptages pour la visualisation
png("pca_cond.png")
vsdta<-vst(dds, blind = F)
plotPCA(vsdta, intgroup="condition")
dev.off()
