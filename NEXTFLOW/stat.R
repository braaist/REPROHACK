library(DESeq2)

##nettoyage de l'environnement
rm(list = objects())
graphics.off()
args=commandArgs(trailingOnly = TRUE)
count_data<- read.csv(args, sep = "\t", header=T, comment.char = '#')
##met data
#count_data<- read.csv("out_count_tab.txt", sep = "\t", header=T, comment.char = '#')
count_data<-count_data[-c(2:6)]
rownames(count_data) <- count_data[,1]
count_data <- count_data[,-1]
col2<-c("MUT","MUT","MUT","WT","WT","WT","WT","MUT")
metadata<-data.frame(row.names = colnames(count_data) ,col2)

##Analyse exploratoires
##Suprression des gènes non exprimés 30452 génes non exprimées, il reste 30159 gènes
##0.05 pour les suites
count_data <- count_data[-which(rowMeans(count_data) < 5),]
##running analysis
condition<-factor((metadata$col2))
data_seq<-DESeqDataSetFromMatrix(count_data, data.frame(condition) , ~condition)
dds<-DESeq(data_seq)
res<-results(dds)
summary(res)

##correction de benjamini alpha=0.1, on a 325 genes differentiellement exprimés
sign_gene<-na.omit(res)
sign_gene<-sign_gene[sign_gene$padj<0.1,]
sign_gene_df<-as.data.frame(sign_gene)
summary(sign_gene)

##les best hits
bestreg<-sign_gene_df[order(sign_gene_df$padj),]
bestreg<-tail(bestreg, n=10)
write.table(bestreg,"bestreg.txt" )


#Controle qualité du modèle

##On regarde l'histogramme des pvalues pour voir leur distribution et si ils ont une distribution uniforme avec un pic à 0
pval<-res$pvalue
png("hist-pval.png")
hist(res$pvalue)
dev.off()

##visualisation of the dispersion

png("Maplot.png")
plotMA(dds)
dev.off()


png("dispersionnplot.png")
plotDispEsts(dds)
dev.off()

##le gène avec la plus petite pvalue et celle qui a la plus grande pour voir les differences d'expression dans les 2 conditionss
plotCounts(dds,gene=which.min(res$padj), intgroup = "condition")
plotCounts(dds,gene=which.max(res$padj), intgroup = "condition")

##Analyse en composante principal
##transormation des données de comptages pour la visualisation
png("pca_cond.png")
vsdta<-vst(dds, blind = F)
plotPCA(vsdta, intgroup="condition")
dev.off()


