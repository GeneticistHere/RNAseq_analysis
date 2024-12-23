library(DESeq2)
library(pheatmap)
library(dplyr)
library(RColorBrewer)
library(ggplot2)
library(ggrepel)
library(apeglm)

count_data<-read.csv('count_matrix.csv', header = TRUE, row.names = 1)
head(count_data)
sample_info=read.csv("design.csv", header = TRUE, row.names = 1)
head(sample_info)
#######
sample_info$Treatment<-factor(sample_info$Treatment)
sample_info$Sequencing<-factor(sample_info$Sequencing)

dds<- DESeqDataSetFromMatrix(countData = count_data, colData = sample_info, design = ~Sequencing + Treatment)
colnames(colData(dds))
dds$Treatment <- factor(dds$Treatment, levels = c("untreated", "treated"))
#Filter the genes
keep <- rowSums(counts(dds)) >= 5
dds <- dds[keep, ]
#Perforn the statistical test(s) to ldentify differentially expressed genes dds <- DESeq(dds)
dds <- DESeq(dds)
deseq_result <- results(dds)
deseq_result

#Change DESeg Object to R object(datafrane)
deseq_result <- as.data.frame (deseq_result)
class(deseq_result)
head (deseq_result)
#Order the result table by increasing p value
deseq_result_ordered<-deseq_result[order(deseq_result$pvalue), ]
head (deseq_result_ordered)

#Select genes with a significant change in gene expression (adjusted p-value below 0.05)
#And logfold change <1 and >1
#Step 1: filter based on p adjusted value
filtered <- deseq_result %>% filter(deseq_result$padj < 0.05)
#Step 2: filter based on fold changes. here we will use a threshold of 1
filtered <- filtered %>% filter (abs(filtered$log2FoldChange) > 1 )
dim(deseq_result)
dim(filtered)

#Save the deseg result. We will save the both the original data(res) and the filtered one(hits)
write.csv(deseq_result,'de_result.all.csv')
write.csv(filtered,'de_result.ftltered.csv')
#Save the normalized counts
normalized_counts <- counts(dds, normalized=TRUE)
head (normalized_counts)
write.csv(normalized_counts, 'normalized_counts.csv' )

#############visualization
plotDispEsts(dds)

############varience stabalizing transformation
vsd <- vst(dds,blind=FALSE)
#use transformed values to generate a pca plot
plotPCA(vsd, intgroup=c("Sequencing", "Treatment"))

#Heatmaps
#Heatmap of sample-to-sample distance matrix (with clustering) based on the normalized counts.
#generate the distance matrix
sampleDists <- dist(t(assay(vsd)))
sampleDistMatrix <- as.matrix(sampleDists)
colnames (sampleDistMatrix)
#set a color scheme 
colors <- colorRampPalette(rev(brewer.pal(9, "Greens")))(255)
#generate the heatmap 
pheatmap(sampleDistMatrix, clustering_distance_rows=sampleDists,
clustering_distance_cols=sampleDists, col=colors)

#Heatmap of log transformed normalized counts.
top_hits <- deseq_result[order(deseq_result$padj), ][1:50,]
top_hits <- row.names(top_hits)
top_hits
rld <- rlog(dds, blind=FALSE)
pheatmap(assay(rld) [top_hits,], cluster_rows=FALSE, show_rownames=TRUE,cluster_cols=FALSE)
pheatmap (assay(rld)[top_hits,])


annot_info <- as.data.frame(colData(dds) [,c('Sequencing', 'Treatment' )])
pheatmap(assay(rld)[top_hits,],cluster_rows=FALSE, show_rownames=TRUE,cluster_cols=FALSE, 
         annotation_col=annot_info)

#Heatmap of Z scores.
cal_z_score <- function(x) {
  (x - mean(x)) / sd(x)
}
zscore_all <- t(apply(normalized_counts, 1, cal_z_score))
zscore_subset <- zscore_all[top_hits,]
pheatmap (zscore_subset)


#####MA Plot
png("MA_plot.png", width=3000, height=2400, res=300)
plotMA(dds,ylim=c(-2,2))
dev.off()
#remove the noise
png("MA_plot_WN.png", width=3000, height=2400, res=300)
resLFC <- lfcShrink(dds, coef="Treatment_treated_vs_untreated", type="apeglm")
plotMA(resLFC,ylim=c(-2,2))
dev.off()

#####volcano plot
#change resLFC to a dataframe
resLFC <- as.data.frame(resLFC)
#Label the genes
resLFC$diffexpressed <- "NO"
resLFC$diffexpressed[resLFC$log2FoldChange>0.1 & resLFC$padj<0.05] <- "UP"
resLFC$diffexpressed[resLFC$log2FoldChange<0.1 & resLFC$padj<0.05] <- "DOWN"
resLFC$delabel<-NA
ggplot(data=resLFC,aes (x=log2FoldChange,y=-log10(pvalue),col=diffexpressed, label=delabel))+
  geom_point()+
  theme_minimal()+
  geom_text_repel()+
  scale_color_manual(values=c('skyblue', 'grey', 'darkred'))+
  theme(text=element_text(size=20))







############addiitonal commands
#####saving all the graphs
# 1. Save PCA plot
pdf("PCA_plot.pdf", width=10, height=8)
plotPCA(vsd, intgroup=c("Sequencing", "Treatment"))
dev.off()

# 2. Save sample distance heatmap
png("sample_distance_heatmap.png", width=3000, height=3000, res=300)
pheatmap(sampleDistMatrix, 
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists, 
         col=colors)
dev.off()

# 3. Save normalized counts heatmap
png("normalized_counts_heatmap.png", width=3000, height=3000, res=300)
pheatmap(assay(rld)[top_hits,], 
         cluster_rows=FALSE, 
         show_rownames=TRUE,
         cluster_cols=FALSE, 
         annotation_col=annot_info)
dev.off()

# 4. Save Z-score heatmap
png("zscore_heatmap.png", width=3000, height=3000, res=300)
pheatmap(zscore_subset)
dev.off()
# 5. Save Volcano plot
png("volcano_plot.png", width=3000, height=2400, res=300)
ggplot(data=resLFC, aes(x=log2FoldChange, y=-log10(pvalue), col=diffexpressed, label=delabel)) +
  geom_point() +
  theme_minimal() +
  geom_text_repel() +
  scale_color_manual(values=c('skyblue', 'grey', 'darkred')) +
  theme(text=element_text(size=20))
dev.off()

# 6. Save dispersion plot
png("dispersion_plot.png", width=3000, height=2400, res=300)
plotDispEsts(dds)
dev.off()
