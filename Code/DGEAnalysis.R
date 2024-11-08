BiocManager::install("org.Mm.eg.db")
BiocManager::install("enrichplot")
BiocManager::install("clusterProfiler")

library(readr)
library(ggplot2)
library(ggrepel)
library(pheatmap)̥
library(dplyr)
library("org.Mm.eg.db")

###########################  Data wrangling ########################## 
#reading and Refining counts data table
cts <- read.table(file.choose(), header = TRUE, sep = "\t")
cts <- cts [,c(1,4,8,12,16,20,24,28,32)]
colnames(cts) <- c("Genes","Heart_0_1","Heart_0_2","Heart_12_1","Heart_12_2","Liver_0_1","Liver_0_2","Liver_12_1","Liver_12_2")
cts <- cts[-c(1:3),]
rownames(cts) <- cts[,1]
cts <- cts[,-c(1)]

#reading metadata 
coldata <- read.csv(file.choose(), header = TRUE, sep = "\t")
coldata <- coldata[,c("Samples","tissue","time")]
head(cts,5)
coldata
coldata$tissue <- factor(coldata$tissue)
coldata$time <- factor(coldata$time)

#check if row names of counts matrix is same as column names in colData
all(rownames(coldata) %in% colnames(cts))
row.names(coldata) <- coldata[,1]
coldata <- coldata[,-c(1)]
head(coldata)

all(rownames(coldata) %in% colnames(cts))
all(rownames(coldata) == colnames(cts))

write.csv(as.data.frame(cts), file="counts_cleaned.csv")

#DESeq data set made for DGE and visualisation
library(DESeq2)
dds <- DESeqDataSetFromMatrix(countData = cts, colData = coldata, design = ~tissue+time+tissue:time)
dds

#Prefiltering 
smallestGroupSize <- 2
keep <- rowSums(counts(dds) >= 10) >= smallestGroupSize
dds <- dds[keep,]
nrow(dds)
dds

########################### DATA REPRODUCIBILTY AND PATTERN OF VARIATION ########################### 

###########################  EDA and visualisation ###########################  

#Running log transformations on dds
vsd <- vst(dds, blind=FALSE)
head(assay(vsd), 3)
rld <- rlog(dds, blind=FALSE)

#Heatmap for counts matrix
top500variableGenes <- head(order(rowVars(assay(vsd)), decreasing = TRUE),500)
pheatmap(assay(vsd)[top500variableGenes,], cluster_rows=FALSE, show_rownames=FALSE, cluster_cols=TRUE, annotation_col=coldata)

#Heatmap for sample wise distances
sampleDists <- dist(t(assay(vsd)))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(vsd$tissue, colnames(vsd), sep=" - ")
colnames(sampleDistMatrix) <- paste(vsd$tissue, colnames(vsd), sep=" - ")
condition <- colData(dds[,c(1:8)])[,c("tissue")]
dfPCA <- as.data.frame(condition)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists, 
         annotation_col=dfPCA)

#PCA plot
plotPCA(vsd, intgroup=c("tissue")) + geom_label(aes(label = name))
plotPCA(vsd, intgroup=c("tissue")) + geom_label_repel(aes(label = name)) + theme(aspect.ratio = 1)

#setting Heart as the reference for statistical tests before running DESeq
dds$tissue <- relevel(dds$tissue, ref = "Heart")

#Differential gene expression assay with DESeqDataSet
dds <- DESeq(dds)
sizeFactors(dds)
colSums(counts(dds, normalized=T))
resultsNames(dds)

normalized_counts <- counts(dds, normalized=TRUE)
normalized_counts

###########################  DIFFERENTIAL EXPRESSION ANALYSIS ##################################

###########################  Differential Expression analysis  ##################################

# setting contrast for comparison of differential expression and setting FDR as 0.05
resultsNames(dds)

# Tissue-specific DEGs (e.g., Heart vs Liver)
res_tissue <- results(dds, contrast = c("tissue", "Liver", "Heart"), alpha = 0.05)
res_tissue
# Time-specific DEGs (e.g., ZT0 vs ZT12)
res_time <- results(dds, contrast = c("time", "0", "12"), alpha = 0.05)
res_time
# Interaction effects DEGs (tissue:time interaction)
res_interaction <- results(dds, name = "tissueLiver.time12", alpha = 0.05)
res_interaction

#summarising DGEs
sum(res_tissue$padj < 0.05, na.rm=TRUE)
summary(res_tissue)
sum(res_time$padj < 0.05, na.rm=TRUE)
summary(res_time)
sum(res_interaction$padj < 0.05, na.rm=TRUE)
summary(res_interaction)

#Filtering the significantly regulated genes
sigRes_tissue <- data.frame(res_tissue) %>% filter(padj < 0.05 & abs(log2FoldChange) >= 1)
sigResOrdered_tissue <- sigRes_tissue[order(sigRes_tissue$padj),]
top_sigResOrdered_tissue <- head(sigResOrdered_tissue,20)
write.csv(as.data.frame(sigResOrdered_tissue), file="tissue_LiverVSHeart.csv")
write.csv(as.data.frame(top_sigResOrdered_tissue), file="Top20_tissue_LiverVSHeart.csv")

sigRes_time <- data.frame(res_time) %>% filter(padj < 0.05 & abs(log2FoldChange) >= 1)
sigResOrdered_time <- sigRes_time[order(sigRes_time$padj),]
top_sigResOrdered_time <- head(sigResOrdered_time,20)
write.csv(as.data.frame(sigResOrdered_time), file="tissue_12VS0.csv")
write.csv(as.data.frame(top_sigResOrdered_time), file="Top20_time_12vs0.csv")

sigRes_interaction <- data.frame(res_interaction) %>% filter(padj < 0.05 & abs(log2FoldChange) >= 1)
sigResOrdered_interaction <- sigRes_interaction[order(sigRes_interaction$padj),]
top_sigResOrdered_interaction <- head(sigResOrdered_interaction,20)
write.csv(as.data.frame(sigResOrdered_interaction), file="tissue_time_interaction.csv")
write.csv(as.data.frame(top_sigResOrdered_interaction), file="Top20_interaction.csv")

# Subset normalized counts for these DEGs
all_sig_genes <- unique(c(rownames(sigRes_tissue),rownames(sigRes_time),rownames(sigRes_interaction)))
normalized_counts_all <- assay(rlog(dds))[all_sig_genes, ]

normalized_counts_all_df <- as.data.frame(normalized_counts_all)
rownames(normalized_counts_all_df) <- sub('\\.[0-9]*$', '', rownames(normalized_counts_all_df))
normalized_counts_sig_tissue_df <- as.data.frame(sigRes_tissue)
rownames(normalized_counts_sig_tissue_df) <- sub('\\.[0-9]*$', '', rownames(normalized_counts_sig_tissue_df))
normalized_counts_sig_time_df <- as.data.frame(sigRes_time)
rownames(normalized_counts_sig_time_df) <- sub('\\.[0-9]*$', '', rownames(normalized_counts_sig_time_df))
normalized_counts_sig_interaction_df <- as.data.frame(sigRes_interaction)
rownames(normalized_counts_sig_interaction_df) <- sub('\\.[0-9]*$', '', rownames(normalized_counts_sig_interaction_df))

########################### Analysing and Interpreting results with visualisations ##################################

# Volcano plots
library(EnhancedVolcano)

EnhancedVolcano(res_tissue,
                lab = rownames(res_tissue),
                x = 'log2FoldChange',
                y = 'pvalue', labSize = 3, axisLabSize = 10, pCutoff = 0.05, ylim = c(0, -log10(10e-15)), titleLabSize = 10, subtitleLabSize = 8,
                captionLabSize = 8, legendLabSize = 8)
EnhancedVolcano(res_time,
                lab = rownames(res_time),
                x = 'log2FoldChange',
                y = 'pvalue', labSize = 3, axisLabSize = 10, pCutoff = 0.05, ylim = c(0, -log10(10e-15)), titleLabSize = 10, subtitleLabSize = 8,
                captionLabSize = 8, legendLabSize = 8)
EnhancedVolcano(res_interaction,
                lab = rownames(res_interaction),
                x = 'log2FoldChange',
                y = 'pvalue', labSize = 3, axisLabSize = 10, pCutoff = 0.05, ylim = c(0, -log10(10e-15)), titleLabSize = 10, subtitleLabSize = 8,
                captionLabSize = 8, legendLabSize = 8)

# Heatmap of all significant DEGs from the three groups of comparison
pheatmap(normalized_counts_all,
         scale = "row",
         clustering_distance_rows = "correlation",
         clustering_distance_cols = "euclidean",
         main = "Clustering of all significant DEGs",
         show_rownames = FALSE)

#Subset normalized counts for each group of top DEGs separately
normalized_counts_top_tissue <- assay(rlog(dds))[rownames(top_sigResOrdered_tissue), ]
normalized_counts_top_time <- assay(rlog(dds))[rownames(top_sigResOrdered_time), ]
normalized_counts_top_interaction <- assay(rlog(dds))[rownames(top_sigResOrdered_interaction), ]

#annotating and analysing the clustering of Top DEGs with heatmap
normalized_counts_top_tissue_df <- as.data.frame(normalized_counts_top_tissue)
rownames(normalized_counts_top_tissue_df) <- sub('\\.[0-9]*$', '', rownames(normalized_counts_top_tissue_df))
normalized_counts_top_tissue_df$symbols <- mapIds(org.Mm.eg.db, keys = rownames(normalized_counts_top_tissue_df), column = c("SYMBOL"), keytype = 'ENSEMBL')
rownames(normalized_counts_top_tissue) <- normalized_counts_top_tissue_df$symbols
pheatmap(normalized_counts_top_tissue,
         scale = "row",
         clustering_distance_rows = "correlation",
         clustering_distance_cols = "euclidean",
         main = "Expression Pattern of Top 20 DEGs between Liver and Heart")

normalized_counts_top_time_df <- as.data.frame(normalized_counts_top_time)
rownames(normalized_counts_top_time_df) <- sub('\\.[0-9]*$', '', rownames(normalized_counts_top_time_df))
normalized_counts_top_time_df$symbols <- mapIds(org.Mm.eg.db, keys = rownames(normalized_counts_top_time_df), column = c("SYMBOL"), keytype = 'ENSEMBL')
rownames(normalized_counts_top_time) <- normalized_counts_top_time_df$symbols
pheatmap(normalized_counts_top_time,
         scale = "row",
         clustering_distance_rows = "correlation",
         clustering_distance_cols = "euclidean",
         main = "Expression Pattern of Top 20 DEGs between time T0 and T12")

normalized_counts_top_interaction_df <- as.data.frame(normalized_counts_top_interaction)
rownames(normalized_counts_top_interaction_df) <- sub('\\.[0-9]*$', '', rownames(normalized_counts_top_interaction_df))
normalized_counts_top_interaction_df$symbols <- mapIds(org.Mm.eg.db, keys = rownames(normalized_counts_top_interaction_df), column = c("SYMBOL"), keytype = 'ENSEMBL')
rownames(normalized_counts_top_interaction) <- normalized_counts_top_interaction_df$symbols
pheatmap(normalized_counts_top_interaction,
         scale = "row",
         clustering_distance_rows = "correlation",
         clustering_distance_cols = "euclidean",
         main = "Expression Pattern of Top 20 DEGs between Liver and Heart across the sampling time")

#listing all top 20 genes from all three comparisons
DEGs <- data.frame(c("tissue",normalized_counts_top_tissue_df$symbols),c("time",normalized_counts_top_time_df$symbols),c("interaction",normalized_counts_top_interaction_df$symbols))
write.csv(x = DEGs, file ="D:/harsha_hp/career mtech/Biostate AI/Data/output/DEGs.csv")

########################### FUNCTIONAL ENRICHMENT ANALYSIS ########################### 

### Pathway Enrichment Analysis using clusterprofiler ###

library(enrichplot)
library(clusterProfiler)
enrichment <- function(x,y){
  plot=enrichGO(
    x,
    org.Mm.eg.db,
    keyType = "ENSEMBL",
    ont = y,
    pvalueCutoff = 0.05,
    pAdjustMethod = "BH",
    qvalueCutoff = 0.2,
    minGSSize = 10
  )
  dotplot(plot)
}
data(rownames(normalized_counts_top_interaction_df))

gene_entrez_sig_tissue <- AnnotationDbi::select(org.Mm.eg.db, keys = rownames(normalized_counts_sig_tissue_df), columns = c("ENTREZID"), keytype = "ENSEMBL")
gene_entrez_sig_time <- AnnotationDbi::select(org.Mm.eg.db, keys = rownames(normalized_counts_sig_time_df), columns = c("ENTREZID"), keytype = "ENSEMBL")
gene_entrez_sig_interaction <- AnnotationDbi::select(org.Mm.eg.db, keys = rownames(normalized_counts_sig_interaction_df), columns = c("ENTREZID"), keytype = "ENSEMBL")

enrichment(rownames(normalized_counts_sig_tissue_df), "MF")
enrichment(rownames(normalized_counts_sig_time_df), "MF")
enrichment(rownames(normalized_counts_sig_interaction_df), "MF")
kegg_tissue <- enrichKEGG(gene = gene_entrez_sig_tissue$ENTREZID, organism = 'mmu', pvalueCutoff = 0.05, keyType = "kegg")
kegg_time <- enrichKEGG(gene = gene_entrez_sig_time$ENTREZID, organism = 'mmu', pvalueCutoff = 0.05, keyType = "kegg")
kegg_interaction <- enrichKEGG(gene = gene_entrez_sig_interaction$ENTREZID, organism = 'mmu', pvalueCutoff = 0.05, keyType = "kegg")
dotplot(kegg_tissue)
dotplot(kegg_time)
dotplot(kegg_interaction)
