#install required packages and loading them

install.packages("pacman")
library("biomaRt")
require(pacman)  
library(reshape2)
library(grid)
library(gridExtra)
library(Seurat)
library(tximport)
library(SeqGSEA)
#loads the package needed to import alevin data
library("tximport", lib.loc="~/R/x86_64-pc-linux-gnu-library/3.6")

#run pacman for downloading the necessary libraries
pacman::p_load(pacman, dplyr, GGally, ggplot2, ggthemes, 
               ggvis, httr, lubridate, plotly, rio, rmarkdown, shiny, 
               stringr, tidyr, ggpubr, grid, gridExtra, Seurat, tximport, biomaRt, SeqGSEA) 

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("fishpond")
BiocManager::install("SeqGSEA")
BiocManager::install('grimbough/biomaRt')
BiocManager::install('org.Hs.eg.db')
BiocManager::install("EnsDb.Hsapiens.v79")

library(AnnotationDbi)
library(EnsDb.Hsapiens.v79)
library(org.Hs.eg.db)

###########################
#                         #  
#Project-4 - Programmer   #
#                         #
###########################

#===========================================================================================================================================================================

#importing the UMI matrix

file_path = file.path("/projectnb/bf528/users/hedgehog_2022/project_4/Data_curator/scripts/salmon_newindex/alevin/quants_mat.gz")
file_avelian <- tximport(file_path, type="alevin")
#===========================================================================================================================================================================

#Creating the SUERAT object

pbmc <- CreateSeuratObject(counts = file_avelian$counts, project = "panc", min.cells = 3, min.features = 200)

#===========================================================================================================================================================================

#Assigning gene names to detect mitochondrial genes

gene_names <- sub("[.][0-9]*$", "", pbmc@assays$RNA@counts@Dimnames[[1]])
geneIDs <- ensembldb::select(EnsDb.Hsapiens.v79, keys= gene_names, keytype = "GENEID", columns = c("SYMBOL","GENEID"))
gene_original_names <- geneIDs$SYMBOL
pbmc@assays$RNA@counts@Dimnames[[1]]->names
pbmc@assays$RNA@counts@Dimnames[[1]] <- gene_original_names
pbmc@assays$RNA@data@Dimnames[[1]] <- gene_original_names

#Adding the mito percentage to the Seural object of RNA data
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")

#===========================================================================================================================================================================

#Visualizing plots to view how the data is

#head(pbmc@meta.data, 5)
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
plot1 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1+plot2

#===========================================================================================================================================================================

#Filtering the data

pbmc<- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

#===========================================================================================================================================================================

#Normalization
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)

#===========================================================================================================================================================================

#Variance filtering
pbmc[["RNA"]]@meta.features <- data.frame(row.names = rownames(pbmc[["RNA"]]))
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
top10 <- head(VariableFeatures(pbmc), 10)
plot1 <- VariableFeaturePlot(pbmc)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2

#===========================================================================================================================================================================

#Dimensionality detection using PCA, elbowplot, Jackstraw

all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
print(pbmc[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(pbmc, dims = 1:2, reduction = "pca")
DimPlot(pbmc, reduction = "pca")
ElbowPlot(pbmc)
pbmc <- JackStraw(pbmc, num.replicate = 100)
pbmc <- ScoreJackStraw(pbmc, dims = 1:20)
JackStrawPlot(pbmc, dims = 1:15)

#===========================================================================================================================================================================

#Clustrering, UMAP
pbmc <- FindNeighbors(pbmc, dims = 1:10)
pbmc <- FindClusters(pbmc, resolution = 0.5)
head(Idents(pbmc), 5)
pbmc <- RunUMAP(pbmc, dims = 1:10)
DimPlot(pbmc, reduction = "umap")
cluster_centers <- Idents(pbmc)

#===========================================================================================================================================================================

#Frequency of cells in each cluster
df <- as.data.frame(table(cluster_centers))
#Clusters histogram
p<-ggplot(data=df, aes(x=cluster_centers, y=Freq)) + ggtitle("Frequency of cells in each cluster") + theme(plot.title = element_text(hjust = 0.5)) + 
  geom_bar(stat="identity")
p + labs(x = "Clusters", y = "Frequency", size = 16)

pie(df$Freq, labels=paste0(round((df$Freq/sum(df$Freq)),2)*100, "", "%"),  main = "Frequency of cells in each cluster")

#saveRDS(pbmc, file = "PBMS.rds")

#===========================================================================================================================================================================

###########################
#                         #  
#Project-4 - Analyst      #
#                         #
###########################

#===========================================================================================================================================================================

pbmc <- readRDS("/projectnb/bf528/users/hedgehog_2022/project_4/Programmer/panc.rds")
pmbc<-NormalizeData(pbmc)

#===========================================================================================================================================================================

#mapping ensembl gene ids
ensembl <- useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl")
x<-as_tibble(listAttributes(ensembl))

gene_map <- as_tibble(
  getBM(
    attributes=c("ensembl_gene_id", "hgnc_symbol"),
    mart=ensembl
  )
)

#write.csv(gene_map, "map.csv")

#===========================================================================================================================================================================

#marker genes

markers <- FindAllMarkers(pbmc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

#===========================================================================================================================================================================

#renaming the ids with names
markers<- markers %>% mutate (genes = sub("[.][0-9]*$","",markers$gene)) 

MAPPER<- markers %>% dplyr::inner_join(
  gene_map,
  by=c("genes" = "ensembl_gene_id")
)

#===========================================================================================================================================================================

#Biologist role writing to directory - not required
#write.csv(MAPPER,"markers.csv")
#===========================================================================================================================================================================

#top markers

top2<- markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC)

top10 <- markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)

#===========================================================================================================================================================================

#heatmap
#Testing if heatmap is okay
#DoHeatmap(pbmc, features = top10$gene) + NoLegend()

#===========================================================================================================================================================================

#annotating ids

new.cluster.ids <- c("Beta", "Epsilon", "Alpha", "Transcription factor", "Ductal", "Delta", "Cytotoxic T", "Gamma", 
                     "Mast", "Stellate", "Acinar", "Vascular")



names(new.cluster.ids) <- levels(pbmc)
pbmc <- RenameIdents(pbmc, new.cluster.ids)

#===========================================================================================================================================================================
#UMAP
DimPlot(pbmc, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
#===========================================================================================================================================================================

#heatmap of differential expression for top 10 genes in all custers
top10 <- markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
top10<- top10 %>% mutate (genes = sub("[.][0-9]*$","",top10$gene)) 
top10<- top10 %>% dplyr::inner_join(
  gene_map,
  by=c("genes" = "ensembl_gene_id")
)

top10 <- markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
x<-DoHeatmap(pbmc, features = top10$gene) + NoLegend()
