library(Seurat)
library(ggplot2)
library(dplyr)
library(tidyr)
library(viridis)
library(Signac)

load_gene_list <- function(file_path) {
  if (!file.exists(file_path)) {
    stop("File not found: ", file_path)
  }
  gene_list <- readLines(file_path)
  if (length(gene_list) == 0) {
    stop("Gene list is empty in file: ", file_path)
  }
  gene_list
}

perform_qc <- function(seurat_obj, min_genes = 200, max_genes = 6000, max_mito_percent = 5) {
  seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^MT-")
  seurat_obj <- subset(
    seurat_obj, 
    subset = nFeature_RNA > min_genes & nFeature_RNA < max_genes & percent.mt < max_mito_percent
  )
  seurat_obj
}

perform_dimensional_reduction <- function(seurat_obj, n_pcs = 30) {
  seurat_obj <- FindVariableFeatures(seurat_obj, selection.method = "vst")
  seurat_obj <- ScaleData(seurat_obj)
  seurat_obj <- RunPCA(seurat_obj, npcs = n_pcs)
  seurat_obj <- RunUMAP(seurat_obj, dims = 1:n_pcs)
  seurat_obj
}

perform_clustering <- function(seurat_obj, resolution = 0.5) {
  seurat_obj <- FindNeighbors(seurat_obj, dims = 1:30)
  seurat_obj <- FindClusters(seurat_obj, resolution = resolution)
  seurat_obj
}

identify_markers <- function(seurat_obj, ident_1, ident_2 = NULL, assay_use = "RNA") {
  markers <- FindMarkers(
    seurat_obj, 
    ident.1 = ident_1, 
    ident.2 = ident_2, 
    assay = assay_use,
    min.pct = 0.25, 
    logfc.threshold = 0.25
  )
  markers
}

plot_violin <- function(seurat_obj, gene_list_path) {
  gene_list <- load_gene_list(gene_list_path)
  valid_genes <- intersect(gene_list, rownames(seurat_obj))
  if (length(valid_genes) == 0) {
    stop("None of the genes in the provided list exist in the Seurat object.")
  }
  rna_data <- GetAssayData(seurat_obj, assay = "integrated", slot = "data")[valid_genes, ]
  rna_data <- log2(rna_data + 1)
  rna_data[is.na(rna_data)] <- 0
  plot_data <- as.data.frame(t(rna_data))
  plot_data$Cluster <- seurat_obj$seurat_clusters
  plot_data_long <- gather(plot_data, key = "Gene", value = "Expression", -Cluster)
  plot_data_long <- filter(plot_data_long, Expression > 0)
  plot_data_long$Cluster <- factor(plot_data_long$Cluster, levels = sort(unique(plot_data_long$Cluster)))
  print(
    ggplot(plot_data_long, aes(x = Cluster, y = Expression, fill = Gene)) +
      geom_violin(trim = TRUE) +
      ggtitle("Violin Plot of Log2(Normalized RNA Expression + 1)") +
      xlab("Cluster") +
      ylab("Log2(Normalized RNA Expression + 1)") +
      theme_minimal() +
      theme(legend.title = element_blank())
  )
}

plot_umap_average <- function(seurat_obj, gene_list_path) {
  gene_list <- load_gene_list(gene_list_path)
  valid_genes <- intersect(gene_list, rownames(seurat_obj))
  if (length(valid_genes) == 0) {
    stop("None of the genes in the provided list exist in the Seurat object.")
  }
  rna_data <- GetAssayData(seurat_obj, assay = "integrated", slot = "data")[valid_genes, ]
  rna_data <- log2(rna_data + 1)
  rna_data[is.na(rna_data)] <- 0
  average_expression <- colMeans(rna_data)
  seurat_obj[["Average_Expression"]] <- average_expression
  print(
    FeaturePlot(seurat_obj, features = "Average_Expression", reduction = "umap") +
      ggtitle("UMAP Colored by Average Expression of Gene List") +
      scale_color_viridis(option = "D")
  )
}

plot_signac_gene <- function(seurat_obj, gene_name) {
  if (!gene_name %in% rownames(seurat_obj)) {
    stop("Gene not found in the Seurat object: ", gene_name)
  }
  print(
    FeaturePlot(seurat_obj, features = gene_name, reduction = "umap") +
      ggtitle(paste("Signac Feature Plot for", gene_name)) +
      scale_color_viridis(option = "A")
  )
}

seurat_obj <- readRDS("all.rds")
DefaultAssay(seurat_obj) <- "integrated"
seurat_obj <- perform_qc(seurat_obj)
seurat_obj <- perform_dimensional_reduction(seurat_obj)
seurat_obj <- perform_clustering(seurat_obj)
