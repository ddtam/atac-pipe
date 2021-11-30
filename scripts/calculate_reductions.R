#!/usr/bin/env Rscript

if (!require(BiocManager)) install.packages("BiocManager")
if (!require(tidyverse)) BiocManager::install("tidyverse")
if (!require(Seurat)) BiocManager::install("Seurat")
if (!require(Signac)) BiocManager::install("Signac")

args <- commandArgs(trailingOnly = TRUE)

seurat_obj <- read_rds(args[1])

seurat_obj <- seurat_obj %>%
    RunTFIDF() %>%
    FindTopFeatures(min.cutoff = 10) %>%
    RunSVD() %>%
    RunUMAP(reduction = "lsi", dims = 2:30) %>%
    FindNeighbors(reduction = "lsi", dims = 2:30) %>%
    FindClusters(resolution = 0.5, algorithm = 3)

seurat_obj <- AddMetaData(
    seurat_obj,
    FetchData(seurat_obj, "seurat_clusters"),
    "clusters_from_scatac"
)

write_rds(
    seurat_obj,
    args[2]
)