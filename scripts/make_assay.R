#!/usr/bin/env Rscript

if (!require(Seurat)) BiocManager::install("Seurat")
if (!require(Signac)) BiocManager::install("Signac")
if (!require(EnsDb.Hsapiens.v75)) BiocManager::install("EnsDb.Hsapiens.v75")
if (!require(tidyverse)) BiocManager::install("tidyverse")

args <- commandArgs(trailingOnly = TRUE)

annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v75)

seqlevelsStyle(annotations) <- "UCSC"
genome(annotations) <- "hg19"

counts <-
    Read10X_h5(arguments[1])

metadata <- read.csv(arguments[2],
    header = T, row.names = 1
)

atac_assay <-
    CreateChromatinAssay(
        counts = counts,
        sep = c(":", "-"),
        annotation = annotations,
        min.cells = 10,
        genome = "hg19",
        fragments = arguments[3]
    )

seurat_obj <-
    CreateSeuratObject(
        counts = atac_assay,
        assay = "peaks",
        meta.data = metadata
    )

write_rds(
    seurat_obj,
    arguments[4]
)