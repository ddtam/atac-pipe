#!/usr/bin/env Rscript

if (!require(BiocManager)) install.packages("BiocManager")
if (!require(tidyverse)) BiocManager::install("tidyverse")
if (!require(Seurat)) BiocManager::install("Seurat")
if (!require(Signac)) BiocManager::install("Signac")

args <- commandArgs(trailingOnly = TRUE)

seurat_obj <-
    read_rds(args[1])

classification_tib <-
    read_tsv(args[2]) %>%
    dplyr::select(-"...1") %>%
    bind_cols("cell_id" = names(seurat_obj$cell_id), .) %>%
    column_to_rownames("cell_id") %>%
    mutate(
        "prediction_confidence" =
            pmap_dbl(subset(., select = -c(prediction)), ~ max(c(...)))
    )

print(colnames(seurat_obj) %>% length())
print(dim(classification_tib))

seurat_obj[["classification"]] <-
    CreateAssayObject(
        counts =
            classification_tib %>%
                dplyr::select(-prediction, -prediction_confidence) %>%
                t()
    )

seurat_obj$predicted_celltype <-
    classification_tib %>%
    dplyr::pull(prediction)

seurat_obj$prediction_confidence <-
    classification_tib %>%
    dplyr::pull(prediction_confidence)

write_rds(seurat_obj, args[3])
