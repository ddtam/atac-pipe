#!/usr/bin/env Rscript
if (!require(BiocManager)) install.packages("BiocManager")
if (!require(tidyverse)) BiocManager::install("tidyverse")
if (!require(Seurat)) BiocManager::install("Seurat")
if (!require(Signac)) BiocManager::install("Signac")
if (!require(motifmatchr)) BiocManager::install("motifmatchr")
if (!require(chromVAR)) BiocManager::install("chromVAR")

if (!require(chromVARmotifs)) {
    if (!require(devtools)) install.packages("devtools")
    devtools::install_github("GreenleafLab/chromVARmotifs")
}

if (!require(BSgenome.Hsapiens.UCSC.hg19)) {
    BiocManager::install("BSgenome.Hsapiens.UCSC.hg19")
}

args <- commandArgs(trailingOnly = TRUE)

data("human_pwms_v2")

seurat_obj <- read_rds(args[1])

seurat_obj <-
    AddMotifs(
        seurat_obj,
        genome = BSgenome.Hsapiens.UCSC.hg19,
        pfm = human_pwms_v2
    ) %>%
    RunChromVAR(
        genome = BSgenome.Hsapiens.UCSC.hg19
    )

write_rds(
    seurat_obj,
    args[2]
)

export_tib <-
    FetchData(seurat_obj,
        vars = paste0(
            "chromvar_",
            str_replace_all(names(human_pwms_v2), "_", "-")
        )
    ) %>%
    set_names(names(.) %>%
        str_remove("chromvar_") %>%
        str_replace_all("-", "_"))

write_tsv(
    export_tib,
    args[3]
)