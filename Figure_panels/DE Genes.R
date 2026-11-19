##-----------------------------------
## DE gene list
## Author: Peng Yang
##-----------------------------------
library(ggplot2)
library(dplyr)
library(Seurat)

## Run from the repository root; processed objects are in Data/
load('Data/merged_data_Harmony_new.RData')
DefaultAssay(merged_combined) <- 'RNA'
Status <- rep('Baseline', length(merged_combined$orig.ident))
Status[which(merged_combined$orig.ident %in% c('PA1', 'PA3 #1', 'PA3 #2',
                                               'PA126', 'PA130', 'PA139', 'PA99', 'PA165')) ] = 'Late Progressor'
Status[which(merged_combined$orig.ident %in% c('PA11', 'PA46', 'PA144')) ] = 'Early Progressor'
merged_combined$Status_new <- as.factor(Status)

## Cutoff for positive versus negative (raw counts > 0)
cutoff_postive <- 0
merged_combined$ESR1_Status <- as.factor(ifelse(merged_combined@assays$RNA@counts['ESR1', ] > cutoff_postive,
                                                'ESR1+', 'ESR1-'))
merged_combined$FOSL1_Status <- as.factor(ifelse(merged_combined@assays$RNA@counts['FOSL1', ] > cutoff_postive,
                                                 'FOSL1+', 'FOSL1-'))

dir.create('DEG', showWarnings = FALSE)

##----------------------------------------------------
## Differentially expressed genes (Wilcoxon, FindMarkers)
## ESR1+ versus ESR1- and FOSL1- versus FOSL1+ cells,
## within each sample group. DEGs are defined downstream
## as |avg_log2FC| >= 1 and p_val_adj < 0.05.
## DEG/DEG_Tumor_cells_ESR1_Progressors.csv is the input
## of the Hallmark GSEA in Figure 1F.
##----------------------------------------------------

groups <- list()
groups[[1]] <- c('Baseline')
groups[[2]] <- c('Early Progressor', 'Late Progressor')
groups[[3]] <- c('Early Progressor')
groups[[4]] <- c('Late Progressor')

groups_names <- c('Basline', 'Progressors', 'Early Progressor', 'Late Progressor')

##----------------------------
## All cells
##----------------------------
for(g in 1:4){

  ## ESR1
  Idents(merged_combined) <- "ESR1_Status"
  de_genes <- FindMarkers(
    subset(merged_combined, subset = Status_new %in% groups[[g]]),
    ident.1 = "ESR1+",
    ident.2 = "ESR1-",
    assay = "RNA",
    slot = "data",
    test.use = "wilcox"
  )
  write.csv(de_genes, file = paste0('DEG/DEG_All_cells_ESR1_', groups_names[g],'.csv'))

  ## FOSL1
  Idents(merged_combined) <- "FOSL1_Status"
  de_genes <- FindMarkers(
    subset(merged_combined, subset = Status_new %in% groups[[g]]),
    ident.1 = "FOSL1-",
    ident.2 = "FOSL1+",
    assay = "RNA",
    slot = "data",
    test.use = "wilcox"
  )
  write.csv(de_genes, file = paste0('DEG/DEG_All_cells_FOSL1_', groups_names[g],'.csv'))
}

##----------------------------
## Tumor cells
##----------------------------
merged_combined_tumor <- subset(merged_combined, subset = Tumor_cells == "tumor cells")

for(g in 1:4){

  ## ESR1
  Idents(merged_combined_tumor) <- "ESR1_Status"
  de_genes <- FindMarkers(
    subset(merged_combined_tumor, subset = Status_new %in% groups[[g]]),
    ident.1 = "ESR1+",
    ident.2 = "ESR1-",
    assay = "RNA",
    slot = "data",
    test.use = "wilcox"
  )
  write.csv(de_genes, file = paste0('DEG/DEG_Tumor_cells_ESR1_', groups_names[g],'.csv'))

  ## FOSL1
  Idents(merged_combined_tumor) <- "FOSL1_Status"
  de_genes <- FindMarkers(
    subset(merged_combined_tumor, subset = Status_new %in% groups[[g]]),
    ident.1 = "FOSL1-",
    ident.2 = "FOSL1+",
    assay = "RNA",
    slot = "data",
    test.use = "wilcox"
  )
  write.csv(de_genes, file = paste0('DEG/DEG_Tumor_cells_FOSL1_', groups_names[g],'.csv'))
}

## All tumor cells, all samples
Idents(merged_combined_tumor) <- "ESR1_Status"
de_genes <- FindMarkers(
  merged_combined_tumor,
  ident.1 = "ESR1+",
  ident.2 = "ESR1-",
  assay = "RNA",
  slot = "data",
  test.use = "wilcox"
)
write.csv(de_genes, file = 'DEG/DEG_Tumor_cells_ESR1.csv')

Idents(merged_combined_tumor) <- "FOSL1_Status"
de_genes <- FindMarkers(
  merged_combined_tumor,
  ident.1 = "FOSL1-",
  ident.2 = "FOSL1+",
  assay = "RNA",
  slot = "data",
  test.use = "wilcox"
)
write.csv(de_genes, file = 'DEG/DEG_Tumor_cells_FOSL1.csv')

##----------------------------
## Non-tumor cells
##----------------------------
merged_combined_nontumor <- subset(merged_combined, subset = Tumor_cells == "non-tumor cells")

for(g in 1:4){

  ## ESR1
  Idents(merged_combined_nontumor) <- "ESR1_Status"
  de_genes <- FindMarkers(
    subset(merged_combined_nontumor, subset = Status_new %in% groups[[g]]),
    ident.1 = "ESR1+",
    ident.2 = "ESR1-",
    assay = "RNA",
    slot = "data",
    test.use = "wilcox"
  )
  write.csv(de_genes, file = paste0('DEG/DEG_nonTumor_cells_ESR1_', groups_names[g],'.csv'))

  ## FOSL1
  Idents(merged_combined_nontumor) <- "FOSL1_Status"
  de_genes <- FindMarkers(
    subset(merged_combined_nontumor, subset = Status_new %in% groups[[g]]),
    ident.1 = "FOSL1-",
    ident.2 = "FOSL1+",
    assay = "RNA",
    slot = "data",
    test.use = "wilcox"
  )
  write.csv(de_genes, file = paste0('DEG/DEG_nonTumor_cells_FOSL1_', groups_names[g],'.csv'))
}
