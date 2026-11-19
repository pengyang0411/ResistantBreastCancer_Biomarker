##------------------------------------
## pySCENIC step 1/3: export tumor-cell
## counts, TF list and cell information
## Author: Peng Yang
##------------------------------------
library(Seurat)
library(Matrix)
library(dplyr)
library(readxl)

##------------------------------------------
## Load the data
##------------------------------------------

## Run from the repository root; processed objects are in Data/
load('Data/merged_data_Harmony_new.RData')
DefaultAssay(merged_combined) <- 'RNA'
Status <- rep('Baseline', length(merged_combined$orig.ident))
Status[which(merged_combined$orig.ident %in% c('PA1', 'PA3 #1', 'PA3 #2',
                                               'PA126', 'PA130', 'PA139', 'PA99', 'PA165')) ] = 'Late Progressor'
Status[which(merged_combined$orig.ident %in% c('PA11', 'PA46', 'PA144')) ] = 'Early Progressor'
merged_combined$Status_new <- factor(Status, levels = c('Baseline', 'Early Progressor', 'Late Progressor'))

## Tumor cells only
tum <- subset(merged_combined, subset = Tumor_cells == 'tumor cells')
rm(merged_combined); gc()

##------------------------------------------
## Downsample to at most 600 tumor cells per
## sample (16 samples -> 8,633 cells); GRNBoost2
## run time and memory scale with cell number
##------------------------------------------
CELLS_PER_SAMPLE <- 600
set.seed(777)
Idents(tum) <- tum$orig.ident
tum <- subset(tum, downsample = CELLS_PER_SAMPLE)
ncol(tum)
table(tum$orig.ident)

output <- 'TF_activity/pyscenic'
dir.create(output, showWarnings = FALSE, recursive = TRUE)

## Raw counts, genes x cells
cnt <- tryCatch(GetAssayData(tum, layer = 'counts'),
                error = function(e) GetAssayData(tum, slot = 'counts'))
Matrix::writeMM(cnt, file.path(output, 'counts.mtx'))
writeLines(rownames(cnt), file.path(output, 'genes.txt'))
writeLines(colnames(cnt), file.path(output, 'barcodes.txt'))

## Candidate TFs: human TFs (Lambert et al., 2018) detected in the data
hmtf <- read_xlsx('Data/human TFs v1.01.xlsx', col_names = 'gene')
writeLines(intersect(hmtf$gene, rownames(cnt)), file.path(output, 'tfs.txt'))

## Cell information for 'pySCENIC_readback.R' and 'Figure S1.R'
grp <- ifelse(tum$Status_new == 'Baseline', 'Baseline', 'Progressor')
write.csv(data.frame(CellID = colnames(tum), Group = grp,
                     Sample = as.character(tum$orig.ident)),
          file.path(output, 'cellinfo.csv'), row.names = FALSE)
