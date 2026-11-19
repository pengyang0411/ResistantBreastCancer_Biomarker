##------------------------------------
## TF activity of tumor cells with
## decoupleR (ULM) and DoRothEA regulons
## Author: Peng Yang
##------------------------------------
library(Seurat)
library(dplyr)
library(tidyr)
library(decoupleR)
library(dorothea)
library(fgsea)

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
tum$ESR1_status <- ifelse(tum@assays$RNA@counts['ESR1', ] > 0, 'ESR1+', 'ESR1-')
tum$Group2 <- ifelse(tum$Status_new == 'Baseline', 'Baseline', 'Progressor')
tum <- NormalizeData(tum, verbose = FALSE)
rm(merged_combined); gc()

output <- 'TF_activity/decoupleR_output'
dir.create(output, showWarnings = FALSE, recursive = TRUE)

##------------------------------------------
## DoRothEA regulons (confidence A, B and C)
##------------------------------------------
data('dorothea_hs', package = 'dorothea')
net <- dorothea_hs %>%
  filter(confidence %in% c('A', 'B', 'C')) %>%
  transmute(source = tf, target = target, mor = mor)
nrow(net)

##------------------------------------------
## Univariate linear model (ULM) activity
## for each TF in each tumor cell
##------------------------------------------
## Seurat v5 uses layer =; older versions use slot =
mat <- as.matrix(tryCatch(GetAssayData(tum, layer = 'data'),
                          error = function(e) GetAssayData(tum, slot = 'data')))

acts <- run_ulm(mat = mat, net = net,
                .source = 'source', .target = 'target',
                .mor = 'mor', minsize = 5)
rm(mat); gc()

## TF x cell activity matrix
act_wide <- acts %>%
  filter(statistic == 'ulm') %>%
  pivot_wider(id_cols = 'source', names_from = 'condition',
              values_from = 'score') %>%
  as.data.frame()
rownames(act_wide) <- act_wide$source
act_wide$source <- NULL
act_mat <- as.matrix(act_wide)[, colnames(tum), drop = FALSE]

saveRDS(act_mat, file = file.path(output, 'decoupleR_ULM_activity_tumor_cells.rds'))

##------------------------------------------
## ER-program activity per cell
## (HALLMARK estrogen response module scores)
##------------------------------------------

## ER program used for the per-cell comparison with FOSL1 activity
## (Supplementary Figure 1O): early + late estrogen response genes
er_gmts <- c('Data/HALLMARK_ESTROGEN_RESPONSE_EARLY.v2025.1.Hs.gmt',
             'Data/HALLMARK_ESTROGEN_RESPONSE_LATE.v2025.1.Hs.gmt')
er_genes <- unique(unlist(lapply(er_gmts, function(f) unname(unlist(gmtPathways(f))))))
er_genes <- intersect(er_genes, rownames(tum))
tum <- AddModuleScore(tum, features = list(er_genes), name = 'ER_program')

## ER programs used for the baseline versus progressor comparison
## (Supplementary Figure 1R): early, late and combined
hall <- gmtPathways('Data/msigdb_v2022.1.Hs_GMTs/h.all.v2022.1.Hs.symbols.gmt')
er_sets <- list(
  early = intersect(hall[['HALLMARK_ESTROGEN_RESPONSE_EARLY']], rownames(tum)),
  late  = intersect(hall[['HALLMARK_ESTROGEN_RESPONSE_LATE']],  rownames(tum)))
er_sets$combined <- union(er_sets$early, er_sets$late)
tum <- AddModuleScore(tum, features = er_sets, name = 'ERset_')   # _1 early, _2 late, _3 combined

##------------------------------------------
## Save the per-cell table used in
## 'Figure_panels/Figure S1.R'
##------------------------------------------
tf_df <- data.frame(
  cell        = colnames(tum),
  Sample      = as.character(tum$orig.ident),
  Status      = as.character(tum$Status_new),
  Group       = as.character(tum$Group2),
  ESR1_status = as.character(tum$ESR1_status),
  FOSL1_act   = act_mat['FOSL1', ],
  ESR1_act    = act_mat['ESR1', ],
  ER_program  = tum$ER_program1,
  ER_early    = tum$ERset_1,
  ER_late     = tum$ERset_2,
  ER_combined = tum$ERset_3
)
write.csv(tf_df, file = file.path(output, 'decoupleR_tumor_cells.csv'), row.names = FALSE)

## Per-sample mean activities (unit of the baseline versus progressor test)
tf_sample <- tf_df %>%
  group_by(Sample, Group) %>%
  summarise(n_cells     = n(),
            FOSL1_act   = mean(FOSL1_act),
            ESR1_act    = mean(ESR1_act),
            ER_combined = mean(ER_combined),
            .groups = 'drop')
tf_sample
write.csv(tf_sample, file = file.path(output, 'decoupleR_per_sample_means.csv'), row.names = FALSE)
