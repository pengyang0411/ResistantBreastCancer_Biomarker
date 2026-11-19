##------------------------------------
## pySCENIC step 3/3: regulon activity
## (AUCell), baseline versus progressor
## Author: Peng Yang
##------------------------------------
library(dplyr)
library(tidyr)

##------------------------------------------
## Load the AUCell matrix (cells x regulons,
## e.g. "FOSL1(+)") and the cell information
##------------------------------------------

## Run from the repository root, after 'TF_activity/run_pyscenic.sh'
py_dir <- 'TF_activity/pyscenic'
output <- 'TF_activity/pyscenic_output'
dir.create(output, showWarnings = FALSE, recursive = TRUE)

auc <- read.csv(file.path(py_dir, 'auc.csv'), row.names = 1, check.names = FALSE)
ci  <- read.csv(file.path(py_dir, 'cellinfo.csv'), stringsAsFactors = FALSE)
dim(auc)

ci <- ci[match(rownames(auc), ci$CellID), ]
stopifnot(all(ci$CellID == rownames(auc)))
ci$Group <- factor(ci$Group, levels = c('Baseline', 'Progressor'))

## Regulons recovered for FOSL1 and ESR1
grep('^FOSL1', colnames(auc), value = TRUE)
grep('^ESR1',  colnames(auc), value = TRUE)   # no ESR1 regulon was recovered

##------------------------------------------
## All regulons: per-sample mean activity,
## Wilcoxon rank-sum test (baseline versus
## progressor) with BH adjustment
##------------------------------------------
all_res <- lapply(colnames(auc), function(rg) {
  d  <- data.frame(Sample = ci$Sample, Group = ci$Group, act = auc[[rg]])
  ps <- d %>% group_by(Sample, Group) %>%
    summarise(m = mean(act, na.rm = TRUE), .groups = 'drop')
  if (length(unique(ps$Group)) < 2) return(NULL)
  wt <- suppressWarnings(wilcox.test(m ~ Group, data = ps, exact = FALSE))
  data.frame(regulon = rg,
             median_Baseline   = median(ps$m[ps$Group == 'Baseline']),
             median_Progressor = median(ps$m[ps$Group == 'Progressor']),
             diff = median(ps$m[ps$Group == 'Progressor']) - median(ps$m[ps$Group == 'Baseline']),
             p = wt$p.value)
}) %>% bind_rows()
all_res$padj <- p.adjust(all_res$p, method = 'BH')
all_res <- all_res %>% arrange(desc(diff))

head(all_res, 15)
all_res[grepl('^FOSL1', all_res$regulon), ]

write.csv(all_res, file = file.path(output, 'SCENIC_allRegulons_baseline_vs_progressor.csv'), row.names = FALSE)

## Per-sample FOSL1 regulon activity (Supplementary Figure 1Q)
fos_col <- grep('^FOSL1', colnames(auc), value = TRUE)[1]
fos_sample <- data.frame(Sample = ci$Sample, Group = ci$Group, act = auc[[fos_col]]) %>%
  group_by(Sample, Group) %>%
  summarise(mean_activity = mean(act, na.rm = TRUE), n_cells = n(), .groups = 'drop')
fos_sample
write.csv(fos_sample, file = file.path(output, 'SCENIC_FOSL1_regulon_per_sample.csv'), row.names = FALSE)
