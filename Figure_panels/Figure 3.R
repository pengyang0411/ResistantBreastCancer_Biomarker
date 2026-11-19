##-------------------------------
## Figure Panel 3
## Author: Peng Yang
##-------------------------------
library(ggplot2)
library(dplyr)
library(Seurat)
library(patchwork)

## Run from the repository root; processed objects are in Data/
load('Data/merged_data_Harmony_new.RData')
DefaultAssay(merged_combined) <- 'RNA'
Status <- rep('Baseline', length(merged_combined$orig.ident))
Status[which(merged_combined$orig.ident %in% c('PA1', 'PA3 #1', 'PA3 #2',
                                               'PA126', 'PA130', 'PA139', 'PA99', 'PA165')) ] = 'Late Progressor'
Status[which(merged_combined$orig.ident %in% c('PA11', 'PA46', 'PA144')) ] = 'Early Progressor'
merged_combined$Status_new <- as.factor(Status)

merged_combined_tumor <- subset(merged_combined, subset = Tumor_cells == 'tumor cells')
rm(merged_combined); gc()

dir.create('Figure3', showWarnings = FALSE)

##-------------------------------
## Figure 3V: UMAP of FOSL1,
## CCNA2, CCND1 and ESR1 in tumor
## cells, with sample-level
## Pearson correlation
##-------------------------------

## Sample x gene percentage of positive tumor cells,
## generated in 'Figure_panels/Figure 1.R' (Figure 1G)
load('Data/tumor_gene_expr.RData')   # -> pct_mat

## Gene pairs (first gene, second gene, colors)
pairs <- list(
  c('CCNA2', 'FOSL1', 'orange',  '#EF4444'),
  c('CCND1', 'FOSL1', 'purple',  '#EF4444'),
  c('ESR1',  'CCNA2', '#3B82F6', 'orange'),
  c('ESR1',  'CCND1', '#3B82F6', 'purple')
)

res_cor_3v <- c()
pt3v <- list()
for(i in 1:length(pairs)){

  g1 <- pairs[[i]][1]
  g2 <- pairs[[i]][2]

  ## Pearson correlation between per-sample positive fractions (n = 16 samples)
  res_cor <- cor.test(pct_mat[, g1], pct_mat[, g2], method = 'pearson')
  res_cor_3v <- rbind(res_cor_3v, data.frame(gene_1 = g1, gene_2 = g2, n_samples = nrow(pct_mat),
                                             pearson_r = unname(res_cor$estimate), p_value = res_cor$p.value))

  ## Blended UMAP of the two genes
  pt_blend <- FeaturePlot(
    merged_combined_tumor, features = c(g1, g2), reduction = 'umap.pca',
    blend = T, order = T, blend.threshold = 0.1, max.cutoff = 'q95',
    pt.size = 0.1, cols = c('grey90', pairs[[i]][3], pairs[[i]][4])
  )

  pt3v[[i]] <- pt_blend[[3]] + xlab(' ') + ylab(' ') + scale_x_reverse() +
    ggtitle(paste0(g1, '+ / ', g2, '+\nPearson correlation = ', round(res_cor$estimate, 4),
                   '\np=', round(res_cor$p.value, 4)))
}
res_cor_3v
write.csv(res_cor_3v, file = 'Figure3/Figure3v_pearson.csv', row.names = F)

pt3v_all <- (pt3v[[1]] + pt3v[[2]]) / (pt3v[[3]] + pt3v[[4]])

pdf('Figure3/Figure3v.pdf', width = 8, height = 8.5)
print(pt3v_all)
dev.off()
