##-------------------------------
## Figure Panel S1
## Author: Peng Yang
##-------------------------------
library(ggplot2)
library(dplyr)
library(tidyr)
library(Seurat)
library(patchwork)
library(openxlsx)

## Run from the repository root; processed objects are in Data/
load('Data/merged_data_Harmony_new.RData')
DefaultAssay(merged_combined) <- 'RNA'
Status <- rep('Baseline', length(merged_combined$orig.ident))
Status[which(merged_combined$orig.ident %in% c('PA1', 'PA3 #1', 'PA3 #2',
                                               'PA126', 'PA130', 'PA139', 'PA99', 'PA165')) ] = 'Late Progressor'
Status[which(merged_combined$orig.ident %in% c('PA11', 'PA46', 'PA144')) ] = 'Early Progressor'
merged_combined$Status_new <- as.factor(Status)

umapColor <- c("#6DCCDD","#9ACBDE","#C8CADF","#EDC6DD","#F092B1","#F27FA5",
               "#F47892", "#F6A395","#F8AD77","#E7B080","#C9AFA2","#ABADC5",
               "#AEB9AC","#B9C984", "#C4DA5D")
SampleNames <- c("PA3",   "PA5",   "PA45",  "PA95",  "PA120", "PA110", "PA131", "PA153",
                 'PA11', 'PA46', 'PA144',
                 'PA1', 'PA3 #1', 'PA3 #2', 'PA126',
                 'PA139', 'PA99', 'PA165')
merged_combined$orig.ident.num <- factor(match(as.character(merged_combined$orig.ident), SampleNames),
                                         levels = 1:length(SampleNames))
global_labels <- paste(1:length(SampleNames), SampleNames)
global_values <- colorRampPalette(umapColor)(length(SampleNames))

grp_cols <- c('Baseline' = 'darkblue', 'Progressor' = 'red')

dir.create('FigureS1', showWarnings = FALSE)

##-------------------------------
## Figure S1A: Tumor versus
## non-tumor cells across clinical
## groups (Fisher's exact test)
##-------------------------------

tab_s1a <- table(Status = merged_combined$Status_new,
                 Condition = merged_combined$Tumor_cells)[, c('tumor cells', 'non-tumor cells')]
tab_s1a
round(100 * prop.table(tab_s1a, 1), 1)

## p-values: one-sided (greater); odds ratios and 95% CI: two-sided
fisher.test(tab_s1a[c('Baseline', 'Early Progressor'), ], alternative = 'greater')$p.value
fisher.test(tab_s1a[c('Baseline', 'Early Progressor'), ])[c('estimate', 'conf.int')]
fisher.test(tab_s1a[c('Baseline', 'Late Progressor'), ], alternative = 'greater')$p.value
fisher.test(tab_s1a[c('Baseline', 'Late Progressor'), ])[c('estimate', 'conf.int')]

df_tumor <- as.data.frame(tab_s1a)
colnames(df_tumor) <- c('Status', 'Condition', 'Count')
df_tumor$Condition <- factor(ifelse(df_tumor$Condition == 'tumor cells', 'Tumor Cells', 'Non-Tumor Cells'),
                             levels = c('Tumor Cells', 'Non-Tumor Cells'))

ptS1a <- ggplot(df_tumor, aes(fill = Condition, y = Count, x = Status)) +
  geom_bar(position = "fill", stat = "identity") +
  scale_fill_manual(values = c('#FD9AA0', '#6DCCFD')) +
  xlab(" ") + ylab(" ")  +
  theme_classic() +
  theme(
    text = element_text(size = 15),
    axis.text.x = element_text(angle = 60, vjust = 1, hjust = 1)
  )

pdf('FigureS1/FigureS1a.pdf', width = 4, height = 3.5)
print(ptS1a)
dev.off()

##-------------------------------
## Figure S1B-G and S1L-N: UMAP
## of samples (B-D), ESR1 (E-G) and
## FOSL1 (L-N) in each clinical group
##-------------------------------

panel_letters <- list('Baseline'         = c('b', 'e', 'l'),
                      'Early Progressor' = c('c', 'f', 'm'),
                      'Late Progressor'  = c('d', 'g', 'n'))

for(status in names(panel_letters)){

  combined_status <- subset(merged_combined, subset = Status_new == status)
  combined_status$orig.ident.num <- droplevels(combined_status$orig.ident.num)
  ids <- as.numeric(levels(combined_status$orig.ident.num))

  ## Samples
  ptS1_sample <- DimPlot(combined_status, reduction = 'umap.pca', group.by = 'orig.ident.num',
                         label = T, repel = F, label.size = 5) +
    NoAxes() + ggtitle(status) +
    scale_x_reverse() +
    theme_void() +
    scale_color_manual(labels = global_labels[ids],
                       values = global_values[ids]) +
    guides(color = guide_legend(override.aes = list(shape = 15, size = 6),
                                label.hjust = 0,
                                label.theme = element_text())) +
    theme(plot.title = element_text(face = 'bold', hjust = 0.5))

  ## ESR1
  ptS1_esr1 <- FeaturePlot(
    combined_status, features = 'ESR1', reduction = 'umap.pca',
    order = TRUE, pt.size = 0.1, max.cutoff = 'q95',
    cols = c('grey90', '#3B82F6')
  ) + ggtitle('ESR1') + xlab(' ') + ylab(' ') + scale_x_reverse()

  ## FOSL1
  ptS1_fosl1 <- FeaturePlot(
    combined_status, features = 'FOSL1', reduction = 'umap.pca',
    order = TRUE, pt.size = 0.1, max.cutoff = 'q95',
    cols = c('grey90', '#EF4444')
  ) + ggtitle('FOSL1') + xlab(' ') + ylab(' ') + scale_x_reverse()

  pdf(paste0('FigureS1/FigureS1', panel_letters[[status]][1], '.pdf'), width = 6, height = 5)
  print(ptS1_sample)
  dev.off()

  pdf(paste0('FigureS1/FigureS1', panel_letters[[status]][2], '.pdf'), width = 5, height = 4)
  print(ptS1_esr1)
  dev.off()

  pdf(paste0('FigureS1/FigureS1', panel_letters[[status]][3], '.pdf'), width = 5, height = 4)
  print(ptS1_fosl1)
  dev.off()
}

##-------------------------------
## Figure S1P: Patient-level ESR1+
## and FOSL1+ tumor-cell fractions
## (Wilcoxon rank-sum test)
##-------------------------------

## Samples without tumor cells (PA1, PA46) do not contribute
table(merged_combined$orig.ident, merged_combined$Tumor_cells)

tum <- subset(merged_combined, subset = Tumor_cells == 'tumor cells')
cell_df <- data.frame(
  sample    = as.character(tum$orig.ident),
  group     = ifelse(tum$Status_new == 'Baseline', 'Baseline', 'Progressor'),
  ESR1_pos  = tum@assays$RNA@counts['ESR1', ]  > 0,
  FOSL1_pos = tum@assays$RNA@counts['FOSL1', ] > 0
)

## Per-sample fractions (the biologically independent unit)
frac_df <- cell_df %>%
  group_by(sample, group) %>%
  summarise(n_tumor        = n(),
            ESR1_pos_frac  = mean(ESR1_pos),
            FOSL1_pos_frac = mean(FOSL1_pos),
            .groups = 'drop')
frac_df
write.csv(frac_df, file = 'FigureS1/FigureS1p_per_sample_fractions.csv', row.names = F)

wt_esr1  <- wilcox.test(frac_df$ESR1_pos_frac[frac_df$group == 'Baseline'],
                        frac_df$ESR1_pos_frac[frac_df$group == 'Progressor'], exact = FALSE)
wt_fosl1 <- wilcox.test(frac_df$FOSL1_pos_frac[frac_df$group == 'Baseline'],
                        frac_df$FOSL1_pos_frac[frac_df$group == 'Progressor'], exact = FALSE)
wt_esr1
wt_fosl1

marker_labs <- c(ESR1_pos_frac = 'ESR1+ fraction', FOSL1_pos_frac = 'FOSL1+ fraction')
plot_long <- frac_df %>%
  pivot_longer(c(ESR1_pos_frac, FOSL1_pos_frac), names_to = 'marker', values_to = 'frac') %>%
  mutate(marker = recode(marker, !!!marker_labs),
         group  = factor(group, levels = c('Baseline', 'Progressor')))
pval_lab <- data.frame(marker = c('ESR1+ fraction', 'FOSL1+ fraction'),
                       label  = paste0('Wilcoxon p = ', signif(c(wt_esr1$p.value, wt_fosl1$p.value), 3)))

ptS1p <- ggplot(plot_long, aes(group, frac)) +
  geom_boxplot(aes(fill = group), outlier.shape = NA, width = 0.55, alpha = 0.35) +
  geom_jitter(aes(color = group), width = 0.12, size = 2, alpha = 0.9) +
  geom_text(data = pval_lab, aes(x = 1.5, y = Inf, label = label),
            inherit.aes = FALSE, vjust = 1.6, size = 4.2) +
  facet_wrap(~ marker, scales = 'free_y') +
  scale_fill_manual(values = grp_cols) +
  scale_color_manual(values = grp_cols) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.18))) +
  labs(x = NULL, y = 'Fraction of tumor cells',
       title = 'Patient-level ESR1+ and FOSL1+ fractions') +
  theme_classic(base_size = 13) +
  theme(strip.text  = element_text(face = 'bold'),
        plot.title  = element_text(face = 'bold', hjust = 0.5),
        legend.position = 'none')

pdf('FigureS1/FigureS1p.pdf', width = 7, height = 4)
print(ptS1p)
dev.off()

rm(merged_combined, tum); gc()

##-------------------------------
## Figure S1H-K: Non-tumor
## compartments, marker genes and
## ESR1+ cells
##-------------------------------

## Non-tumor cells with Azimuth labels (CT_L1, CT_L2)
load('Data/combined_nonTumor_final.RData')
DefaultAssay(combined_nontumor) <- 'RNA'

## CT_L1 labels grouped into compartments
compartment_map <- c(
  "CD4 T" = "Immune (T)", "CD8 T" = "Immune (T)", "other T" = "Immune (T)",
  "B" = "Immune (B)", "NK" = "Immune (NK)",
  "Mono" = "Immune (Myeloid)", "DC" = "Immune (Myeloid)",
  "CAFs" = "Fibroblast/Stromal (CAF)",
  "PVL"  = "Perivascular (PVL)",
  "Endothelial" = "Endothelial",
  "other" = "Other"
)
## unname(): a named vector would be matched to cell barcodes by Seurat
combined_nontumor$Compartment <- unname(compartment_map[as.character(combined_nontumor$CT_L1)])
combined_nontumor$Compartment[is.na(combined_nontumor$Compartment)] <- "Other"

combined_nontumor$ESR1_pos <- combined_nontumor@assays$RNA@counts["ESR1", ] > 0
combined_nontumor$ESR1_status <- ifelse(combined_nontumor$ESR1_pos, "ESR1+", "ESR1-")

## Figure S1H: ESR1+ non-tumor cells
ptS1h <- DimPlot(combined_nontumor, reduction = 'umap.harmony', group.by = "ESR1_status",
                 order = "ESR1+", cols = c("ESR1-" = "grey85", "ESR1+" = "red")) +
  ggtitle("ESR1+ non-tumor cells") +
  theme(plot.title = element_text(face = "bold", hjust = 0.5))

pdf('FigureS1/FigureS1h.pdf', width = 7, height = 5.5)
print(ptS1h)
dev.off()

## Figure S1I: Non-tumor compartments
ptS1i <- DimPlot(combined_nontumor, reduction = 'umap.harmony', group.by = "CT_L1",
                 label = TRUE, repel = TRUE) +
  ggtitle("Non-tumor compartments") +
  theme(plot.title = element_text(face = "bold", hjust = 0.5))

pdf('FigureS1/FigureS1i.pdf', width = 7, height = 5.5)
print(ptS1i)
dev.off()

## Figure S1J: Marker genes per compartment (FindAllMarkers)
if (inherits(combined_nontumor[["RNA"]], "Assay5"))
  combined_nontumor <- JoinLayers(combined_nontumor)
combined_nontumor <- NormalizeData(combined_nontumor, verbose = FALSE)
Idents(combined_nontumor) <- factor(as.character(combined_nontumor$Compartment))

all_markers <- FindAllMarkers(combined_nontumor, assay = "RNA",
                              only.pos = FALSE, min.pct = 0.1,
                              logfc.threshold = 0.25, verbose = FALSE)

top10_up <- all_markers %>%
  filter(p_val_adj < 0.05) %>%
  group_by(cluster) %>% slice_max(avg_log2FC, n = 10, with_ties = FALSE) %>%
  ungroup() %>% arrange(cluster, desc(avg_log2FC))
top10_down <- all_markers %>%
  filter(p_val_adj < 0.05) %>%
  group_by(cluster) %>% slice_min(avg_log2FC, n = 10, with_ties = FALSE) %>%
  ungroup() %>% arrange(cluster, avg_log2FC)

wb <- createWorkbook()
addWorksheet(wb, "All_markers");   writeData(wb, "All_markers", all_markers)
addWorksheet(wb, "Top10_up");      writeData(wb, "Top10_up", top10_up)
addWorksheet(wb, "Top10_down");    writeData(wb, "Top10_down", top10_down)
saveWorkbook(wb, 'FigureS1/FigureS1j_nonTumor_compartment_markers.xlsx', overwrite = TRUE)

## Top 5 positive markers per compartment
top5 <- all_markers %>%
  filter(p_val_adj < 0.05, avg_log2FC > 0) %>%
  group_by(cluster) %>% slice_max(avg_log2FC, n = 5, with_ties = FALSE) %>%
  ungroup()
dp_genes <- unique(top5$gene)

ptS1j <- DotPlot(combined_nontumor, features = dp_genes, group.by = "Compartment") +
  RotatedAxis() +
  labs(title = "Top marker genes per non-tumor compartment", x = NULL, y = NULL) +
  theme(plot.title = element_text(face = "bold", hjust = 0.5),
        axis.text.x = element_text(size = 7))

pdf('FigureS1/FigureS1j.pdf', width = max(10, length(dp_genes) * 0.22), height = 5.5)
print(ptS1j)
dev.off()

## Figure S1K: ESR1+ rate per compartment
esr1_by_compartment <- data.frame(
  Compartment = combined_nontumor$Compartment,
  ESR1_pos    = combined_nontumor$ESR1_pos
) %>%
  group_by(Compartment) %>%
  summarise(n_cells            = n(),
            n_ESR1pos          = sum(ESR1_pos),
            pct_ESR1pos_within = round(100 * mean(ESR1_pos), 3),
            .groups = "drop") %>%
  mutate(pct_of_all_ESR1pos = round(100 * n_ESR1pos / sum(n_ESR1pos), 2)) %>%
  arrange(desc(n_ESR1pos))
esr1_by_compartment
write.csv(esr1_by_compartment, file = 'FigureS1/FigureS1k_ESR1pos_by_compartment.csv', row.names = F)

bar_df <- esr1_by_compartment %>%
  mutate(Compartment = factor(Compartment, levels = Compartment[order(pct_ESR1pos_within)]))

ptS1k <- ggplot(bar_df, aes(x = Compartment, y = pct_ESR1pos_within)) +
  geom_col(fill = "#3B6FB6", color = "black", linewidth = 0.3, width = 0.7) +
  geom_text(aes(label = sprintf("%.1f%%", pct_ESR1pos_within)), hjust = -0.15, size = 3.4) +
  coord_flip() +
  scale_y_continuous(expand = expansion(mult = c(0, 0.18))) +
  labs(x = NULL, y = "ESR1+ cells within compartment (%)",
       title = "ESR1+ rate per non-tumor compartment") +
  theme_classic(base_size = 12) +
  theme(plot.title = element_text(face = "bold", hjust = 0.5))

pdf('FigureS1/FigureS1k.pdf', width = 5.5, height = 4.5)
print(ptS1k)
dev.off()

rm(combined_nontumor); gc()

##-------------------------------
## Figure S1O: Per-cell FOSL1 ULM
## activity versus ER-program
## activity
##-------------------------------

## Generated by 'TF_activity/decoupleR_ULM.R'
tf_df <- read.csv('TF_activity/decoupleR_output/decoupleR_tumor_cells.csv', stringsAsFactors = FALSE)
tf_df$Group <- factor(tf_df$Group, levels = c('Baseline', 'Progressor'))

sp <- cor.test(tf_df$FOSL1_act, tf_df$ER_program, method = 'spearman')
sp

ptS1o <- ggplot(tf_df, aes(ER_program, FOSL1_act)) +
  geom_point(aes(color = Group), size = 0.3, alpha = 0.3) +
  geom_smooth(method = "lm", se = FALSE, color = "black", linewidth = 0.7) +
  scale_color_manual(values = grp_cols) +
  annotate("text", x = Inf, y = Inf, hjust = 1.05, vjust = 1.5,
           label = sprintf("Spearman rho = %.2f", unname(sp$estimate)), size = 4) +
  labs(title = "FOSL1 vs. ER program activity (per cell)",
       x = "ER-program activity (HALLMARK estrogen response)",
       y = "FOSL1 ULM activity") +
  theme_classic(base_size = 13) +
  theme(plot.title = element_text(face = "bold", hjust = 0.5))

pdf('FigureS1/FigureS1o.pdf', width = 5.5, height = 4.5)
print(ptS1o)
dev.off()

##-------------------------------
## Figure S1Q-R: FOSL1 regulon
## activity and ER-program activity,
## baseline versus progressor
##-------------------------------

## Violin = per-cell values; dots = per-sample means;
## Wilcoxon rank-sum test on the per-sample means
sample_violin <- function(d, ylab, title_text) {
  ps <- d %>% group_by(Sample, Group) %>%
    summarise(m = mean(val, na.rm = TRUE), .groups = "drop")
  wt <- wilcox.test(m ~ Group, data = ps, exact = FALSE)
  print(ps)
  print(wt)
  ggplot(d, aes(Group, val)) +
    geom_violin(aes(fill = Group), scale = "width", alpha = 0.3, color = NA) +
    geom_boxplot(aes(fill = Group), width = 0.12, outlier.shape = NA, alpha = 0.5) +
    geom_point(data = ps, aes(Group, m), size = 2.6, color = "black") +
    annotate("text", x = 1.5, y = max(d$val, na.rm = TRUE),
             label = paste0("Wilcoxon p = ", signif(wt$p.value, 3)), vjust = 1.3, size = 4) +
    scale_fill_manual(values = grp_cols) +
    labs(title = title_text,
         subtitle = "Violin: per-cell analysis   Dots: per-sample means",
         x = NULL, y = ylab) +
    theme_classic(base_size = 13) +
    theme(legend.position = "none",
          plot.title = element_text(face = "bold", hjust = 0.5),
          plot.subtitle = element_text(hjust = 0.5, size = 9))
}

## Figure S1Q: FOSL1 regulon (pySCENIC/AUCell), from 'TF_activity/run_pyscenic.sh'
auc <- read.csv('TF_activity/pyscenic/auc.csv', row.names = 1, check.names = FALSE)
ci  <- read.csv('TF_activity/pyscenic/cellinfo.csv', stringsAsFactors = FALSE)
ci  <- ci[match(rownames(auc), ci$CellID), ]
stopifnot(all(ci$CellID == rownames(auc)))
fos_col <- grep('^FOSL1', colnames(auc), value = TRUE)[1]
fos_col

d_q <- data.frame(Sample = ci$Sample,
                  Group  = factor(ci$Group, levels = c('Baseline', 'Progressor')),
                  val    = auc[[fos_col]])
ptS1q <- sample_violin(d_q, "Regulon activity (AUCell)",
                       "FOSL1 regulon activity: baseline vs progressor")

pdf('FigureS1/FigureS1q.pdf', width = 6.5, height = 4.4)
print(ptS1q)
dev.off()

## Figure S1R: ER program (HALLMARK estrogen response early + late)
d_r <- data.frame(Sample = tf_df$Sample, Group = tf_df$Group, val = tf_df$ER_combined)
ptS1r <- sample_violin(d_r, "ER program score (combined)",
                       "ER-program activity: Baseline vs Progressor")

pdf('FigureS1/FigureS1r.pdf', width = 5, height = 4.4)
print(ptS1r)
dev.off()

## ESR1 ULM activity (DoRothEA), baseline versus progressor (reported in the Results)
ps_esr1 <- tf_df %>% group_by(Sample, Group) %>%
  summarise(m = mean(ESR1_act, na.rm = TRUE), .groups = "drop")
wilcox.test(m ~ Group, data = ps_esr1, exact = FALSE)
