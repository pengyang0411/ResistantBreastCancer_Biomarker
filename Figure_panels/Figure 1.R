##-------------------------------
## Figure Panel 1
## Author: Peng Yang
##-------------------------------
library(ggplot2)
library(dplyr)
library(tidyr)
library(Seurat)
library(patchwork)
library(readxl)
library(ggalluvial)

## Run from the repository root; processed objects are in Data/
load('Data/merged_data_Harmony_new.RData')
DefaultAssay(merged_combined) <- 'RNA'
Status <- rep('Baseline', length(merged_combined$orig.ident))
Status[which(merged_combined$orig.ident %in% c('PA1', 'PA3 #1', 'PA3 #2',
                                               'PA126', 'PA130', 'PA139', 'PA99', 'PA165')) ] = 'Late Progressor'
Status[which(merged_combined$orig.ident %in% c('PA11', 'PA46', 'PA144')) ] = 'Early Progressor'
merged_combined$Status_new <- as.factor(Status)

## Sample numbering used in Figure 1B (PA130 is not in the final object)
SampleNames <- c("PA3",   "PA5",   "PA45",  "PA95",  "PA120", "PA110", "PA131", "PA153",
                 'PA11', 'PA46', 'PA144',
                 'PA1', 'PA3 #1', 'PA3 #2', 'PA126',
                 'PA139', 'PA99', 'PA165')
merged_combined$orig.ident.num <- factor(match(as.character(merged_combined$orig.ident), SampleNames))
Idents(merged_combined) <- merged_combined$orig.ident.num

## Cutoff for positive versus negative (raw counts > 0)
cutoff_postive <- 0
merged_combined$ESR1_Status <- factor(ifelse(merged_combined@assays$RNA@counts['ESR1', ] > cutoff_postive,
                                             'ESR1+', 'ESR1-'), levels = c('ESR1+', 'ESR1-'))
merged_combined$FOSL1_Status <- factor(ifelse(merged_combined@assays$RNA@counts['FOSL1', ] > cutoff_postive,
                                              'FOSL1+', 'FOSL1-'), levels = c('FOSL1+', 'FOSL1-'))

merged_combined_tumor    <- subset(merged_combined, subset = Tumor_cells == 'tumor cells')
merged_combined_nontumor <- subset(merged_combined, subset = Tumor_cells == 'non-tumor cells')

dir.create('Figure1', showWarnings = FALSE)

##-------------------------------
## Figure 1B: UMAP of ESR1
## expression in tumor and
## non-tumor cells
##-------------------------------

pt1b_tumor <- FeaturePlot(
  merged_combined_tumor, features = 'ESR1', reduction = 'umap.pca',
  order = TRUE, pt.size = 0.1, max.cutoff = 'q95',
  cols = c('grey90', '#3B82F6')
) + ggtitle('Tumor cells') + xlab(' ') + ylab(' ') + scale_x_reverse() + NoLegend()
pt1b_tumor <- LabelClusters(pt1b_tumor, id = 'ident', repel = TRUE, size = 3.5)

pt1b_nontumor <- FeaturePlot(
  merged_combined_nontumor, features = 'ESR1', reduction = 'umap.pca',
  order = TRUE, pt.size = 0.1, max.cutoff = 'q95',
  cols = c('grey90', '#3B82F6')
) + ggtitle('Non-tumor cells') + xlab(' ') + ylab(' ') + scale_x_reverse() + NoLegend()
pt1b_nontumor <- LabelClusters(pt1b_nontumor, id = 'ident', repel = TRUE, size = 3.5)

pdf('Figure1/Figure1b_tumor.pdf', height = 4, width = 4)
print(pt1b_tumor)
dev.off()

pdf('Figure1/Figure1b_nontumor.pdf', height = 4, width = 4)
print(pt1b_nontumor)
dev.off()

##-------------------------------
## Figure 1C: Fraction of ESR1+
## cells in each sample, all
## cells and tumor cells
##-------------------------------

## All the cells
df_ESR1 <- c()
for(sample in unique(as.character(merged_combined$orig.ident))){

  df_ESR1 <- rbind(df_ESR1, c(sample, as.character(merged_combined$Status_new[merged_combined$orig.ident == sample][1]),
                              sum(merged_combined$ESR1_Status[merged_combined$orig.ident == sample] == 'ESR1+'), 'ESR1+'))
  df_ESR1 <- rbind(df_ESR1, c(sample, as.character(merged_combined$Status_new[merged_combined$orig.ident == sample][1]),
                              sum(merged_combined$ESR1_Status[merged_combined$orig.ident == sample] == 'ESR1-'), 'ESR1-'))
}

colnames(df_ESR1) <- c('Sample', 'Status', 'Count', 'Condition')
df_ESR1 <- as.data.frame(df_ESR1)
df_ESR1$Status <- as.factor(df_ESR1$Status)
df_ESR1$Sample <- as.factor(df_ESR1$Sample)
df_ESR1$Count <- as.numeric(df_ESR1$Count)
df_ESR1$Condition <- factor(df_ESR1$Condition, levels = c('ESR1+', 'ESR1-'))

pt1c_all <- ggplot(df_ESR1, aes(fill = Condition, y = Count, x = Sample)) +
  geom_bar(position = "fill", stat = "identity") +
  facet_grid(~ Status, drop = T, scales = 'free', space = 'free') +
  scale_fill_manual(values = c('#FD9AA0', '#6DCCFD')) +
  xlab(" ") + ylab("Fraction in all cells")  +
  theme_classic() +
  theme(
    text = element_text(size = 15),
    axis.text.x = element_text(angle = 60, vjust = 1, hjust = 1),
    strip.background = element_blank(),
    strip.text.x = element_blank()
  )

pdf('Figure1/Figure1c_all_cells.pdf', width = 6.5, height = 2.5)
print(pt1c_all)
dev.off()

## Tumor cells
df_ESR1 <- c()
for(sample in unique(as.character(merged_combined$orig.ident))){

  idx <- merged_combined$orig.ident == sample & merged_combined$Tumor_cells == 'tumor cells'
  df_ESR1 <- rbind(df_ESR1, c(sample, as.character(merged_combined$Status_new[merged_combined$orig.ident == sample][1]),
                              sum(merged_combined$ESR1_Status[idx] == 'ESR1+'), 'ESR1+'))
  df_ESR1 <- rbind(df_ESR1, c(sample, as.character(merged_combined$Status_new[merged_combined$orig.ident == sample][1]),
                              sum(merged_combined$ESR1_Status[idx] == 'ESR1-'), 'ESR1-'))
}

colnames(df_ESR1) <- c('Sample', 'Status', 'Count', 'Condition')
df_ESR1 <- as.data.frame(df_ESR1)
df_ESR1$Status <- as.factor(df_ESR1$Status)
df_ESR1$Sample <- as.factor(df_ESR1$Sample)
df_ESR1$Count <- as.numeric(df_ESR1$Count)
df_ESR1$Condition <- factor(df_ESR1$Condition, levels = c('ESR1+', 'ESR1-'))

pt1c_tumor <- ggplot(df_ESR1, aes(fill = Condition, y = Count, x = Sample)) +
  geom_bar(position = "fill", stat = "identity") +
  facet_grid(~ Status, drop = T, scales = 'free', space = 'free') +
  scale_fill_manual(values = c('#FD9AA0', '#6DCCFD')) +
  xlab(" ") + ylab("Fraction in tumor cells")  +
  theme_classic() +
  theme(
    text = element_text(size = 15),
    axis.text.x = element_text(angle = 60, vjust = 1, hjust = 1),
    strip.background = element_blank(),
    strip.text.x = element_blank()
  )

pdf('Figure1/Figure1c_tumor_cells.pdf', width = 6.5, height = 2.5)
print(pt1c_tumor)
dev.off()

## Per-sample ESR1+ tumor-cell fraction (range reported in the Results)
df_ESR1_sample <- df_ESR1 %>%
  group_by(Sample, Status) %>%
  summarise(n_tumor = sum(Count),
            frac_ESR1_pos = sum(Count[Condition == 'ESR1+']) / sum(Count),
            .groups = 'drop') %>%
  filter(n_tumor > 0)
print(df_ESR1_sample)
range(df_ESR1_sample$frac_ESR1_pos)

##-------------------------------
## Figure 1D: ESR1+ versus ESR1-
## tumor cells across clinical
## groups (Fisher's exact test)
##-------------------------------

tab_1d <- table(Status = merged_combined_tumor$Status_new,
                Condition = merged_combined_tumor$ESR1_Status)
tab_1d
prop.table(tab_1d, 1)

## Fisher exact test between baseline and early progressor
fisher.test(tab_1d[c('Baseline', 'Early Progressor'), c('ESR1+', 'ESR1-')])
## Fisher exact test between baseline and late progressor
fisher.test(tab_1d[c('Baseline', 'Late Progressor'), c('ESR1+', 'ESR1-')])
## Fisher exact test between early and late progressor
fisher.test(tab_1d[c('Early Progressor', 'Late Progressor'), c('ESR1+', 'ESR1-')])

df_1d <- as.data.frame(tab_1d)
colnames(df_1d) <- c('Status', 'Condition', 'Count')
df_1d$Condition <- factor(df_1d$Condition, levels = c('ESR1-', 'ESR1+'))

pt1d <- ggplot(df_1d, aes(fill = Condition, y = Count, x = Status)) +
  geom_bar(position = "fill", stat = "identity") +
  scale_fill_manual(values = c('ESR1-' = 'red', 'ESR1+' = 'darkblue')) +
  xlab(" ") + ylab("Fraction in tumor cells")  +
  theme_classic() +
  theme(
    text = element_text(size = 15),
    axis.text.x = element_text(angle = 60, vjust = 1, hjust = 1)
  )

pdf('Figure1/Figure1d.pdf', width = 4, height = 4)
print(pt1d)
dev.off()

##-------------------------------
## Figure 1E: ESR1+ versus ESR1-
## tumor cells in the longitudinal
## biopsies of patient #3
##-------------------------------

pa3_levels <- c('PA3', 'PA3 #1', 'PA3 #2')
merged_combined_pa3 <- subset(merged_combined_tumor, subset = orig.ident %in% pa3_levels)
merged_combined_pa3$timepoint <- factor(as.character(merged_combined_pa3$orig.ident), levels = pa3_levels)

tab_1e <- table(Sample = merged_combined_pa3$timepoint,
                Condition = merged_combined_pa3$ESR1_Status)
tab_1e
prop.table(tab_1e, 1)

## Fisher exact test between PA3 (baseline) and each progression biopsy
fisher.test(tab_1e[c('PA3', 'PA3 #1'), c('ESR1+', 'ESR1-')])
fisher.test(tab_1e[c('PA3', 'PA3 #2'), c('ESR1+', 'ESR1-')])

df_1e <- as.data.frame(tab_1e)
colnames(df_1e) <- c('Sample', 'Condition', 'Count')
df_1e$Condition <- factor(df_1e$Condition, levels = c('ESR1-', 'ESR1+'))

pt1e <- ggplot(df_1e, aes(fill = Condition, y = Count, x = Sample)) +
  geom_bar(position = "fill", stat = "identity") +
  scale_fill_manual(values = c('ESR1-' = 'red', 'ESR1+' = 'darkblue')) +
  scale_x_discrete(labels = c('PA3', 'PA3#1', 'PA3#2')) +
  xlab(" ") + ylab("Fraction in tumor cells")  +
  theme_classic() +
  theme(
    text = element_text(size = 15),
    axis.text.x = element_text(angle = 60, vjust = 1, hjust = 1)
  )

pdf('Figure1/Figure1e.pdf', width = 3.5, height = 4)
print(pt1e)
dev.off()

##-------------------------------
## Figure 1F: Hallmark GSEA of
## ESR1- versus ESR1+ tumor cells
## in CDK4/6i progressors
##-------------------------------

## Input ranked list: DEG/DEG_Tumor_cells_ESR1_Progressors.csv, generated by
## 'Figure_panels/DE Genes.R' (FindMarkers, ESR1+ vs ESR1- tumor cells,
## early + late progressors). Pre-ranked Hallmark GSEA was run on this table.

##-------------------------------
## Figure 1G: Transcription
## factors negatively correlated
## with ESR1 (Supplementary Table 2)
##-------------------------------

## Percentage of positive tumor cells for each gene in each sample
expr <- merged_combined_tumor@assays$RNA@data
expr_bin <- expr > 0
samples <- as.character(merged_combined_tumor$orig.ident)
sample_levels <- unique(samples)

pct_mat <- matrix(
  0,
  nrow = length(sample_levels),
  ncol = nrow(expr),
  dimnames = list(sample_levels, rownames(expr))
)
for (s in sample_levels) {
  cells_s <- which(samples == s)
  pct_mat[s, ] <- Matrix::rowMeans(expr_bin[, cells_s, drop = FALSE]) * 100
}
rm(expr, expr_bin); gc()

## Sample x gene matrix, reused for Figure 1Q-R and Figure 3V
save(pct_mat, file = 'Data/tumor_gene_expr.RData')
write.csv(round(pct_mat, 4), file = 'Data/tumor_gene_positive_percentage.csv', row.names = T)

## Pearson correlation of every gene with ESR1 across samples
esr1_vec <- pct_mat[, 'ESR1']
r <- apply(pct_mat, 2, function(g) cor(esr1_vec, g, use = "pairwise.complete.obs"))
n <- apply(pct_mat, 2, function(g) sum(complete.cases(esr1_vec, g)))
tval <- r * sqrt((n - 2) / pmax(1e-12, 1 - r^2))
pval <- 2 * pt(-abs(tval), df = n - 2)
padj_all <- p.adjust(pval, method = "BH")

## Human transcription factors (Lambert et al., 2018)
hmtf <- read_xlsx('Data/human TFs v1.01.xlsx', col_names = 'gene')
genes_tf <- hmtf$gene[hmtf$gene %in% colnames(pct_mat)]

res_tf <- data.frame(
  TF             = genes_tf,
  cor            = as.numeric(r[genes_tf]),
  n              = as.integer(n[genes_tf]),
  p              = as.numeric(pval[genes_tf]),
  padj_all_genes = as.numeric(padj_all[genes_tf])   # BH across all detected genes
)
res_tf$padj_TF <- p.adjust(res_tf$p, method = "BH")  # BH across the transcription factors
res_tf <- res_tf[order(res_tf$cor), ]
length(genes_tf)
head(res_tf, 12)

write.csv(res_tf, file = 'Figure1/Supplementary_Table2_ESR1_TF_correlation.csv', row.names = F)

## Top 10 transcription factors negatively correlated with ESR1
top_tf <- res_tf %>% filter(TF != 'ESR1') %>% head(10)
top_tf$TF <- factor(top_tf$TF, levels = rev(top_tf$TF))

pt1g <- ggplot(top_tf, aes(x = 1, y = TF, fill = cor)) +
  geom_tile(color = 'white') +
  scale_fill_gradient(low = 'darkblue', high = '#AFB6E0',
                      name = 'Pearson correlation\ncoefficient') +
  xlab(' ') + ylab(' ') +
  ggtitle('Top 10 TFs with negative\ncorrelation to ESR1') +
  theme_classic() +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.line = element_blank(),
    axis.text.y = element_text(face = 'italic'),
    plot.title = element_text(hjust = 0.5)
  )

pdf('Figure1/Figure1g.pdf', width = 3, height = 4)
print(pt1g)
dev.off()

##-------------------------------
## Figure 1H-K: UMAP of FOSL1
## overlaid with ESR1 in tumor
## cells, with Fisher's exact test
## for ESR1-/FOSL1+ co-distribution
##-------------------------------

groups <- list()
groups[['All samples']] <- c('Baseline', 'Early Progressor', 'Late Progressor')
groups[['Baseline']] <- c('Baseline')
groups[['Early Progressor']] <- c('Early Progressor')
groups[['Late Progressor']] <- c('Late Progressor')

## Fisher exact test between ESR1 and FOSL1 status in tumor cells
res_fisher_1h <- c()
for(g in names(groups)){

  idx <- merged_combined_tumor$Status_new %in% groups[[g]]
  res_fisher <- fisher.test(table(ESR1 = merged_combined_tumor$ESR1_Status[idx],
                                  FOSL1 = merged_combined_tumor$FOSL1_Status[idx]))
  res_fisher_1h <- rbind(res_fisher_1h,
                         data.frame(Group = g,
                                    n_samples = length(unique(merged_combined_tumor$orig.ident[idx])),
                                    n_cells = sum(idx),
                                    odds_ratio = unname(res_fisher$estimate),
                                    p_value = res_fisher$p.value))
}
res_fisher_1h
write.csv(res_fisher_1h, file = 'Figure1/Figure1h_k_fisher.csv', row.names = F)

p_label <- function(p) ifelse(p < 1e-4, 'p<0.0001', paste0('p=', round(p, 4)))

pt1hk <- list()
for(g in names(groups)){

  combined_status_tumor <- subset(merged_combined_tumor, subset = Status_new %in% groups[[g]])

  pt_blend <- FeaturePlot(
    combined_status_tumor, features = c('ESR1', 'FOSL1'), reduction = 'umap.pca',
    blend = T, order = T, blend.threshold = 0.1, max.cutoff = 'q95',
    pt.size = 0.1, cols = c('grey90', '#3B82F6', '#EF4444')  # blue for ESR1, red for FOSL1
  )

  pt1hk[[g]] <- pt_blend[[3]] + xlab(' ') + ylab(' ') + scale_x_reverse() +
    ggtitle(paste0(g, '\n', p_label(res_fisher_1h$p_value[res_fisher_1h$Group == g])))
}

pdf('Figure1/Figure1h.pdf', height = 4, width = 4)
print(pt1hk[['All samples']])
dev.off()

pdf('Figure1/Figure1i.pdf', height = 4, width = 4)
print(pt1hk[['Baseline']])
dev.off()

pdf('Figure1/Figure1j.pdf', height = 4, width = 4)
print(pt1hk[['Early Progressor']])
dev.off()

pdf('Figure1/Figure1k.pdf', height = 4, width = 4)
print(pt1hk[['Late Progressor']])
dev.off()

##-------------------------------
## Figure 1N: FOSL1+ versus FOSL1-
## tumor cells in the longitudinal
## biopsies of patient #3
##-------------------------------

tab_1n <- table(Sample = merged_combined_pa3$timepoint,
                Condition = merged_combined_pa3$FOSL1_Status)
tab_1n
round(100 * prop.table(tab_1n, 1), 2)
round(100 * prop.table(tab_1e, 1), 2)

## Fisher exact test between PA3 (baseline) and each progression biopsy
fisher.test(tab_1n[c('PA3', 'PA3 #1'), c('FOSL1+', 'FOSL1-')])
fisher.test(tab_1n[c('PA3', 'PA3 #2'), c('FOSL1+', 'FOSL1-')])

df_1n <- as.data.frame(tab_1n)
colnames(df_1n) <- c('Sample', 'Condition', 'Count')
df_1n$Condition <- factor(df_1n$Condition, levels = c('FOSL1+', 'FOSL1-'))

pt1n <- ggplot(df_1n, aes(fill = Condition, y = Count, x = Sample)) +
  geom_bar(position = "fill", stat = "identity") +
  scale_fill_manual(values = c('FOSL1+' = 'red', 'FOSL1-' = 'darkblue')) +
  scale_x_discrete(labels = c('PA3', 'PA3#1', 'PA3#2')) +
  xlab(" ") + ylab("Fraction in tumor cells")  +
  theme_classic() +
  theme(
    text = element_text(size = 15),
    axis.text.x = element_text(angle = 60, vjust = 1, hjust = 1)
  )

pdf('Figure1/Figure1n.pdf', width = 3.5, height = 4)
print(pt1n)
dev.off()

##-------------------------------
## Figure 1O: Alluvial plot of the
## four FOSL1/ESR1 states across
## the patient #3 biopsies
##-------------------------------

## The three biopsies are independent samples; ribbons connect shared
## categories for readability and do not represent lineage tracing.
state_levels <- c("FOSL1+/ESR1-", "FOSL1+/ESR1+", "FOSL1-/ESR1+", "FOSL1-/ESR1-")
state_colors <- c("FOSL1+/ESR1-" = "#ff9966",
                  "FOSL1+/ESR1+" = "#ff99cc",
                  "FOSL1-/ESR1+" = "#66cccc",
                  "FOSL1-/ESR1-" = "#9999cc")

cell_df <- data.frame(timepoint = merged_combined_pa3$timepoint,
                      state = factor(paste0(as.character(merged_combined_pa3$FOSL1_Status), "/",
                                            as.character(merged_combined_pa3$ESR1_Status)),
                                     levels = state_levels))

state_df <- cell_df %>%
  dplyr::count(timepoint, state, name = "n_cells") %>%
  tidyr::complete(timepoint, state, fill = list(n_cells = 0)) %>%
  group_by(timepoint) %>%
  mutate(n_tumor  = sum(n_cells),
         Fraction = n_cells / n_tumor,
         Pct      = 100 * Fraction) %>%
  ungroup() %>%
  mutate(subject = as.numeric(state))
print(as.data.frame(state_df))

stratum_width <- 0.6
pt1o <- ggplot(state_df,
               aes(x = timepoint, stratum = state, alluvium = subject,
                   y = Fraction, fill = state, label = state)) +
  scale_x_discrete(expand = expansion(mult = .12, add = stratum_width / 2)) +
  geom_flow(width = stratum_width) +
  geom_stratum(alpha = .5, width = stratum_width) +
  geom_text(stat = "stratum", size = 2.8, min.y = 0.04) +
  scale_fill_manual(values = state_colors) +
  ylab("Fraction of tumor cells") + xlab(" ") +
  theme_void() +
  theme(legend.position = "none",
        axis.title.x    = element_text(),
        axis.text.x     = element_text())

## Call-out for the FOSL1+/ESR1+ stratum at PA3 #2 (below the min.y label cutoff)
callout_tp    <- "PA3 #2"
callout_state <- "FOSL1+/ESR1+"
a     <- state_df %>% dplyr::filter(timepoint == callout_tp)
a     <- a[match(state_levels, as.character(a$state)), ]
k     <- match(callout_state, state_levels)
y_mid <- 1 - sum(a$Fraction[seq_len(k - 1)]) - a$Fraction[k] / 2
x_pos <- match(callout_tp, levels(state_df$timepoint))
bar_half <- stratum_width / 2

pt1o <- pt1o +
  annotate("segment",
           x = x_pos + bar_half + 0.13, xend = x_pos + bar_half + 0.02,
           y = y_mid, yend = y_mid,
           arrow = grid::arrow(length = grid::unit(0.14, "cm"), type = "closed"),
           linewidth = .35) +
  annotate("text", x = x_pos + bar_half + 0.17, y = y_mid, hjust = 0,
           label = callout_state, size = 2.9, fontface = "bold") +
  coord_cartesian(clip = "off") +
  theme(plot.margin = margin(8, 74, 6, 10))

pdf('Figure1/Figure1o.pdf', width = 6, height = 4.5)
print(pt1o)
dev.off()

##-------------------------------
## Figure 1P: Ki67+ tumor-cell
## fractions by ESR1 or FOSL1
## status
##-------------------------------

meta_df <- merged_combined_tumor@meta.data
meta_df$ESR1_Status <- factor(as.character(meta_df$ESR1_Status), levels = c('ESR1-', 'ESR1+'))
meta_df$FOSL1_Status <- factor(as.character(meta_df$FOSL1_Status), levels = c('FOSL1-', 'FOSL1+'))
meta_df$Ki67_Positive <- merged_combined_tumor@assays$RNA@counts['MKI67', ] > cutoff_postive
meta_df$Ki67_Status <- factor(ifelse(meta_df$Ki67_Positive, 'Ki67+', 'Ki67-'), levels = c('Ki67-', 'Ki67+'))
meta_df$ProgressorGroup <- factor(ifelse(meta_df$Status_new == 'Baseline', 'Baseline', 'Progressor'),
                                  levels = c('Baseline', 'Progressor'))

## Fisher exact test p-value
get_fisher_p <- function(df, group_var, ki67_var = "Ki67_Status") {
  tab <- table(df[[group_var]], df[[ki67_var]])
  if (all(dim(tab) == c(2, 2))) {
    p <- fisher.test(tab)$p.value
  } else {
    p <- NA
  }
  p
}

## Ki67+ percentage
make_ki67_percent_df <- function(df, group_var) {
  df %>%
    group_by(.data[[group_var]]) %>%
    summarise(
      n_total = n(),
      n_ki67_pos = sum(Ki67_Positive, na.rm = TRUE),
      pct_ki67_pos = 100 * n_ki67_pos / n_total,
      .groups = "drop"
    ) %>%
    rename(Group = 1)
}

plot_ki67_percent <- function(df, group_var, title_text, fill_color = "darkred") {
  plot_df <- make_ki67_percent_df(df, group_var)
  pval <- get_fisher_p(df, group_var)

  p_label_ki67 <- ifelse(is.na(pval),
                         "Fisher's exact test: NA",
                         paste0("Fisher's exact test\np = ", signif(pval, 3)))

  ggplot(plot_df, aes(x = Group, y = pct_ki67_pos)) +
    geom_col(width = 0.55, fill = fill_color, color = "black", linewidth = 0.3) +
    ylim(0, 30) +
    labs(title = title_text, x = NULL, y = "Ki67+ cells (%)") +
    annotate("text", x = mean(seq_len(nrow(plot_df))), y = 28.5, label = p_label_ki67, size = 4) +
    theme_classic(base_size = 14) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5),
      axis.text.x = element_text(angle = 45, hjust = 1)
    )
}

df_baseline   <- meta_df %>% filter(ProgressorGroup == "Baseline")
df_progressor <- meta_df %>% filter(ProgressorGroup == "Progressor")
df_all        <- meta_df

pt1p_all        <- plot_ki67_percent(df_all,        "ESR1_Status",  "All samples")
pt1p_baseline   <- plot_ki67_percent(df_baseline,   "ESR1_Status",  "Baseline samples")
pt1p_prog_esr1  <- plot_ki67_percent(df_progressor, "ESR1_Status",  "Progressor samples")
pt1p_prog_fosl1 <- plot_ki67_percent(df_progressor, "FOSL1_Status", "Progressor samples")

pt1p <- (pt1p_all + pt1p_baseline) / (pt1p_prog_esr1 + pt1p_prog_fosl1)

pdf('Figure1/Figure1p.pdf', width = 6, height = 8)
print(pt1p)
dev.off()

## Ki67+ percentages shown in Figure 1P
bind_rows(
  make_ki67_percent_df(df_all,        "ESR1_Status")  %>% mutate(Setting = 'All samples'),
  make_ki67_percent_df(df_baseline,   "ESR1_Status")  %>% mutate(Setting = 'Baseline samples'),
  make_ki67_percent_df(df_progressor, "ESR1_Status")  %>% mutate(Setting = 'Progressor samples'),
  make_ki67_percent_df(df_progressor, "FOSL1_Status") %>% mutate(Setting = 'Progressor samples')
)

##-------------------------------
## Figure 1Q-R: Correlation of
## PCNA+ fraction with ESR1+
## (baseline) or FOSL1+ (late
## progressor) fractions
##-------------------------------

## Per-sample percentages of positive tumor cells (pct_mat, Figure 1G)
df_pct <- data.frame(Sample = rownames(pct_mat),
                     ESR1   = pct_mat[, 'ESR1'],
                     FOSL1  = pct_mat[, 'FOSL1'],
                     PCNA   = pct_mat[, 'PCNA'])
df_pct$Status <- as.character(merged_combined_tumor$Status_new)[match(df_pct$Sample,
                                                                      as.character(merged_combined_tumor$orig.ident))]
df_pct

## Pearson correlation. The p-values shown in Figure 1Q-R are one-sided
## (alternative = 'greater'); two-sided p-values are printed for reference.
df_1q <- df_pct[df_pct$Status == 'Baseline', ]
res_1q <- cor.test(df_1q$PCNA, df_1q$ESR1, method = 'pearson', alternative = 'greater')
res_1q
cor.test(df_1q$PCNA, df_1q$ESR1, method = 'pearson')$p.value

df_1r <- df_pct[df_pct$Status == 'Late Progressor', ]
res_1r <- cor.test(df_1r$PCNA, df_1r$FOSL1, method = 'pearson', alternative = 'greater')
res_1r
cor.test(df_1r$PCNA, df_1r$FOSL1, method = 'pearson')$p.value

pt1q <- ggplot(df_1q, aes(x = PCNA, y = ESR1)) +
  geom_point(size = 2) +
  geom_smooth(method = 'lm', se = FALSE, color = '#1F77B4') +
  annotate('text', x = 20, y = 78, hjust = 0, size = 3.5,
           label = paste0('Pearson r = ', round(res_1q$estimate, 2), '\np = ', round(res_1q$p.value, 2))) +
  scale_x_continuous(limits = c(0, 80)) +
  scale_y_continuous(limits = c(0, 80)) +
  labs(title = 'Baseline samples', x = 'PCNA+ (%)', y = 'ESR1+ (%)') +
  theme_classic(base_size = 13) +
  theme(plot.title = element_text(face = "bold", hjust = 0.5))

pt1r <- ggplot(df_1r, aes(x = PCNA, y = FOSL1)) +
  geom_point(size = 2) +
  geom_smooth(method = 'lm', se = FALSE, color = '#D62728') +
  annotate('text', x = 20, y = 19.5, hjust = 0, size = 3.5,
           label = paste0('Pearson r = ', round(res_1r$estimate, 2), '\np = ', round(res_1r$p.value, 2))) +
  scale_x_continuous(limits = c(0, 80)) +
  scale_y_continuous(limits = c(0, 20)) +
  labs(title = 'Late progressor samples', x = 'PCNA+ (%)', y = 'FOSL1+ (%)') +
  theme_classic(base_size = 13) +
  theme(plot.title = element_text(face = "bold", hjust = 0.5))

pdf('Figure1/Figure1q.pdf', width = 3.5, height = 3.5)
print(pt1q)
dev.off()

pdf('Figure1/Figure1r.pdf', width = 3.5, height = 3.5)
print(pt1r)
dev.off()
