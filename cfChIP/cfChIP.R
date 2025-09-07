---
  title: "cfChIP manuscript analyses"
date: "2024-06-04"
output:
  html_document: default
pdf_document: default
---
  
  This document contains the analyses for the cfChIP-seq manuscript. To run it, install the required libraries and change the root dir below as needed
```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
library(rtracklayer)
library(GenomicRanges)
library(ggplot2)
library(dplyr)
library(tidyr)
library(stringr)
library(data.table)
library(readxl)
library(GGally)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(org.Hs.eg.db)
library(ggpubr)
library(ROCR)
options(bitmapType='cairo')
png("xzvf.png")
source('./common_functions_3groups.R')

```

#Create an oncoprint-style plot of the plasma cohort
``` {r cohort_plot}
if(TRUE) {
set.seed(1)

library(tidyverse)
library(ComplexHeatmap)

metadata_file = "tRCC_oncoprint.xlsx"
meta = read_xlsx(metadata_file, na = 'NA')
meta = meta  %>% 
  dplyr::select('study_name', 'antibody', 'cancer_type', 'sample_id', 'individual_id', 'sex', 'fusion_partner') 

tRCC_pts = meta[meta$cancer_type == "tRCC",]
ccRCC_pts = meta[meta$cancer_type == "ccRCC",]
healthy = meta[meta$cancer_type == "HP",]

ann = as.data.frame(meta)
rownames(ann) = ann$study_name

pdata = xtabs(~antibody + sample_id, meta) %>% as.matrix()
pdata = ifelse(pdata > 0, "yes", "")

pdata = pdata[,c(24:13,34:51,1:12,25:33)]

tmp = xtabs(~antibody, ann) # number of samples for each antibody
rownames(pdata) = paste0(names(tmp), ' (N=', tmp, ')')
  
heatmap_legend_param1 = list(title = "Oncoprint", at = c("yes"), 
                            labels = c("yes"))

col = c("yes" = "black")
alter_fun = list(
  background = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0, "pt"), h-unit(2, "pt"), 
              gp = gpar(fill = "white", col = "grey70"))
  },
  yes = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0, "pt"), h-unit(2, "pt"), 
              gp = gpar(fill = col["yes"], col = "grey70"))
  }
)

type_colors <- c("tRCC" = "#DC000099",
                 "ccRCC" = "#B09C8599",
                 "HP" = "#8491B488")
sex_colors <- c("Male" = "#1f78b4", "Female" = "pink")
fusion_colors <- c(
  "Absent" = "#FFFFFF",       # White
  "ASPL-TFE3" = "#8DCCF6",    # Blue cyan
  "PRCC-TFE3" = "#00CC99",    # Green cyan
  "NONO-TFE3" = "#FF9933",    # Dull orange
  "RBM10-TFE3" = "#FFA07A",   # Light Salmon
  "ASPSCR1-TFE3" = "#FF6347", # Tomato
  "MED15-TFE3" = "#9370DB",
  "Unkown" = "grey80" # Medium Purple
)

  p = oncoPrint(pdata,
                  show_column_names = FALSE, 
                  show_pct = FALSE,
                  heatmap_width = unit(14, "cm"), 
                  heatmap_height = unit(8, "cm"),
                  row_order = rownames(pdata)[c(1,2,4,3)],
                  column_order = c(unique(tRCC_pts$sample_id),
                                   unique(ccRCC_pts$sample_id),
                                   unique(healthy$sample_id)), 
                column_split = factor(c(rep("tRCC", length(unique(tRCC_pts$sample_id))), 
                                    
                                            rep("ccRCC", length(unique(ccRCC_pts$sample_id))), 
                                        rep("HP", length(unique(healthy$sample_id)))),
                                      levels = c("tRCC", "ccRCC", "HP")),
                  alter_fun = alter_fun, col = col, 
                  remove_empty_columns = TRUE, remove_empty_rows = FALSE,
                  show_heatmap_legend = FALSE,
                  top_annotation = HeatmapAnnotation("Fusion" = meta$fusion_partner[match(colnames(pdata), meta$sample_id)],
                                                     Type = case_when(
                                                       colnames(pdata) %in% tRCC_pts$sample_id ~ 'tRCC',
                                                       colnames(pdata) %in% ccRCC_pts$sample_id ~ 'ccRCC',
                                                       colnames(pdata) %in% healthy$sample_id ~ 'HP'),
                                                     Sex = meta$sex[match(colnames(pdata), meta$sample_id)],  
                                                     col = list(Type = type_colors,
                                                                Sex = sex_colors,
                                                                "Fusion" = fusion_colors),
                                                     show_legend = TRUE, 
                                                     simple_anno_size = unit(0.35, "cm"), # Keep annotation unified
                                                     annotation_name_side = "right"
                  ),
                  column_title = "cell-free epigenomic datasets", heatmap_legend_param = heatmap_legend_param1)

  print(p)
  
}

```

#MeDIP cell-free analysis
```{r TF_binding_site_signal_analysis}

# load metadata:
meta = read.csv2("data_plasmas_RCC.csv",header=TRUE, sep = ';')
meta <- subset(meta, antibody == 'H3K4me3/H3K27ac/MeDIP') #select one
meta_sub = meta

# prepare sites
# TODO: trim these - not using all of them
sites_file_list = list(   
  
  "/Users/simongarinet/Documents/Travail/Figures_ctDNA/scripts/data/sites/top_common_DHS.bed",
  "path/to/tRCCup.bed"
)

exclude_regions_file = "/Users/simongarinet/Documents/Travail/Figures_ctDNA/scripts/data/sites/hg19-blacklist.v2.bed"
exclude_regions = import(exclude_regions_file, format = 'bed')

# prepare sites
message('tiling sites')
sites_list = tile_sites(sites_file_list, exclude_regions = exclude_regions)

sites_names = str_replace(sites_file_list, '.*/', '') %>%
  str_replace('.bed', '') 
#TODO: should incorporate the below calls to signal_at_sites into a loop to avoid rewriting code

dir.create('/Users/simongarinet/Documents/output/out/')

# --- Breast cancer ER analyses:
if(TRUE) {
  
  keep_sites = grepl('MeDIP|H3K27Ac|H3K4me3|DHS|hg19', sites_names)
  toplot = NULL
  count = 1
  tot = nrow(meta_sub)
  for(rds in meta_sub$rds_file) {
    message('processing ', rds, ' --- ', count, ' of ', tot)
    counts = signal_at_sites(
      frags_file = rds,
      sites_list = sites_list[keep_sites],
      sites_names = sites_names[keep_sites],
      remove_peaks_wider_than = 4000
    )
    counts$cancer_type = meta_sub$cancer_type[meta_sub$rds_file == rds]
    counts$cancer_subtype = meta_sub$cancer_subtype[meta_sub$rds_file == rds]
    counts$study_name = meta_sub$study_name[meta_sub$rds_file == rds]
    toplot = rbind(toplot, counts)
    count = count + 1
  }
  
  saveRDS(toplot, paste0('.output/profile.RDS'))
  
}


#TODO: remove any of these that aren't used
plot_params_list = list(
  list('type', 'Antibody',"site", 'top_common_DHS')
)

dir.create('/Volumes/sviswanathan/users/sgarinet/Liquid_Biopsies_RCC/ctDNA/Results/MeDIP_Baca/out/H3K4me3/scores')
for(plot_params in plot_params_list){
  cancer_type = plot_params[[1]]
  mark = plot_params[[2]]
  sites_to_plot = plot_params[[3]]
  normalize_to = plot_params[[4]]
  
    toplot = readRDS('.output/profile.RDS')

  toplot = merge(toplot[, !colnames(toplot) %in% c('cancer_type', 'cancer_subtype')], # don't include cancer_type and cancer_subtype because these may have been modified in toplot for the comparisons above.,
                 meta_sub,
                 by = 'study_name')
  
  # since comparisons below are made based on subtype, use the cancer type in place of subtype if subtype is NA
  #toplot$cancer_subtype = ifelse(is.na(toplot$cancer_subtype), toplot$cancer_type, toplot$cancer_subtype)
  for (boxplot_arg in c(FALSE, TRUE)) {
    p = plot_signal_at_sites(
      subset(toplot,
             antibody %in% mark & sites %in% c(sites_to_plot, normalize_to)),
      normalize_to_these_sites = normalize_to,
      subtract_shoulder = TRUE,
      group_by_this = 'cancer_subtype',
      main_group = 'tRCC',
      auc_boxplot = boxplot_arg,
      #out_file = paste0('/Volumes/sviswanathan/users/sgarinet/Liquid_Biopsies_RCC/ctDNA/Results/MeDIP_Baca/out/H3K4me3/scores', cancer_type, '_', paste0(mark, collapse = ''), '_', sites_to_plot[1], '_', normalize_to, '.tsv')
      out_file = paste0('/Users/simongarinet/Documents/', cancer_type, '_', paste0(mark, collapse = ''), '_', sites_to_plot[1], '_', normalize_to, '.tsv')
      
      
    )
    
    p = p + theme(axis.text.y = element_text(size = 14),
                  axis.text.x = element_text(size = 9))
    #pdf(file="test.pdf")
    # if (boxplot_arg == TRUE) {
    p = p  +
      stat_compare_means(comparisons = list(c("tRCC", "ccRCC"), c("tRCC", "HP"), c("ccRCC", "HP")), #not sure about the names here
                         method = "wilcox.test")
    #}
    suppressWarnings(print(p))
    
  }
  
}


```

#ROC curves for tRCC vs Healthy and ccRCC
```{r ROC curves}

to_keep = "ccRCC" #either keep HP or ccRCC (insert c("ccRCC", "HP) to keep both)

ctDNA_cut_off = 0 # to filter samples with ctDNA fraction < ctDNA_cut_off 

# Calculate the prediction of H3K27ac signal at TFE3 binding sites alone

scores_file = '/Users/simongarinet/Documents/RCC_H3K27Ac_GTRD_common3CL_hg19_lessDEpeaks_top_common_DHS.tsv'

scores_TFE3 = read.table(scores_file, sep = '\t', header = TRUE) %>%
  subset(sites != "top_common_DHS") %>% #filter out sites used for normalization
  dplyr::filter(cancer_subtype %in% c("tRCC", to_keep)) %>% #keep tRCC samples with either healthy or ccRCC
  dplyr::filter(ctDNA_estimate >= ctDNA_cut_off | cancer_subtype == "HP")

pred = prediction(scores_TFE3$auc, scores_TFE3$cancer_subtype == 'tRCC')  
perf = performance(pred, 'tpr', 'fpr')
auc_TFE3 = unlist(performance(pred, 'auc')@y.values) %>% round(digits = 2)

perf_df = data.frame(FPR = unlist(perf@x.values), 
                     TPR = unlist(perf@y.values),
                     comparison = paste0('tRCC vs.', to_keep, ' using TFE3 at TFE3 bindings sites'))

# Calculate the prediction of H3K27ac signal at tRCC-up sites alone

scores_file = '/Users/simongarinet/Documents/RCC_H3K27Ac_H3K27Ac_tRCC_up_NEW_1FC_hg19_top_common_DHS.tsv'

scores_H3K27ac = read.table(scores_file, sep = '\t', header = TRUE) %>%
  subset(sites != "top_common_DHS") %>% #filter out sites used for normalization
  dplyr::filter(cancer_subtype %in% c("tRCC", to_keep)) %>% #keep tRCC samples with either healthy or ccRCC
  dplyr::filter(ctDNA_estimate >= ctDNA_cut_off | cancer_subtype == "HP")

pred = prediction(scores_H3K27ac$auc, scores_H3K27ac$cancer_subtype == 'tRCC')  
perf = performance(pred, 'tpr', 'fpr')
auc_H3K27ac = unlist(performance(pred, 'auc')@y.values) %>% round(digits = 2)

perf_df = rbind(perf_df, data.frame(FPR = unlist(perf@x.values),
                                    TPR = unlist(perf@y.values),
                                    comparison = paste0('tRCC vs.', to_keep, ' using H3K27ac tRCC-up')))

# Calculate the prediction of H3K4me3 signal at tRCC-up sites alone

scores_file = '/Users/simongarinet/Documents/RCC_H3K4me3_ccRCC_tRCC_CL_H3K4me3_tRCCup_2FC_hg19_top_common_DHS.tsv'

scores = read.table(scores_file, sep = '\t', header = TRUE)

scores = subset(scores, sites == "ccRCC_tRCC_CL_H3K4me3_tRCCup_2FC_hg19") %>%
  subset(sites != "top_common_DHS") %>% #filter out sites used for normalization
  dplyr::filter(cancer_subtype %in% c("tRCC", to_keep)) %>% #keep tRCC samples with either healthy or ccRCC
  dplyr::filter(ctDNA_estimate >= ctDNA_cut_off | cancer_subtype == "HP")

pred = prediction(scores$auc, scores$cancer_subtype == 'tRCC')  
perf = performance(pred, 'tpr', 'fpr')
auc_H3K4me3 = unlist(performance(pred, 'auc')@y.values) %>% round(digits = 2)

perf_df = rbind(perf_df, data.frame(FPR = unlist(perf@x.values), 
                                    TPR = unlist(perf@y.values),
                                    comparison = paste0('tRCC vs.', to_keep, ' using H3K4me3 tRCC-up')))

### Integrate the 3 scores together

scores_file = '/Users/simongarinet/Documents/RCC_H3K27Ac_GTRD_common3CL_hg19_lessDEpeaks_top_common_DHS.tsv'
scores_TFE3 = read.table(scores_file, sep = '\t', header = TRUE) %>%
  subset(sites != "top_common_DHS") %>% #filter out sites used for normalization
  dplyr::filter(cancer_subtype %in% c("tRCC", to_keep)) %>% #keep tRCC samples with either healthy or ccRCC
  dplyr::filter(ctDNA_estimate >= ctDNA_cut_off | cancer_subtype == "HP")

scores_file = '/Users/simongarinet/Documents/RCC_H3K4me3_ccRCC_tRCC_CL_H3K4me3_tRCCup_2FC_hg19_top_common_DHS.tsv'
scores_H3K4me3 = read.table(scores_file, sep = '\t', header = TRUE) %>%
  subset(sites != "top_common_DHS") %>% #filter out sites used for normalization
  dplyr::filter(cancer_subtype %in% c("tRCC", to_keep)) %>% #keep tRCC samples with either healthy or ccRCC
  dplyr::filter(ctDNA_estimate >= ctDNA_cut_off | cancer_subtype == "HP")

scores_file = '/Users/simongarinet/Documents/RCC_H3K27Ac_H3K27Ac_tRCC_up_NEW_1FC_hg19_top_common_DHS.tsv'
scores_H3K27ac = read.table(scores_file, sep = '\t', header = TRUE) %>%
  subset(sites != "top_common_DHS") %>% #filter out sites used for normalization
  dplyr::filter(cancer_subtype %in% c("tRCC", to_keep)) %>% #keep tRCC samples with either healthy or ccRCC
  dplyr::filter(ctDNA_estimate >= ctDNA_cut_off | cancer_subtype == "HP")

scores_final = scores_H3K27ac
scores_final$auc = scores_H3K27ac$auc + scores_TFE3$auc + scores_H3K4me3$auc
scores_final = scores_final %>% 
  dplyr::filter(ctDNA_estimate >= ctDNA_cut_off | cancer_subtype == "HP")

pred = prediction(scores_final$auc, scores_final$cancer_subtype == 'tRCC')  
perf = performance(pred, 'tpr', 'fpr')
auc_integration = unlist(performance(pred, 'auc')@y.values) %>% round(digits = 2)


perf_df = rbind(perf_df, data.frame(FPR = unlist(perf@x.values), 
                                    TPR = unlist(perf@y.values),
                                    comparison = paste0('tRCC vs. ', to_keep, ' integration')))

my_theme <- function() {
  theme_minimal() +
    theme(
      panel.grid.minor = element_blank(),
      panel.background = element_rect(fill = 'white'),
      axis.line = element_line(color = 'black'),
      legend.position = 'right',
      legend.key = element_rect(fill = 'white')
    )
}

# Plot
p <- ggplot(perf_df, aes(y = TPR, x = FPR, col = comparison)) +
  geom_path(size = ifelse(perf_df$comparison == paste0("tRCC vs. ", to_keep, " integration"), 1.5, 0.5), alpha = 0.8) +  # Increase line size and transparency
  ggtitle('tRCC vs. ccRCC - all samples') +
  geom_abline(slope = 1, intercept = 0, linetype = 'dashed', color = 'gray') +
  scale_color_manual(values = c('darkgreen', 'darkblue', 'darkred', 'darkorange'),
                     name = '',
                     labels = c(paste0('Integration (AUC = ', auc_integration, ')'),
                                paste0('H3K27ac tRCC up sites (AUC=', auc_H3K27ac, ')'),
                                paste0('H3K4me3 tRCC up sites (AUC=', auc_H3K4me3, ')'),
                                paste0('H3K27ac at TFE3 binding sites (AUC=', auc_TFE3, ')'))) +
  my_theme() +
  theme(legend.margin = margin(t = 10, unit = 'pt'),
        plot.margin = unit(c(1, 1, 1, 1), "lines"),
        legend.position = "right",
        panel.grid = element_blank())

print(p)

library(ggpubr)
# Plot with reordered x-axis and custom colors
boxplot <- ggplot(scores_final, aes(x = factor(cancer_subtype, levels = c("tRCC", "ccRCC", "HP")), 
                             y = log2(auc + 1), 
                             col = ctDNA_estimate)) +
  geom_boxplot(width = 0.7, outlier.shape = NA) +
  geom_jitter(width = 0.3, size = 1) +  # Apply coloring to jitter points based on 'Patient_color'
#  scale_color_manual(values = patient_color_mapping) +
#  scale_color_manual(values = subtype_colors, name = "Type", guide = "none") +
  scale_color_gradient(low = "darkblue", high = "yellow", name = "ctDNA fraction") + 
  scale_alpha(range = c(0.8, 1), guide = "none") +
  theme_classic() +
  guides(color = guide_legend(title = NULL)) +
  labs(title = "",
       x = "",
       y = "Integrated signal (arbitrary value)") +
  theme(plot.margin = unit(c(1, 5, 1, 5), "cm"))
 # geom_hline(yintercept = 5.365, linetype = "dashed", color = "#DC000099")

print(boxplot)


```

``` {r Leave-one-out-cross-validation tRCC vs HP}

library(ROCR)

prediction_df <- NULL
cut_offs_summary_final <- NULL

scores_final$cancer_subtype = ifelse(scores_final$cancer_subtype == "tRCC", "tRCC", "Other") #Other could be either HP or ccRCC, or all together

for (sample_to_test in scores_final$study_name) {
  
  sample_to_test_score = scores_final$auc[scores$study_name==sample_to_test]
  sample_to_test_subtype = scores_final$cancer_subtype[scores$study_name==sample_to_test]
  
  scores_final_subset = scores_final %>% 
    dplyr::filter(study_name != sample_to_test)

# Step 1: Generate the ROC curve
pred = prediction(scores_final_subset$auc, scores_final_subset$cancer_subtype == 'tRCC')  
perf = performance(pred, 'tpr', 'fpr')

# Step 2: Calculate AUC (Area Under the Curve)
auc_value = unlist(performance(pred, 'auc')@y.values) %>% round(digits = 2)

# Step 3: Find the optimal cutoff (maximizing Youden's index)
# Youden's index = sensitivity + specificity - 1, where sensitivity = TPR and specificity = 1 - FPR
fpr_values = unlist(perf@x.values)
tpr_values = unlist(perf@y.values)

# Calculate Youden's index for each threshold
youden_index = tpr_values - fpr_values

# Find the threshold corresponding to the maximum Youden's index
optimal_cutoff_index = which.max(youden_index)
optimal_cutoff = pred@cutoffs[[1]][optimal_cutoff_index]

cut_offs_summary = data.frame(cut_off = optimal_cutoff)

cut_offs_summary_final = rbind(cut_offs_summary_final, cut_offs_summary)

# Step 4: Apply the cutoff to classify the samples
scores_final_subset$classified = ifelse(scores_final_subset$auc >= optimal_cutoff, 'tRCC', 'Other')

row_to_add <- data.frame(sample_id = sample_to_test,
                         cancer_subtype = sample_to_test_subtype,
                         prediction = case_when(
                           sample_to_test_score >= optimal_cutoff ~ "tRCC",
                           sample_to_test_score < optimal_cutoff ~ "Other"
                           ))

prediction_df = rbind(prediction_df, row_to_add)

}

conf_mat <- table(prediction_df$cancer_subtype, prediction_df$prediction)

TP <- conf_mat["tRCC", "tRCC"]
FN <- conf_mat["tRCC", "Other"]
TN <- conf_mat["Other", "Other"]
FP <- conf_mat["Other", "tRCC"]

sensitivity <- TP / (TP + FN)
specificity <- TN / (TN + FP)

precision <- TP / (TP + FP)
recall <- TP / (TP + FN)  # same as sensitivity

cat("Precision: ", round(precision, 3), "\n")
cat("Recall: ", round(recall, 3), "\n")

cat("Sensitivity: ", round(sensitivity, 3), "\n")
cat("Specificity: ", round(specificity, 3), "\n")

LOOCV_cut_off = log2(mean(cut_offs_summary_final$cut_off)+1)
cat("Average cutoff: ", round(LOOCV_cut_off,2))
```
