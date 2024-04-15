# setwd("/ghds/pharma/customers/Ambrx/Infinity/ARX_SOW01_Infinity/ARX_SOW01_Infinity_03/data/subtype")
# dataDir <- "/ghds/pharma/customers/Ambrx/Infinity/ARX_SOW01_Infinity/ARX_SOW01_Infinity_03/data/subtype/METHYLATION_DATA"
# catalin_dataDir <- "/ghds/bioinformatics/02_DEVELOPMENT/230420_CRC_MOLECULAR_SUBTYPES/results/230420_METHYLATION_DATA"
# saveDir <- "/ghds/pharma/customers/Ambrx/Infinity/ARX_SOW01_Infinity/ARX_SOW01_Infinity_03/data/subtype/CLUSTER_SAMPLES"
# if(!file.exists(saveDir)) dir.create(saveDir, recursive = T)
# catalin_saveDir <- "/ghds/bioinformatics/02_DEVELOPMENT/230420_CRC_MOLECULAR_SUBTYPES/results/230421_CLUSTER_CRC_SAMPLES"

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("ComplexHeatmap")

library(ghcnv)
library(ggplot2)
library(reshape)
library(quantmod)
library(plotly)
library(cluster)
library(ComplexHeatmap)
library(GenomicRanges)
library(RColorBrewer)
library(MASS)
library(parallel)
library(viridis)

############################
##  Load sample info
############################

load(file = file.path(dataDir, "/sample_sirius.RData"))
anno = read.table("/ghds/pharma/customers/Ambrx/Infinity/ARX_SOW01_Infinity/ARX_SOW01_Infinity_03/data/subtype/explore_2/prostate_sample_sirius_ann.csv",
                  sep=",", header=T, stringsAsFactors = F)

############################
##  Load panels annotation
############################

repo <- "/ghds/bioinformatics/02_DEVELOPMENT/230325_S3_APPLICATION_PDL1/results"
# epi_panel, epi_probes, genomic_panel, genomic_probes
load(file = file.path(repo, "230403_PDL1_PREDICTION/panel_annotation.RData"))

###################################################
##  CpG strata pos ctrl regions normalized data
###################################################

# load(file = file.path(dataDir, "counts_scaled_hyper_epi_strata.RData"))
load(file = file.path(dataDir, "counts_scaled_hyper_epi_strata.230424.RData"))
counts_strata_norm <- counts_scaled_hyper_epi_strata[, match(anno$run_sample_id, colnames(counts_scaled_hyper_epi_strata))]

exp_num_clusters = 2

# counts_scaled_hyper_epi_strata

##############################
##  STRATA normalized data
##  clustering
##############################

##  Select regions with high variation
#i_samples <- which(anno$pred_frac > 0.1)
#sd_cohort <- apply(counts_strata_norm[,i_samples], 1, mad, na.rm = T)
#m_cohort <- apply(counts_strata_norm[,i_samples], 1, median, na.rm = T)
#i_targets <- which(m_cohort > -1)
#i_plot <- i_targets[order(sd_cohort[i_targets], decreasing = T)][1:200]
#plot(density(sd_cohort))
#range(sd_cohort[i_plot])

cc <- apply(counts_strata_norm, 1, function(x) cor(x, anno$pred_frac, use = "complete"))
i_samples <- which(anno$pred_frac > 0.01)
sd_cohort <- apply(counts_strata_norm[,i_samples], 1, mad, na.rm = T)
m_cohort <- apply(counts_strata_norm[,i_samples], 1, median, na.rm = T)
plot(density(m_cohort)); plot(density(sd_cohort)); plot(density(cc, na.rm = T))
i_targets <- which(m_cohort > -1 & 
                     sd_cohort > 0.7 &
                     abs(cc) < 0.1)
i_plot_targets <- i_targets[order(sd_cohort[i_targets], decreasing = T)][1:pmin(200, length(i_targets))]
i_plot_samples <- which(anno$pred_frac > 0.01 )


pdf(file = file.path(saveDir, "HEATMAP_COUNTS_STRATA_NORM_20230607.pdf"), height = 6, width = 12)
##   HEATMAP 
col <- list(TP53_deletion = c("No" = "green", "loh_deletion" = "red", "homozygous_deletion" = "darkred"),
            RB1_deletion = c("No" = "green", "loh_deletion" = "red", "homozygous_deletion" = "darkred"),
            TP53_mutation = c("No" = "lightblue", "Yes" = "purple"),
            AR_mutation = c("No" = "lightblue", "Yes" = "purple"),
            TF = plasma(11))
names(col$TF) <- as.character(seq(0, 1, by = 0.1))
ha <- HeatmapAnnotation(TP53_deletion = anno$TP53_deletion[i_plot_samples], 
                        RB1_deletion = anno$RB1_deletion[i_plot_samples], 
                        TP53_mutation = anno$TP53_mutation[i_plot_samples],
                        AR_mutation = anno$AR_mutation[i_plot_samples],
                        TF = round(pmin(anno$pred_frac[i_plot_samples], 1), 1),
                        source = anno$source[i_plot_samples],
                        col = col)

data_plot <- counts_strata_norm[i_plot_targets, i_plot_samples]
heatmap.4 <- ComplexHeatmap::Heatmap(data_plot, 
                                     row_title = "targeted regions",
                                     column_title = "Prostate Cancer samples", column_title_gp = gpar(fontsize = 18),
                                     name = "scores", column_km = exp_num_clusters,
                                     top_annotation = ha, 
                                     row_names_gp = gpar(fontsize = 6), column_names_gp = gpar(fontsize = 6))
print(heatmap.4)


##  TF vs cluster boxplot
clusterlist <- column_order(heatmap.4)
clu_df <- lapply(names(clusterlist), function(i){
  out <- data.frame(run_sample_id = colnames(data_plot[,clusterlist[[i]]]),
                    Cluster = paste0("cluster", i),
                    stringsAsFactors = FALSE)
  return(out)
}) %>%  
  do.call(rbind, .)

clu_df <- clu_df[match(anno$run_sample_id[i_plot_samples], clu_df$run_sample_id),]
dd <- data.frame(run_sample_id = anno$run_sample_id[i_plot_samples],
                 cluster = clu_df$Cluster,
                 TF = pmin(anno$pred_frac[i_plot_samples], 1),
                 TP53_deletion = anno$TP53_deletion[i_plot_samples],
                 RB1_deletion = anno$RB1_deletion[i_plot_samples],
                 TP53_mutation = anno$TP53_mutation[i_plot_samples],
                 AR_mutation = anno$AR_mutation[i_plot_samples])
g <- ggplot(dd, aes(x = cluster, y = TF)) +
  geom_boxplot() +
  ggtitle("Clusters TF") +
  theme(plot.title = element_text(size = 20, face = "bold"), # family = "Tahoma",
        text = element_text(size = 14), # , family = "Tahoma"
        axis.title = element_text(face="bold"),
        axis.text.x=element_text(size = 18, angle = 0),
        axis.text.y=element_text(size = 18))
print(g)

dev.off()

write.table(dd, file = file.path(saveDir, "METHYLATION_STRATA_NORM_CLUSTERS_20230607.tsv"),
            sep = "\t", col.names = T, row.names = F, quote = F)



#########################################################
# Beltran 2016 regions
#########################################################
sel_epi_region = read.table("/ghds/pharma/customers/Ambrx/Infinity/ARX_SOW01_Infinity/ARX_SOW01_Infinity_03/data/subtype/explore_1/sel_epi_regions_beltran_2016.csv",
                            sep=",", header=T, stringsAsFactors = F)

sel_epi_region_hypo = sel_epi_region[sel_epi_region["score"]=="hypo",]
sel_epi_region_hyper = sel_epi_region[sel_epi_region["score"]=="hyper",]

##################### hypo
cc <- apply(counts_strata_norm, 1, function(x) cor(x, anno$pred_frac, use = "complete"))
i_samples <- which(anno$pred_frac > 0.01)
sd_cohort <- apply(counts_strata_norm[,i_samples], 1, mad, na.rm = T)
m_cohort <- apply(counts_strata_norm[,i_samples], 1, median, na.rm = T)
plot(density(m_cohort)); plot(density(sd_cohort)); plot(density(cc, na.rm = T))
i_targets <- match(sel_epi_region_hypo$region_id, names(m_cohort))

# i_plot_targets <- i_targets[order(sd_cohort[i_targets], decreasing = T)][1:pmin(200, length(i_targets))]
i_plot_targets <- i_targets[order(sd_cohort[i_targets], decreasing = T)][1:length(i_targets)]
i_plot_samples <- which(anno$pred_frac > 0.01 )


pdf(file = file.path(saveDir, "HEATMAP_COUNTS_STRATA_NORM_beltran_hypo_20230607.pdf"), height = 6, width = 12)
##   HEATMAP 
col <- list(TP53_deletion = c("No" = "green", "loh_deletion" = "red", "homozygous_deletion" = "darkred"),
            RB1_deletion = c("No" = "green", "loh_deletion" = "red", "homozygous_deletion" = "darkred"),
            TP53_mutation = c("No" = "lightblue", "Yes" = "purple"),
            AR_mutation = c("No" = "lightblue", "Yes" = "purple"),
            TF = plasma(11))
names(col$TF) <- as.character(seq(0, 1, by = 0.1))
ha <- HeatmapAnnotation(TP53_deletion = anno$TP53_deletion[i_plot_samples], 
                        RB1_deletion = anno$RB1_deletion[i_plot_samples], 
                        TP53_mutation = anno$TP53_mutation[i_plot_samples],
                        AR_mutation = anno$AR_mutation[i_plot_samples],
                        TF = round(pmin(anno$pred_frac[i_plot_samples], 1), 1),
                        source = anno$source[i_plot_samples],
                        col = col)

data_plot <- counts_strata_norm[i_plot_targets, i_plot_samples]
heatmap.4 <- ComplexHeatmap::Heatmap(data_plot, 
                                     row_title = "targeted regions",
                                     column_title = "Prostate Cancer samples", column_title_gp = gpar(fontsize = 18),
                                     name = "scores", column_km = exp_num_clusters,
                                     top_annotation = ha, 
                                     row_names_gp = gpar(fontsize = 6), column_names_gp = gpar(fontsize = 6))
print(heatmap.4)


##  TF vs cluster boxplot
clusterlist <- column_order(heatmap.4)
clu_df <- lapply(names(clusterlist), function(i){
  out <- data.frame(run_sample_id = colnames(data_plot[,clusterlist[[i]]]),
                    Cluster = paste0("cluster", i),
                    stringsAsFactors = FALSE)
  return(out)
}) %>%  
  do.call(rbind, .)

clu_df <- clu_df[match(anno$run_sample_id[i_plot_samples], clu_df$run_sample_id),]
dd <- data.frame(run_sample_id = anno$run_sample_id[i_plot_samples],
                 cluster = clu_df$Cluster,
                 TF = pmin(anno$pred_frac[i_plot_samples], 1),
                 TP53_deletion = anno$TP53_deletion[i_plot_samples],
                 RB1_deletion = anno$RB1_deletion[i_plot_samples],
                 TP53_mutation = anno$TP53_mutation[i_plot_samples],
                 AR_mutation = anno$AR_mutation[i_plot_samples])
g <- ggplot(dd, aes(x = cluster, y = TF)) +
  geom_boxplot() +
  ggtitle("Clusters TF") +
  theme(plot.title = element_text(size = 20, face = "bold"), # family = "Tahoma",
        text = element_text(size = 14), # , family = "Tahoma"
        axis.title = element_text(face="bold"),
        axis.text.x=element_text(size = 18, angle = 0),
        axis.text.y=element_text(size = 18))
print(g)
dev.off()

write.table(dd, file = file.path(saveDir, "METHYLATION_STRATA_NORM_CLUSTERS_beltran_hypo_20230607.tsv"),
            sep = "\t", col.names = T, row.names = F, quote = F)


##################### hyper
cc <- apply(counts_strata_norm, 1, function(x) cor(x, anno$pred_frac, use = "complete"))
i_samples <- which(anno$pred_frac > 0.01)
sd_cohort <- apply(counts_strata_norm[,i_samples], 1, mad, na.rm = T)
m_cohort <- apply(counts_strata_norm[,i_samples], 1, median, na.rm = T)
plot(density(m_cohort)); plot(density(sd_cohort)); plot(density(cc, na.rm = T))
i_targets <- match(sel_epi_region_hyper$region_id, names(m_cohort))

# i_plot_targets <- i_targets[order(sd_cohort[i_targets], decreasing = T)][1:pmin(200, length(i_targets))]
i_plot_targets <- i_targets[order(sd_cohort[i_targets], decreasing = T)][1:length(i_targets)]
i_plot_samples <- which(anno$pred_frac > 0.01 )


pdf(file = file.path(saveDir, "HEATMAP_COUNTS_STRATA_NORM_beltran_hyper_20230607.pdf"), height = 6, width = 12)
##   HEATMAP 
col <- list(TP53_deletion = c("No" = "green", "loh_deletion" = "red", "homozygous_deletion" = "darkred"),
            RB1_deletion = c("No" = "green", "loh_deletion" = "red", "homozygous_deletion" = "darkred"),
            TP53_mutation = c("No" = "lightblue", "Yes" = "purple"),
            AR_mutation = c("No" = "lightblue", "Yes" = "purple"),
            TF = plasma(11))
names(col$TF) <- as.character(seq(0, 1, by = 0.1))
ha <- HeatmapAnnotation(TP53_deletion = anno$TP53_deletion[i_plot_samples], 
                        RB1_deletion = anno$RB1_deletion[i_plot_samples], 
                        TP53_mutation = anno$TP53_mutation[i_plot_samples],
                        AR_mutation = anno$AR_mutation[i_plot_samples],
                        TF = round(pmin(anno$pred_frac[i_plot_samples], 1), 1),
                        source = anno$source[i_plot_samples],
                        col = col)

data_plot <- counts_strata_norm[i_plot_targets, i_plot_samples]
heatmap.4 <- ComplexHeatmap::Heatmap(data_plot, 
                                     row_title = "targeted regions",
                                     column_title = "Prostate Cancer samples", column_title_gp = gpar(fontsize = 18),
                                     name = "scores", column_km = exp_num_clusters,
                                     top_annotation = ha, 
                                     row_names_gp = gpar(fontsize = 6), column_names_gp = gpar(fontsize = 6))
print(heatmap.4)


##  TF vs cluster boxplot
clusterlist <- column_order(heatmap.4)
clu_df <- lapply(names(clusterlist), function(i){
  out <- data.frame(run_sample_id = colnames(data_plot[,clusterlist[[i]]]),
                    Cluster = paste0("cluster", i),
                    stringsAsFactors = FALSE)
  return(out)
}) %>%  
  do.call(rbind, .)

clu_df <- clu_df[match(anno$run_sample_id[i_plot_samples], clu_df$run_sample_id),]
dd <- data.frame(run_sample_id = anno$run_sample_id[i_plot_samples],
                 cluster = clu_df$Cluster,
                 TF = pmin(anno$pred_frac[i_plot_samples], 1),
                 TP53_deletion = anno$TP53_deletion[i_plot_samples],
                 RB1_deletion = anno$RB1_deletion[i_plot_samples],
                 TP53_mutation = anno$TP53_mutation[i_plot_samples],
                 AR_mutation = anno$AR_mutation[i_plot_samples])
g <- ggplot(dd, aes(x = cluster, y = TF)) +
  geom_boxplot() +
  ggtitle("Clusters TF") +
  theme(plot.title = element_text(size = 20, face = "bold"), # family = "Tahoma",
        text = element_text(size = 14), # , family = "Tahoma"
        axis.title = element_text(face="bold"),
        axis.text.x=element_text(size = 18, angle = 0),
        axis.text.y=element_text(size = 18))
print(g)
dev.off()

write.table(dd, file = file.path(saveDir, "METHYLATION_STRATA_NORM_CLUSTERS_beltran_hyper_20230607.tsv"),
            sep = "\t", col.names = T, row.names = F, quote = F)






library(ComplexHeatmap)
library(dplyr)
library(tidyr)
library(tibble)

# Set the working directory to the location of the CSV file
setwd("/ghsfa/projects/pharma/projects/sirius_pharma/hazhang_projects/On_treatment_effect_on_methylation_03172024")

# Load the dataset
df_consolidated <- read.csv("df_consolidated_short_small.csv")

# Use dplyr's select to avoid namespace issues

# Remove the "dplyr::" if encounter error
heatmap_data <- df_consolidated %>%
  dplyr::select(sample_id, region_annotation, norm_peak_count) %>%
  pivot_wider(names_from = sample_id, values_from = norm_peak_count) %>%
  as_tibble() %>%
  column_to_rownames("region_annotation")

# Convert the data to a matrix for the heatmap
heatmap_matrix <- as.matrix(heatmap_data)

# Prepare the stratification data
stratification_data <- df_consolidated %>%
  dplyr::select(sample_id, MB_Lung_v4_call, Visit_name, epiTF_TFv2beta) %>%
  distinct() %>%
  mutate(
    MB_Lung_v4_call = as.numeric(MB_Lung_v4_call == "Detected"),
    Visit_category = case_when(
      Visit_name %in% c('CCRTC1D1', 'CCRTC2D1') ~ "Baseline",
      Visit_name %in% c('CCRTC3D1', 'MAINTC1D1', 'MAINTC2D1', 'MAINTC3D1', 'MAINTC8D1') ~ "Treatment",
      Visit_name == 'UPPROG' ~ "Progression",
      TRUE ~ as.character(NA)  # handle any other unexpected cases
    ),
    epiTF_category = cut(epiTF_TFv2beta,
                         breaks = c(-Inf, 0.25, 0.5, 0.75, Inf),
                         labels = c("Low", "Medium", "High", "Very High"),
                         include.lowest = TRUE)
  ) %>%
  arrange(sample_id)

# Ensure the order of stratification matches the columns in the heatmap_matrix
stratification_vector <- stratification_data$MB_Lung_v4_call
names(stratification_vector) <- stratification_data$sample_id

# Create the heatmap annotation object with the correct labels and colors
annotation <- HeatmapAnnotation(
  MB_Lung_v4_call = stratification_vector,
  Visit_category = stratification_data$Visit_category[match(colnames(heatmap_matrix), stratification_data$sample_id)],
  epiTF_category = stratification_data$epiTF_category[match(colnames(heatmap_matrix), stratification_data$sample_id)],
  col = list(
    MB_Lung_v4_call = c("0" = "green", "1" = "red"),
    Visit_category = c("Baseline" = "blue", "Treatment" = "yellow", "Progression" = "purple"),
    epiTF_category = c("Low" = "gray", "Medium" = "orange", "High" = "pink", "Very High" = "red")
  )
)

# Generate the heatmap
pdf("HEATMAP_CONSOLIDATED_NORM.pdf", width = 12, height = 10)
Heatmap(heatmap_matrix, 
        name = "norm_peak_count",
        top_annotation = annotation,
        show_row_names = TRUE,
        show_column_names = TRUE,
        clustering_distance_rows = "euclidean",
        clustering_distance_columns = "euclidean",
        clustering_method_rows = "complete",
        clustering_method_columns = "complete",
        row_title = "Region Annotation",
        column_title = "Sample ID",
        na_col = "white",
        row_names_gp = gpar(fontsize = 3)  # Smaller font size for row names
)
dev.off()

cat("The heatmap has been saved to 'HEATMAP_CONSOLIDATED_NORM.pdf'")