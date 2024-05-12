setwd("/Users/hazhang/Desktop/R_heatmap")

# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")

# BiocManager::install("ComplexHeatmap")

# library(ghcnv) # issue
library(ggplot2)

# install.packages("reshape")
library(reshape)

# install.packages("quantmod")
library(quantmod)

# install.packages("plotly")
library(plotly)

# install.packages("cluster")
library(cluster)

library(ComplexHeatmap)


# BiocManager::install("GenomicRanges")
library(GenomicRanges)

# library(RColorBrewer)
# library(MASS)
# library(parallel)

# install.packages("viridis")
# install.packages("viridisLite")
# library(viridis)

library(ComplexHeatmap)
library(dplyr)
library(tidyr)
library(tibble)

setwd("/Users/hazhang/Desktop/R_heatmap/individuals")

#Process patient 0058-00877
df_0171_01551 <- read.csv("df_0058_00877.csv")

heatmap_data_0171_01551 <- df_0171_01551 %>%
  mutate(Visit_name = factor(Visit_name, levels = c('CCRTC1D1', 'CCRTC2D1', 'CCRTC3D1', 'MAINTC1D1', 'MAINTC2D1', 'MAINTC3D1', 'MAINTC8D1', 'UPPROG'))) %>%
  arrange(Visit_name, sample_id) %>%
  select(sample_id, region_annotation, norm_region_count) %>%
  pivot_wider(names_from = sample_id, values_from = norm_region_count) %>%
  as_tibble() %>%
  column_to_rownames("region_annotation")

heatmap_matrix_0171_01551 <- as.matrix(heatmap_data_0171_01551)

stratification_data_0171_01551 <- df_0171_01551 %>%
  select(sample_id, MB_Lung_v4_call, Visit_name, epiTF_TFv2beta) %>%
  distinct() %>%
  mutate(
    MB_Lung_v4_call = as.factor(as.numeric(MB_Lung_v4_call == "Detected")),
    Visit_name = factor(Visit_name, levels = c('CCRTC1D1', 'CCRTC2D1', 'CCRTC3D1', 'MAINTC1D1', 'MAINTC2D1', 'MAINTC3D1', 'MAINTC8D1', 'UPPROG')),
    Visit_category = case_when(
      Visit_name %in% c('CCRTC1D1', 'CCRTC2D1') ~ "Baseline",
      Visit_name %in% c('CCRTC3D1', 'MAINTC1D1', 'MAINTC2D1', 'MAINTC3D1', 'MAINTC8D1') ~ "Treatment",
      Visit_name == 'UPPROG' ~ "Progression"
    ),
    epiTF_category = cut(epiTF_TFv2beta,
                         breaks = c(-Inf, 0.0001, 0.0005, 0.001, 0.005, 0.01, 0.05, 0.1, 0.5, Inf),
                         labels = c("0-0.0001", "0.0001-0.0005", "0.0005-0.001", "0.001-0.005", "0.005-0.01", "0.01-0.05", "0.05-0.1", "0.1-0.5", ">0.5"),
                         include.lowest = TRUE)
  ) %>%
  arrange(match(Visit_name, levels(Visit_name)))

ordered_sample_ids <- stratification_data_0171_01551$sample_id
heatmap_matrix_0171_01551 <- as.matrix(heatmap_data_0171_01551[, ordered_sample_ids])

annotation_colors <- list(
  MB_Lung_v4_call = c("0" = "green", "1" = "red"),
  Visit_category = c("Baseline" = "blue", "Treatment" = "yellow", "Progression" = "purple"),
  epiTF_category = setNames(colorRampPalette(c("blue", "red"))(9), 
                            c("0-0.0001", "0.0001-0.0005", "0.0005-0.001", "0.001-0.005", "0.005-0.01", "0.01-0.05", "0.05-0.1", "0.1-0.5", ">0.5"))
  # Not setting epiTF_TFv2beta here, using default
)

annotation <- HeatmapAnnotation(
  MB_Lung_v4_call = stratification_data_0171_01551$MB_Lung_v4_call,
  Visit_category = stratification_data_0171_01551$Visit_category,
  epiTF_category = stratification_data_0171_01551$epiTF_category,
  epiTF_TFv2beta = stratification_data_0171_01551$epiTF_TFv2beta
  # col parameter is not needed for epiTF_TFv2beta
)

# Generate the heatmap with disabled column clustering
pdf("HEATMAP_0058_00877_new.pdf", width = 12, height = 10)
Heatmap(heatmap_matrix_0171_01551, 
        name = "norm_region_count",
        top_annotation = annotation,
        show_row_names = TRUE,
        show_column_names = TRUE,
        clustering_distance_rows = "euclidean",
        clustering_distance_columns = NA,  # Disable column clustering
        clustering_method_rows = "complete",
        column_order = order(factor(colnames(heatmap_matrix_0171_01551), levels = ordered_sample_ids)),  # Force column order
        row_title = "Region Annotation",
        column_title = "Sample Visit",
        na_col = "white",
        row_names_side = "left",
        row_names_gp = gpar(fontsize = 3),
)
dev.off()

cat("The heatmap has been saved to 'HEATMAP_0058_00877_new.pdf'")
