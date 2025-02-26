# Load necessary libraries
library(ComplexHeatmap)
library(dplyr)
library(tidyr)
library(tibble)

# Set the working directory and read the data
setwd("/Users/hazhang/Desktop/R_heatmap")
df_consolidated <- read.csv("df_consolidated_inform_DMR_top30LR_01152025.csv")

# Sort the data frame based on sample_type
df_sorted <- df_consolidated %>%
  arrange(sample_type, request_id)

# Prepare the data for the heatmap
heatmap_data <- df_sorted %>%
  select(region_id, request_id, adjusted_norm_peak_count) %>%
  pivot_wider(names_from = request_id, values_from = adjusted_norm_peak_count) %>%
  replace(is.na(.), 0) %>%
  as_tibble()

# Convert to matrix, with region_id as row names
heatmap_matrix <- as.matrix(heatmap_data[-1])
rownames(heatmap_matrix) <- heatmap_data$region_id

# Create the top annotation based on sorted dataframe
annotations <- df_sorted %>%
  select(request_id, sample_type) %>%
  distinct() %>%
  `$`(sample_type) %>%
  factor(levels = c('ILD_Positive', "Prednisone_Positive", "IO_Control", "Osimertinib_Control"))

top_annotation <- HeatmapAnnotation(
  sample_type = annotations,
  col = list(sample_type = c("ILD_Positive" = "blue", "Prednisone_Positive" = "gray", "IO_Control" = "red", "Osimertinib_Control" = "green")),
  which = "column"
)


# Open PDF device
pdf("inform_LRtop30_dmr.pdf", width = 12, height = 10)

# Generate the heatmap
Heatmap(heatmap_matrix[, order(annotations)],  # Sort columns based on the order of annotations
        name = "adjusted_norm_peak_count",
        top_annotation = top_annotation,
        show_row_names = TRUE,
        show_column_names = TRUE,
        clustering_distance_rows = "euclidean",
        clustering_distance_columns = "euclidean",
        clustering_method_rows = "complete",
        clustering_method_columns = "complete",
        row_title = "Region ID",
        column_title = "Sample ID",
        na_col = "white",
        row_names_side = "left",
        row_names_gp = gpar(fontsize = 10),
        column_names_gp = gpar(fontsize = 10)
)

# Close PDF device
dev.off()



