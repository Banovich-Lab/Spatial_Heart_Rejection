################################################################################
# Angela Oill
# Spatial Heart Rejection Analysis code
# Analysis: Expression clustered heatmap
################################################################################
# Correlation Heatmap
# Pseudo-bulk per patient across all cell types. Then do a correlation heat map and add 
# annotations on the top

# Load libraries ----
library(ComplexHeatmap)
library(circlize)  


# Read in object ---- 
clustered_obj <- readRDS("/scratch/aoill/projects/heart_transplant/GEO/heart_spatial_obj.rds")

## Subset to pre-treatment patients only ----
clustered_obj_pre <- subset(clustered_obj, subset = biopsy_timing == "pre")

# Pseudo-bulk
# Extract the normalized expression data (data is log-transformed counts)
expr_matrix <- as.data.frame(as.matrix(clustered_obj_pre@assays$RNA@data))
expr_matrix_t <- t(expr_matrix)
expr_matrix_t <- as.data.frame(expr_matrix_t)

# Add sample information to the data
expr_matrix_t$patient_id <- clustered_obj_pre@meta.data$patient_id

# Summarize by Sample to calculate the mean expression for each gene and sample
mean_expression_per_sample <- expr_matrix_t %>%
  group_by(patient_id) %>%
  summarise(across(everything(), mean, na.rm = TRUE))

# Convert back to a data frame 
mean_expression_per_sample <- as.data.frame(mean_expression_per_sample)
rownames(mean_expression_per_sample) <- mean_expression_per_sample$patient_id
mean_expression_per_sample$patient_id <- NULL
#head(t(mean_expression_per_sample))
avg_expression_df <- t(mean_expression_per_sample)

# Add in sample info for plotting
# make a col called age_cat for pediatric or adult
clustered_obj_pre@meta.data$age_cat <- ifelse(clustered_obj_pre@meta.data$transplant_age < 18, "Pediatric", "Adult")

sample_data <- clustered_obj_pre@meta.data %>% 
  dplyr::select(patient_id, patient_rejection_type, patient_cellular_grading, patient_antibody_grading,
                age_cat) %>%
  unique()
#outcome_data$Sample <- gsub("_", "-", outcome_data$Sample)
head(sample_data)
rownames(sample_data) <- NULL

# add to expression matrix
mean_expression_per_sample$patient_id <- rownames(mean_expression_per_sample)
head(mean_expression_per_sample)
mean_expression_per_sample_info <- left_join(mean_expression_per_sample, sample_data,
                                             by = "patient_id")

# add rownames back in
rownames(mean_expression_per_sample_info) <- mean_expression_per_sample_info$patient_id
head(mean_expression_per_sample_info)

# run correlation
cor_matrix <- cor(avg_expression_df, method = "pearson")

# Set sample IDs as row names if not already
sample_data_tmp <- sample_data %>%
  column_to_rownames(var = "patient_id")  # Ensure patient_id colnames in cor_matrix

# Make sure the order of cor_matrix is the same as annotation_df
annotation_df <- annotation_df[colnames(cor_matrix), , drop = FALSE]

colnames(annotation_df) <- c("Rejection Type", "ACR Grade", "AMR Grade", "Age Category")

# Define colors for annotations
anno_colors <- list(
  `Age Category` = c("Pediatric" = "blue", "Adult" = "red"),
  `ACR Grade` = c("0" = "#28b463", "1" = "#ca6f1e", "2" = "#cc0000", 
                  "3" = "#0b5394", "NA" = "grey"),
  #`AMR Grade` = c("pAMR1-h" = "#DAF7A6", "pAMR1-i" = "#FFC300", 
  #                    "pAMR0" = "#FF5733", "pAMR2" = "#C70039", "NA" = "grey"),
  `AMR Grade` = c("pAMR0" = "#28b463", "pAMR1-h" = "#ca6f1e", 
                  "pAMR1-i" = "#cc0000", "pAMR2" = "#0b5394", "NA" = "grey"),
  
  `Rejection Type` = c("cellular_rejection" = "green", 
                       "mixed_rejection" = "purple", 
                       "antibody_rejection" = "orange")
)

# Create a HeatmapAnnotation object
col_anno <- HeatmapAnnotation(df = annotation_df, 
                              col = anno_colors)

ht <- Heatmap(cor_matrix, 
              name = "Correlation",
              top_annotation = col_anno,  # Add sample metadata annotations
              #col = colorRamp2(c(-1, 0, 1), c("blue", "white", "red")),  # Correlation color scale
              col = colorRamp2(c(min(cor_matrix), max(cor_matrix)), c("white", "blue")),  # Correlation color scale
              clustering_distance_rows = "euclidean",
              clustering_distance_columns = "euclidean",
              clustering_method_rows = "complete",
              clustering_method_columns = "complete",
              show_column_names = F,
              show_row_names = F)


pdf("/home/aoill/plots/00_heart_manuscript/pre_treatment_expression_correlation_heatmap_small_20250225.pdf",
    width = 8.25, height = 6)
draw(ht)
dev.off()