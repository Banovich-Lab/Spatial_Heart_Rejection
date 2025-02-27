################################################################################
# Spatial Heart Rejection Analysis code
# Analysis: Expression PCA
################################################################################


# Load libraries ----
library(factoextra)
library(FactoMineR)
library(ggpubr)

# Read in object ---- 
clustered_obj <- readRDS("/scratch/aoill/projects/heart_transplant/GEO/heart_spatial_obj.rds")

# Pseudo-bulk ----
# subset to only pre treatment so it matches the cell type pca
clustered_obj_pre <- subset(clustered_obj, subset = biopsy_timing == "pre")

# Extract the normalized expression data (data is log-transformed counts)
expr_matrix <- as.data.frame(as.matrix(clustered_obj_pre@assays$RNA@data))
expr_matrix_t <- t(expr_matrix)
expr_matrix_t <- as.data.frame(expr_matrix_t)

#all_genes <- colnames(expr_matrix_t)

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

# mean_expression_per_sample will be used for PCA but need to add rejection type
# info for plotting
mean_expression_per_sample$patient_id <- rownames(mean_expression_per_sample)

tmp <- clustered_obj_pre@meta.data %>%
  dplyr::select(patient_id, patient_biopsy, biopsy_timing, #rejection_type, 
                sex, donor_sex, #cellular_grading, antibody_grading, 
                treatment,
                patient_cellular_grading, patient_antibody_grading, patient_rejection_type) %>%
  unique()
mean_expression_per_sample_info <- left_join(mean_expression_per_sample, tmp)
rownames(mean_expression_per_sample_info) <- rownames(mean_expression_per_sample)


smpl.pca <- PCA(mean_expression_per_sample_info[,1:477], graph = F)

saveRDS(smpl.pca, "/home/aoill/projects/heart_transplant/00_final/smpl.pca_CR_pre_pseudo_all_together.rds")


# PCA loadings

# By cell type ----
cts_to_analyze <- levels(as.factor(clustered_obj_pre@meta.data$ct_second_pass))
pca_plots <- list()  # Empty list to store plots


custom_colors <- c("antibody_rejection" = "#F8766D", 
                   "cellular_rejection" = "#00BFC4",  
                   "mixed_rejection" = "#7CAE00")
# Define custom shapes 
custom_shapes <- c("antibody_rejection" = 16,
                   "cellular_rejection" = 17,
                   "mixed_rejection" = 15)

for (ctn in 1:length(cts_to_analyze)) {
  ct <- cts_to_analyze[ctn]
  print(ct)
  
  # subset object to cell type
  obj_CT <- subset(clustered_obj_pre, subset = ct_second_pass == ct)
  
  # Extract the normalized expression data (data is log-transformed counts)
  expr_matrix <- as.data.frame(as.matrix(obj_CT@assays$RNA@data))
  expr_matrix_t <- t(expr_matrix)
  expr_matrix_t <- as.data.frame(expr_matrix_t)
  
  #all_genes <- colnames(expr_matrix_t)
  
  # Add sample information to the data
  expr_matrix_t$patient_id <- obj_CT@meta.data$patient_id
  
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
  
  # mean_expression_per_sample will be used for PCA but need to add rejection type
  # info for plotting
  mean_expression_per_sample$patient_id <- rownames(mean_expression_per_sample)
  
  tmp <- obj_CT@meta.data %>%
    dplyr::select(patient_id, patient_biopsy, biopsy_timing, #rejection_type, 
                  sex, donor_sex, #cellular_grading, antibody_grading, 
                  treatment,
                  patient_cellular_grading, patient_antibody_grading, patient_rejection_type) %>%
    unique()
  mean_expression_per_sample_info <- left_join(mean_expression_per_sample, tmp)
  rownames(mean_expression_per_sample_info) <- rownames(mean_expression_per_sample)
  
  
  ## All samples, pesudo-bulked by patient and biopsy timing across all cell types together ----
  
  smpl.pca <- PCA(mean_expression_per_sample_info[,1:477], graph = F)
  
  
  # Generate PCA plot
  pca_plots[[ct]] <- fviz_pca_ind(smpl.pca,
                                  label = "none", # hide individual labels
                                  habillage = as.factor(mean_expression_per_sample_info$patient_rejection_type), # color by groups
                                  addEllipses = F, # Concentration ellipses
                                  title = ct
  ) +
    labs(color = "Rejection type", shape = "Rejection type") +
    scale_color_manual(values = custom_colors) +
    scale_shape_manual(values = custom_shapes) +
    theme(plot.title = element_text(hjust = 0.5))
  
  
}

# Arrange all plots into a grid
pdf("/home/aoill/plots/00_heart_manuscript/pca_exp_cell_type.pdf",
    width = 17,
    height = 19)
ggarrange(plotlist = pca_plots, ncol = 5, nrow = ceiling(length(pca_plots) / 5),
          common.legend = TRUE, legend = "bottom")
dev.off()


