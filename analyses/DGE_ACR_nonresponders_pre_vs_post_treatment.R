################################################################################
# Spatial Heart Rejection Analysis code
# Analysis: Differential gene expression ACR nonresponders pre- vs post-treatment 
################################################################################

# Load libraries ----
library(Seurat)
library(limma)
library(edgeR)

# Read in object ---- 
clustered_obj <- readRDS("/scratch/aoill/projects/heart_transplant/GEO/heart_spatial_obj.rds")


# Select CR patients (patient_rejection_type)
# start with clustered_obj and patient_rejection_type
cellular_rejection_smpls <- clustered_obj@meta.data %>% 
  filter(patient_rejection_type == "cellular_rejection") %>%
  filter(resolution_broad == "persistence") %>%
  dplyr::pull(patient_id) %>% unique()

clustered_obj_pre_post_persistence <- subset(clustered_obj, subset = patient_id %in% cellular_rejection_smpls)
cts_to_analyze <- levels(as.factor(clustered_obj_pre_post_persistence@meta.data$ct_second_pass))
table(clustered_obj_pre_post_persistence$biopsy_timing)

# Loop through each cell type and perform DGE analysis on all cell types
all_results <- c()
all_results_efit <- c()

for (ctn in 1:length(cts_to_analyze)) {
  ct <- cts_to_analyze[ctn]
  print(ct)
  
  # subset object to cell type
  obj_CT <- subset(clustered_obj_pre_post_persistence, subset = ct_second_pass == ct)
  
  # Extract the normalized expression data (data is log-transformed counts)
  expr_matrix <- as.data.frame(as.matrix(obj_CT@assays$RNA@data))
  expr_matrix_t <- t(expr_matrix)
  expr_matrix_t <- as.data.frame(expr_matrix_t)
  
  #all_genes <- colnames(expr_matrix_t)
  
  # Add sample information to the data
  expr_matrix_t$patient_biopsy <- obj_CT@meta.data$patient_biopsy
  
  # Summarize by Sample to calculate the mean expression for each gene and sample
  mean_expression_per_sample <- expr_matrix_t %>%
    group_by(patient_biopsy) %>%
    summarise(across(everything(), mean, na.rm = TRUE))
  
  # Convert back to a data frame 
  mean_expression_per_sample <- as.data.frame(mean_expression_per_sample)
  rownames(mean_expression_per_sample) <- mean_expression_per_sample$patient_biopsy
  mean_expression_per_sample$patient_biopsy <- NULL
  #head(t(mean_expression_per_sample))
  avg_expression_df <- t(mean_expression_per_sample)
  
  # Make design matrix. Need to make sure the order of the categories is the same
  # as the order of the columns in the average expression data frame
  # Get the correct order for the rejection type
  rejection_cat <- c()
  for (i in colnames(avg_expression_df)) {
    #print(i)
    # get just the sample ID (this is the sample ID plus cluster number)
    # can use the sub function for this
    #tmp <- sub("-[^-]+$", "", i) # might need to replace dashes with underscores
    #sample_id <- sub("-(?=[^\\-]*$)", "_", i, perl = TRUE)
    #sample_id <- gsub("-", "_", i)
    #print(sample_id)
    sample_id <- i
    
    # get rejection type for this sample
    i_rt <- obj_CT@meta.data %>%
      filter(patient_biopsy == sample_id) %>%
      dplyr::select(patient_biopsy, biopsy_timing) %>%
      unique() %>% dplyr::pull(biopsy_timing)
    
    rejection_cat <- c(rejection_cat, i_rt)
  }
  rejection_cat <- as.factor(rejection_cat)
  
  # check the number of levels 
  #length(levels(rejection_cat))
  if (length(rejection_cat) <3) {
    print(paste("Skipping, ", ct, " because too few samples. Number of samples = ", length(rejection_cat), sep = ""))
  } else {
    # Set design matrix
    design <- model.matrix(~0+rejection_cat)
    
    # remove rejection_cat from column name to make cleaner
    colnames(design) <- gsub("rejection_cat", "", colnames(design))
    colnames(design) <- gsub(" ", "", colnames(design))
    
    # Set contrasts
    contr.matrix <- makeContrasts(
      PrevsPost = post-pre,
      #AntibodyvsMixed = antibody_rejection-mixed_rejection,
      #CellularvsMixed = cellular_rejection-mixed_rejection,
      levels = colnames(design))
    #contr.matrix
    
    d0 <- DGEList(avg_expression_df)
    
    # Filtering low-exp genes 
    #keep <- rowSums(avg_expression_df) > 0
    #d <- d0[keep,] 
    #dim(d0)
    #dim(d)
    
    # I want to set a filter to keep genes that have expression in at least 3 samples
    keep <- rownames(avg_expression_df)[rowSums(avg_expression_df > 0) >= 3]
    d <- d0[keep,] 
    
    dim(d0)
    dim(d)
    
    
    # NEW
    y <- voom(d, design, plot = T)
    fit <- lmFit(y, design)
    tmp <- contrasts.fit(fit, contr.matrix)
    tmp <- eBayes(tmp)
    # tmp is efit
    all_results[[ct]] <- topTable(tmp, number = 477)
    all_results[[ct]]$gene <- rownames(all_results[[ct]])
    all_results[[ct]]$cell_type <- ct
    all_results_efit[[ct]] <- tmp
  }
  
}


# Merge the list of data frames into one data frame
all_results_merged <- do.call(rbind, all_results)

all_results_merged_sig <- all_results_merged %>% filter(adj.P.Val <=0.05)

length(unique(all_results_merged_sig$cell_type))
length(unique(all_results_merged_sig$gene))
# p 0.05
# Corrected analysis fixed expression and more stringent filtering -
#     14 genes across 3 cell types 

# Save results
#saveRDS(all_results, "/scratch/aoill/projects/heart_transplant/00_final/DGE_pre_vs_post_CR_persistence_NEW_CORRECTED_no_calcnormfactor_exp_fixed.rds")
#all_results <- readRDS("/scratch/aoill/projects/heart_transplant/DGE_pre_vs_post_CR_persistence_NEW_CORRECTED_no_calcnormfactor_exp_fixed.rds")


# Define the groups based on cell types
endothelial_cells <- c("Activated endothelial", "BMX+ Activated endothelial", 
                       "Endothelial", "Lymphatic endothelial", "Proliferating endothelial")

mesenchymal_cells <- c("Adipocytes", "Cardiomyocytes", "Fibroblasts", 
                       "Myofibroblasts", "Pericytes", "POSTN+ Fibroblasts", 
                       "Proliferating pericytes", "vSMCs")

immune_cells <- c("B cells", "CD4+ T cells", "CD8+ T cells", "cDC1", "cDC2", 
                  "Macrophages", "Mast", "mDC", "NK", "pDC", "Plasma", 
                  "Proliferating DCs", "Proliferating T cells", 
                  "SPP1+ Macrophages", "Treg")

# Combine all data frames into one
combined_results <- bind_rows(all_results, .id = "cell_type")

# Create a new column to classify each cell type into its lineage
combined_results <- combined_results %>%
  mutate(lineage = case_when(
    cell_type %in% endothelial_cells ~ "Endothelial",
    cell_type %in% mesenchymal_cells ~ "Mesenchymal",
    cell_type %in% immune_cells ~ "Immune",
    TRUE ~ "Other"
  ))

# Create a new column to categorize significance (adj.P.Val < 0.1)
combined_results <- combined_results %>%
  mutate(significance = ifelse(adj.P.Val < 0.05, "Significant", "Non-significant"))


# save merged as a csv 
# THIS IS SUPPLEMENTARY TABLE 5
write.csv(combined_results, "/home/aoill/projects/heart_transplant/00_final/DGE_pre_vs_post_CR_persistence_exp_fixed.csv",
          row.names = F, quote = F)
