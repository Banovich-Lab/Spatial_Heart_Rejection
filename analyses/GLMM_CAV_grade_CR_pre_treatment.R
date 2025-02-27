# Load required libraries
library(Seurat)
library(dplyr)
library(glmmTMB)
library(stats)

# Read in object ---- 
clustered_obj <- readRDS("/scratch/aoill/projects/heart_transplant/GEO/heart_spatial_obj.rds")

# CR pre-treatment - CAV score ----
## Subset CR pre-treatment samples ----
pre_CR_smpls <-  clustered_obj@meta.data %>% 
  filter(biopsy_rejection_type == "cellular_rejection") %>%
  filter(biopsy_timing == "pre") %>%
  dplyr::pull(Sample) %>% unique()


clustered_obj_linear_regression <- subset(clustered_obj, 
                                          subset = Sample %in% pre_CR_smpls)


## Get cell types to analyze ----
cts_to_analyze <- levels(as.factor(clustered_obj_linear_regression@meta.data$ct_second_pass))


## Run linear regression for each gene in each cell type ----
all_CT_results <- list()

for (ctn in 1:length(cts_to_analyze)) {
  ct <- cts_to_analyze[ctn]
  print(ct)
  
  # Step: get average expression for each biopsy (Sample) in the seurat object 
  # for the cell type
  obj_CT <- subset(clustered_obj_linear_regression, subset = ct_second_pass == ct)
  
  # Extract the normalized expression data (data is log-transformed counts)
  expr_matrix <- as.data.frame(as.matrix(obj_CT@assays$RNA@data))
  expr_matrix_t <- t(expr_matrix)
  expr_matrix_t <- as.data.frame(expr_matrix_t)
  
  all_genes <- colnames(expr_matrix_t)
  
  # Add sample information to the data
  expr_matrix_t$Sample <- obj_CT@meta.data$Sample
  
  # Summarize by Sample to calculate the mean expression for each gene and sample
  mean_expression_per_sample <- expr_matrix_t %>%
    group_by(Sample) %>%
    summarise(across(everything(), mean, na.rm = TRUE))
  
  # Convert back to a data frame 
  mean_expression_per_sample <- as.data.frame(mean_expression_per_sample)
  
  
  # Step: For each gene run linear regression including patient ID in the model 
  lm_results <- c()
  
  # Initialize an empty list to store failed gene results
  failed_genes <- list()

  for (gene in 1:length(all_genes)) {
    genei <- all_genes[gene]
    print(genei)
    
    # select gene
    genei_expression_counts <- mean_expression_per_sample %>%
      select(Sample, genei)
    
    # add sample info
    outcome_data <- obj_CT@meta.data %>% 
      dplyr::select(Sample, patient_id, cav_grade, age_at_biopsy) %>%
      unique()
    
    genei_expression_counts_info <- left_join(genei_expression_counts, outcome_data)
    colnames(genei_expression_counts_info) <- c("Sample", "expression", "patient_id", 
                                                "cav_grade", "age_at_biopsy")
    genei_df <- genei_expression_counts_info

    
    # Fit the full model using tryCatch
    tryCatch({
      # Fit the model
      LM_result <- glmmTMB(expression ~ cav_grade + age_at_biopsy + (1 | patient_id), 
                           data = genei_df, family = gaussian())
      
      # Check for convergence issues in the summary
      model_summary <- summary(LM_result)
      if (any(grepl("Model convergence problem", model_summary$warnings))) {
        # If there is a convergence problem, skip this gene
        warning(paste("Skipping gene", genei, "due to convergence issues"))
        
        # Optionally store the failed gene information
        failed_genes[[genei]] <- list("gene" = genei, "error" = "Model convergence problem")
        next # Skip to the next iteration of the loop
      }
      
      # Extract coefficients (fixed effects - cav_grade)
      coef_summary <- model_summary$coefficients$cond
      beta <- coef_summary["cav_grade", "Estimate"]
      std_error <- coef_summary["cav_grade", "Std. Error"]
      p_value <- coef_summary["cav_grade", "Pr(>|z|)"]
      
      # Add results to a data frame (you can customize this as needed)
      genei_df_no_0 <- genei_df %>% filter(expression > 0)
      n_samples <- nrow(genei_df)
      n_samples_no0 <- nrow(genei_df_no_0)
      row_to_add <- cbind(ct, genei, n_samples, n_samples_no0, beta, std_error, p_value)
      lm_results <- rbind(lm_results, row_to_add)
      
    }, warning = function(w) {
      # If a warning occurs, catch it and check if it's a convergence warning
      if (grepl("Model convergence problem", w$message)) {
        message(paste("Skipping gene", genei, "due to convergence issues"))
        # Optionally log the failed gene
        failed_genes[[genei]] <- list("gene" = genei, "warning" = w$message)
      } else {
        # For other warnings, proceed as usual
        warning(w)
      }
    }, error = function(e) {
      # If an error occurs, skip the gene and log the error
      message(paste("Skipping gene", genei, "due to an error:", e$message))
      failed_genes[[genei]] <- list("gene" = genei, "error" = e$message)
    })
    
  }
  
  
  lm_results_df <- as.data.frame(lm_results)
  # Expression in at least 5 samples
  lm_results_df_clean <- lm_results_df %>% filter(as.numeric(n_samples_no0) >4) 
  # Correct for multiple testing
  lm_results_df_clean$padj <- p.adjust(lm_results_df_clean$p_value, method = "fdr")
  
  #Add results to a list of list
  all_CT_results <- append(all_CT_results, list(lm_results_df_clean))
  
  #saveRDS(all_CT_results, 
  #        "/scratch/aoill/projects/heart_transplant/GLMM_results_CAV_CR_all_rm_model_errors_fix_exp_20241126_in_progress_results.rds")
  
}


names(all_CT_results) <- cts_to_analyze
cr_pre_cav_grade_all_CT_results <- all_CT_results

# Filter out empty lists
cr_pre_cav_grade_all_CT_results_filter <- cr_pre_cav_grade_all_CT_results[sapply(cr_pre_cav_grade_all_CT_results, length) != 0]
cr_pre_cav_grade_all_CT_results_filter <- cr_pre_cav_grade_all_CT_results_filter[sapply(cr_pre_cav_grade_all_CT_results_filter, nrow) != 0]

# Add cell type name as a column to each data frame in all_CT_results
#for (ct in names(cr_pre_cav_grade_all_CT_results_filter)) {
#  cr_pre_cav_grade_all_CT_results_filter[[ct]]$ct <- ct
#}
#saveRDS(cr_pre_cav_grade_all_CT_results_filter, 
#        "/scratch/aoill/projects/heart_transplant/GLMM_results_CAV_CR_all_rm_model_errors_fix_exp_20241126_FIX.rds")



# Merge the list of data frames into one data frame
cr_pre_cav_grade_all_CT_results_filter_merged <- do.call(rbind, cr_pre_cav_grade_all_CT_results_filter)
unique(cr_pre_cav_grade_all_CT_results_filter_merged$ct)
table(cr_pre_cav_grade_all_CT_results_filter_merged$ct)


# Keep only significant genes
cr_pre_cav_grade_all_CT_results_filter_merged_sig <- cr_pre_cav_grade_all_CT_results_filter_merged %>%
  filter(padj <= 0.05)
cr_pre_cav_grade_all_CT_results_filter_merged_sig$n_samples_no0 <- as.numeric(cr_pre_cav_grade_all_CT_results_filter_merged_sig$n_samples_no0
                                                                              )
nrow(cr_pre_cav_grade_all_CT_results_filter_merged_sig)
length(unique(cr_pre_cav_grade_all_CT_results_filter_merged_sig$genei))
length(unique(cr_pre_cav_grade_all_CT_results_filter_merged_sig$ct))
#   padj 0.05 - 116 significant genes across 14 cell types



# save merged as a csv 
# THIS IS SUPPLEMENTARY TABLE 8
write.csv(cr_pre_cav_grade_all_CT_results_filter_merged, "/home/aoill/projects/heart_transplant/00_final/GLMM_CAV_grade_CR_pre_treatment_biopsies_NEW.csv",
          row.names = F, quote = F)
