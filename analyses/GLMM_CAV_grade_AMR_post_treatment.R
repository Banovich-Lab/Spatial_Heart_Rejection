# Load required libraries
library(Seurat)
library(dplyr)
library(glmmTMB)
library(stats)
library(tidyverse)


# Read in object ---- 
clustered_obj <- readRDS("/scratch/aoill/projects/heart_transplant/GEO/heart_spatial_obj.rds")

# Get list of genes to analyze from each cell type group ----
# Define the groups based on cell types
endothelial_cells <- c("Activated endothelial", "BMX+ Activated endothelial", 
                       "Endothelial", "Proliferating endothelial")

lymphatic_endothelial <- c("Lymphatic endothelial")

cms <- c("Cardiomyocytes")

adipocyte <- c("Adipocytes")

fibroblast <- c("Fibroblasts", "Myofibroblasts", "POSTN+ Fibroblasts")

mural <- c("Pericytes", "Proliferating pericytes", "vSMCs")

myeloid <- c("cDC1", "cDC2", "Macrophages", "mDC", 
             "Proliferating DCs", 
             "SPP1+ Macrophages")

mast <- c("Mast")

lymphoid <- c("B cells", "CD4+ T cells", "CD8+ T cells", "NK", "pDC", "Plasma", 
              "Proliferating T cells", "Treg")

clustered_obj_meta_rns <- rownames(clustered_obj@meta.data)
clustered_obj_meta <- clustered_obj@meta.data
# Create a new column to classify each cell type into its Lineage
clustered_obj@meta.data <- clustered_obj@meta.data %>%
  mutate(Lineage = case_when(
    ct_second_pass %in% endothelial_cells ~ "Endothelial",
    ct_second_pass %in% lymphatic_endothelial ~ "Lymphatic_Endothelial",
    ct_second_pass %in% cms ~ "Cardiomyocyte",
    ct_second_pass %in% adipocyte ~ "Adipocyte",
    ct_second_pass %in% fibroblast ~ "Fibroblast",
    ct_second_pass %in% mural ~ "Mural",
    ct_second_pass %in% mast ~ "Mast",
    ct_second_pass %in% myeloid ~ "Myeloid",
    ct_second_pass %in% lymphoid ~ "Lymphoid",
    TRUE ~ "Other"
  ))
rownames(clustered_obj@meta.data) <- clustered_obj_meta_rns
unique(clustered_obj@meta.data$Lineage)

# read in results file from heart atlas
# this has the proportion of cells in a group that express each gene in our panel
gene_prop_info_atlas <- read.csv("/scratch/aoill/projects/heart_transplant/scrna_atlas/genes_to_keep_proportions_new_groups.csv")

# gene genes with at least x percent expression for the different lineages
# we do not have Neural or Mesothelial cells
gene_prop_info_atlas_filter <- gene_prop_info_atlas %>%
  filter(Lineage != "Neural") %>%
  filter(Lineage != "Mesothelial") %>%
  filter(Proprtion_Expressed >= 0.01)

table(gene_prop_info_atlas_filter$Lineage)

# AMR post-treatment - CAV score ----
## Subset AMR post-treatment samples ----
pre_AMR_smpls <-  clustered_obj@meta.data %>% 
  filter(biopsy_rejection_type == "antibody_rejection ") %>%
  filter(biopsy_timing == "post") %>%
  dplyr::pull(Sample) %>% unique()


clustered_obj_linear_regression <- subset(clustered_obj, 
                                          subset = Sample %in% pre_AMR_smpls)


## Get cell types to analyze ----
cts_to_analyze <- levels(as.factor(clustered_obj_linear_regression@meta.data$ct_second_pass))



## Run linear regression for each gene in each cell type ----
all_CT_results <- list()
all_CT_results_with_issues <- list()

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
  
  # get genes to analyze for this particular cell type
  #all_genes <- colnames(expr_matrix_t)
  # grab what lineage the cell type is from
  lin <- unique(obj_CT@meta.data$Lineage)
  # get which genes to include in analysis
  genes_to_analyze <- gene_prop_info_atlas_filter %>% 
    filter(Lineage == lin) %>% 
    dplyr::pull(Gene)
  
  # Add sample information to the data
  expr_matrix_t$Sample <- obj_CT@meta.data$Sample
  
  # Summarize by Sample to calculate the mean expression for each gene and sample
  mean_expression_per_sample <- expr_matrix_t %>%
    dplyr::group_by(Sample) %>%
    dplyr::summarise(across(everything(), mean, na.rm = TRUE))
  
  # Convert back to a data frame 
  mean_expression_per_sample <- as.data.frame(mean_expression_per_sample)
  
  
  # Remove cell types with less than 5 samples. Will also remove genes without
  # expression in at least 5 samples
  n_samples <- nrow(mean_expression_per_sample)
  
  if (n_samples < 5) {
    print(paste("Not enough samples to do test. Cell type: ", ct, sep = ""))
    lm_results <- c()
    
  } else {
    # Step: For each gene run linear regression including patient ID in the model 
    lm_results <- c()
    for (gene in 1:length(genes_to_analyze)) {
      genei <- genes_to_analyze[gene]
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
      
      # Only analyze genes where there is expression in at least 5
      # samples
      n_samples_no0 <- nrow(genei_df %>% filter(expression > 0))
      if (n_samples_no0 < 5) {
        print(paste("Not enough samples with expression to do test. Gene: ", genei, sep = ""))
      } else {
        # Fit the model
        LM_result <- glmmTMB(expression ~ cav_grade + age_at_biopsy + (1 | patient_id), 
                             data = genei_df, family = gaussian())
        # Check for convergence issues in the summary
        model_summary <- summary(LM_result)
        
        ## get variance from the test
        patient_id_variance <- model_summary$varcor$cond$patient_id
        patient_id_std_dev_numeric <- as.numeric(sqrt(patient_id_variance))
        
        convergence_val <- LM_result$fit$convergence
        convergence_message <- LM_result$fit$message
        
        # Extract coefficients (fixed effects - cav_grade)
        coef_summary <- model_summary$coefficients$cond
        beta <- coef_summary["cav_grade", "Estimate"]
        std_error <- coef_summary["cav_grade", "Std. Error"]
        p_value <- coef_summary["cav_grade", "Pr(>|z|)"]
        
        # Add results to a data frame (you can customize this as needed)
        genei_df_no_0 <- genei_df %>% filter(expression > 0)
        n_samples <- nrow(genei_df)
        n_samples_no0 <- nrow(genei_df_no_0)
        row_to_add <- cbind(ct, genei, n_samples, n_samples_no0, beta, std_error, p_value, convergence_val, convergence_message, patient_id_std_dev_numeric)
        lm_results <- rbind(lm_results, row_to_add)
        
      } # END OF GENE FILTER IF ELSE STATEMENT
      
    } # GENE FOR LOOP ENDS HERE
    
    
    lm_results_df <- as.data.frame(lm_results)
    
    # Remove gene cell type pairs with convergence issues and NAs (likely issue 
    # with patient_id_std_dev being close to 0 meaning the estimated variance 
    # for the patient_id random effect was found to be effectively zero...it's 
    # possible we just do a simpler model without the random effect...but not 
    # sure if that is what I want to do)
    lm_results_df$p_value <- as.numeric(lm_results_df$p_value)
    lm_results_df_no_conv_iss <- lm_results_df %>% filter(convergence_val == 0) %>%
      drop_na(p_value)
    # Correct for multiple testing
    lm_results_df_no_conv_iss$padj <- p.adjust(lm_results_df_no_conv_iss$p_value, method = "fdr")
    
    # Add results to a list of list
    all_CT_results[[ct]] <- lm_results_df_no_conv_iss
    all_CT_results_with_issues[[ct]] <- lm_results_df
    
    
  } # END OF n_samples IF ELSE STATEMENT 
  
  
} # END OF CELL TYPE FOR LOOP


# Prep and output results
amr_post_cav_grade_all_CT_results <- all_CT_results

saveRDS(amr_post_cav_grade_all_CT_results, 
        "/scratch/aoill/projects/heart_transplant/00_final/GLMM_CAV_grade_AMR_post_treatment_biopsies_revisions.rds")
saveRDS(all_CT_results_with_issues, 
        "/scratch/aoill/projects/heart_transplant/00_final/GLMM_CAV_grade_AMR_post_treatment_biopsies_revisions_with_issues.rds")



amr_post_cav_grade_all_CT_results <- readRDS("/scratch/aoill/projects/heart_transplant/00_final/GLMM_CAV_grade_AMR_post_treatment_biopsies_revisions.rds")

for (i in names(amr_post_cav_grade_all_CT_results)) {
  print(i)
  amr_post_cav_grade_all_CT_results[[i]]$p_value <- as.numeric(amr_post_cav_grade_all_CT_results[[i]]$p_value)
  amr_post_cav_grade_all_CT_results[[i]] <- amr_post_cav_grade_all_CT_results[[i]] %>%
    drop_na(p_value)
  amr_post_cav_grade_all_CT_results[[i]]$padj <- p.adjust(amr_post_cav_grade_all_CT_results[[i]]$p_value, method = "fdr")
}

saveRDS(amr_post_cav_grade_all_CT_results, 
        "/scratch/aoill/projects/heart_transplant/00_final/GLMM_CAV_grade_AMR_post_treatment_biopsies_revisions_fixed.rds")

all_CT_results_with_issues_merged <- do.call(rbind, all_CT_results_with_issues)

# Merge the list of data frames into one data frame
amr_post_cav_grade_all_CT_results_merged <- do.call(rbind, amr_post_cav_grade_all_CT_results)
unique(amr_post_cav_grade_all_CT_results_merged$ct)
table(amr_post_cav_grade_all_CT_results_merged$ct)

# Keep only significant genes
amr_post_cav_grade_all_CT_results_merged_sig <- amr_post_cav_grade_all_CT_results_merged %>%
  filter(padj <= 0.05)
amr_post_cav_grade_all_CT_results_merged_sig$n_samples_no0 <- as.numeric(amr_post_cav_grade_all_CT_results_merged_sig$n_samples_no0
)
nrow(amr_post_cav_grade_all_CT_results_merged_sig)
length(unique(amr_post_cav_grade_all_CT_results_merged_sig$genei))
length(unique(amr_post_cav_grade_all_CT_results_merged_sig$ct))
#   padj 0.05 - 59 significant genes across 12 cell types, 66 gene cell type pairs
#   padj 0.05 - 91 significant genes across 17 cell types, 119 gene cell type pairs



# save merged as a csv 
# THIS IS SUPPLEMENTARY TABLE 8
write.csv(amr_post_cav_grade_all_CT_results_merged, "/home/aoill/projects/heart_transplant/00_final/revisions/GLMM_CAV_grade_AMR_post_treatment_biopsies_revisions_fixed.csv",
          row.names = F, quote = F)

write.csv(all_CT_results_with_issues_merged, "/home/aoill/projects/heart_transplant/00_final/revisions/GLMM_CAV_grade_AMR_post_treatment_biopsies_revisions_with_issues.csv",
          row.names = F, quote = F)


# Make a table with which gene cell type pairs had convergenc issues or small 
# variance issues
no_iss_results <- amr_post_cav_grade_all_CT_results_merged
iss_results <- all_CT_results_with_issues_merged
no_iss_results$padj <- NULL

nrow(setdiff(iss_results, no_iss_results))

glmms_with_issues <- setdiff(iss_results, no_iss_results)

glmms_with_issues_tbl <- as.data.frame(table(glmms_with_issues$convergence_message))
glmms_with_issues_tbl$test <- "AMR post-treatment across CAV"
colnames(glmms_with_issues_tbl) <- c("Convergence message", "Freq", "Test")

glmms_with_issues_tbl <- glmms_with_issues_tbl %>% dplyr::select(Test, `Convergence message`, Freq)

# output table
write.csv(glmms_with_issues_tbl, 
          "/home/aoill/projects/heart_transplant/00_final/revisions/GLMM_CAV_grade_AMR_post_treatment_biopsies_convergence_issues.csv",
          row.names = F, quote = F)


