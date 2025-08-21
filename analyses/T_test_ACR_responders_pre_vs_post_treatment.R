################################################################################
# Angela Oill
# Analysis: T-test comparing ACR responders pre- vs post-treatment 
################################################################################

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


clustered_obj <- readRDS("/scratch/aoill/projects/heart_transplant/GSE290577_heart_spatial_obj.rds")
clustered_obj_meta_rns <- rownames(clustered_obj@meta.data)
clustered_obj_meta <- clustered_obj@meta.data
# Create a new column to classify each cell type into its Lineage
clustered_obj@meta.data <- clustered_obj@meta.data %>%
  mutate(Lineage = case_when(
    ct_second_pass %in% endothelial_cells ~ "Endothelial",
    ct_second_pass %in% lymphatic_endothelial ~ "Lymphatic_Endothelial",
    ct_second_pass %in% cms ~ "Cardiomyocyte",
    ct_second_pass %in% adipocyte ~ "Adipocytes",
    ct_second_pass %in% fibroblast ~ "Fibroblasts",
    ct_second_pass %in% mural ~ "Mural",
    ct_second_pass %in% mast ~ "Mast",
    ct_second_pass %in% myeloid ~ "Myeloid",
    ct_second_pass %in% lymphoid ~ "Lymphoid",
    TRUE ~ "Other"
  ))
rownames(clustered_obj@meta.data) <- clustered_obj_meta_rns
unique(clustered_obj@meta.data$Lineage)


# CR Patients that had resolution pre vs post treatment
# Select CR patients (patient_rejection_type)
# start with clustered_obj and patient_rejection_type
cellular_rejection_smpls <- clustered_obj@meta.data %>% 
  filter(patient_rejection_type == "cellular_rejection") %>%
  filter(resolution_broad == "resolved") %>%
  filter(patient_id != "Patient_17") %>%
  dplyr::pull(patient_id) %>% unique()

clustered_obj_pre_post_resolved <- subset(clustered_obj, subset = patient_id %in% cellular_rejection_smpls)
cts_to_analyze <- levels(as.factor(clustered_obj_pre_post_resolved@meta.data$ct_second_pass))

# sample numbers
clustered_obj_pre_post_resolved@meta.data$patient_id
clustered_obj_pre_post_resolved@meta.data$biopsy_timing

tmp <- clustered_obj_pre_post_resolved@meta.data %>% dplyr::select(patient_id, biopsy_timing) %>% unique()
table(tmp$biopsy_timing)



# read in results file
# see: scRNA_heart_expression_comparison.ipynb
gene_prop_info_atlas <- read.csv("/scratch/aoill/projects/heart_transplant/scrna_atlas/genes_to_keep_proportions_new_groups.csv")

# gene genes with at least x percent expression for the different lineages
gene_prop_info_atlas_filter <- gene_prop_info_atlas %>%
  filter(Lineage != "Neural") %>%
  filter(Lineage != "Mesothelial") %>%
  filter(Proprtion_Expressed >= 0.01)


table(gene_prop_info_atlas_filter$Lineage)


### Run test ----
t_test_results <- c()

for (ctn in 1:length(cts_to_analyze)) {
  ct <- cts_to_analyze[ctn]
  print(ct)
  
  # subset object to cell type
  obj_CT <- subset(clustered_obj_pre_post_resolved, subset = ct_second_pass == ct)
  
  # grab what lineage the cell type is from
  lin <- unique(obj_CT@meta.data$Lineage)
  # get which genes to include in analysis
  genes_to_analyze <- gene_prop_info_atlas_filter %>% 
    filter(Lineage == lin) %>% 
    dplyr::pull(Gene)
  
  # Extract the normalized expression data (data is log-transformed counts)
  expr_matrix <- as.data.frame(as.matrix(obj_CT@assays$RNA@data))
  expr_matrix_t <- t(expr_matrix)
  expr_matrix_t <- as.data.frame(expr_matrix_t)
  

  # Add sample information to the data
  expr_matrix_t$patient_biopsy <- obj_CT@meta.data$patient_biopsy
  
  # Summarize by patient_biopsy to calculate the mean expression for each gene and sample
  mean_expression_per_sample <- expr_matrix_t %>%
    group_by(patient_biopsy) %>%
    summarise(across(everything(), mean, na.rm = TRUE))
  
  # Convert back to a data frame 
  mean_expression_per_sample <- as.data.frame(mean_expression_per_sample)
  # add sample info
  sample_info <- obj_CT@meta.data %>% dplyr::select(patient_biopsy, biopsy_timing) %>% unique()
  mean_expression_per_sample_info <- left_join(mean_expression_per_sample, sample_info)
  
  
  if (min(table(mean_expression_per_sample_info$biopsy_timing)) < 3 || (length(table(mean_expression_per_sample_info$biopsy_timing)) < 2) ) {
    print(paste("Not enough samples or groups to do test. Cell type: ", ct, sep = ""))
    
    print(table(mean_expression_per_sample_info$biopsy_timing))
    
    # add empyt list so this gets added to list of results instead of the previous cell types results
    t_test_results_all_genes <- c()
    
  } else {
    
    t_test_results_all_genes <- c()
    for (genei in genes_to_analyze) {
      #print(genei)
      # Expression filter
      #sum(genei_mean_expression_per_sample$expression > 0) >= 3
      genei_mean_expression_per_sample <- mean_expression_per_sample_info %>% dplyr::select(patient_biopsy, biopsy_timing, genei)
      colnames(genei_mean_expression_per_sample) <- c("patient_biopsy", "biopsy_timing", "expression")
      
      if (sum(genei_mean_expression_per_sample$expression > 0) >= 3) {
        #genei_mean_expression_per_sample <- mean_expression_per_sample_info %>% dplyr::select(patient_biopsy, biopsy_timing, genei)
        #colnames(genei_mean_expression_per_sample) <- c("patient_biopsy", "biopsy_timing", "expression")
        
        # do T test
        ttest_result_ct_gene <- t.test(
          genei_mean_expression_per_sample%>% filter(biopsy_timing == "pre") %>% select("expression"),
          genei_mean_expression_per_sample%>% filter(biopsy_timing == "post") %>% select("expression"),
        )
        
        
        # Calculate logFC
        # Define a pseudocount to handle cases where mean expression might be zero.
        pseudocount <- .1
        
        # 1. Calculate mean expression for pretreatment
        mean_pre <- mean(genei_mean_expression_per_sample$expression[genei_mean_expression_per_sample$biopsy_timing == "pre"])
        
        # 2. Calculate mean expression for post
        mean_post <- mean(genei_mean_expression_per_sample$expression[genei_mean_expression_per_sample$biopsy_timing == "post"])
        
        # 3. Add pseudocount to the means
        mean_pre_adjusted <- mean_pre + pseudocount
        mean_post_adjusted <- mean_post + pseudocount
        
        # 4. Calculate the log2 fold change
        # Formula: log2(Mean Expression in pre / Mean Expression in post)
        # A positive value means higher expression in pre
        # A negative value means higher expression in post
        log2_fold_change <- log2(mean_pre_adjusted / mean_post_adjusted)
        #log2(mean_cellular_rejection / mean_antibody_rejection)  ``
        
        
        # Add all results to larger dataframe/list
        t_result_genei <- as.data.frame(cbind(ct, genei, ttest_result_ct_gene$statistic, log2_fold_change, ttest_result_ct_gene$p.value))
        colnames(t_result_genei) <- c("ct", "gene", "t_statistic", "Log2FC", "p_value")
        rownames(t_result_genei) <- NULL
        
        t_test_results_all_genes <- rbind(t_test_results_all_genes, t_result_genei)
        
      } else {
        print(paste("Not enough samples with expression in gene: ", genei))
      }
      
      
    }
    
    # Add multiple testing correction
    t_test_results_all_genes$padj <- p.adjust(t_test_results_all_genes$p_value, method = "fdr")
    
  }
  
  t_test_results[[ct]] <- t_test_results_all_genes
}

# Merge the list of data frames into one data frame
all_results_merged <- do.call(rbind, t_test_results)
all_results_merged$t_statistic <- as.numeric(all_results_merged$t_statistic)
all_results_merged$Log2FC <- as.numeric(all_results_merged$Log2FC)
all_results_merged$p_value <- as.numeric(all_results_merged$p_value)
all_results_merged$padj <- as.numeric(all_results_merged$padj)

all_results_merged_sig <- all_results_merged %>% filter(padj <=0.05)
length(unique(all_results_merged_sig$ct))
length(unique(all_results_merged_sig$gene))
(nrow(all_results_merged_sig))
# 8 sig genes across 5 cell types, 8 gene cell type pairs



# Define the groups based on cell types
endothelial_cells <- c("Activated endothelial", "BMX+ Activated endothelial", 
                       "Endothelial", "Lymphatic endothelial", "Proliferating endothelial")

mesenchymal_cells <- c("Adipocytes", "Fibroblasts", 
                       "Myofibroblasts", "Pericytes", "POSTN+ Fibroblasts", 
                       "Proliferating pericytes", "vSMCs")
cm <- c("Cardiomyocytes")

immune_cells <- c("B cells", "CD4+ T cells", "CD8+ T cells", "cDC1", "cDC2", 
                  "Macrophages", "Mast", "mDC", "NK", "pDC", "Plasma", 
                  "Proliferating DCs", "Proliferating T cells", 
                  "SPP1+ Macrophages", "Treg")

# Combine all data frames into one
combined_results <- bind_rows(t_test_results, .id = "ct")

# Create a new column to classify each cell type into its lineage
combined_results <- combined_results %>%
  mutate(lineage = case_when(
    ct %in% endothelial_cells ~ "Endothelial",
    ct %in% mesenchymal_cells ~ "Mesenchymal",
    ct %in% cm ~ "Cardiomyocytes",
    ct %in% immune_cells ~ "Immune",
    TRUE ~ "Other"
  ))

# Create a new column to categorize significance (padj < 0.1)
combined_results <- combined_results %>%
  mutate(significance = ifelse(padj < 0.05, "Significant", "Non-significant"))




### Save output ----
out_table <- combined_results %>% dplyr::select(t_statistic, Log2FC, p_value, padj, gene, ct, lineage, significance)
colnames(out_table) <- c("t_statistic", "Log2FC", "p_value", "padj", "gene", "cell_type", "lineage", "significance")
head(out_table)

write_csv(out_table, "/home/aoill/projects/heart_transplant/00_final/revisions/t_test_ACR_responders_pre_vs_post_treatment_results.csv")
#out_table <- read.csv("/home/aoill/projects/heart_transplant/00_final/revisions/t_test_ACR_responders_pre_vs_post_treatment_results.csv")
