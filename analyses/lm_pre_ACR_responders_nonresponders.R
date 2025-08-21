################################################################################
# Angela Oill
# Linear model for resistence vs resolution accounting for biopsy grade
# Analysis: ACR pre-treatment responders vs non-responders
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


# CR Pre-treatment Resolved vs Persistence
# Select CR patients (patient_rejection_type)
# start with clustered_obj and patient_rejection_type
cellular_rejection_smpls <- clustered_obj@meta.data %>% 
  filter(patient_rejection_type == "cellular_rejection") %>%
  filter(biopsy_cellular_grading %in% c(0, 1, 2, 3)) %>% # to make sure only keeping samples with a CR grade
  filter(biopsy_timing == "pre") %>%
  filter(resolution_broad %in% c("resolved", "persistence")) %>% # there were some NAs this will exclude those
  dplyr::pull(Sample) %>% unique()

clustered_obj_pre <- subset(clustered_obj, subset = Sample %in% cellular_rejection_smpls)
cts_to_analyze <- levels(as.factor(clustered_obj_pre@meta.data$ct_second_pass))

# sample numbers
clustered_obj_pre@meta.data$Sample
clustered_obj_pre@meta.data$resolution_broad

tmp <- clustered_obj_pre@meta.data %>% dplyr::select(Sample, resolution_broad) %>% unique()
table(tmp$resolution_broad)

# get genes to analyze
# read in results file
# see: scRNA_heart_expression_comparison.ipynb
gene_prop_info_atlas <- read.csv("/scratch/aoill/projects/heart_transplant/scrna_atlas/genes_to_keep_proportions_new_groups.csv")

# gene genes with at least x percent expression for the different lineages
gene_prop_info_atlas_filter <- gene_prop_info_atlas %>%
  filter(Lineage != "Neural") %>%
  filter(Lineage != "Mesothelial") %>%
  filter(Proprtion_Expressed >= 0.01)

table(gene_prop_info_atlas_filter$Lineage)


# All cell types and genes ----
### Run test ----
lm_test_results <- c()

for (ctn in 1:length(cts_to_analyze)) {
  ct <- cts_to_analyze[ctn]
  print(ct)
  
  # subset object to cell type
  obj_CT <- subset(clustered_obj_pre, subset = ct_second_pass == ct)
  
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
  expr_matrix_t$Sample <- obj_CT@meta.data$Sample
  
  # Summarize by Sample to calculate the mean expression for each gene and sample
  mean_expression_per_sample <- expr_matrix_t %>%
    group_by(Sample) %>%
    summarise(across(everything(), mean, na.rm = TRUE))
  
  # Convert back to a data frame 
  mean_expression_per_sample <- as.data.frame(mean_expression_per_sample)
  # add sample info
  sample_info <- obj_CT@meta.data %>% dplyr::select(Sample, resolution_broad, biopsy_cellular_grading) %>% unique()
  mean_expression_per_sample_info <- left_join(mean_expression_per_sample, sample_info)
  
  
  if (min(table(mean_expression_per_sample_info$resolution_broad)) < 3 || (length(table(mean_expression_per_sample_info$resolution_broad)) < 2) ) {
    print(paste("Not enough samples or groups to do test. Cell type: ", ct, sep = ""))
    
    print(table(mean_expression_per_sample_info$resolution_broad))
    
    # add empyt list so this gets added to list of results instead of the previous cell types results
    lm_test_results_all_genes <- c()
    
  } else {
    
    lm_test_results_all_genes <- c()
    for (genei in genes_to_analyze) {
      #print(genei)
      # Expression filter
      #sum(genei_mean_expression_per_sample$expression > 0) >= 3
      genei_mean_expression_per_sample <- mean_expression_per_sample_info %>% dplyr::select(Sample, resolution_broad, biopsy_cellular_grading, genei)
      colnames(genei_mean_expression_per_sample) <- c("Sample", "resolution_broad", "biopsy_cellular_grading", "expression")
      
      genei_mean_expression_per_sample$resolution_broad <- factor(genei_mean_expression_per_sample$resolution_broad, levels = c("resolved", "persistence"))
      genei_mean_expression_per_sample$resolution_broad <- as.factor(genei_mean_expression_per_sample$resolution_broad)
      genei_mean_expression_per_sample$biopsy_cellular_grading <- as.numeric(genei_mean_expression_per_sample$biopsy_cellular_grading)
      
      
      if (sum(genei_mean_expression_per_sample$expression > 0) >= 3) {
        #genei_mean_expression_per_sample <- mean_expression_per_sample_info %>% dplyr::select(Sample, resolution_broad, genei)
        #colnames(genei_mean_expression_per_sample) <- c("Sample", "resolution_broad", "expression")
        
        # do lm
        model <- lm(expression ~ resolution_broad + biopsy_cellular_grading, data = genei_mean_expression_per_sample)
        model_summary <- summary(model)
        coefficients <- as.data.frame(model_summary$coefficients)
        
        coefficients
        
        # Add all results to larger dataframe/list
        response_coeff <- coefficients["resolution_broadpersistence", ]
        lm_result_genei <- data.frame(
          cell_type = ct,
          gene = genei,
          estimate = response_coeff["Estimate"],
          std_Error = response_coeff["Std. Error"],
          t_value = response_coeff["t value"],
          p_value = response_coeff["Pr(>|t|)"]
        )
        colnames(lm_result_genei) <- c("cell_type", "gene", "estimate", "std_Error", "t_value", "p_value")
        rownames(lm_result_genei) <- NULL
        
        lm_test_results_all_genes <- rbind(lm_test_results_all_genes, lm_result_genei)
        
        
      } else {
        print(paste("Not enough samples with expression in gene: ", genei))
      }
      
      
    }
    
    # Add multiple testing correction
    lm_test_results_all_genes$padj <- p.adjust(lm_test_results_all_genes$p_value, method = "fdr")
    
  }
  
  lm_test_results[[ct]] <- lm_test_results_all_genes
}


# Merge the list of data frames into one data frame
all_results_merged <- do.call(rbind, lm_test_results)
all_results_merged$estimate <- as.numeric(all_results_merged$estimate)
all_results_merged$p_value <- as.numeric(all_results_merged$p_value)
all_results_merged$padj <- as.numeric(all_results_merged$padj)

all_results_merged_sig <- all_results_merged %>% filter(padj <=0.05)
length(unique(all_results_merged_sig$cell_type))
length(unique(all_results_merged_sig$gene))
(nrow(all_results_merged_sig))
# 21 sig genes across 6 cell types, 23 gene cell type pairs



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
combined_results <- bind_rows(lm_test_results, .id = "ct")

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
out_table <- combined_results %>% dplyr::select(estimate, std_Error, t_value, p_value, padj, gene, ct, lineage, significance)
colnames(out_table) <- c("estimate", "std_Error", "t_value", "p_value", "padj", "gene", "cell_type", "lineage", "significance")
head(out_table)

write_csv(out_table, "/home/aoill/projects/heart_transplant/00_final/revisions/lm_ACR_resp_vs_non_resp_results.csv")

