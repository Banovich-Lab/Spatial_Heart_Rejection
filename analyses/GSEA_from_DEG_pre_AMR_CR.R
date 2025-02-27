################################################################################
# Spatial Heart Rejection Analysis code
# Analysis: GSEA from Differential gene expression pre-treatment AMR vs ACR
################################################################################
# I want to perform a GSEA on this. I will target a minsize of 5 genes, since 
# we're doing a targeted 477 gene panel. The 15-500 cutoff is for "defaults that 
# are appropriate for datasets with 10,000 to 20,000 features" 
# (https://www.gsea-msigdb.org/gsea/doc/GSEAUserGuideFrame.html).
# NES > 0 means that pathway is enriched in ACR, while NES < 0 means that 
# pathway is enriched in AMR.


# Load libraries ----
library(tidyverse)
library(ggpubr)
library(msigdbr)
library(fgsea)
library(forcats)

# Read in data ----
acr_v_amr <- read_csv("/home/aoill/projects/heart_transplant/00_final/DGE_CR_vs_AMR_fixed_exp_pre_only_smpl_fixed.csv")


# Load MSigDB gene sets
msigdb_hallmark <- msigdbr(species = "Homo sapiens", category = "H")
gene_sets <- split(msigdb_hallmark$gene_symbol, msigdb_hallmark$gs_name)

perform_gsea_by_cell_type <- function(df, gene_sets) {
  cell_type_list <- split(df, df$cell_type)
  
  gsea_results <- lapply(cell_type_list, function(cell_df) {
    cell_df <- cell_df %>% filter(!is.na(logFC) & !is.na(P.Value) & !is.na(gene))
    if (nrow(cell_df) < 5) {
      return(data.frame())
    }
    cell_df$GSEA_score <- sign(cell_df$logFC) * (-log10(cell_df$P.Value))
    ranked_genes <- setNames(cell_df$GSEA_score, cell_df$gene)
    ranked_genes <- ranked_genes[is.finite(ranked_genes)]
    gsea_res <- fgsea(pathways = gene_sets, stats = ranked_genes, minSize = 5, maxSize = 500, nperm = 1000)
    gsea_res <- as.data.frame(gsea_res) %>%
      mutate(cell_type = unique(cell_df$cell_type))
    
    return(gsea_res)
  })
  
  
  gsea_results_combined <- bind_rows(gsea_results)
  return(gsea_results_combined)
}

# Set seed
set.seed(2)

# Run GSEA on a per-cell-type basis
gsea_results.acr_v_amr <- perform_gsea_by_cell_type(acr_v_amr, gene_sets)


# Save csv
write_csv(gsea_results.acr_v_amr, "/home/aoill/projects/heart_transplant/00_final/20250225_gsea_ACR_v_AMR.csv")



