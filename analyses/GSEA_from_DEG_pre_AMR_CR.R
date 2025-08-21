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
acr_v_amr <- read_csv("/home/aoill/projects/heart_transplant/00_final/revisions/t_test_AMR_vs_CR_results.csv")


# Load MSigDB gene sets
msigdb_hallmark <- msigdbr(species = "Homo sapiens", category = "H")
gene_sets <- split(msigdb_hallmark$gene_symbol, msigdb_hallmark$gs_name)

perform_gsea_by_cell_type <- function(df, gene_sets) {
  cell_type_list <- split(df, df$cell_type)
  
  gsea_results <- lapply(cell_type_list, function(cell_df) {
    cell_df <- cell_df %>% filter(!is.na(Log2FC) & !is.na(p_value) & !is.na(gene))
    if (nrow(cell_df) < 5) {
      return(data.frame())
    }
    cell_df$GSEA_score <- sign(cell_df$Log2FC) * (-log10(cell_df$p_value))
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

# Remove the HALLMARK_ portion of the pathway
gsea_results.acr_v_amr <- gsea_results.acr_v_amr %>%
  mutate(pathway = str_remove(pathway, "^HALLMARK_"))

# Add significance stars for pathways with padj < 0.05
gsea_results.acr_v_amr <- gsea_results.acr_v_amr %>%
  mutate(significance = ifelse(padj < 0.05, "***", ""))
table(gsea_results.acr_v_amr$pathway)

# Rename pathways to make them easier to read
gsea_results.acr_v_amr <- gsea_results.acr_v_amr %>%
  mutate(
    pathway = case_when(
      pathway == "EPITHELIAL_MESENCHYMAL_TRANSITION" ~ "EMT",
      pathway == "APICAL_JUNCTION" ~ "Apical junction",
      pathway == "KRAS_SIGNALING_UP" ~ "KRAS signaling up",
      pathway == "IL2_STAT5_SIGNALING" ~ "IL2-STAT5 signaling",
      pathway == "INFLAMMATORY_RESPONSE" ~ "Inflammatory response",
      pathway == "TNFA_SIGNALING_VIA_NFKB" ~ "TNF-alpha signaling via NFKβ",
      pathway == "COMPLEMENT" ~ "Complement",
      pathway == "INTERFERON_GAMMA_RESPONSE" ~ "IFNy-gamma response",
      pathway == "ALLOGRAFT_REJECTION" ~ "Allograft rejection",
      pathway == "ADIPOGENESIS" ~ "Adipogenesis",
      pathway == "ANDROGEN_RESPONSE" ~ "Androgen response",
      pathway == "ANGIOGENESIS" ~ "Angiogenesis",
      pathway == "APICAL_SURFACE" ~ "Apical surface",
      pathway == "APOPTOSIS" ~ "Apoptosis",
      pathway == "COAGULATION" ~ "Coagulation",
      pathway == "E2F_TARGETS" ~ "E2F targets",
      pathway == "ESTROGEN_RESPONSE_EARLY" ~ "Estrogen response (early)",
      pathway == "ESTROGEN_RESPONSE_LATE" ~ "Estrogen response (late)",
      pathway == "FATTY_ACID_METABOLISM" ~ "Fatty acid metabolism",
      pathway == "G2M_CHECKPOINT" ~ "G2M checkpoint",
      pathway == "GLYCOLYSIS" ~ "Glycolysis",
      pathway == "HEME_METABOLISM" ~ "Heme metabolism",
      pathway == "HYPOXIA" ~ "Hypoxia",
      pathway == "IL6_JAK_STAT3_SIGNALING" ~ "IL6-JAK-STAT3 signaling",
      pathway == "INTERFERON_ALPHA_RESPONSE" ~ "IFN-alpha response",
      pathway == "KRAS_SIGNALING_DN" ~ "KRAS signaling down",
      pathway == "MITOTIC_SPINDLE" ~ "Mitotic spindle",
      pathway == "MYOGENESIS" ~ "Myogenesis",
      pathway == "P53_PATHWAY" ~ "p53 pathway",
      pathway == "PANCREAS_BETA_CELLS" ~ "Pancreas β cells",
      pathway == "UV_RESPONSE_DN" ~ "UV response (down)",
      pathway == "UV_RESPONSE_UP" ~ "UV response (up)",
      pathway == "XENOBIOTIC_METABOLISM" ~ "Xenobiotic metabolism",
      TRUE ~ pathway
    )
  )

# View them for easier writing
gsea_results.acr_v_amr %>%
  filter(padj < 0.05) %>%
  group_by(pathway) %>%
  print()

# Save csv
write_csv(gsea_results.acr_v_amr, "/home/aoill/projects/heart_transplant/00_final/20250805_gsea_ACR_v_AMR.csv")
