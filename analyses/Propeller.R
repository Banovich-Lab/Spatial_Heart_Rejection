################################################################################
# Spatial Heart Rejection Analysis code
# Analysis: Cell type proportion differences test with Propeller
################################################################################

# Load libraries ----
library(speckle)
library(Seurat)

# Read in object ---- 
clustered_obj <- readRDS("/scratch/aoill/projects/heart_transplant/GEO/heart_spatial_obj.rds")


# Subset out pre-treatment biopsies ----
filtered_smpls <- clustered_obj@meta.data %>% 
  filter(biopsy_timing == "pre") %>%
  dplyr::pull(Sample) %>% unique()

clustered_obj_filtered_smpls <- subset(clustered_obj, 
                                       subset =  Sample %in% filtered_smpls)

# Run propeller ----
propeller_results <- propeller(
  clusters = clustered_obj_filtered_smpls$ct_second_pass, 
  sample = clustered_obj_filtered_smpls$Sample, 
  group = clustered_obj_filtered_smpls$rejection_type)


# Save results ----
write.csv(propeller_results, 
          "/home/aoill/projects/heart_transplant/propeller_pre-treatment_by_biopsy_rejection_type_correct.csv",
          row.names = T, quote = F)
