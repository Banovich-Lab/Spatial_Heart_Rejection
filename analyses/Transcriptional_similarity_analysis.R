################################################################################
# Angela Oill
# Spatial Heart Rejection Analysis code
# Analysis: Transcriptional similarity analysis 
################################################################################
# Nearest neighbor analysis to compare transcriptional similarity among biopsies

# Load libraries ----
library(Seurat)
library(circlize)
library(RColorBrewer)
library(ComplexHeatmap)

# Load functions ----
check_all_infinite <- function(matrix) {
  
  for (i in 1:nrow(matrix)) {
    
    for (j in 1:ncol(matrix)) {
      
      if (!is.infinite(matrix[i, j])) {
        
        return(FALSE)  
        
      }
      
    }
    
  }
  
  return(TRUE)
  
} 
 

# Load object ----
clustered_obj <- readRDS("/scratch/aoill/projects/heart_transplant/GEO/heart_spatial_obj.rds")


# Merged results - CR, pre-treatment, by biopsy level grading ----
# Save fe_results for all cell types, make a col called cell_cell_type and that 
# will be on the y axis of the heat map and do not do any clustering
cts_to_analyze <- levels(as.factor(clustered_obj@meta.data$ct_second_pass))

fe_results_all <- c()

for (ctn in 1:length(cts_to_analyze)) {
  ct <- cts_to_analyze[ctn]
  print(ct)
  
  
  DefaultAssay(clustered_obj) <- "RNA"
  obj_CT <- subset(clustered_obj, subset = ct_second_pass == ct)
  # Need to select CR patients so all grade biopsies are selected 
  obj_CT_CR <- subset(obj_CT, subset = patient_rejection_type == "cellular_rejection") 
  obj_CT_CR_pre <- subset(obj_CT_CR, subset = biopsy_timing == "pre")
  
  # Calculate nearest neighbors
  gene_features <- rownames(obj_CT_CR_pre)
  
  # Perform PCA 
  obj_CT_CR_pre <- RunPCA(obj_CT_CR_pre, npcs = 20, features = gene_features)  
  
  # Find nearest neighbors using the first 20 PCs 
  obj_CT_CR_pre <- FindNeighbors(obj_CT_CR_pre, dims = 1:20)
  
  
  # Extract the non-zero elements of the sparse matrix
  nn_graph <- obj_CT_CR_pre@graphs[["RNA_nn"]]
  
  # Excluding self as nearest neighbor
  # Initialize results
  results <- data.frame(
    CellID = rownames(nn_graph),
    ClosestNeighborID = character(nrow(nn_graph)),
    Distance = numeric(nrow(nn_graph))
  )
  
  # Iterate over rows in the sparse matrix
  for (i in seq_len(nrow(nn_graph))) {
    #print(i)
    # Extract the row as a sparse vector
    row <- nn_graph[i, ]
    
    # Remove self-neighbor (set to -Inf or NA for exclusion)
    row[i] <- -Inf  # Set the self-distance to a minimum value
    
    # Find the index and value of the maximum non-self element
    max_index <- which.max(row)
    max_distance <- row[max_index]
    
    # Record the results
    results$ClosestNeighborID[i] <- colnames(nn_graph)[max_index]
    results$Distance[i] <- max_distance
  }
  
  #length(unique(results$ClosestNeighborID))
  
  colnames(results) <- c("X", "Nearest_Neighbor_ID", "Distance")
  #length(unique(results$Nearest_Neighbor_ID))
  
  # Join results with metadata
  meta_cobj_CT_CR_pre <- obj_CT_CR_pre@meta.data
  
  meta_cobj_CT_CR_pre_join <- left_join(meta_cobj_CT_CR_pre, results)
  
  
  # Get nearest neighbor's info into meta
  cell_info <- meta_cobj_CT_CR_pre_join %>% 
    dplyr::select(X, biopsy_timing, biopsy_rejection_type, patient_rejection_type,
                  biopsy_cellular_grading, patient_cellular_grading)
  # relabel col names so we know this will be the nearest neighbor results
  colnames(cell_info) <- c("Nearest_Neighbor_ID", "Nearest_Neighbor_biopsy_timing", 
                           "Nearest_Neighbor_rejection_type", "Nearest_Neighbor_OLD_rejection_type",
                           "Nearest_Neighbor_cellular_grading", "Nearest_Neighbor_OLD_cellular_grading")
  
  meta_cobj_CT_CR_pre_join_neighbor_info <- left_join(meta_cobj_CT_CR_pre_join, cell_info)
  
  cell_neighbor_biopsy_level_info <- meta_cobj_CT_CR_pre_join_neighbor_info %>% 
    dplyr::select(biopsy_cellular_grading, Nearest_Neighbor_cellular_grading)
  #head(cell_neighbor_biopsy_level_info)
  
  
  # Remove NAs
  cell_neighbor_biopsy_level_info <- na.omit(cell_neighbor_biopsy_level_info)
  #nrow(cell_neighbor_biopsy_level_info)
  
  
  fe_results <- c()
  for (i in 0:3) {
    #print(i)
    for (j in 0:3) {
      #print(j)
      
      # a is cell is grade i, neighbor is grade j
      a <- cell_neighbor_biopsy_level_info %>% 
        filter(biopsy_cellular_grading == i) %>% 
        filter(Nearest_Neighbor_cellular_grading == j) %>% nrow()
      
      # b is cell is grade i, neighbor is NOT grade j
      b <- cell_neighbor_biopsy_level_info %>% 
        filter(biopsy_cellular_grading == i) %>% 
        filter(Nearest_Neighbor_cellular_grading != j) %>% nrow()
      
      # c is cell is NOT grade i, neighbor is grade j
      c <- cell_neighbor_biopsy_level_info %>% 
        filter(biopsy_cellular_grading != i) %>% 
        filter(Nearest_Neighbor_cellular_grading == j) %>% nrow()
      
      # d is cell is NOT grade i, neighbor is NOT grade j
      d <- cell_neighbor_biopsy_level_info %>% 
        filter(biopsy_cellular_grading != i) %>% 
        filter(Nearest_Neighbor_cellular_grading != j) %>% nrow()
      
      (matrix(c(a, c, b, d), nrow=2))
      
      tmp <- fisher.test(matrix(c(a, c, b, d), nrow=2))
      
      tmp2 <- data.frame(cell_type = ct, cell = i, cell_cell_type = paste(i, ct, sep = "_"), 
                         neighbor = j,
                         pval = tmp$p.value, or = tmp$estimate)
      fe_results <- rbind(fe_results, tmp2)
      
    }
    
  }
  
  
  # adjust p value and add log OR
  fe_results <- fe_results %>% 
    group_by(cell) %>% 
    mutate(fdr = p.adjust(pval, method = 'fdr'))
  fe_results$logOR <- log(fe_results$or)
  
  
  fe_results_all <- rbind(fe_results_all, fe_results)
  
}

# Save fe_results_all
saveRDS(fe_results_all, "/scratch/aoill/projects/heart_transplant/transcriptional_similarity_results.rds")


