################################################################################
# Angela Oill
# Spatial Heart Rejection Analysis code
# Analysis: Cell-type proximity
# Code 2/2 for this analysis
################################################################################

# Load libraries ----
library(data.table)
library(dplyr)
library(circlize)
library(RColorBrewer)
library(ComplexHeatmap)

'%!in%' <- function(x,y)!('%in%'(x,y))
filter <- dplyr::filter
select <- dplyr::select
#### load data and define code paramters
#inDir <- '/scratch/aoill/projects/heart_transplant/prox_analysis/input/'
outDir <- '/scratch/aoill/projects/heart_transplant/prox_analysis_all_cell_types/output/prox_enrichment'
options(contrasts=c("contr.sum","contr.poly"))


# Read in object ----
clustered_obj <- readRDS("/scratch/aoill/projects/heart_transplant/GEO/heart_spatial_obj.rds")

## CR biopsies pre-treatment ----
conf_clustered_obj_pre <- subset(clustered_obj, subset = biopsy_timing == "pre" )
conf_clustered_obj_pre_CR <- subset(conf_clustered_obj_pre, subset = biopsy_rejection_type == "cellular_rejection")

# get metadata 
prox_analysis_meta <- conf_clustered_obj_pre_CR@meta.data

meta <- prox_analysis_meta
colnames(meta)
meta <- meta[, grep('^leiden', colnames(meta), invert = T)]

cell_types <- unique(conf_clustered_obj_pre_CR@meta.data$ct_second_pass)


enrich_score <- c()

for(i in cell_types){
  
  celltype_id <- i
  print(celltype_id)
  #cell_type_id_out <- gsub("\\/", "\\_", celltype_id)
  cell_type_id_out <- gsub(" ", "_", gsub("\\+", "", celltype_id))
  
  inpath <-  '/scratch/aoill/projects/heart_transplant/prox_analysis_all_cell_types/output/prox_files'
  f <- paste0(cell_type_id_out, '_cell_dist_adjdegree_CR_biopsies.rds')
  df <- readRDS(file.path(inpath, f))
  df <- df[order(df$V1, decreasing = F ),]

  df$celltypeA <- meta$ct_second_pass[match(df$cell.a, rownames(meta))]
  df$celltypeB <- meta$ct_second_pass[match(df$cell.b, rownames(meta))]
  
  # set cell type threshold
  d_threshold <- df %>% group_by(cell.a) %>%
    mutate(idx = row_number(cell.a)) %>%
    group_by(idx) %>%
    summarise(d_threshold = mean(V1) + 2*sd(V1)) %>%
    filter(idx == 1) %>%
    pull(d_threshold)
  d_threshold
  
  ## assign each cell to a range category
  degrees <- df$degree
  range <- seq(30, 360, 30) - 15
  
  a=data.table(degrees=degrees)
  a[,merge:=degrees]
  b=data.table(range=range)
  b[,merge:=range]
  setkeyv(a,c('merge'))
  setkeyv(b,c('merge'))
  Merge_a_b=b[a,roll='nearest']
  df$angle <- Merge_a_b$range[match(df$degree, Merge_a_b$degrees)]
  df <- df %>% group_by(cell.a, angle) %>%
    mutate(idx = row_number(cell.a))
  
  tmp0 <- df %>% 
    filter(V1 < d_threshold & idx == 1)
  
  tmp0 <- tmp0[!duplicated(tmp0),]
  tmp0 <- tmp0[,c('cell.b', 'celltypeB', 'sid')] ## rm V1, cell.a, degree, angle, idx
  tmp0 <- tmp0[!duplicated(tmp0),]
  
  # proximal cells by cell type
  proximal_cell_counts <- table(tmp0$celltypeB) %>% as.data.frame()
  proximal_cell_counts$ct_total = sum(proximal_cell_counts$Freq)
  
  all_cell_counts <- table(meta$ct_second_pass) %>% as.data.frame()
  all_cell_counts$cell_total = sum(all_cell_counts$Freq)
  
  tmp1 <- merge(proximal_cell_counts, all_cell_counts, by = 'Var1', all = T)
  colnames(tmp1) <- c('Celltype', 'proximal_ct', 'proximal_cell', 'total_ct', 'total_cell')
  tmp1[is.na(tmp1)] = 0
  enrich_score_ct <- c()
  for(j in 1:nrow(tmp1)){
    ct <- tmp1[j, 1]
    total <- tmp1[j, 5]
    list1 <- tmp1[j, 4]
    list2 <- tmp1[j, 3]
    overlap <- tmp1[j, 2]
    res <- fisher.test(matrix(c((total - (list1 + list2) + overlap),  (list1-overlap),  (list2 - overlap), overlap), nrow=2))
    tmp2 <- data.frame(prox_to = celltype_id, prox_ct = ct,
                       pval = res$p.value, or = res$estimate)
    enrich_score_ct <- rbind(enrich_score_ct, tmp2)
  }
  
  enrich_score <- rbind(enrich_score_ct, enrich_score)
}

enrich_score <- enrich_score %>% 
  group_by(prox_to) %>% 
  mutate(fdr = p.adjust(pval, method = 'fdr'))
enrich_score$logOR <- log(enrich_score$or)

# Save for plotting
saveRDS(enrich_score, "/home/aoill/projects/heart_transplant/00_final/cell_prox_CR_biopsies_pre.rds")
#enrich_score <- readRDS("/home/aoill/projects/heart_transplant/00_final/cell_prox_CR_biopsies_pre.rds")



## AMR biopsies pre-treatment ----
conf_clustered_obj_pre <- subset(clustered_obj, subset = biopsy_timing == "pre" )
conf_clustered_obj_pre_CR <- subset(conf_clustered_obj_pre, subset = biopsy_rejection_type == "antibody_rejection ")

# get metadata 
prox_analysis_meta <- conf_clustered_obj_pre_CR@meta.data


meta <- prox_analysis_meta
colnames(meta)
meta <- meta[, grep('^leiden', colnames(meta), invert = T)]

#cell_types <- names(table(meta$ct_endo_vcam1))
cell_types <- unique(conf_clustered_obj_pre@meta.data$ct_second_pass)

# Some cell types didn't have adiquate data so removing them from plotting
cell_types <- c("Endothelial", "CD8+ T cells", "NK", "Cardiomyocytes",
                "pDC", "CD4+ T cells", "Treg", "Macrophages", "B cells", 
                "Activated endothelial", "Fibroblasts", "Mast", "cDC1",
                "Proliferating T cells", "Plasma",
                "cDC2", "Pericytes", "POSTN+ Fibroblasts", "Adipocytes",                
                "Proliferating endothelial", "BMX+ Activated endothelial",
                "SPP1+ Macrophages", "vSMCs", "Myofibroblasts", 
                "Proliferating pericytes")

enrich_score <- c()


for(i in cell_types){
  
  celltype_id <- i
  print(celltype_id)
  #cell_type_id_out <- gsub("\\/", "\\_", celltype_id)
  cell_type_id_out <- gsub(" ", "_", gsub("\\+", "", celltype_id))
  
  inpath <-  '/scratch/aoill/projects/heart_transplant/prox_analysis_all_cell_types/output/prox_files'
  f <- paste0(cell_type_id_out, '_cell_dist_adjdegree_AMR_biopsies.rds')
  df <- readRDS(file.path(inpath, f))
  df <- df[order(df$V1, decreasing = F ),]

  
  df$celltypeA <- meta$ct_second_pass[match(df$cell.a, rownames(meta))]
  df$celltypeB <- meta$ct_second_pass[match(df$cell.b, rownames(meta))]
  
  # set cell type threshold
  d_threshold <- df %>% group_by(cell.a) %>%
    mutate(idx = row_number(cell.a)) %>%
    group_by(idx) %>%
    summarise(d_threshold = mean(V1) + 2*sd(V1)) %>%
    filter(idx == 1) %>%
    pull(d_threshold)
  d_threshold
  
  ## assign each cell to a range category
  degrees <- df$degree
  range <- seq(30, 360, 30) - 15
  
  a=data.table(degrees=degrees)
  a[,merge:=degrees]
  b=data.table(range=range)
  b[,merge:=range]
  setkeyv(a,c('merge'))
  setkeyv(b,c('merge'))
  Merge_a_b=b[a,roll='nearest']
  df$angle <- Merge_a_b$range[match(df$degree, Merge_a_b$degrees)]
  df <- df %>% group_by(cell.a, angle) %>%
    mutate(idx = row_number(cell.a))
  
  tmp0 <- df %>% 
    filter(V1 < d_threshold & idx == 1)
  
  tmp0 <- tmp0[!duplicated(tmp0),]
  tmp0 <- tmp0[,c('cell.b', 'celltypeB', 'sid')] ## rm V1, cell.a, degree, angle, idx
  tmp0 <- tmp0[!duplicated(tmp0),]
  
  # proximal cells by cell type
  proximal_cell_counts <- table(tmp0$celltypeB) %>% as.data.frame()
  proximal_cell_counts$ct_total = sum(proximal_cell_counts$Freq)
  
  all_cell_counts <- table(meta$ct_second_pass) %>% as.data.frame()
  all_cell_counts$cell_total = sum(all_cell_counts$Freq)
  
  tmp1 <- merge(proximal_cell_counts, all_cell_counts, by = 'Var1', all = T)
  colnames(tmp1) <- c('Celltype', 'proximal_ct', 'proximal_cell', 'total_ct', 'total_cell')
  tmp1[is.na(tmp1)] = 0
  enrich_score_ct <- c()
  for(j in 1:nrow(tmp1)){
    ct <- tmp1[j, 1]
    total <- tmp1[j, 5]
    list1 <- tmp1[j, 4]
    list2 <- tmp1[j, 3]
    overlap <- tmp1[j, 2]
    res <- fisher.test(matrix(c((total - (list1 + list2) + overlap),  (list1-overlap),  (list2 - overlap), overlap), nrow=2))
    tmp2 <- data.frame(prox_to = celltype_id, prox_ct = ct,
                       pval = res$p.value, or = res$estimate)
    enrich_score_ct <- rbind(enrich_score_ct, tmp2)
  }
  
  enrich_score <- rbind(enrich_score_ct, enrich_score)
}

enrich_score <- enrich_score %>% 
  group_by(prox_to) %>% 
  mutate(fdr = p.adjust(pval, method = 'fdr'))
enrich_score$logOR <- log(enrich_score$or)


# Save for plotting
saveRDS(enrich_score, "/home/aoill/projects/heart_transplant/00_final/cell_prox_AMR_biopsies_pre.rds")
#enrich_score <- readRDS("/home/aoill/projects/heart_transplant/00_final/cell_prox_AMR_biopsies_pre.rds")

