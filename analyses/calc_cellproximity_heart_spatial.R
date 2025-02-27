################################################################################
# Angela Oill
# Spatial Heart Rejection Analysis code
# Analysis: Cell-type proximity
# Code 1/2 for this analysis
################################################################################

# Functions and libraries ----
library(sf)
library(dplyr)

filter <- dplyr::filter
select <- dplyr::select

## load functions 
circleFun <- function(center = c(0,0),diameter = 1, npoints = 100){
  r = diameter / 2
  tt <- seq(0,2*pi,length.out = npoints)
  xx <- center[1] + r * cos(tt)
  yy <- center[2] + r * sin(tt)
  return(data.frame(x = xx, y = yy))
}

## calculate distance function
calc_d <- function(pt1, pt2) {
  d <- sqrt( (( pt2[2] - pt1[2] )^2 ) + (( pt2[1] - pt1[1] )^2 ) )
  return(d)
}

getDir <- function(origin, pt2){
  pt2 <- pt2 - origin
  deg <- abs(atan(pt2[2]/pt2[1]) * (180/pi))
  if(pt2[1] > 0 & pt2[2] > 0){deg <- deg
  } else if(pt2[1] < 0 & pt2[2] > 0){ deg <- 180-deg
  } else if(pt2[1] < 0 & pt2[2] < 0){ deg <- 180 + deg
  } else if(pt2[1] > 0 & pt2[2] < 0){ deg <- 360 - deg}
  return(deg)
}

assign.angle <- function(degrees, range){
  require(data.table)
  a=data.table(degrees=degrees)
  a[,merge:=degrees]
  
  b=data.table(range=range)
  b[,merge:=range]
  
  setkeyv(a,c('merge'))
  setkeyv(b,c('merge'))
  Merge_a_b=b[a,roll='nearest']
  return(Merge_a_b)
}

anchorPoint <- function(degrees, anchor){
  object <- (degrees + (360 - anchor)) 
  object <- ifelse(object > 360, object - 360, object)
  return(object)
}


# Read in object ----
clustered_obj <- readRDS("/scratch/aoill/projects/heart_transplant/GEO/heart_spatial_obj.rds")

# Analysis ----
## CR biopsies pre-treatment ----
conf_clustered_obj_pre <- subset(clustered_obj, subset = biopsy_timing == "pre" )
conf_clustered_obj_pre_CR <- subset(conf_clustered_obj_pre, subset = biopsy_rejection_type == "cellular_rejection")

# get metadata 
prox_analysis_meta <- conf_clustered_obj_pre_CR@meta.data

# Get cell types and proximal cells. This runs for each cell type 
cts_to_analyze <- levels(as.factor(conf_clustered_obj_pre_CR@meta.data$ct_second_pass))


for (ctn in 1:length(cts_to_analyze)) {
  ct <- cts_to_analyze[ctn]
  print(ct) 
  
  
  # Do for all cells
  cell_type_id <- ct
  # remove spaces with underscore and + with nothing
  cell_type_id_out <- gsub(" ", "_", gsub("\\+", "", cell_type_id))
  
  print(cell_type_id)
  print(cell_type_id_out)
  
  
  
  ## load data and define code paramters
  inDir <- '/scratch/aoill/projects/heart_transplant/prox_analysis_all_cell_types/input/'
  outDir <- '/scratch/aoill/projects/heart_transplant/prox_analysis_all_cell_types/output/'
  
  x <- prox_analysis_meta
  #x[1:5,1:5]
  
  # I think this would be every biopsy (can try also doing by patient)
  sample_ids <- x %>%
    pull(Sample) %>% as.character %>% unique() 
  
  nCores <- 24
  library(doParallel)
  registerDoParallel(cores=nCores)
  proximal_cell_pop <- c()
  for (i in 1:length(sample_ids)) {
    sid <- sample_ids[i]; print(sid)
    obj.sid <- subset(x, Sample == sid)
    
    # 1. get cell coordinates
    meta <- obj.sid
    cells_coord <- meta %>% 
      filter(ct_second_pass == cell_type_id)
    print(nrow(cells_coord))
    if(nrow(cells_coord) <= 1){
      print("skip is Sample does not have the right cells")
    } else {
      # 2. define radius
      cell_counts <- c()
      r <- 60 ## set radius ~30 is 3 cells 
      print(r)
      for(j in 1:nrow(cells_coord)){
        #print(j)
        cell_id <- rownames(cells_coord)[j]
        split_cell.id <- strsplit(cell_id, ";") %>% unlist()
        cell_x <- cells_coord$x_centroid[j]
        cell_y <- cells_coord$y_centroid[j]
        
        circ_coord <- circleFun(c(cell_x, cell_y), r, 1000)
        circ_coord <- circ_coord[!duplicated(circ_coord),]
        
        # 3. get cell composition within each ring
        ## 3a. subset meta data to include all points with ring boundary 
        all_cell_coord <- data.frame(x = meta$x_centroid, y = meta$y_centroid)
        tmp <- circ_coord[,c(1:2)]
        io <- splancs::inout(all_cell_coord, tmp)
        
        if(mean(io) == 0){next}
        extract_cells <- as.data.frame(meta)[which(io == TRUE),] ## get cell within radius
        if(nrow(extract_cells) == 1){next}
        
        # calculate distance and degrees
        ## Calculate distance between cells
        tmp <- extract_cells[,c('x_centroid', 'y_centroid')]
        tmp <- rbind(tmp[which(rownames(tmp) %in% cell_id),], tmp)
        tmp <- tmp[!duplicated(tmp),]
        
        toCheck <- combn(rownames(tmp), 2, simplify = FALSE)
        names(toCheck) <-
          sapply(toCheck, paste, collapse = " - ")
        toCheck <- toCheck[grep(pattern = cell_id, toCheck)]
        
        ## calculate distance
        tmp.d <- sapply(toCheck, function(j){
          calc_d(tmp[j[1],c(1,2)], tmp[j[2],c(1,2)]) })
        names(tmp.d) <- gsub('\\.y_centroid', '', names(tmp.d))
        tmp.d <- as.data.frame(do.call(rbind, tmp.d))
        tmp.d$cell_id <- rownames(tmp.d)
        tmp.d <- tidyr::separate(tmp.d, cell_id, sep = "\\ - ", c('cell.a', 'cell.b'))
        rownames(tmp.d) <- NULL
        
        ## calculate degrees
        tmp.degree <- sapply(toCheck, function(j){
          getDir(c(cell_x, cell_y), tmp[j[2],c(1,2)]) })
        names(tmp.degree) <- gsub('\\.y_centroid', '', names(tmp.degree))
        tmp.degree <- as.data.frame(do.call(rbind, tmp.degree))
        tmp.degree$cell_id <- rownames(tmp.degree)
        tmp.degree <- tidyr::separate(tmp.degree, cell_id, sep = "\\ - ", c('cell.a', 'cell.b'))
        rownames(tmp.degree) <- NULL
        
        ## add degrees to distance table
        tmp.d$degree <- tmp.degree$V1[match(tmp.d$cell.b, tmp.degree$cell.b)]
        tmp.d$celltypeA <- extract_cells$ct_second_pass[match(tmp.d$cell.a, rownames(extract_cells))]
        tmp.d$celltypeB <- extract_cells$ct_second_pass[match(tmp.d$cell.b, rownames(extract_cells))]
        
        
        ## assign each cell to a range category - this will be calculated in the logistic regression code
        #range <- seq(30, 360, 30) - 15
        range <- seq(60, 360, 60) -30
        Merge_a_b <- assign.angle(degrees = tmp.d$degree, range = range)
        tmp.d$angle <- Merge_a_b$range[match(tmp.d$degree, Merge_a_b$degrees)]
        tmp.d$cell.a <- cell_id
        tmp.d <- tmp.d[!duplicated(tmp.d),]
        # append
        cell_counts <- rbind(cell_counts, tmp.d)
      }
      cell_counts$sid <- sid
      
      # print(paste("Ncol cell_counts ", ncol(cell_counts), " sample ", i, sep = ""))
      
      proximal_cell_pop <- rbind(cell_counts, proximal_cell_pop)
      #return(cell_counts)
    }
  }
  
  stopImplicitCluster() 
  out_path <- '/scratch/aoill/projects/heart_transplant/prox_analysis_all_cell_types/output/prox_files'
  
  out_fname <- paste0(cell_type_id_out, '_cell_dist_adjdegree_CR_biopsies.rds')
  saveRDS(proximal_cell_pop, file.path(out_path, out_fname))  
  
  
  
  
  
}


## AMR biopsies pre-treatment ----
conf_clustered_obj_pre <- subset(clustered_obj, subset = biopsy_timing == "pre" )
conf_clustered_obj_pre_CR <- subset(conf_clustered_obj_pre, subset = biopsy_rejection_type == "antibody_rejection ")

# get metadata 
prox_analysis_meta <- conf_clustered_obj_pre_CR@meta.data

# Get cell types and proximal cells
cts_to_analyze <- levels(as.factor(conf_clustered_obj_pre_CR@meta.data$ct_second_pass))


for (ctn in 1:length(cts_to_analyze)) {
  ct <- cts_to_analyze[ctn]
  print(ct) 
  
  
  # Do for all cells
  cell_type_id <- ct
  # remove spaces with underscore and + with nothing
  cell_type_id_out <- gsub(" ", "_", gsub("\\+", "", cell_type_id))
  
  print(cell_type_id)
  print(cell_type_id_out)
  
  
  
  ## load data and define code paramters
  inDir <- '/scratch/aoill/projects/heart_transplant/prox_analysis_all_cell_types/input/'
  outDir <- '/scratch/aoill/projects/heart_transplant/prox_analysis_all_cell_types/output/'
  
  x <- prox_analysis_meta
  #x[1:5,1:5]
  
  # I think this would be every biopsy (can try also doing by patient)
  sample_ids <- x %>%
    pull(Sample) %>% as.character %>% unique() 
  
  nCores <- 24
  library(doParallel)
  registerDoParallel(cores=nCores)
  proximal_cell_pop <- c()
  for (i in 1:length(sample_ids)) {
    sid <- sample_ids[i]; print(sid)
    obj.sid <- subset(x, Sample == sid)
    
    # 1. get cell coordinates
    meta <- obj.sid
    cells_coord <- meta %>% 
      filter(ct_second_pass == cell_type_id)
    print(nrow(cells_coord))
    if(nrow(cells_coord) <= 1){
      print("skip is Sample does not have the right cells")
    } else {
      # 2. define radius
      cell_counts <- c()
      r <- 60 ## set radius ~30 is 3 cells 
      print(r)
      for(j in 1:nrow(cells_coord)){
        #print(j)
        cell_id <- rownames(cells_coord)[j]
        split_cell.id <- strsplit(cell_id, ";") %>% unlist()
        cell_x <- cells_coord$x_centroid[j]
        cell_y <- cells_coord$y_centroid[j]
        
        circ_coord <- circleFun(c(cell_x, cell_y), r, 1000)
        circ_coord <- circ_coord[!duplicated(circ_coord),]
        
        # 3. get cell composition within each ring
        ## 3a. subset meta data to include all points with ring boundary 
        all_cell_coord <- data.frame(x = meta$x_centroid, y = meta$y_centroid)
        tmp <- circ_coord[,c(1:2)]
        io <- splancs::inout(all_cell_coord, tmp)
        
        if(mean(io) == 0){next}
        extract_cells <- as.data.frame(meta)[which(io == TRUE),] ## get cell within radius
        if(nrow(extract_cells) == 1){next}
        
        # calculate distance and degrees
        ## Calculate distance between cells
        tmp <- extract_cells[,c('x_centroid', 'y_centroid')]
        tmp <- rbind(tmp[which(rownames(tmp) %in% cell_id),], tmp)
        tmp <- tmp[!duplicated(tmp),]
        
        toCheck <- combn(rownames(tmp), 2, simplify = FALSE)
        names(toCheck) <-
          sapply(toCheck, paste, collapse = " - ")
        toCheck <- toCheck[grep(pattern = cell_id, toCheck)]
        
        ## calculate distance
        tmp.d <- sapply(toCheck, function(j){
          calc_d(tmp[j[1],c(1,2)], tmp[j[2],c(1,2)]) })
        names(tmp.d) <- gsub('\\.y_centroid', '', names(tmp.d))
        tmp.d <- as.data.frame(do.call(rbind, tmp.d))
        tmp.d$cell_id <- rownames(tmp.d)
        tmp.d <- tidyr::separate(tmp.d, cell_id, sep = "\\ - ", c('cell.a', 'cell.b'))
        rownames(tmp.d) <- NULL
        
        ## calculate degrees
        tmp.degree <- sapply(toCheck, function(j){
          getDir(c(cell_x, cell_y), tmp[j[2],c(1,2)]) })
        names(tmp.degree) <- gsub('\\.y_centroid', '', names(tmp.degree))
        tmp.degree <- as.data.frame(do.call(rbind, tmp.degree))
        tmp.degree$cell_id <- rownames(tmp.degree)
        tmp.degree <- tidyr::separate(tmp.degree, cell_id, sep = "\\ - ", c('cell.a', 'cell.b'))
        rownames(tmp.degree) <- NULL
        
        ## add degrees to distance table
        tmp.d$degree <- tmp.degree$V1[match(tmp.d$cell.b, tmp.degree$cell.b)]
        tmp.d$celltypeA <- extract_cells$ct_second_pass[match(tmp.d$cell.a, rownames(extract_cells))]
        tmp.d$celltypeB <- extract_cells$ct_second_pass[match(tmp.d$cell.b, rownames(extract_cells))]
        
        
        ## assign each cell to a range category - this will be calculated in the logistic regression code
        #range <- seq(30, 360, 30) - 15
        range <- seq(60, 360, 60) -30
        Merge_a_b <- assign.angle(degrees = tmp.d$degree, range = range)
        tmp.d$angle <- Merge_a_b$range[match(tmp.d$degree, Merge_a_b$degrees)]
        tmp.d$cell.a <- cell_id
        tmp.d <- tmp.d[!duplicated(tmp.d),]
        # append
        cell_counts <- rbind(cell_counts, tmp.d)
      }
      cell_counts$sid <- sid
      
      # print(paste("Ncol cell_counts ", ncol(cell_counts), " sample ", i, sep = ""))
      
      proximal_cell_pop <- rbind(cell_counts, proximal_cell_pop)
      #return(cell_counts)
    }
  }
  
  stopImplicitCluster() 
  out_path <- '/scratch/aoill/projects/heart_transplant/prox_analysis_all_cell_types/output/prox_files'
  
  out_fname <- paste0(cell_type_id_out, '_cell_dist_adjdegree_AMR_biopsies.rds')
  saveRDS(proximal_cell_pop, file.path(out_path, out_fname))  
  
  
  
  
  
}

