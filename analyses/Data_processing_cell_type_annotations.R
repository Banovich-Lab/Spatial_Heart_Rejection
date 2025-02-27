################################################################################
# Angela Oill
# Heart transplant data processing and cell type annotation
################################################################################

# Load libraries ----
library(Seurat)
library(SeuratObject)
library(SeuratDisk)
library(tidyverse)
library(gplots)
library(tibble)
library(ggpubr)
library(ggrepel)
library(RColorBrewer)
library(presto)


# Load in data ----
# Set path where Xenium data is located on the cluster
path <- "/tgen_labs/banovich/xenium_run_folders/20240517__213503__KA_HrtRej032024/"
rb1 <- "output-XETG00048__0022213__KA_032024A__20240517__214455"
rb2 <- "output-XETG00048__0022215__KA_032024B__20240517__214455"
# This loads in the high quality transcripts at the cell segmented level. Cells
# are segmented using three methods and which method was used for each cell is
# stored somewhere that I will need to find
# 477 genes
xenium_obj_1 <- LoadXenium(paste(path, rb1, sep = ""), fov = "fov")
# total number of cells: 97262
xenium_obj_2 <- LoadXenium(paste(path, rb2, sep = ""), fov = "fov")
# total number of cells: 94842

## Load in relevant metadata ---- 
### RB1 ----
# Read in cell metadata
cell_meta <- read.csv(paste(path, rb1, "/cells.csv.gz", sep = ""),
                      header = T)
# Will want to join on cell id (cell_id) for this file

# Which cells correspond to which sample
path2 <- "/scratch/aoill/projects/heart_transplant/RB1/"
# This will need to be done in a for loop because each sample has their own file
# ls *_cells_stats.csv | cut -f1,2 -d'_' | awk '{ print "\""$1"\"\," }'
smpls <- smpl_list_RB1

smpl_cell_ids_all <- c()
for (smpl in smpls) {
  print(smpl)
  smpl_cell_ids <- read.csv(paste(path2, smpl, "_cells_stats.csv", sep = ""),
                            header = T, comment.char = "#")
  # add sample 
  smpl_cell_ids$Sample <- smpl
  smpl_cell_ids_sub <- smpl_cell_ids %>% 
    dplyr::select(Cell.ID, Sample)
  smpl_cell_ids_all <- rbind(smpl_cell_ids_all, smpl_cell_ids_sub)
}


# Join the two files
cell_meta_merge_rb1 <- merge(cell_meta, smpl_cell_ids_all, 
                             by.x = "cell_id",
                             by.y = "Cell.ID")


### RB2 ----
# Read in cell metadata
cell_meta <- read.csv(paste(path, rb2, "/cells.csv.gz", sep = ""),
                      header = T)
# Will want to join on cell id (cell_id) for this file

# Which cells correspond to which sample
path2 <- "/scratch/aoill/projects/heart_transplant/RB2/"
# This will need to be done in a for loop because each sample has their own file
# ls *_cells_stats.csv | cut -f1,2 -d'_' | awk '{ print "\""$1"\"\," }'
smpls <- smpl_list_RB2

smpl_cell_ids_all <- c()
for (smpl in smpls) {
  print(smpl)
  smpl_cell_ids <- read.csv(paste(path2, smpl, "_cells_stats.csv", sep = ""),
                            header = T, comment.char = "#")
  # add sample 
  smpl_cell_ids$Sample <- smpl
  smpl_cell_ids_sub <- smpl_cell_ids %>% 
    dplyr::select(Cell.ID, Sample)
  smpl_cell_ids_all <- rbind(smpl_cell_ids_all, smpl_cell_ids_sub)
}


# Join the two files
cell_meta_merge_rb2 <- merge(cell_meta, smpl_cell_ids_all, 
                             by.x = "cell_id",
                             by.y = "Cell.ID")



# Filter out cells that were not assigned to a sample ----
## RB1 ----
# Since some of the samples overlap they will not be included
xenium_obj_1@meta.data$cell_id <- rownames(xenium_obj_1@meta.data)
cells_to_keep <- cell_meta_merge_rb1$cell_id
xenium_obj_1_sub <- subset(xenium_obj_1, subset = cell_id %in% cells_to_keep)

## RB2 ----
xenium_obj_2@meta.data$cell_id <- rownames(xenium_obj_2@meta.data)
cells_to_keep <- cell_meta_merge_rb2$cell_id
xenium_obj_2_sub <- subset(xenium_obj_2, subset = cell_id %in% cells_to_keep)

# Add metadata to filtered object ----
## RB1 ----
xenium_obj_1_sub_meta <- xenium_obj_1_sub@meta.data
xenium_obj_1_sub_meta_cell_meta <- full_join(xenium_obj_1_sub_meta, cell_meta_merge_rb1)
rownames(xenium_obj_1_sub_meta_cell_meta) <- rownames(xenium_obj_1_sub@meta.data)
xenium_obj_1_sub@meta.data <- xenium_obj_1_sub_meta_cell_meta
xenium_obj_1_sub@meta.data$slide <- "slide1"

## RB2 ----
xenium_obj_2_sub_meta <- xenium_obj_2_sub@meta.data
xenium_obj_2_sub_meta_cell_meta <- full_join(xenium_obj_2_sub_meta, cell_meta_merge_rb2)
rownames(xenium_obj_2_sub_meta_cell_meta) <- rownames(xenium_obj_2_sub@meta.data)
xenium_obj_2_sub@meta.data <- xenium_obj_2_sub_meta_cell_meta
xenium_obj_2_sub@meta.data$slide <- "slide2"


# total number of cells, RB1: 95807
# total number of cells, RB2: 93814


# Make SP reduction for spatial dim plot ----
## RB1 ----
# Add spatial dimension reduction object separately
position_xy <- cbind((xenium_obj_1_sub@meta.data$x_centroid), (xenium_obj_1_sub@meta.data$y_centroid)*-1)
rownames(position_xy) <- rownames(xenium_obj_1_sub@meta.data)
colnames(position_xy) <- c("SP_1", "SP_2")
xenium_obj_1_sub[["sp"]] <- CreateDimReducObject(
  embeddings = position_xy, key = "SP_", assay = DefaultAssay(xenium_obj_1_sub))
DimPlot(xenium_obj_1_sub, reduction = "sp", group.by = "Sample", 
        label = TRUE,
        raster=FALSE) + NoLegend()

## RB2 ----
# Add spatial dimension reduction object separately
position_xy <- cbind((xenium_obj_2_sub@meta.data$x_centroid), (xenium_obj_2_sub@meta.data$y_centroid)*-1)
rownames(position_xy) <- rownames(xenium_obj_2_sub@meta.data)
colnames(position_xy) <- c("SP_1", "SP_2")
xenium_obj_2_sub[["sp"]] <- CreateDimReducObject(
  embeddings = position_xy, key = "SP_", assay = DefaultAssay(xenium_obj_2_sub))
DimPlot(xenium_obj_2_sub, reduction = "sp", group.by = "Sample", 
        #label = TRUE,
        raster=FALSE) + NoLegend()

### Adjust x and y centroid ----
xenium_obj_1_sub@meta.data$adj_x_centroid <- xenium_obj_1_sub@meta.data$x_centroid
xenium_obj_2_sub@meta.data$adj_x_centroid <- xenium_obj_2_sub@meta.data$x_centroid+8000

xenium_obj_1_sub@meta.data$adj_y_centroid <- xenium_obj_1_sub@meta.data$y_centroid
xenium_obj_2_sub@meta.data$adj_y_centroid <- xenium_obj_2_sub@meta.data$y_centroid

# Merge objects ----
# Merge objects (cannot do spatial DimPlots for this)
merged_spatial_unfiltered <- merge(x = xenium_obj_1_sub, y = xenium_obj_2_sub)
# Warning: Some cell names are duplicated across objects provided. Renaming to enforce unique cell names.
# need to check what these were renamed to and whether I should rename them myself
# Warning message:
# Key ‘Xenium_’ taken, using ‘fov2_’ instead 

## NEW - update cell id column so that there aren't any duplicated cell ids
merged_spatial_unfiltered@meta.data$cell_id <- rownames(merged_spatial_unfiltered@meta.data)
#View(table(merged_spatial_unfiltered@meta.data$cell_id))


# Add spatial dimension reduction object separately
position_xy <- cbind(merged_spatial_unfiltered$adj_x_centroid,
                     (merged_spatial_unfiltered$adj_y_centroid)*-1)
row.names(position_xy) <- row.names(merged_spatial_unfiltered@meta.data)
colnames(position_xy) <- c("SP_1", "SP_2")
merged_spatial_unfiltered[["sp"]] <- CreateDimReducObject(
  embeddings = position_xy, key = "SP_", assay = DefaultAssay(merged_spatial_unfiltered))


DimPlot(merged_spatial_unfiltered, reduction = "sp", 
        group.by = "Sample", label = F,
        raster = F) + NoLegend()


DimPlot(merged_spatial_unfiltered, reduction = "sp", 
        group.by = "slide", 
        label = F,
        raster = F) + NoLegend()

# Save merged, unfiltered data ----
saveRDS(merged_spatial_unfiltered, "/scratch/aoill/projects/heart_transplant/xenium_spatial_cells_only_unfiltered.rds")
ncol(merged_spatial_unfiltered)
# total number of cells: 189621
#merged_spatial_unfiltered <- readRDS("/scratch/aoill/projects/heart_transplant/xenium_spatial_cells_only_unfiltered.rds")


# QC ----
# remove cells with 0 counts
merged_spatial_filtered <- subset(merged_spatial_unfiltered, subset = nCount_Xenium > 0)
ncol(merged_spatial_filtered)
# total number of cells: 189109


VlnPlot(merged_spatial_filtered, features = c("nFeature_Xenium", "nCount_Xenium"
                                              #"nCount_BlankCodeword", "nFeature_BlankCodeword",
                                              #"nCount_ControlCodeword", "nFeature_ControlCodeword",
                                              #"nCount_ControlProbe", "nFeature_ControlProbe"
), 
ncol = 2, pt.size = 0)


## Calculating percent blanks ----
merged_spatial_filtered@meta.data$percent_BlankCodeword <- merged_spatial_filtered@meta.data$nCount_BlankCodeword/merged_spatial_filtered@meta.data$total_counts
merged_spatial_filtered@meta.data$percent_ControlCodeword <- merged_spatial_filtered@meta.data$nCount_ControlCodeword/merged_spatial_filtered@meta.data$total_counts
merged_spatial_filtered@meta.data$percent_ControlProbe <- merged_spatial_filtered@meta.data$nCount_ControlProbe/merged_spatial_filtered@meta.data$total_counts

# Plot percent blanks...
VlnPlot(merged_spatial_filtered, features = c("nFeature_Xenium", "nCount_Xenium",
                                              "percent_BlankCodeword", "percent_ControlCodeword",
                                              "percent_ControlProbe"
), 
ncol = 3, pt.size = 0,
raster = F)


## Plot QC metrics ----
summary(as.factor(merged_spatial_filtered$Sample))

pdf("/home/aoill/plots/n_cells_heart_transplant.pdf",
    height = 6, width =  22)
merged_spatial_filtered@meta.data %>%
  ggplot(aes(y = fct_rev(fct_infreq(Sample)))) +
  geom_bar() + xlab("# Cells") +
  ylab("Sample") +
  ggtitle("Unfiltered") + 
  coord_flip() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()



# Smoothscatter plots
# nCount_Xenium vs. nFeature_Xenium
smoothScatter(log(merged_spatial_filtered$nCount_Xenium),
              log(merged_spatial_filtered$nFeature_Xenium),
              cex = 0.5, pch = 16,
              xlab = "log(nCount_Xenium)",
              ylab = "log(nFeature_Xenium)")
abline(v = log(12), h = log(10), lty = "dashed", col = "black")
text(2.5, 5, col = "black", adj = c(0, -.1),
     "nCount_Xenium >= 12 & nFeature_Xenium >= 10")


# nCount vs. nucleus_area
smoothScatter(merged_spatial_filtered$nucleus_area,
              log(merged_spatial_filtered$nCount_Xenium),
              cex = 0.5, pch = 16)
abline(v = c(6, 70), h = log(12), lty = "dashed", col = "black")
text(90, 0.7, col = "black", adj = c(0, -.1),
     "nCount_Xenium >= 12 & nucleus_area between 6-70")

# nFeature vs. nucleus_area
smoothScatter(merged_spatial_filtered$nucleus_area,
              log(merged_spatial_filtered$nFeature_Xenium),
              cex = 0.5, pch = 16)
abline(v = c(6, 70), h = log(10), lty = "dashed", col = "black")
text(90, 0.7, col = "black", adj = c(0, -.1),
     "nFeature_Xenium >= 10 & nucleus_area between 6-70")


# Filter merged and individual data
merged_spatial_filtered_filtered <- subset(merged_spatial_filtered,
                                           subset = nCount_Xenium >= 12 & nFeature_Xenium >= 10 &
                                             percent_BlankCodeword <= 5 & percent_ControlCodeword <= 5 & 
                                             percent_ControlProbe <= 5 & 
                                             nucleus_area >= 6 & nucleus_area <= 70
)

# Number of nuclei before and after filtering
bf_cells <- table(merged_spatial_filtered$Sample)
aft_cells <- table(merged_spatial_filtered_filtered$Sample)
diff_cells <- bf_cells - aft_cells
prop_kept_cells <- round(aft_cells/bf_cells*100, 2)
prop_kept_cells

merged_spatial_filtered_filtered@meta.data %>%
  ggplot(aes(y = fct_rev(fct_infreq(Sample)))) +
  geom_bar() + xlab("# Cells") +
  ylab("Sample") +
  ggtitle("Filtered") + 
  coord_flip() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

bf_cells <- ncol(merged_spatial_filtered)
aft_cells <- ncol(merged_spatial_filtered_filtered)
# DimPlots of before and after for each sample
a <- DimPlot(merged_spatial_filtered, reduction = "sp", group.by = "Sample", 
             label = F,
             raster=FALSE) + NoLegend() +
  labs(title = paste0("Unfiltered, ", bf_cells, " cells"))
b <- DimPlot(merged_spatial_filtered_filtered, reduction = "sp", group.by = "Sample", 
             label = F,
             raster=FALSE) + NoLegend() +
  labs(title = paste0("Filtered, ", aft_cells, " cells"))

pdf("/home/aoill/plots/sp_dim_heart_transplant_filtering.pdf",
    width = 40, height = 40)
ggarrange(a,b)
dev.off()

# Output filtered data for clustering through the rapids pipeline ----
saveRDS(merged_spatial_filtered_filtered, "/scratch/aoill/projects/heart_transplant/xenium_spatial_cells_only_filtered.rds")
#merged_spatial_filtered_filtered <- readRDS("/scratch/aoill/projects/heart_transplant/xenium_spatial_cells_only_filtered.rds")


# Read in clustered data ----
clustered_obj <- read_rds("/scratch/aoill/projects/heart_transplant/xenium_spatial_cells_only_filtered_clustered_NC100_NN20_PC20_2024_05_22.rds")
#saveRDS(clustered_obj, "/scratch/aoill/projects/heart_transplant/xenium_spatial_cells_only_filtered_clustered_NC100_NN20_PC20_2024_05_22.rds")

## add sample metadata ----
sample_metadata <- read.csv("/home/aoill/projects/heart_transplant/deidentified_spatial_information_edited.csv")
sample_metadata$Sample <- sample_metadata$Donor_Block_ID_uniq
sample_metadata_duprm <- unique(sample_metadata)

clustered_obj_metadata <- clustered_obj@meta.data

clustered_obj_metadata_join <- left_join(clustered_obj_metadata, sample_metadata_duprm)
clustered_obj@meta.data <- clustered_obj_metadata_join


clustered_obj@meta.data <- clustered_obj@meta.data %>% 
  mutate(rejection_type = case_when(
    cellular_grading >= 2 & (antibody_grading == "pAMR0" | antibody_grading == "pAMR1-h" | is.na(antibody_grading) ) ~ "cellular_rejection",
    cellular_grading < 2 & (antibody_grading == "pAMR2" | antibody_grading == "pAMR1-i" ) ~ "antibody_rejection",
    cellular_grading >= 2 & (antibody_grading == "pAMR2" | antibody_grading == "pAMR1-i" ) ~ "mixed_rejection",
    cellular_grading < 2 & (antibody_grading == "pAMR0" | antibody_grading == "pAMR1-h" | is.na(antibody_grading) ) ~ "no_rejection",
    TRUE ~ "other"))
rownames(clustered_obj@meta.data) <- rownames(clustered_obj_metadata)

# Fix cluster types to character
clustered_obj$leiden_res0.1 <- as.character(clustered_obj$leiden_res0.1)
clustered_obj$leiden_res0.2 <- as.character(clustered_obj$leiden_res0.2)
clustered_obj$leiden_res0.3 <- as.character(clustered_obj$leiden_res0.3)
clustered_obj$leiden_res0.5 <- as.character(clustered_obj$leiden_res0.5)
clustered_obj$leiden_res1.0 <- as.character(clustered_obj$leiden_res1.0)
clustered_obj$leiden_res1.5 <- as.character(clustered_obj$leiden_res1.5)
clustered_obj$leiden_res2.0 <- as.character(clustered_obj$leiden_res2.0)


# Define the number of colors you want
nb.cols <- length(unique(clustered_obj_smpl@meta.data$leiden_res0.5))
mycolors <- colorRampPalette(brewer.pal(8, "Paired"))(nb.cols)

randomcoloR::distinctColorPalette(nb.cols)
mycolors <- randomcoloR::distinctColorPalette(nb.cols)
DimPlot(clustered_obj_smpl, 
        reduction = "SP",
        raster = F,
        group.by = "leiden_res0.5",
        label = F) + 
  scale_color_manual(values = mycolors) 
#+ NoLegend()


FeaturePlot(clustered_obj, 
            features = c(
              "PTPRC", # Immune
              "EPCAM", # Epithelial
              "PECAM1", # Endothelial
              "TTN", "FHL2" # Cardiomyocyte
            ),
            reduction = "umap",
            raster=FALSE,
            ncol = 3
)

VlnPlot(clustered_obj, 
        features = c("PTPRC", # Immune
                     "EPCAM", # Epithelial
                     "PECAM1", # Endothelial
                     "TTN", "FHL2" # Cardiomyocyte
        ),
        group.by = "leiden_res0.5",
        pt.size = 0, raster=FALSE,
        ncol = 3)


FeaturePlot(clustered_obj, 
            features = c(
              #"PTPRC", # Immune
              #"CD163", "MRC1", # Macrophage
              #"CD14", "FCGR3A" # Monocyte
              #"GNLY", "NKG7", #NK
              #"CD3E" # T,
              "MS4A2", # mast
              "MS4A1" # B
            ),
            reduction = "umap",
            raster=FALSE,
            ncol = 2
)


FeaturePlot(clustered_obj, 
            features = c(
              "ACTA2" # smooth muscle
            ),
            reduction = "umap",
            raster=FALSE,
            ncol = 1
)

# Cardiomyocytes
FeaturePlot(clustered_obj, 
            features = c(
              # Cardiomyocytes
              "TTN", "FHL2"
              # "PROX1", "PTGDS", "S100A1"
              #"PEBP4", "PLIN4", "SMYD2"
              #"DES", # Smooth muscle cells, Heart - Cardiomyocytes
              
            ),
            reduction = "umap",
            raster=FALSE,
            ncol = 2
)

VlnPlot(clustered_obj, 
        features = c(
          "TTN", "FHL2" # Cardiomyocyte
        ),
        group.by = "leiden_res0.5",
        pt.size = 0, raster=FALSE,
        ncol = 2)

# Vascular smooth muscle cells (VSMCs) 
FeaturePlot(clustered_obj, 
            features = c(
              "ACTA2", "MYH11"
            ),
            reduction = "umap",
            raster=FALSE,
            ncol = 2
)

# Fibroblasts:
FeaturePlot(clustered_obj, 
            features = c(
              "DCN", "C7", "FBLN1", "LTBP2", "OGN", "PDGFRA"
              # "CTSK", 
              # "DPT", 
              # "ASPN", # Fibroblasts and smooth muscle
              # "FBN1", # also smooth muscle
              # "DST" # also cardiomyocytes
              
              #"COL5A2", "CRISPLD2" # could also be marking smooth muscle
              #"CD34", # also endothelial 
            ),
            reduction = "umap",
            raster=FALSE,
            ncol = 3
)

# Pericytes
FeaturePlot(clustered_obj, 
            features = c(
              "PDGFRB", "ACTA2"
            ),
            reduction = "umap",
            raster=FALSE,
            ncol = 2
)


# Adipocytes
FeaturePlot(clustered_obj, 
            features = c(
              "ADIPOQ"
            ),
            reduction = "umap",
            raster=FALSE)

# Endocardial cells
FeaturePlot(clustered_obj, 
            features = c(
              "VWF", "PECAM1", "BMX", "NPR3"
            ),
            reduction = "umap",
            raster=FALSE,
            ncol = 2
)



# proliferating
FeaturePlot(clustered_obj, 
            features = c(
              "MKI67"
            ),
            reduction = "umap",
            raster=FALSE
)

# Neuronal (from internet)
FeaturePlot(clustered_obj, 
            features = c(
              "LGI4"
            ),
            reduction = "umap",
            raster=FALSE
)

# Larger/rough cell types
#Immune: 0, 2, 5, 6, 10, 11, 12
#Endothelial: 8
#Epithelial: none?
#Mesenchymal: 7, 9
#Cardiomyocytes: 1, 4
#?: 3


# All lineages 
FeaturePlot(clustered_obj, 
            features = c(
              "PTPRC", # Immune
              "EPCAM", # Epithelial
              "PECAM1", # Endothelial
              "ACTA2", "PDGFRB", # VSMCs, pericytes
              "DCN", # Fibroblasts
              "TTN", "FHL2" # Cardiomyocyte
            ),
            reduction = "umap",
            raster=FALSE, 
            ncol = 4)


VlnPlot(clustered_obj, 
        features = c("PTPRC", # Immune
                     "EPCAM", # Epithelial
                     "PECAM1", # Endothelial
                     "ACTA2", "PDGFRB", # VSMCs, pericytes
                     "DCN", # Fibroblasts
                     "TTN", "FHL2" # Cardiomyocyte
        ),
        group.by = "leiden_res0.5",
        pt.size = 0, raster=FALSE,
        ncol = 4)


# Finding cluster markers for a heatmap/dotplot
markers <- presto::wilcoxauc(clustered_obj, group_by = "leiden_res0.5", 
                             assay = "data", seurat_assay = "RNA")
topmarkers <- top_markers(markers)

topmarkers

all_markers <- topmarkers %>%
  select(-rank) %>% 
  unclass() %>% 
  stack() %>%
  pull(values) %>%
  unique() %>%
  .[!is.na(.)]


# Pass 1 - assign cells to the correct lineage ----
## Immune ----
cell_id_imm <- clustered_obj@meta.data %>% 
  filter(leiden_res0.5 %in% c("0", "2", "3", "5", "6", "10", "11", "12")) %>% 
  dplyr::pull(cell_id)

merged_spatial_filtered_imm <- subset(merged_spatial_filtered_filtered, 
                                      subset = cell_id %in% cell_id_imm)


saveRDS(merged_spatial_filtered_imm, "/scratch/aoill/projects/heart_transplant/new/merged_spatial_filtered_imm_pass1.rds")

clustered_obj_imm <- read_rds("/scratch/aoill/projects/heart_transplant/new/clustered_obj_imm_pass1_NC100_NN20_PC20_2024_06_10.rds")

# Fix cluster types to character
clustered_obj_imm$leiden_res0.1 <- as.character(clustered_obj_imm$leiden_res0.1)
clustered_obj_imm$leiden_res0.2 <- as.character(clustered_obj_imm$leiden_res0.2)
clustered_obj_imm$leiden_res0.3 <- as.character(clustered_obj_imm$leiden_res0.3)
clustered_obj_imm$leiden_res0.7 <- as.character(clustered_obj_imm$leiden_res0.7)
clustered_obj_imm$leiden_res0.7 <- as.character(clustered_obj_imm$leiden_res0.7)
clustered_obj_imm$leiden_res0.75 <- as.character(clustered_obj_imm$leiden_res0.75)
clustered_obj_imm$leiden_res1.0 <- as.character(clustered_obj_imm$leiden_res1.0)
clustered_obj_imm$leiden_res1.5 <- as.character(clustered_obj_imm$leiden_res1.5)
clustered_obj_imm$leiden_res2.0 <- as.character(clustered_obj_imm$leiden_res2.0)


DimPlot(clustered_obj_imm, 
        reduction = "umap",
        raster = F,
        group.by = "leiden_res0.7",
        label = T
) + NoLegend()
#LabelClusters(p, id = "ident",  fontface = "bold", color = "red")

DimPlot(clustered_obj_imm, 
        reduction = "umap",
        raster = F,
        group.by = "leiden_res0.7",
        split.by = "leiden_res0.7",
        ncol = 4)

FeaturePlot(clustered_obj_imm, 
            features = c(
              "nCount_Xenium",
              "nFeature_Xenium"
            ),
            reduction = "umap",
            raster=FALSE,
            ncol = 2
)

VlnPlot(clustered_obj_imm, 
        features = c( "nCount_Xenium",
                      "nFeature_Xenium"
        ),
        group.by = "leiden_res0.7",
        pt.size = 0,
        ncol = 2, raster=FALSE)

FeaturePlot(clustered_obj_imm, 
            features = c(
              "PTPRC"#, # Immune
              #"EPCAM", # Epithelial
              #"PECAM1" # Endothelial
            ),
            reduction = "umap",
            raster=FALSE,
            ncol = 1
)
VlnPlot(clustered_obj_imm, 
        features = c("PTPRC"#, # Immune
                     #"EPCAM", # Epithelial
                     #"PECAM1" # Endothelial
        ),
        group.by = "leiden_res0.7",
        pt.size = 0,
        ncol = 1, raster=FALSE)


FeaturePlot(clustered_obj_imm, 
            features = c(
              "TTN", "DES", "FHL2", "S100A1"
            ),
            reduction = "umap",
            raster=FALSE,
            ncol = 2
)
VlnPlot(clustered_obj_imm, 
        features = c("TTN", "DES", "FHL2", "S100A1"
        ),
        group.by = "leiden_res0.7",
        pt.size = 0,
        ncol = 2, raster=FALSE)

DefaultAssay(clustered_obj_imm) <- "RNA" 
clustered_obj_imm <- NormalizeData(clustered_obj_imm)
VariableFeatures(clustered_obj_imm) <- rownames(clustered_obj_imm)
clustered_obj_imm <- ScaleData(clustered_obj_imm)

DotPlot(clustered_obj_imm, 
        features = c(#"PECAM1",
          "PTPRC",
          "CD163", "MRC1", "MARCO", "FCGR1A", # Macrophage
          "LYZ", "CD14", "FCGR3A",  "MS4A7", # Monocyte
          "FCER1A", #mo-DCs, cDC1, or pDC
          "CD1A", "CD1C", "CD209", "ITGAM", # monocyte derived DCs
          "ITGAX", "XCR1", # cDC1
          "GNLY", "NKG7", "KLRB1", #"KLRC1", #NK 
          "KLRD1", #NK
          "CD3E", "CD3D", "TRAC", "CD8A", "CD4", "FOXP3", # T
          "MS4A2", # mast
          "MS4A1", # B
          "CD79A", # B and plasma
          "TNFRSF17", "DERL3", # plasma
          "LILRA4", "CCR7", # pDCs
          "MKI67"
        ), 
        group.by = "leiden_res0.7", 
        scale = TRUE) + coord_flip()

DotPlot(clustered_obj_imm, 
        features = c(
          "TSPY3", "PSG2", "RTKN2", "ELF5", "SOX2", "CELF4",
          "EPCAM", "PECAM1",
          "VWF", "BMX", # Ednocardial cells
          "MYH11",  # Vascular smooth muscle cells (VSMCs) 
          "PDGFRB", "ACTA2", # Pericytes
          "DCN", "C7", "FBLN1", "LTBP2", 
          "OGN", "PDGFRA", # Fibroblasts
          "ADIPOQ", # Adipocytes
          "TTN", "DES", "FHL2", "S100A1", # cardiomyocytes
          "PTPRC",
          "CD163", "MRC1", "MARCO", "FCGR1A", # Macrophage
          "LYZ", "CD14", "FCGR3A",  "MS4A7", # Monocyte
          "FCER1A", #mo-DCs, cDC1, or pDC
          "CD1A", "CD1C", "CD209", "ITGAM", # monocyte derived DCs
          "ITGAX", "XCR1", # cDC1
          "GNLY", "NKG7", "KLRB1", #"KLRC1", #NK 
          "KLRD1", #NK
          "CD3E", "CD3D", "TRAC", "CD8A", "CD4", "FOXP3", # T
          "MS4A2", # mast
          "MS4A1", # B
          "CD79A", # B and plasma
          "TNFRSF17", "DERL3", # plasma
          "LILRA4", "CCR7", # pDCs
          "MKI67"
        ), 
        group.by = "leiden_res0.7", 
        scale = TRUE) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

# T cells/NK
FeaturePlot(clustered_obj_imm, 
            features = c(
              "GNLY", "NKG7", "KLRB1", #"KLRC1", #NK 
              "KLRD1", #NK
              "CD3E", "CD3D", "TRAC", "CD8A", "CD4", "FOXP3" # T
              # "IL7R", "CCR7", # T
            ),
            reduction = "umap",
            raster=FALSE,
            ncol = 5
)
VlnPlot(clustered_obj_imm, 
        features = c(
          "GNLY", "NKG7", "KLRB1", #"KLRC1", #NK 
          "KLRD1", #NK
          "CD3E", "CD3D", "TRAC", "CD8A", "CD4", "FOXP3" # T
          # "IL7R", "CCR7", # T
          
        ),
        group.by = "leiden_res0.7",
        pt.size = 0,
        ncol = 4, raster=FALSE)



# Mono/mac
FeaturePlot(clustered_obj_imm, 
            features = c(
              "CD163", "MRC1", "MARCO", "FCGR1A", # Macrophage
              "LYZ", "CD14", "FCGR3A",  "MS4A7" # Monocyte
            ),
            reduction = "umap",
            raster=FALSE,
            ncol = 4
)

VlnPlot(clustered_obj_imm, 
        features = c(
          "CD163", "MRC1", "MARCO", "FCGR1A", # Macrophage
          "LYZ", "CD14", "FCGR3A",  "MS4A7" # Monocyte
          
        ),
        group.by = "leiden_res0.7",
        pt.size = 0,
        ncol = 4, raster=FALSE)



# moDCs
FeaturePlot(clustered_obj_imm, 
            features = c(
              "FCER1A", #mo-DCs, cDC1, or pDC
              "CD1A", "CD1C",	"MRC1", "CD209", "ITGAM" # monocyte derived DCs
              
            ),
            reduction = "umap",
            raster=FALSE,
            ncol = 4 
)
VlnPlot(clustered_obj_imm, 
        features = c(
          "FCER1A", #mo-DCs, cDC1, or pDC
          "CD1A", "CD1C",	"MRC1", "CD209", "ITGAM" # monocyte derived DCs
          
        ),
        group.by = "leiden_res0.7",
        pt.size = 0,
        ncol = 3, raster=FALSE)

# DC1 and 2
FeaturePlot(clustered_obj_imm, 
            features = c(
              "CD8A", "ITGAX", "XCR1", # cDC1
              "CD1C", "ITGAM" # cDC2 or mo-DC
            ),
            reduction = "umap",
            raster=FALSE,
            ncol = 3 
)
VlnPlot(clustered_obj_imm, 
        features = c(
          "CD8A", "ITGAX", "XCR1", # cDC1
          "CD1C", "ITGAM" # cDC2 or mo-DC
          
        ),
        group.by = "leiden_res0.7",
        pt.size = 0,
        ncol = 3, raster=FALSE)


# pDCs
FeaturePlot(clustered_obj_imm, 
            features = c(
              "LILRA4", "CCR7" # pDCs
            ),
            reduction = "umap",
            raster=FALSE,
            ncol = 2 
)
VlnPlot(clustered_obj_imm, 
        features = c(
          "LILRA4", "CCR7" # pDCs
        ),
        group.by = "leiden_res0.7",
        pt.size = 0,
        ncol = 2, raster=FALSE)


# B and plasma
FeaturePlot(clustered_obj_imm, 
            features = c(
              "MS4A1", # B
              "CD79A", # B and plasma
              "TNFRSF17", "DERL3" # plasma
            ),
            reduction = "umap",
            raster=FALSE,
            ncol = 2 
)
VlnPlot(clustered_obj_imm, 
        features = c(
          "MS4A1", # B
          "CD79A", # B and plasma
          "TNFRSF17", "DERL3" # plasma       
        ),
        group.by = "leiden_res0.7",
        pt.size = 0,
        ncol = 2, raster=FALSE)

FeaturePlot(clustered_obj_imm, 
            features = c(
              "BANK1"
            ),
            reduction = "umap",
            raster=FALSE,
            ncol = 2 
)

# Mast
FeaturePlot(clustered_obj_imm, 
            features = c(
              "MS4A2" # mast
            ),
            reduction = "umap",
            raster=FALSE 
)
VlnPlot(clustered_obj_imm, 
        features = c(
          "MS4A2" # mast      
        ),
        group.by = "leiden_res0.7",
        pt.size = 0, raster=FALSE)


FeaturePlot(clustered_obj_imm, 
            features = c(
              "MS4A2", # mast
              "SLC18A2"#, "VWA5A" # granulocyte
            ),
            reduction = "umap",
            raster=FALSE 
)

# Proliferating
FeaturePlot(clustered_obj_imm, 
            features = c(
              "MKI67"
            ),
            reduction = "umap",
            raster=FALSE
)

FeaturePlot(clustered_obj_imm, 
            features = c(
              "TTN", "FHL2"
            ),
            reduction = "umap",
            raster=FALSE
)

# Top markers
markers <- presto::wilcoxauc(clustered_obj_imm, group_by = "leiden_res0.7", 
                             assay = "data", seurat_assay = "RNA")
topmarkers <- top_markers(markers)
topmarkers <- top_markers(markers, n = 10, auc_min = 0.6, pct_in_min = 50, 
                          pct_out_max = 100)


topmarkers



### Clusters that are not immune ----
# Mesenchymal 9, 7, 0
# Cardiomyocyte: 4

# Putative lineages
# Myeloid: 3, 5, 12, 13, 15
# Lymphoid: 1, 2, 6, 8, 10, 11, 14  
# But will re-cluster and umap in pass 2 then go down to these groups


## Mesenchymal ----
# This does not include cardiomyocytes 
# 7, 9
cell_id_mes <- clustered_obj@meta.data %>% 
  filter(leiden_res0.5 %in% c("7", "9")) %>% 
  dplyr::pull(cell_id)
merged_spatial_filtered_mes <- subset(merged_spatial_filtered_filtered, 
                                      subset = cell_id %in% cell_id_mes)


saveRDS(merged_spatial_filtered_mes, "/scratch/aoill/projects/heart_transplant/new/merged_spatial_filtered_mes_pass1.rds")
clustered_obj_mes <- read_rds("/scratch/aoill/projects/heart_transplant/new/clustered_obj_mes_pass1_NC100_NN20_PC20_2024_06_10.rds")

# Fix cluster types to character
clustered_obj_mes$leiden_res0.1 <- as.character(clustered_obj_mes$leiden_res0.1)
clustered_obj_mes$leiden_res0.2 <- as.character(clustered_obj_mes$leiden_res0.2)
clustered_obj_mes$leiden_res0.3 <- as.character(clustered_obj_mes$leiden_res0.3)
clustered_obj_mes$leiden_res0.5 <- as.character(clustered_obj_mes$leiden_res0.5)
clustered_obj_mes$leiden_res1.0 <- as.character(clustered_obj_mes$leiden_res1.0)
clustered_obj_mes$leiden_res1.5 <- as.character(clustered_obj_mes$leiden_res1.5)
clustered_obj_mes$leiden_res2.0 <- as.character(clustered_obj_mes$leiden_res2.0)


DimPlot(clustered_obj_mes, 
        reduction = "umap",
        raster = F,
        group.by = "leiden_res0.7",
        label = T
) + NoLegend()

DimPlot(clustered_obj_mes, 
        reduction = "umap",
        raster = F,
        group.by = "leiden_res0.5",
        split.by = "leiden_res0.5",
        ncol = 4)

FeaturePlot(clustered_obj_mes, 
            features= c(
              "PTPRC", "LYZ", "CD3E" # Immune
              #"EPCAM", # Epithelial
              #"PECAM1" # Endothelial
            ),
            reduction = "umap",
            raster=FALSE,
            ncol = 3
)
VlnPlot(clustered_obj_mes, 
        features = c("PTPRC", "LYZ", "CD3E", # Immune
                     #"EPCAM", # Epithelial
                     #"PECAM1" # Endothelial
                     "COL1A1", "COL1A2"
        ),
        group.by = "leiden_res0.5",
        pt.size = 0,
        ncol = 3, raster=FALSE)




markers <- presto::wilcoxauc(clustered_obj_mes, group_by = "leiden_res0.5", 
                             assay = "data", seurat_assay = "RNA")
topmarkers <- top_markers(markers, n = 10, auc_min = 0.6, pct_in_min = 50, 
                          pct_out_max = 100)


topmarkers

topmarkers <- top_markers(markers)
topmarkers



DefaultAssay(clustered_obj_mes) <- "RNA" 
clustered_obj_mes <- NormalizeData(clustered_obj_mes)
VariableFeatures(clustered_obj_mes) <- rownames(clustered_obj_mes)
clustered_obj_mes <- ScaleData(clustered_obj_mes)


DotPlot(clustered_obj_mes, 
        features = c(
          "EPCAM", "PECAM1",
          "VWF", "BMX", # Ednocardial cells
          "MYH11",  # Vascular smooth muscle cells (VSMCs) 
          "PDGFRB", "ACTA2", # Pericytes
          "DCN", "C7", "FBLN1", "LTBP2", 
          "OGN", "PDGFRA", # Fibroblasts
          "ADIPOQ", # Adipocytes
          "TTN", "DES", "FHL2", "S100A1", # cardiomyocytes
          "PTPRC",
          "CD163", "MRC1", "MARCO", "FCGR1A", # Macrophage
          "LYZ", "CD14", "FCGR3A",  "MS4A7", # Monocyte
          "FCER1A", #mo-DCs, cDC1, or pDC
          "CD1A", "CD1C", "CD209", "ITGAM", # monocyte derived DCs
          "ITGAX", "XCR1", # cDC1
          "GNLY", "NKG7", "KLRB1", #"KLRC1", #NK 
          "KLRD1", #NK
          "CD3E", "CD3D", "TRAC", "CD8A", "CD4", "FOXP3", # T
          "MS4A2", # mast
          "MS4A1", # B
          "CD79A", # B and plasma
          "TNFRSF17", "DERL3", # plasma
          "LILRA4", "CCR7", # pDCs
          "MKI67"
        ), 
        group.by = "leiden_res0.5", 
        scale = TRUE) + coord_flip()



# Fibroblasts: 
FeaturePlot(clustered_obj_mes, 
            features = c(
              #"SEMA3C", # myofibroblast
              #"PMP22", 
              #"FAP", 
              "DCN", "C7", "FBLN1", "LTBP2", 
              "OGN", "PDGFRA" # Fibroblasts
              #"SFRP4", "MFAP5", # Adventicial Fibroblasts
              #"SFRP2", "PTGDS", # activated
              #"ACTA2" # could be expressed in myofibroblast
              
            ),
            reduction = "umap",
            raster=FALSE,
            ncol = 3
)

FeaturePlot(clustered_obj_mes, 
            features = c(
              "SEMA3C", # myofibroblast
              "SFRP4", "MFAP5", # Adventicial Fibroblasts
              "SFRP2", "PTGDS" # activated
              #"ACTA2" # could be expressed in myofibroblast
              
            ),
            reduction = "umap",
            raster=FALSE,
            ncol = 3
)

FeaturePlot(clustered_obj_mes, 
            features = c(
              "ACTA2", "MYH11",  # Vascular smooth muscle cells (VSMCs) 
              "PDGFRB", "ACTA2"#, # Pericytes
              #"ADIPOQ" # Adipocytes
              #"VWF", "BMX" # Ednocardial cells
              
            ),
            reduction = "umap",
            raster=FALSE,
            ncol = 2
)


# Cardiomyocytes
#"ASCL3", "MTRNR2L11", "MYBPC1", "PPP1R1A", "SH2D3C", # these don't mark cardiomyocytes well
# "DES" # Smooth muscle cells, Heart - Cardiomyocytes
FeaturePlot(clustered_obj_mes, 
            features = c(
              # Cardiomyocytes
              "TTN", "FHL2",
              # Below are good, featureplot just gets fussy when plotting a lot
              # of genes
              "PROX1", "PTGDS", "S100A1",
              "PEBP4", "PLIN4", "SMYD2"
              #"DES", # Smooth muscle cells, Heart - Cardiomyocytes
              
            ),
            reduction = "umap",
            raster=FALSE,
            ncol = 4
)


# Activated fibroblasts/myofibroblasts
# ACTA2, FSP1/S100A4 TNC
FeaturePlot(clustered_obj_mes, 
            features = c(
              "ACTA2", "FSP1", "S100A4", "TNC"
            ),
            reduction = "umap",
            raster=FALSE,
            ncol = 2
)

# Cardiac Pericytes 
FeaturePlot(clustered_obj_mes, 
            features = c(
              "PDGFRB", "NG2", "CSPG4", "TNC"
            ),
            reduction = "umap",
            raster=FALSE,
            ncol = 2
)

# Cardiac fibroblasts:
FeaturePlot(clustered_obj_mes, 
            features = c(
              "COL1A1", "COL3A1", "FN1", "POSTN", "DCN", "SPARC", "TGFB1"
            ),
            reduction = "umap",
            raster=FALSE,
            ncol = 3
)


# Myofibroblasts (activated fibroblasts):
FeaturePlot(clustered_obj_mes2, 
            features = c(
              "ACTA2", "TAGLN", "VIM", "FAP", "DDR2", "CTGF"
            ),
            reduction = "umap",
            raster=FALSE,
            ncol = 2
)


#Perivascular fibroblasts:
FeaturePlot(clustered_obj_mes, 
            features = c(
              "PDGFRB", "NG2", "RGS5", "DES", "RERG"
            ),
            reduction = "umap",
            raster=FALSE,
            ncol = 2
)


# Endocardial fibroblasts:
FeaturePlot(clustered_obj_mes, 
            features = c(
              "NPPA", "NPPB", "CAV1", "TNNI3", "CHNG5"
            ),
            reduction = "umap",
            raster=FALSE,
            ncol = 2
)

# Interstitial fibroblasts:
FeaturePlot(clustered_obj_mes, 
            features = c(
              "PDGFRA", "DDR1", "CD90", "CD29", "CD105"
            ),
            reduction = "umap",
            raster=FALSE,
            ncol = 2
)


### Clusters that are not mesenchymal ----
# Cardiomyocyte: 4
# Immune: 7 

# Putative cell types
# vSMCs:0
# Proliferating pericytes: 1
# Fibroblast: 2, 5
# Pericytes: 3
# Myofibroblasts: 6
# But will re-cluster and umap in pass 2 then go down to these groups



## Endothelial ----
# 8
cell_id_endo <- clustered_obj@meta.data %>% 
  filter(leiden_res0.5 %in% c("8")) %>% 
  dplyr::pull(cell_id)
merged_spatial_filtered_endo <- subset(merged_spatial_filtered_filtered, 
                                       subset = cell_id %in% cell_id_endo)

saveRDS(merged_spatial_filtered_endo, "/scratch/aoill/projects/heart_transplant/new/merged_spatial_filtered_endo_pass1.rds")


clustered_obj_endo <- read_rds("/scratch/aoill/projects/heart_transplant/new/clustered_obj_endo_pass1_NC100_NN20_PC20_2024_06_10.rds")


# Fix cluster types to character
clustered_obj_endo$subcluster <- as.character(clustered_obj_endo$subcluster)
clustered_obj_endo$subcluster2 <- as.character(clustered_obj_endo$subcluster2)
clustered_obj_endo$leiden_res0.1 <- as.character(clustered_obj_endo$leiden_res0.1)
clustered_obj_endo$leiden_res0.2 <- as.character(clustered_obj_endo$leiden_res0.2)
clustered_obj_endo$leiden_res0.3 <- as.character(clustered_obj_endo$leiden_res0.3)
clustered_obj_endo$leiden_res0.5 <- as.character(clustered_obj_endo$leiden_res0.5)
clustered_obj_endo$leiden_res1.0 <- as.character(clustered_obj_endo$leiden_res1.0)
clustered_obj_endo$leiden_res1.5 <- as.character(clustered_obj_endo$leiden_res1.5)
clustered_obj_endo$leiden_res2.0 <- as.character(clustered_obj_endo$leiden_res2.0)


DimPlot(clustered_obj_endo, 
        reduction = "umap",
        raster = F,
        group.by = "subcluster2",
        label = T
) + NoLegend()
#LabelClusters(p, id = "ident",  fontface = "bold", color = "red")

DimPlot(clustered_obj_endo, 
        reduction = "umap",
        raster = F,
        group.by = "subcluster",
        split.by = "subcluster",
        ncol = 4)


markers <- presto::wilcoxauc(clustered_obj_endo, group_by = "subcluster2", 
                             assay = "data", seurat_assay = "RNA")
topmarkers <- top_markers(markers, n = 10, auc_min = 0.6, pct_in_min = 50, 
                          pct_out_max = 100)


topmarkers


topmarkers <- top_markers(markers)
topmarkers

DefaultAssay(clustered_obj_endo) <- "RNA" 
clustered_obj_endo <- NormalizeData(clustered_obj_endo)
VariableFeatures(clustered_obj_endo) <- rownames(clustered_obj_endo)
clustered_obj_endo <- ScaleData(clustered_obj_endo)

DotPlot(clustered_obj_endo, 
        features = c(
          "EPCAM", "PECAM1",
          "VWF", "BMX", # Ednocardial cells
          "MYH11",  # Vascular smooth muscle cells (VSMCs) 
          "PDGFRB", "ACTA2", # Pericytes
          "DCN", "C7", "FBLN1", "LTBP2", 
          "OGN", "PDGFRA", # Fibroblasts
          "ADIPOQ", # Adipocytes
          "TTN", "DES", "FHL2", "S100A1", # cardiomyocytes
          "PTPRC",
          "CD163", "MRC1", "MARCO", "FCGR1A", # Macrophage
          "LYZ", "CD14", "FCGR3A",  "MS4A7", # Monocyte
          "FCER1A", #mo-DCs, cDC1, or pDC
          "CD1A", "CD1C", "CD209", "ITGAM", # monocyte derived DCs
          "ITGAX", "XCR1", # cDC1
          "GNLY", "NKG7", "KLRB1", #"KLRC1", #NK 
          "KLRD1", #NK
          "CD3E", "CD3D", "TRAC", "CD8A", "CD4", "FOXP3", # T
          "MS4A2", # mast
          "MS4A1", # B
          "CD79A", # B and plasma
          "TNFRSF17", "DERL3", # plasma
          "LILRA4", "CCR7", # pDCs
          "MKI67"
        ), 
        group.by = "subcluster2", 
        scale = TRUE) + coord_flip()



FeaturePlot(clustered_obj_endo, 
            features = c(
              "DCN", "C7", "FBLN1", "LTBP2", 
              "OGN", "PDGFRA", # Fibroblasts
              "CHST15"
            ),
            reduction = "umap",
            raster=FALSE,
            ncol = 3
)

FeaturePlot(clustered_obj_endo, 
            features = c(
              "PECAM1", "PTPRC"
            ),
            reduction = "umap",
            raster=FALSE,
            ncol = 1
)
VlnPlot(clustered_obj_endo, 
        features = c(
          "PECAM1"
        ),
        group.by = "subcluster",
        pt.size = 0,
        ncol = 1, raster=FALSE)

FeaturePlot(clustered_obj_endo, 
            features = c(
              "PECAM1", "PTPRC", "LYZ"
            ),
            reduction = "umap",
            raster=FALSE,
            ncol = 2
)
VlnPlot(clustered_obj_endo, 
        features = c(
          "PECAM1", "PTPRC", "LYZ"
        ),
        group.by = "subcluster",
        pt.size = 0,
        ncol = 3, raster=FALSE)


FeaturePlot(clustered_obj_endo, 
            features = c(
              "DLL4", "NOTCH4", "EPHB4", "NRP1",
              "CDH5", "EPHB2", "EFNB2", "HEY2",
              "DES", "TTN", "FHL2", "S100A1",
              "MYBPC1", # cardiomyocytes
              "ABCC11"
            ),
            reduction = "umap",
            raster=FALSE,
            ncol = 4
)

# fibroblasts
#POSTN
FeaturePlot(clustered_obj_endo, 
            features = c(
              "POSTN", "TFPI", "SELP", "TMEM100"
            ),
            reduction = "umap",
            raster=FALSE,
            ncol = 4
)

# cluster 2 top
FeaturePlot(clustered_obj_endo, 
            features = c(
              "CXCL11","CXCL9",  "FKBP1A", "UBD", 
              "CTSS",   "VCAM1",  "PLA1A",  "TAP1",   "ICAM1",  "CXCL10"
            ),
            reduction = "umap",
            raster=FALSE,
            ncol = 5
)
# Venous Endothelial Cells	
# COUP-TFII (NR2F2), VEGFR3 (Vascular endothelial growth factor receptor 3), CLDN5 (Claudin-5), NRP2 (Neuropilin-2)
FeaturePlot(clustered_obj_endo, 
            features = c(
              "NR2F2", "VEGFR3", "CLDN5", "NRP2"
            ),
            reduction = "umap",
            raster=FALSE,
            ncol = 4
)



FeaturePlot(clustered_obj_endo, 
            features = c(
              "PROX1", "LYVE1", "PDPN", "VEGFR3"
            ),
            reduction = "umap",
            raster=FALSE,
            ncol = 3
)

FeaturePlot(clustered_obj_endo, 
            features = c(
              "CD34", "CDH5", "CLDN5", "VWF"
            ),
            reduction = "umap",
            raster=FALSE,
            ncol = 3
)


FeaturePlot(clustered_obj_endo, 
            features = c(
              "VCAM1", "CD31", "CDH5", "NRG1",
              "VWF", "PECAM1", "BMX", "NPR3", "TIE1", "FLK1", "ESM1", "THBD"
            ),
            reduction = "umap",
            raster=FALSE,
            ncol = 3
)


# Endocardial cells:
#VWF, PECAM1, BMX, NPR3

FeaturePlot(clustered_obj_endo, 
            features = c(
              "VWF", "PECAM1", "BMX", "NPR3"
            ),
            reduction = "umap",
            raster=FALSE,
            ncol = 3
)


### Clusters that are not endothelial ----
# Cardiomyocytes: 0
# Mesenchymal: 1,0 (vSMCc), 1,1 (fibroblasts)
# Immune: 2,2, 2,3, 2,4

# Putative cell types
# Endocardial cells: 3,0 
# Lymphatic endothelial cells: 3,1
# Proliferating endothelial cells: 2,1
# Endothelial cells: 1,2, 1,3, 4, 2,0, 
# But will re-cluster and umap in pass 2 then go down to these groups

## Cardiomyocytes ----
# 1, 4
cell_id_card <- clustered_obj@meta.data %>% 
  filter(leiden_res0.5 %in% c("1", "4")) %>% 
  dplyr::pull(cell_id)
merged_spatial_filtered_card <- subset(merged_spatial_filtered_filtered, 
                                       subset = cell_id %in% cell_id_card)

saveRDS(merged_spatial_filtered_card, "/scratch/aoill/projects/heart_transplant/new/merged_spatial_filtered_card_pass1.rds")


clustered_obj_card <- read_rds("/scratch/aoill/projects/heart_transplant/new/clustered_obj_card_pass1_NC100_NN20_PC20_2024_06_10.rds")


# Fix cluster types to character
clustered_obj_card$leiden_res0.1 <- as.character(clustered_obj_card$leiden_res0.1)
clustered_obj_card$leiden_res0.2 <- as.character(clustered_obj_card$leiden_res0.2)
clustered_obj_card$leiden_res0.3 <- as.character(clustered_obj_card$leiden_res0.3)
clustered_obj_card$leiden_res0.4 <- as.character(clustered_obj_card$leiden_res0.4)
clustered_obj_card$leiden_res0.5 <- as.character(clustered_obj_card$leiden_res0.5)
clustered_obj_card$leiden_res1.0 <- as.character(clustered_obj_card$leiden_res1.0)
clustered_obj_card$leiden_res1.5 <- as.character(clustered_obj_card$leiden_res1.5)
clustered_obj_card$leiden_res2.0 <- as.character(clustered_obj_card$leiden_res2.0)


DimPlot(clustered_obj_card, 
        reduction = "umap",
        raster = F,
        group.by = "leiden_res0.4",
        label = T
) + NoLegend()


markers <- presto::wilcoxauc(clustered_obj_card, group_by = "leiden_res0.4", 
                             assay = "data", seurat_assay = "RNA")
topmarkers <- top_markers(markers, n = 10, auc_min = 0.6, pct_in_min = 50, 
                          pct_out_max = 100)


topmarkers


topmarkers <- top_markers(markers)
topmarkers

DefaultAssay(clustered_obj_card) <- "RNA" 
clustered_obj_card <- NormalizeData(clustered_obj_card)
VariableFeatures(clustered_obj_card) <- rownames(clustered_obj_card)
clustered_obj_card <- ScaleData(clustered_obj_card)

DotPlot(clustered_obj_card, 
        features = c(
          "EPCAM", "PECAM1",
          "VWF", "BMX", # Ednocardial cells
          "MYH11",  # Vascular smooth muscle cells (VSMCs) 
          "PDGFRB", "ACTA2", # Pericytes
          "DCN", "C7", "FBLN1", "LTBP2", 
          "OGN", "PDGFRA", # Fibroblasts
          "ADIPOQ", # Adipocytes
          "TTN", "DES", "FHL2", "S100A1", # cardiomyocytes
          "PTPRC",
          "CD163", "MRC1", "MARCO", "FCGR1A", # Macrophage
          "LYZ", "CD14", "FCGR3A",  "MS4A7", # Monocyte
          "FCER1A", #mo-DCs, cDC1, or pDC
          "CD1A", "CD1C", "CD209", "ITGAM", # monocyte derived DCs
          "ITGAX", "XCR1", # cDC1
          "GNLY", "NKG7", "KLRB1", #"KLRC1", #NK 
          "KLRD1", #NK
          "CD3E", "CD3D", "TRAC", "CD8A", "CD4", "FOXP3", # T
          "MS4A2", # mast
          "MS4A1", # B
          "CD79A", # B and plasma
          "TNFRSF17", "DERL3", # plasma
          "LILRA4", "CCR7", # pDCs
          "MKI67"
        ), 
        group.by = "leiden_res0.4", 
        scale = TRUE) + coord_flip()


FeaturePlot(clustered_obj_card, 
            features = c(
              "PTPRC", "TTN", "DES", "FHL2", "S100A1"
            ),
            reduction = "umap",
            raster=FALSE,
            ncol = 3
)


FeaturePlot(clustered_obj_card, 
            features = c(
              "PTPRC", "TTN"
            ),
            reduction = "umap",
            raster=FALSE,
            blend = T
)

### Clusters that are not cardiomyocytes ----
# All are cardiomyocytes


# Pass 2 - annotate each lineage ----
## Immune ----
### Get cells from all objects ----

# Original immune object
# Putative lineages
# Myeloid: 3, 5, 12, 13, 15
# Lymphoid: 1, 2, 6, 8, 10, 11, 14 
#DimPlot(clustered_obj_imm, reduction = "umap",
#        group.by = "leiden_res0.7")
cell_id_imm1 <- clustered_obj_imm@meta.data %>% 
  filter(leiden_res0.7 %in% c("3", "5", "12", "13", "15", "1", "2", "6", "8", 
                              "10", "11", "14")) %>% 
  dplyr::pull(cell_id)


# Mesenchymal object
#DimPlot(clustered_obj_mes, reduction = "umap",
#        group.by = "leiden_res0.5")
cell_id_imm2 <- clustered_obj_mes@meta.data %>% 
  filter(leiden_res0.5 %in% c("7")) %>% 
  dplyr::pull(cell_id)

# Endothelial object
#DimPlot(clustered_obj_endo, reduction = "umap",
#        group.by = "subcluster2", label = T)
cell_id_imm3 <- clustered_obj_endo@meta.data %>% 
  filter(subcluster2 %in% c("2,2", "2,3", "2,4")) %>% 
  dplyr::pull(cell_id)


length(cell_id_imm1)
length(cell_id_imm2)
length(cell_id_imm3)

cell_id_imm_pass2 <- c(cell_id_imm1, cell_id_imm2, cell_id_imm3)
length(cell_id_imm_pass2)

merged_spatial_filtered_imm_pass2 <- subset(merged_spatial_filtered_filtered, 
                                            subset = cell_id %in% cell_id_imm_pass2)
dim(merged_spatial_filtered_imm_pass2)

saveRDS(merged_spatial_filtered_imm_pass2, "/scratch/aoill/projects/heart_transplant/new/merged_spatial_filtered_imm_pass2.rds")


clustered_obj_imm_pass2 <- read_rds("/scratch/aoill/projects/heart_transplant/new/clustered_obj_imm_pass2_NC100_NN20_PC20_2024_06_13.rds")

# Fix cluster types to character
clustered_obj_imm_pass2$subcluster2 <- as.character(clustered_obj_imm_pass2$subcluster2)
clustered_obj_imm_pass2$subcluster <- as.character(clustered_obj_imm_pass2$subcluster)
clustered_obj_imm_pass2$leiden_res0.1 <- as.character(clustered_obj_imm_pass2$leiden_res0.1)
clustered_obj_imm_pass2$leiden_res0.2 <- as.character(clustered_obj_imm_pass2$leiden_res0.2)
clustered_obj_imm_pass2$leiden_res0.3 <- as.character(clustered_obj_imm_pass2$leiden_res0.3)
clustered_obj_imm_pass2$leiden_res0.5 <- as.character(clustered_obj_imm_pass2$leiden_res0.5)
clustered_obj_imm_pass2$leiden_res0.7 <- as.character(clustered_obj_imm_pass2$leiden_res0.7)
clustered_obj_imm_pass2$leiden_res0.75 <- as.character(clustered_obj_imm_pass2$leiden_res0.75)
clustered_obj_imm_pass2$leiden_res1.0 <- as.character(clustered_obj_imm_pass2$leiden_res1.0)
clustered_obj_imm_pass2$leiden_res1.5 <- as.character(clustered_obj_imm_pass2$leiden_res1.5)
clustered_obj_imm_pass2$leiden_res2.0 <- as.character(clustered_obj_imm_pass2$leiden_res2.0)




DimPlot(clustered_obj_imm_pass2, 
        reduction = "umap",
        raster = F,
        group.by = "subcluster2",
        label = T
) + NoLegend()

DimPlot(clustered_obj_imm_pass2, 
        reduction = "umap",
        raster = F,
        group.by = "subcluster",
        split.by = "subcluster",
        ncol = 4)

FeaturePlot(clustered_obj_imm_pass2, 
            features = c(
              "nCount_Xenium",
              "nFeature_Xenium"
            ),
            reduction = "umap",
            raster=FALSE,
            ncol = 2
)

VlnPlot(clustered_obj_imm_pass2, 
        features = c( "nCount_Xenium",
                      "nFeature_Xenium"
        ),
        group.by = "subcluster",
        pt.size = 0,
        ncol = 2, raster=FALSE)

FeaturePlot(clustered_obj_imm_pass2, 
            features = c(
              "PTPRC"#, # Immune
              #"EPCAM", # Epithelial
              #"PECAM1" # Endothelial
            ),
            reduction = "umap",
            raster=FALSE,
            ncol = 1
)
VlnPlot(clustered_obj_imm_pass2, 
        features = c("PTPRC"#, # Immune
                     #"EPCAM", # Epithelial
                     #"PECAM1" # Endothelial
        ),
        group.by = "subcluster2",
        pt.size = 0,
        ncol = 1, raster=FALSE)


FeaturePlot(clustered_obj_imm_pass2, 
            features = c(
              "TTN", "DES", "FHL2", "S100A1"
            ),
            reduction = "umap",
            raster=FALSE,
            ncol = 2
)
VlnPlot(clustered_obj_imm_pass2, 
        features = c("TTN", "DES", "FHL2", "S100A1"
        ),
        group.by = "subcluster",
        pt.size = 0,
        ncol = 2, raster=FALSE)

DefaultAssay(clustered_obj_imm_pass2) <- "RNA" 
clustered_obj_imm_pass2 <- NormalizeData(clustered_obj_imm_pass2)
VariableFeatures(clustered_obj_imm_pass2) <- rownames(clustered_obj_imm_pass2)
clustered_obj_imm_pass2 <- ScaleData(clustered_obj_imm_pass2)

DotPlot(clustered_obj_imm_pass2, 
        features = c(#"PECAM1",
          "PTPRC",
          "CD163", "MRC1", "MARCO", "FCGR1A", # Macrophage
          "LYZ", "CD14", "FCGR3A",  "MS4A7", # Monocyte
          "FCER1A", #mo-DCs, cDC1, or pDC
          "CD1A", "CD1C", "CD209", "ITGAM", # monocyte derived DCs
          "ITGAX", "XCR1", # cDC1
          "GNLY", "NKG7", "KLRB1", #"KLRC1", #NK 
          "KLRD1", #NK
          "CD3E", "CD3D", "TRAC", "CD8A", "CD4", "FOXP3", # T
          "MS4A2", # mast
          "MS4A1", # B
          "CD79A", # B and plasma
          "TNFRSF17", "DERL3", # plasma
          "LILRA4", "CCR7", # pDCs
          "MKI67"
        ), 
        group.by = "subcluster2", 
        scale = TRUE) + coord_flip()

DotPlot(clustered_obj_imm_pass2, 
        features = c(
          #"TSPY3", "PSG2", "RTKN2", "ELF5", "SOX2", "CELF4",
          "EPCAM", "PECAM1",
          "VWF", "BMX", # Ednocardial cells
          "MYH11",  # Vascular smooth muscle cells (VSMCs) 
          "PDGFRB", "ACTA2", # Pericytes
          "DCN", "C7", "FBLN1", "LTBP2", 
          "OGN", "PDGFRA", # Fibroblasts
          "ADIPOQ", # Adipocytes
          "TTN", "DES", "FHL2", "S100A1", # cardiomyocytes
          "PTPRC",
          "CD163", "MRC1", "MARCO", "FCGR1A", # Macrophage
          "LYZ", "CD14", "FCGR3A",  "MS4A7", # Monocyte
          "FCER1A", #mo-DCs, cDC1, or pDC
          "CD1A", "CD1C", "CD209", "ITGAM", # monocyte derived DCs
          "ITGAX", "XCR1", # cDC1
          "GNLY", "NKG7", "KLRB1", #"KLRC1", #NK 
          "KLRD1", #NK
          "CD3E", "CD3D", "TRAC", "CD8A", "CD4", "FOXP3", # T
          "MS4A2", # mast
          "MS4A1", # B
          "CD79A", # B and plasma
          "TNFRSF17", "DERL3", # plasma
          "LILRA4", "CCR7", # pDCs
          "MKI67"
        ), 
        group.by = "subcluster2", 
        scale = TRUE) + coord_flip() #theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

# non-immune
FeaturePlot(clustered_obj_imm_pass2, 
            features = c(
              "PECAM1", "VWF", # Endothelial cells
              "MYH11",  # Vascular smooth muscle cells (VSMCs) 
              "PDGFRB", "ACTA2", # Pericytes
              "DCN" # Fibroblasts
            ),
            reduction = "umap",
            raster=FALSE,
            ncol = 3
)
# T cells/NK
FeaturePlot(clustered_obj_imm_pass2, 
            features = c(
              "GNLY", "NKG7", "KLRB1", #"KLRC1", #NK 
              "KLRD1", #NK
              "CD3E", "CD3D", "TRAC", "CD8A", "CD4", "FOXP3" # T
              # "IL7R", "CCR7", # T
            ),
            reduction = "umap",
            raster=FALSE,
            ncol = 5
)
VlnPlot(clustered_obj_imm_pass2, 
        features = c(
          "GNLY", "NKG7", "KLRB1", #"KLRC1", #NK 
          "KLRD1", #NK
          "CD3E", "CD3D", "TRAC", "CD8A", "CD4", "FOXP3" # T
          # "IL7R", "CCR7", # T
          
        ),
        group.by = "subcluster",
        pt.size = 0,
        ncol = 4, raster=FALSE)



# Mono/mac
FeaturePlot(clustered_obj_imm_pass2, 
            features = c(
              "CD163", "MRC1", "MARCO", "FCGR1A", # Macrophage
              "LYZ", "CD14", "FCGR3A",  "MS4A7" # Monocyte
            ),
            reduction = "umap",
            raster=FALSE,
            ncol = 4
)

VlnPlot(clustered_obj_imm_pass2, 
        features = c(
          "CD163", "MRC1", "MARCO", "FCGR1A", # Macrophage
          "LYZ", "CD14", "FCGR3A",  "MS4A7" # Monocyte
          
        ),
        group.by = "subcluster",
        pt.size = 0,
        ncol = 4, raster=FALSE)



# moDCs
FeaturePlot(clustered_obj_imm_pass2, 
            features = c(
              "FCER1A", #mo-DCs, cDC1, or pDC
              "CD1A", "CD1C",	"MRC1", "CD209", "ITGAM" # monocyte derived DCs
              
            ),
            reduction = "umap",
            raster=FALSE,
            ncol = 4 
)
VlnPlot(clustered_obj_imm_pass2, 
        features = c(
          "FCER1A", #mo-DCs, cDC1, or pDC
          "CD1A", "CD1C",	"MRC1", "CD209", "ITGAM" # monocyte derived DCs
          
        ),
        group.by = "subcluster",
        pt.size = 0,
        ncol = 3, raster=FALSE)

# DC1 and 2
FeaturePlot(clustered_obj_imm_pass2, 
            features = c(
              "CD8A", "ITGAX", "XCR1", # cDC1
              "CD1C", "ITGAM" # cDC2 or mo-DC
            ),
            reduction = "umap",
            raster=FALSE,
            ncol = 3 
)
VlnPlot(clustered_obj_imm_pass2, 
        features = c(
          "CD8A", "ITGAX", "XCR1", # cDC1
          "CD1C", "ITGAM" # cDC2 or mo-DC
          
        ),
        group.by = "subcluster",
        pt.size = 0,
        ncol = 3, raster=FALSE)


# pDCs
FeaturePlot(clustered_obj_imm_pass2, 
            features = c(
              "LILRA4", "CCR7" # pDCs
            ),
            reduction = "umap",
            raster=FALSE,
            ncol = 2 
)
VlnPlot(clustered_obj_imm_pass2, 
        features = c(
          "LILRA4", "CCR7" # pDCs
        ),
        group.by = "subcluster",
        pt.size = 0,
        ncol = 2, raster=FALSE)


# B and plasma
FeaturePlot(clustered_obj_imm_pass2, 
            features = c(
              "MS4A1", # B
              "CD79A", # B and plasma
              "TNFRSF17", "DERL3" # plasma
            ),
            reduction = "umap",
            raster=FALSE,
            ncol = 4 
)
VlnPlot(clustered_obj_imm_pass2, 
        features = c(
          "MS4A1", # B
          "CD79A", # B and plasma
          "TNFRSF17", "DERL3" # plasma       
        ),
        group.by = "subcluster",
        pt.size = 0,
        ncol = 2, raster=FALSE)

FeaturePlot(clustered_obj_imm_pass2, 
            features = c(
              "BANK1"
            ),
            reduction = "umap",
            raster=FALSE,
            ncol = 2 
)

# Mast
FeaturePlot(clustered_obj_imm_pass2, 
            features = c(
              "MS4A2" # mast
            ),
            reduction = "umap",
            raster=FALSE 
)
VlnPlot(clustered_obj_imm_pass2, 
        features = c(
          "MS4A2" # mast      
        ),
        group.by = "subcluster",
        pt.size = 0, raster=FALSE)


FeaturePlot(clustered_obj_imm_pass2, 
            features = c(
              "MS4A2", # mast
              "SLC18A2"#, "VWA5A" # granulocyte
            ),
            reduction = "umap",
            raster=FALSE 
)

# Proliferating
FeaturePlot(clustered_obj_imm_pass2, 
            features = c(
              "MKI67"
            ),
            reduction = "umap",
            raster=FALSE
)

FeaturePlot(clustered_obj_imm_pass2, 
            features = c(
              "TTN", "FHL2"
            ),
            reduction = "umap",
            raster=FALSE
)

# Top markers
markers <- presto::wilcoxauc(clustered_obj_imm_pass2, group_by = "subcluster2", 
                             assay = "data", seurat_assay = "RNA")
topmarkers <- top_markers(markers)
topmarkers <- top_markers(markers, n = 10, auc_min = 0.6, pct_in_min = 50, 
                          pct_out_max = 100)


topmarkers


### FINAL Non-immune cell types from immune object ----
# subcluster 2
# Endothelial: 13
# Fibroblasts: 7,0, 7,1
# Pericytes: 7,2
# These are final annotations for these clusters

# Myeloid: 1, 10, 11,0, 11,1, 5, 9, 4
# Lymphoid: 0, 12, 2, 3, 8, 6

### Meyloid ----
# subcluster 2
# Myeloid: 1, 10, 11,0, 11,1, 5, 9, 4
cell_id_meyloid <- clustered_obj_imm_pass2@meta.data %>% 
  filter(subcluster2 %in% c("1", "10", "11,0", "11,1", "5", "9", "4")) %>% 
  dplyr::pull(cell_id)

length(cell_id_meyloid)




merged_spatial_filtered_meyloid <- subset(merged_spatial_filtered_filtered, 
                                          subset = cell_id %in% cell_id_meyloid)
dim(merged_spatial_filtered_meyloid)

saveRDS(merged_spatial_filtered_meyloid, "/scratch/aoill/projects/heart_transplant/new/merged_spatial_filtered_meyloid.rds")


clustered_obj_meyloid <- read_rds("/scratch/aoill/projects/heart_transplant/new/clustered_obj_meyloid_NC100_NN20_PC20_2024_06_13.rds")


# Fix cluster types to character
clustered_obj_meyloid$subcluster <- as.character(clustered_obj_meyloid$subcluster)
clustered_obj_meyloid$leiden_res0.1 <- as.character(clustered_obj_meyloid$leiden_res0.1)
clustered_obj_meyloid$leiden_res0.2 <- as.character(clustered_obj_meyloid$leiden_res0.2)
clustered_obj_meyloid$leiden_res0.3 <- as.character(clustered_obj_meyloid$leiden_res0.3)
clustered_obj_meyloid$leiden_res0.5 <- as.character(clustered_obj_meyloid$leiden_res0.5)
clustered_obj_meyloid$leiden_res0.7 <- as.character(clustered_obj_meyloid$leiden_res0.7)
clustered_obj_meyloid$leiden_res0.75 <- as.character(clustered_obj_meyloid$leiden_res0.75)
clustered_obj_meyloid$leiden_res1.0 <- as.character(clustered_obj_meyloid$leiden_res1.0)
clustered_obj_meyloid$leiden_res1.5 <- as.character(clustered_obj_meyloid$leiden_res1.5)
clustered_obj_meyloid$leiden_res2.0 <- as.character(clustered_obj_meyloid$leiden_res2.0)




DimPlot(clustered_obj_meyloid, 
        reduction = "umap",
        raster = F,
        group.by = "subcluster",
        label = T
) + NoLegend()

DimPlot(clustered_obj_meyloid, 
        reduction = "umap",
        raster = F,
        group.by = "subcluster",
        split.by = "subcluster",
        ncol = 4)

FeaturePlot(clustered_obj_meyloid, 
            features = c(
              "nCount_Xenium",
              "nFeature_Xenium"
            ),
            reduction = "umap",
            raster=FALSE,
            ncol = 2
)

VlnPlot(clustered_obj_meyloid, 
        features = c( "nCount_Xenium",
                      "nFeature_Xenium"
        ),
        group.by = "subcluster",
        pt.size = 0,
        ncol = 2, raster=FALSE)

FeaturePlot(clustered_obj_meyloid, 
            features = c(
              "PTPRC"#, # Immune
              #"EPCAM", # Epithelial
              #"PECAM1" # Endothelial
            ),
            reduction = "umap",
            raster=FALSE,
            ncol = 1
)
VlnPlot(clustered_obj_meyloid, 
        features = c("PTPRC"#, # Immune
                     #"EPCAM", # Epithelial
                     #"PECAM1" # Endothelial
        ),
        group.by = "subcluster2",
        pt.size = 0,
        ncol = 1, raster=FALSE)


FeaturePlot(clustered_obj_meyloid, 
            features = c(
              "TTN", "DES", "FHL2", "S100A1"
            ),
            reduction = "umap",
            raster=FALSE,
            ncol = 2
)
VlnPlot(clustered_obj_meyloid, 
        features = c("TTN", "DES", "FHL2", "S100A1"
        ),
        group.by = "subcluster",
        pt.size = 0,
        ncol = 2, raster=FALSE)

DefaultAssay(clustered_obj_meyloid) <- "RNA" 
clustered_obj_meyloid <- NormalizeData(clustered_obj_meyloid)
VariableFeatures(clustered_obj_meyloid) <- rownames(clustered_obj_meyloid)
clustered_obj_meyloid <- ScaleData(clustered_obj_meyloid)

DotPlot(clustered_obj_meyloid, 
        features = c(#"PECAM1",
          "PTPRC",
          "CD163", "MRC1", "MARCO", "FCGR1A", # Macrophage
          "LYZ", "CD14", "FCGR3A",  "MS4A7", # Monocyte
          "FCER1A", #mo-DCs, cDC1, or pDC
          "CD1A", "CD1C", "CD209", "ITGAM", # monocyte derived DCs
          "ITGAX", "XCR1", # cDC1
          "GNLY", "NKG7", "KLRB1", #"KLRC1", #NK 
          "KLRD1", #NK
          "CD3E", "CD3D", "TRAC", "CD8A", "CD4", "FOXP3", # T
          "MS4A2", # mast
          "MS4A1", # B
          "CD79A", # B and plasma
          "TNFRSF17", "DERL3", # plasma
          "LILRA4", "CCR7", # pDCs
          "MKI67"
        ), 
        group.by = "subcluster", 
        scale = TRUE) + coord_flip()

DotPlot(clustered_obj_meyloid, 
        features = c(
          #"TSPY3", "PSG2", "RTKN2", "ELF5", "SOX2", "CELF4",
          "EPCAM", "PECAM1",
          "VWF", "BMX", # Ednocardial cells
          "MYH11",  # Vascular smooth muscle cells (VSMCs) 
          "PDGFRB", "ACTA2", # Pericytes
          "DCN", "C7", "FBLN1", "LTBP2", 
          "OGN", "PDGFRA", # Fibroblasts
          "ADIPOQ", # Adipocytes
          "TTN", "DES", "FHL2", "S100A1", # cardiomyocytes
          "PTPRC",
          "CD163", "MRC1", "MARCO", "FCGR1A", # Macrophage
          "LYZ", "CD14", "FCGR3A",  "MS4A7", # Monocyte
          "FCER1A", #mo-DCs, cDC1, or pDC
          "CD1A", "CD1C", "CD209", "ITGAM", # monocyte derived DCs
          "ITGAX", "XCR1", # cDC1
          "GNLY", "NKG7", "KLRB1", #"KLRC1", #NK 
          "KLRD1", #NK
          "CD3E", "CD3D", "TRAC", "CD8A", "CD4", "FOXP3", # T
          "MS4A2", # mast
          "MS4A1", # B
          "CD79A", # B and plasma
          "TNFRSF17", "DERL3", # plasma
          "LILRA4", "CCR7", # pDCs
          "MKI67"
        ), 
        group.by = "subcluster", 
        scale = TRUE) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))


markers <- presto::wilcoxauc(clustered_obj_meyloid, group_by = "subcluster", 
                             assay = "data", seurat_assay = "RNA")
topmarkers <- top_markers(markers)
topmarkers <- top_markers(markers, n = 10, auc_min = 0.6, pct_in_min = 50, 
                          pct_out_max = 100)


topmarkers


# Dendritic cell markers
# moDCs
FeaturePlot(clustered_obj_meyloid, 
            features = c(
              "FCER1A", #mo-DCs, cDC1, or pDC
              "CD1A", "CD1C",	"MRC1", "CD209", "ITGAM" # monocyte derived DCs
              
            ),
            reduction = "umap",
            raster=FALSE,
            ncol = 4 
)
VlnPlot(clustered_obj_meyloid, 
        features = c(
          "FCER1A", #mo-DCs, cDC1, or pDC
          "CD1A", "CD1C",	"MRC1", "CD209", "ITGAM" # monocyte derived DCs
          
        ),
        group.by = "subcluster",
        pt.size = 0,
        ncol = 3, raster=FALSE)

# DC1 and 2
FeaturePlot(clustered_obj_meyloid, 
            features = c(
              "CD8A", "ITGAX", "XCR1", # cDC1
              "CD1C", "ITGAM", # cDC2 or mo-DC
              "CLEC10A", "CX3CR1", # cDC2
              "LAMP3", "CD80" # mature DC
            ),
            reduction = "umap",
            raster=FALSE,
            ncol = 3
)
VlnPlot(clustered_obj_meyloid, 
        features = c(
          "CD8A", "ITGAX", "XCR1", # cDC1
          "CD1C", "ITGAM", # cDC2 or mo-DC
          "CLEC10A", "CX3CR1", # cDC2
          "LAMP3", "CD80" # mature DC
          
        ),
        group.by = "subcluster",
        pt.size = 0,
        ncol = 5, raster=FALSE)


FeaturePlot(clustered_obj_meyloid, 
            features = c(
              "CCR7", "CCL22"
            ),
            reduction = "umap",
            raster=FALSE,
            ncol = 1
)


# Macrophage
FeaturePlot(clustered_obj_meyloid, 
            features = c(
              "CD68", "CD163", # macrophage
              # CD14 could be expressed in macrophages but is a monocyte marker. 
              # If expresses other macrophage markers like CD68 and CD163, label 
              # as macrophagee 
              "CD14" 
            ),
            reduction = "umap",
            raster=FALSE,
            ncol = 3
)
VlnPlot(clustered_obj_meyloid, 
        features = c(
          "CD68", "CD163", # macrophage
          # CD14 could be expressed in macrophages but is a monocyte marker. 
          # If expresses other macrophage markers like CD68 and CD163, label 
          # as macrophagee 
          "CD14" 
        ),
        group.by = "subcluster",
        pt.size = 0,
        ncol = 3, raster=FALSE)

# More macrophage markers
FeaturePlot(clustered_obj_meyloid, 
            features = c(
              "MRC1", "MARCO", "FCGR1A", # Macrophage
              "LYZ", "FCGR3A",  "MS4A7" # Monocyte
            ),
            reduction = "umap",
            raster=FALSE,
            ncol = 3
)


# Mast
FeaturePlot(clustered_obj_meyloid, 
            features = c(
              "MS4A2"#, # mast
              #"SLC18A2", "VWA5A" # granulocyte
            ),
            reduction = "umap",
            raster=FALSE 
)
VlnPlot(clustered_obj_meyloid, 
        features = c(
          "MS4A2"#, # mast
          #"SLC18A2", "VWA5A" # granulocyte
        ),
        group.by = "subcluster",
        pt.size = 0, raster=FALSE)


# pDCs
FeaturePlot(clustered_obj_meyloid, 
            features = c(
              "LILRA4" # pDCs
            ),
            reduction = "umap",
            raster=FALSE,
            ncol = 1
)
VlnPlot(clustered_obj_meyloid, 
        features = c(
          "LILRA4"
        ),
        group.by = "subcluster",
        pt.size = 0,
        ncol = 1, raster=FALSE)


FeaturePlot(clustered_obj_meyloid, 
            features = c(
              "PECAM1", "VWF"
            ),
            reduction = "umap",
            raster=FALSE,
            ncol = 2
)


#### FINAL Myeloid object annotations ----
# subcluster
#Cardiomyocytes: 3

#Mast: 9
#pDC: 8
#Proliferating myeloid: 6,2
#Mature DC (mDCs): 6,1 
#cDC2: 5
#cDC1: 6,0

#Macrophages: 0, 1, 2, 4, 7, 10


### Lymphoid ----
# subcluster 2
# Lymphoid: 0, 12, 2, 3, 8, 6

cell_id_lymphoid <- clustered_obj_imm_pass2@meta.data %>% 
  filter(subcluster2 %in% c("0", "12", "2", "3", "8", "6")) %>% 
  dplyr::pull(cell_id)

length(cell_id_lymphoid)


merged_spatial_filtered_lymphoid <- subset(merged_spatial_filtered_filtered, 
                                           subset = cell_id %in% cell_id_lymphoid)
dim(merged_spatial_filtered_lymphoid)

saveRDS(merged_spatial_filtered_lymphoid, "/scratch/aoill/projects/heart_transplant/new/merged_spatial_filtered_lymphoid.rds")


clustered_obj_lymphoid <- read_rds("/scratch/aoill/projects/heart_transplant/new/clustered_obj_lymphoid_NC100_NN20_PC20_2024_06_13.rds")


# Fix cluster types to character
clustered_obj_lymphoid$subcluster <- as.character(clustered_obj_lymphoid$subcluster)
clustered_obj_lymphoid$leiden_res0.1 <- as.character(clustered_obj_lymphoid$leiden_res0.1)
clustered_obj_lymphoid$leiden_res0.2 <- as.character(clustered_obj_lymphoid$leiden_res0.2)
clustered_obj_lymphoid$leiden_res0.3 <- as.character(clustered_obj_lymphoid$leiden_res0.3)
clustered_obj_lymphoid$leiden_res0.5 <- as.character(clustered_obj_lymphoid$leiden_res0.5)
clustered_obj_lymphoid$leiden_res0.7 <- as.character(clustered_obj_lymphoid$leiden_res0.7)
clustered_obj_lymphoid$leiden_res0.75 <- as.character(clustered_obj_lymphoid$leiden_res0.75)
clustered_obj_lymphoid$leiden_res1.0 <- as.character(clustered_obj_lymphoid$leiden_res1.0)
clustered_obj_lymphoid$leiden_res1.5 <- as.character(clustered_obj_lymphoid$leiden_res1.5)
clustered_obj_lymphoid$leiden_res2.0 <- as.character(clustered_obj_lymphoid$leiden_res2.0)


DimPlot(clustered_obj_lymphoid, 
        reduction = "umap",
        raster = F,
        group.by = "subcluster",
        label = T
) + NoLegend()



DimPlot(clustered_obj_lymphoid, 
        reduction = "umap",
        raster = F,
        group.by = "subcluster",
        split.by = "subcluster",
        label = F,
        ncol = 4
) + NoLegend()


DefaultAssay(clustered_obj_lymphoid) <- "RNA" 
clustered_obj_lymphoid <- NormalizeData(clustered_obj_lymphoid)
VariableFeatures(clustered_obj_lymphoid) <- rownames(clustered_obj_lymphoid)
clustered_obj_lymphoid <- ScaleData(clustered_obj_lymphoid)

DotPlot(clustered_obj_lymphoid, 
        features = c(
          "PTPRC",
          "GNLY", "NKG7", "KLRB1", #"KLRC1", #NK 
          "KLRD1", #NK
          "CD3E", "CD3D", "CD8A", "CD4",  "FOXP3", # T
          # "IL7R", "CCR7", # T
          "MS4A1", # B
          "CD79A", # B and plasma
          "TNFRSF17", "DERL3", # plasma
          "LILRA4", # pDCs
          "MKI67",
          "CD163", "MRC1", "MARCO", "FCGR1A", # Macrophage
          "LYZ", "CD14", "FCGR3A",  "MS4A7", # Monocyte
          "FCER1A", "CD1A", "CD1C",	# DCs
          "MS4A2" # mast
          
          
        ), 
        group.by = "subcluster", 
        scale = TRUE) + coord_flip()



DotPlot(clustered_obj_lymphoid, 
        features = c(
          #"TSPY3", "PSG2", "RTKN2", "ELF5", "SOX2", "CELF4",
          "EPCAM", "PECAM1",
          "VWF", "BMX", # Ednocardial cells
          "MYH11",  # Vascular smooth muscle cells (VSMCs) 
          "PDGFRB", "ACTA2", # Pericytes
          "DCN", "C7", "FBLN1", "LTBP2", 
          "OGN", "PDGFRA", # Fibroblasts
          "ADIPOQ", # Adipocytes
          "TTN", "DES", "FHL2", "S100A1", # cardiomyocytes
          "PTPRC",
          "CD163", "MRC1", "MARCO", "FCGR1A", # Macrophage
          "LYZ", "CD14", "FCGR3A",  "MS4A7", # Monocyte
          "FCER1A", #mo-DCs, cDC1, or pDC
          "CD1A", "CD1C", "CD209", "ITGAM", # monocyte derived DCs
          "ITGAX", "XCR1", # cDC1
          "GNLY", "NKG7", "KLRB1", #"KLRC1", #NK 
          "KLRD1", #NK
          "CD3E", "CD3D", "TRAC", "CD8A", "CD4", "FOXP3", # T
          "MS4A2", # mast
          "MS4A1", # B
          "CD79A", # B and plasma
          "TNFRSF17", "DERL3", # plasma
          "LILRA4", "CCR7", # pDCs
          "MKI67"
        ), 
        group.by = "subcluster", 
        scale = TRUE) + coord_flip()

FeaturePlot(clustered_obj_lymphoid, 
            features = c(
              "TTN", "DES", "FHL2", "S100A1"
            ),
            reduction = "umap",
            raster=FALSE,
            ncol = 2
)
FeaturePlot(clustered_obj_lymphoid, 
            features = c(
              "CD68", "CD163", # macrophage
              # CD14 could be expressed in macrophages but is a monocyte marker. 
              # If expresses other macrophage markers like CD68 and CD163, label 
              # as macrophagee 
              "CD14" 
            ),
            reduction = "umap",
            raster=FALSE,
            ncol = 3
)

FeaturePlot(clustered_obj_lymphoid, 
            features = c(
              "CD68", "TRAC"
            ),
            reduction = "umap",
            raster=FALSE,
            blend = T
)
FeaturePlot(clustered_obj_lymphoid, 
            features = c(
              "TTN", "TRAC"
            ),
            reduction = "umap",
            raster=FALSE,
            blend = T
)

# T cells/NK
FeaturePlot(clustered_obj_lymphoid, 
            features = c(
              #"GNLY", "NKG7", "KLRB1", #"KLRC1", #NK 
              #"KLRD1", #NK
              "CD3E", "CD3D", "TRAC", "CD8A", "CD4", "IL7R", "CCR7", "FOXP3",  # T
              # "IL7R", "CCR7", # T
              "MKI67"
            ),
            reduction = "umap",
            raster=FALSE,
            ncol = 4
)
FeaturePlot(clustered_obj_lymphoid, 
            features = c(
              "GNLY", "NKG7", "KLRB1", #"KLRC1", #NK 
              "KLRD1", #NK
              "CD3E", "CD3D", "CD8A", "CD4", "IL7R", "CCR7", "FOXP3",  # T
              # "IL7R", "CCR7", # T
              "MKI67"
            ),
            reduction = "umap",
            raster=FALSE,
            ncol = 4
)
VlnPlot(clustered_obj_lymphoid, 
        features = c(
          "GNLY", "NKG7", "KLRB1", #"KLRC1", #NK 
          "KLRD1", #NK
          "CD3E", "CD3D", "CD8A", "CD4", "IL7R", "CCR7", "FOXP3",  # T
          # "IL7R", "CCR7", # T
          "MKI67"
          
        ),
        group.by = "subcluster",
        pt.size = 0,
        ncol = 4, raster=FALSE)

FeaturePlot(clustered_obj_lymphoid, 
            features = c(
              "FOXP3"
            ),
            reduction = "umap",
            raster=FALSE,
            split.by = "subcluster"
) + patchwork::plot_layout(ncol = 5, nrow = 2)

# B and plasma
FeaturePlot(clustered_obj_lymphoid, 
            features = c(
              "MS4A1", # B
              "CD79A", # B and plasma
              "TNFRSF17", "DERL3" # plasma
            ),
            reduction = "umap",
            raster=FALSE,
            ncol = 2 
)
VlnPlot(clustered_obj_lymphoid, 
        features = c(
          "MS4A1", # B
          "CD79A", # B and plasma
          "TNFRSF17", "DERL3" # plasma       
        ),
        group.by = "subcluster",
        pt.size = 0,
        ncol = 2, raster=FALSE)

# Top markers to see whats going on with 5
markers <- presto::wilcoxauc(clustered_obj_lymphoid, group_by = "subcluster", 
                             assay = "data", seurat_assay = "RNA")
topmarkers <- top_markers(markers)

topmarkers

#### FINAL Lymphoid object annotations (minus T cells) ----
# subcluster
#B cells: 6
#Plasma: 0
#NK: 7,0, 7,2


#### T cells ----
# Subset T cells and re-cluster
cell_id_T <- clustered_obj_lymphoid@meta.data %>% 
  filter(subcluster %in% c("1", "2", "3", "4", "5", "7,1")) %>% 
  dplyr::pull(cell_id)

length(cell_id_T)


merged_spatial_filtered_T <- subset(merged_spatial_filtered_filtered, 
                                    subset = cell_id %in% cell_id_T)
dim(merged_spatial_filtered_T)

saveRDS(merged_spatial_filtered_T, "/scratch/aoill/projects/heart_transplant/new/merged_spatial_filtered_T.rds")


clustered_obj_T <- read_rds("/scratch/aoill/projects/heart_transplant/new/clustered_obj_T_NC100_NN20_PC20_2024_06_13.rds")


# Fix cluster types to character
clustered_obj_T$subcluster <- as.character(clustered_obj_T$subcluster)
clustered_obj_T$leiden_res0.1 <- as.character(clustered_obj_T$leiden_res0.1)
clustered_obj_T$leiden_res0.2 <- as.character(clustered_obj_T$leiden_res0.2)
clustered_obj_T$leiden_res0.3 <- as.character(clustered_obj_T$leiden_res0.3)
clustered_obj_T$leiden_res0.5 <- as.character(clustered_obj_T$leiden_res0.5)
clustered_obj_T$leiden_res0.7 <- as.character(clustered_obj_T$leiden_res0.7)
clustered_obj_T$leiden_res0.75 <- as.character(clustered_obj_T$leiden_res0.75)
clustered_obj_T$leiden_res1.0 <- as.character(clustered_obj_T$leiden_res1.0)
clustered_obj_T$leiden_res1.5 <- as.character(clustered_obj_T$leiden_res1.5)
clustered_obj_T$leiden_res2.0 <- as.character(clustered_obj_T$leiden_res2.0)


DimPlot(clustered_obj_T, 
        reduction = "umap",
        raster = F,
        group.by = "subcluster",
        label = T
) + NoLegend()



DimPlot(clustered_obj_T, 
        reduction = "umap",
        raster = F,
        group.by = "subcluster",
        split.by = "subcluster",
        label = F,
        ncol = 4
) + NoLegend()


DefaultAssay(clustered_obj_T) <- "RNA" 
clustered_obj_T <- NormalizeData(clustered_obj_T)
VariableFeatures(clustered_obj_T) <- rownames(clustered_obj_T)
clustered_obj_T <- ScaleData(clustered_obj_T)

DotPlot(clustered_obj_T, 
        features = c(
          "PTPRC",
          "GNLY", "NKG7", "KLRB1", #"KLRC1", #NK 
          "KLRD1", #NK
          "CD3E", "CD3D", "CD8A", "CD4",  "FOXP3", # T
          # "IL7R", "CCR7", # T
          "MS4A1", # B
          "CD79A", # B and plasma
          "TNFRSF17", "DERL3", # plasma
          "LILRA4", # pDCs
          "MKI67",
          "CD163", "MRC1", "MARCO", "FCGR1A", # Macrophage
          "LYZ", "CD14", "FCGR3A",  "MS4A7", # Monocyte
          "FCER1A", "CD1A", "CD1C",	# DCs
          "MS4A2" # mast
          
          
        ), 
        group.by = "subcluster", 
        scale = TRUE) + coord_flip()



DotPlot(clustered_obj_T, 
        features = c(
          #"TSPY3", "PSG2", "RTKN2", "ELF5", "SOX2", "CELF4",
          "EPCAM", "PECAM1",
          "VWF", "BMX", # Ednocardial cells
          "MYH11",  # Vascular smooth muscle cells (VSMCs) 
          "PDGFRB", "ACTA2", # Pericytes
          "DCN", "C7", "FBLN1", "LTBP2", 
          "OGN", "PDGFRA", # Fibroblasts
          "ADIPOQ", # Adipocytes
          "TTN", "DES", "FHL2", "S100A1", # cardiomyocytes
          "PTPRC",
          "CD163", "MRC1", "MARCO", "FCGR1A", # Macrophage
          "LYZ", "CD14", "FCGR3A",  "MS4A7", # Monocyte
          "FCER1A", #mo-DCs, cDC1, or pDC
          "CD1A", "CD1C", "CD209", "ITGAM", # monocyte derived DCs
          "ITGAX", "XCR1", # cDC1
          "GNLY", "NKG7", "KLRB1", #"KLRC1", #NK 
          "KLRD1", #NK
          "CD3E", "CD3D", "TRAC", "CD8A", "CD4", "FOXP3", # T
          "MS4A2", # mast
          "MS4A1", # B
          "CD79A", # B and plasma
          "TNFRSF17", "DERL3", # plasma
          "LILRA4", "CCR7", # pDCs
          "MKI67"
        ), 
        group.by = "subcluster", 
        scale = TRUE) + coord_flip()

FeaturePlot(clustered_obj_T, 
            features = c(
              #"TTN", "DES", "FHL2", "S100A1"
              "CD163", "MRC1", "MARCO", "FCGR1A", # Macrophage
              "LYZ", "CD14", "FCGR3A",  "MS4A7" # Monocyte
            ),
            reduction = "umap",
            raster=FALSE,
            ncol = 4
)
VlnPlot(clustered_obj_T, 
        features = c(
          "CD163", "MRC1", "MARCO", "FCGR1A", # Macrophage
          "LYZ", "CD14", "FCGR3A",  "MS4A7" # Monocyte
          
        ),
        group.by = "subcluster",
        pt.size = 0,
        ncol = 4, raster=FALSE)


# T cells/NK
FeaturePlot(clustered_obj_T, 
            features = c(
              #"GNLY", "NKG7", "KLRB1", #"KLRC1", #NK 
              #"KLRD1", #NK
              "CD3E", "CD3D", "TRAC", "CD8A", "CD4", "IL7R", "CCR7", "FOXP3",  # T
              # "IL7R", "CCR7", # T
              "MKI67"
            ),
            reduction = "umap",
            raster=FALSE,
            ncol = 4
)
VlnPlot(clustered_obj_T, 
        features = c(
          "CD3E", "CD3D", "TRAC", "CD8A", "CD4", "IL7R", "CCR7", "FOXP3",  # T
          # "IL7R", "CCR7", # T
          "MKI67"
          
        ),
        group.by = "subcluster",
        pt.size = 0,
        ncol = 4, raster=FALSE)


FeaturePlot(clustered_obj_T, 
            features = c(
              "GNLY", "NKG7", "KLRB1", #"KLRC1", #NK 
              "KLRD1" #NK
              #"CD3E", "CD3D", "CD8A", "CD4", "IL7R", "CCR7", "FOXP3",  # T
              # "IL7R", "CCR7", # T
              #"MKI67"
            ),
            reduction = "umap",
            raster=FALSE,
            ncol = 4
)
VlnPlot(clustered_obj_T, 
        features = c(
          "GNLY", "NKG7", "KLRB1", #"KLRC1", #NK 
          "KLRD1", #NK
          "CD3E", "CD3D", "CD8A", "CD4", "IL7R", "CCR7", "FOXP3",  # T
          # "IL7R", "CCR7", # T
          "MKI67"
          
        ),
        group.by = "subcluster",
        pt.size = 0,
        ncol = 4, raster=FALSE)


markers <- presto::wilcoxauc(clustered_obj_T, group_by = "subcluster", 
                             assay = "data", seurat_assay = "RNA")
topmarkers <- top_markers(markers)

topmarkers

##### FINAL T cell object annotations ----
# subcluster
#Treg: 5,2
#Proliferating T cells: 2
#CD4+ T cells: 0, 1,0, 1,1, 1,2, 5,0, 5,1 
#CD8+ T cells: 1,3, 1,4, 3, 4, 6


## Mesenchymal ----
### Get cells from all objects ----

# Immune object
# Mesenchymal 9, 7, 0
#DimPlot(clustered_obj_imm, reduction = "umap",
#        group.by = "leiden_res0.7")
cell_id_mes1 <- clustered_obj_imm@meta.data %>% 
  filter(leiden_res0.7 %in% c("0", "7", "9")) %>% 
  dplyr::pull(cell_id)


# Mesenchymal object
# 0, 1, 2, 3, 5, 6
#DimPlot(clustered_obj_mes, reduction = "umap",
#        group.by = "leiden_res0.5")
cell_id_mes2 <- clustered_obj_mes@meta.data %>% 
  filter(leiden_res0.5 %in% c("0", "1", "2", "3", "5", "6")) %>% 
  dplyr::pull(cell_id)

# Endothelial object
# 1,0, 1,1
#DimPlot(clustered_obj_endo, reduction = "umap",
#        group.by = "subcluster2", label = T)
cell_id_mes3 <- clustered_obj_endo@meta.data %>% 
  filter(subcluster2 %in% c("1,0", "1,1")) %>% 
  dplyr::pull(cell_id)


length(cell_id_mes1)
length(cell_id_mes2)
length(cell_id_mes3)

cell_id_mes_pass2 <- c(cell_id_mes1, cell_id_mes2, cell_id_mes3)
length(cell_id_mes_pass2)

merged_spatial_filtered_mes_pass2 <- subset(merged_spatial_filtered_filtered, 
                                            subset = cell_id %in% cell_id_mes_pass2)
dim(merged_spatial_filtered_mes_pass2)

#saveRDS(merged_spatial_filtered_mes_pass2, "/scratch/aoill/projects/heart_transplant/new/merged_spatial_filtered_mes_pass2_ALT.rds")
saveRDS(merged_spatial_filtered_mes_pass2, "/scratch/aoill/projects/heart_transplant/new/merged_spatial_filtered_mes_pass2.rds")


#clustered_obj_mes_pass2 <- read_rds("/scratch/aoill/projects/heart_transplant/new/clustered_obj_mes_pass2_NC100_NN20_PC20_2024_06_13_ALT.rds")
clustered_obj_mes_pass2 <- read_rds("/scratch/aoill/projects/heart_transplant/new/clustered_obj_mes_pass2_NC100_NN20_PC20_2024_06_13.rds")

# Fix cluster types to character
clustered_obj_mes_pass2$subcluster <- as.character(clustered_obj_mes_pass2$subcluster)
clustered_obj_mes_pass2$leiden_res0.1 <- as.character(clustered_obj_mes_pass2$leiden_res0.1)
clustered_obj_mes_pass2$leiden_res0.2 <- as.character(clustered_obj_mes_pass2$leiden_res0.2)
clustered_obj_mes_pass2$leiden_res0.3 <- as.character(clustered_obj_mes_pass2$leiden_res0.3)
clustered_obj_mes_pass2$subcluster <- as.character(clustered_obj_mes_pass2$subcluster)
clustered_obj_mes_pass2$leiden_res1.0 <- as.character(clustered_obj_mes_pass2$leiden_res1.0)
clustered_obj_mes_pass2$leiden_res1.5 <- as.character(clustered_obj_mes_pass2$leiden_res1.5)
clustered_obj_mes_pass2$leiden_res2.0 <- as.character(clustered_obj_mes_pass2$leiden_res2.0)


DimPlot(clustered_obj_mes_pass2, 
        reduction = "umap",
        raster = F,
        group.by = "subcluster",
        label = T
) + NoLegend()

DimPlot(clustered_obj_mes_pass2, 
        reduction = "umap",
        raster = F,
        group.by = "subcluster",
        split.by = "subcluster",
        ncol = 4)

FeaturePlot(clustered_obj_mes_pass2, 
            features= c(
              "PTPRC", "LYZ", "CD3E" # Immune
              #"EPCAM", # Epithelial
              #"PECAM1" # Endothelial
            ),
            reduction = "umap",
            raster=FALSE,
            ncol = 3
)
VlnPlot(clustered_obj_mes_pass2, 
        features = c("PTPRC", "LYZ", "CD3E", # Immune
                     #"EPCAM", # Epithelial
                     #"PECAM1" # Endothelial
                     "COL1A1", "COL1A2"
        ),
        group.by = "subcluster",
        pt.size = 0,
        ncol = 3, raster=FALSE)




markers <- presto::wilcoxauc(clustered_obj_mes_pass2, group_by = "subcluster", 
                             assay = "data", seurat_assay = "RNA")
topmarkers <- top_markers(markers, n = 10, auc_min = 0.6, pct_in_min = 50, 
                          pct_out_max = 100)


topmarkers

topmarkers <- top_markers(markers)
topmarkers



DefaultAssay(clustered_obj_mes_pass2) <- "RNA" 
clustered_obj_mes_pass2 <- NormalizeData(clustered_obj_mes_pass2)
VariableFeatures(clustered_obj_mes_pass2) <- rownames(clustered_obj_mes_pass2)
clustered_obj_mes_pass2 <- ScaleData(clustered_obj_mes_pass2)


DotPlot(clustered_obj_mes_pass2, 
        features = c(
          "EPCAM", "PECAM1",
          "VWF", "BMX", # Ednocardial cells
          "MYH11",  # Vascular smooth muscle cells (VSMCs) 
          "PDGFRB", "ACTA2", # Pericytes
          "DCN", "C7", "FBLN1", "LTBP2", 
          "OGN", "PDGFRA", # Fibroblasts
          "ADIPOQ", # Adipocytes
          "TTN", "DES", "FHL2", "S100A1", # cardiomyocytes
          "PTPRC",
          "CD163", "MRC1", "MARCO", "FCGR1A", # Macrophage
          "LYZ", "CD14", "FCGR3A",  "MS4A7", # Monocyte
          "FCER1A", #mo-DCs, cDC1, or pDC
          "CD1A", "CD1C", "CD209", "ITGAM", # monocyte derived DCs
          "ITGAX", "XCR1", # cDC1
          "GNLY", "NKG7", "KLRB1", #"KLRC1", #NK 
          "KLRD1", #NK
          "CD3E", "CD3D", "TRAC", "CD8A", "CD4", "FOXP3", # T
          "MS4A2", # mast
          "MS4A1", # B
          "CD79A", # B and plasma
          "TNFRSF17", "DERL3", # plasma
          "LILRA4", "CCR7", # pDCs
          "MKI67"
        ), 
        group.by = "subcluster", 
        scale = TRUE) + coord_flip()


FeaturePlot(clustered_obj_mes_pass2, 
            features = c(
              "PTPRC",
              "CD68", "CD163", # macrophage
              # CD14 could be expressed in macrophages but is a monocyte marker. 
              # If expresses other macrophage markers like CD68 and CD163, label 
              # as macrophagee 
              "CD14" 
            ),
            reduction = "umap",
            raster=FALSE,
            ncol = 2
)

FeaturePlot(clustered_obj_mes_pass2, 
            features = c(
              "DES", "TTN", "FHL2", "S100A1"
            ),
            reduction = "umap",
            raster=FALSE,
            ncol = 2
)
FeaturePlot(clustered_obj_mes_pass2, 
            features = c(
              "PECAM1", "VWF"
            ),
            reduction = "umap",
            raster=FALSE,
            ncol = 2
)
FeaturePlot(clustered_obj_mes_pass2, 
            features = c(
              "LYZ"
            ),
            reduction = "umap",
            raster=FALSE,
            ncol = 1
)


# Fibroblasts: 
FeaturePlot(clustered_obj_mes_pass2, 
            features = c(
              "SEMA3C", # myofibroblast
              #"PMP22", 
              #"FAP", 
              "DCN", "C7", "FBLN1", "LTBP2", 
              "OGN", "PDGFRA" # Fibroblasts
              #"SFRP4", "MFAP5", # Adventicial Fibroblasts
              #"SFRP2", "PTGDS", # activated
              #"ACTA2" # could be expressed in myofibroblast
              
            ),
            reduction = "umap",
            raster=FALSE,
            ncol = 4
)
VlnPlot(clustered_obj_mes_pass2, 
        features = c(
          "SEMA3C", # myofibroblast
          #"PMP22", 
          #"FAP", 
          "DCN", "C7", "FBLN1", "LTBP2", 
          "OGN", "PDGFRA" # Fibroblasts
          #"SFRP4", "MFAP5", # Adventicial Fibroblasts
          #"SFRP2", "PTGDS", # activated
          #"ACTA2" # could be expressed in myofibroblast
        ),
        group.by = "subcluster",
        pt.size = 0,
        ncol = 4, raster=FALSE)


FeaturePlot(clustered_obj_mes_pass2, 
            features = c(
              "SEMA3C", # myofibroblast
              "SFRP4", "MFAP5", # Adventicial Fibroblasts
              "SFRP2", "PTGDS" # activated
              #"ACTA2" # could be expressed in myofibroblast
              
            ),
            reduction = "umap",
            raster=FALSE,
            ncol = 3
)


FeaturePlot(clustered_obj_mes_pass2, 
            features = c(
              "ACTA2", "MYH11",  # Vascular smooth muscle cells (VSMCs) 
              "PDGFRB", "ACTA2"#, # Pericytes
              #"ADIPOQ" # Adipocytes
              #"VWF", "BMX" # Ednocardial cells
              
            ),
            reduction = "umap",
            raster=FALSE,
            ncol = 2
)
VlnPlot(clustered_obj_mes_pass2, 
        features = c(
          "ACTA2", "MYH11",  # Vascular smooth muscle cells (VSMCs) 
          "PDGFRB"#, "ACTA2"#, # Pericytes
        ),
        group.by = "subcluster",
        pt.size = 0,
        ncol = 3, raster=FALSE)
VlnPlot(clustered_obj_mes_pass2, 
        features = c(
          "PECAM1", "VWF"
        ),
        group.by = "subcluster",
        pt.size = 0,
        ncol = 2, raster=FALSE)

# Cardiomyocytes
#"ASCL3", "MTRNR2L11", "MYBPC1", "PPP1R1A", "SH2D3C", # these don't mark cardiomyocytes well
# "DES" # Smooth muscle cells, Heart - Cardiomyocytes
FeaturePlot(clustered_obj_mes_pass2, 
            features = c(
              # Cardiomyocytes
              "TTN", "FHL2",
              # Below are good, featureplot just gets fussy when plotting a lot
              # of genes
              "PROX1", "PTGDS", "S100A1",
              "PEBP4", "PLIN4", "SMYD2"
              #"DES", # Smooth muscle cells, Heart - Cardiomyocytes
              
            ),
            reduction = "umap",
            raster=FALSE,
            ncol = 4
)


# Activated fibroblasts/myofibroblasts
# ACTA2, FSP1/S100A4 TNC
FeaturePlot(clustered_obj_mes_pass2, 
            features = c(
              "ACTA2", "FSP1", "S100A4", "TNC"
            ),
            reduction = "umap",
            raster=FALSE,
            ncol = 2
)

# Cardiac Pericytes 
FeaturePlot(clustered_obj_mes_pass2, 
            features = c(
              "PDGFRB", "NG2", "CSPG4", "TNC"
            ),
            reduction = "umap",
            raster=FALSE,
            ncol = 2
)


# Cardiac fibroblasts:
FeaturePlot(clustered_obj_mes_pass2, 
            features = c(
              #"COL1A1", "COL3A1", "FN1", "POSTN", "DCN", "SPARC", "TGFB1"
              "COL1A1", "DCN", "POSTN"
            ),
            reduction = "umap",
            raster=FALSE,
            ncol = 3
)

VlnPlot(clustered_obj_mes_pass2, 
        features = c(
          "COL1A1", "DCN", "POSTN"
        ),
        group.by = "subcluster",
        pt.size = 0,
        ncol = 3, raster=FALSE)


# Myofibroblasts (activated fibroblasts):
FeaturePlot(clustered_obj_mes2, 
            features = c(
              "ACTA2", "TAGLN", "VIM", "FAP", "DDR2", "CTGF"
            ),
            reduction = "umap",
            raster=FALSE,
            ncol = 2
)


#Perivascular fibroblasts:
FeaturePlot(clustered_obj_mes_pass2, 
            features = c(
              "PDGFRB", "NG2", "RGS5", "DES", "RERG"
            ),
            reduction = "umap",
            raster=FALSE,
            ncol = 2
)


# Endocardial fibroblasts:
FeaturePlot(clustered_obj_mes_pass2, 
            features = c(
              "NPPA", "NPPB", "CAV1", "TNNI3", "CHNG5"
            ),
            reduction = "umap",
            raster=FALSE,
            ncol = 2
)

# Interstitial fibroblasts:
FeaturePlot(clustered_obj_mes_pass2, 
            features = c(
              "PDGFRA", "DDR1", "CD90", "CD29", "CD105"
            ),
            reduction = "umap",
            raster=FALSE,
            ncol = 2
)

### FINAL Mesenchymal object annotations ----
# subcluster
#Fibroblast associated CD8+ T cells: 0
#Cardiomyocyte: 4,0, 4,3
#Endothelial: 3,0
#Macrophage: 1

#vSMCs: 5, 3,1
#Fibroblasts: 2,0, 2,2, 2,3
#POSTN+ Fibroblasts: 2,1 
#Proliferating pericytes: 4,1  
#Myofibroblasts: 6
#Pericytes: 4,2, 4,4, 4,5, 4,6


## Endothelial ----
### Get cells from all objects ----

# Immune object (NONE)

# Mesenchymal object (NONE)

# Endothelial object
# Putative cell types
# Endocardial cells: 3,0 
# Lymphatic endothelial cells: 3,1
# Proliferating endothelial cells: 2,1
# Endothelial cells: 1,2, 1,3, 4, 2,0, 
# But will re-cluster and umap in pass 2 then go down to these groups

#DimPlot(clustered_obj_endo, reduction = "umap",
#        group.by = "leiden_res0.5", label = T)
cell_id_endo_pass2 <- clustered_obj_endo@meta.data %>% 
  filter(leiden_res0.5 %in% c("3,0", "3,1", "2,1", "1,2", "1,3", "4", "2,0")) %>% 
  dplyr::pull(cell_id)


length(cell_id_endo_pass2)



merged_spatial_filtered_endo_pass2 <- subset(merged_spatial_filtered_filtered, 
                                             subset = cell_id %in% cell_id_endo_pass2)
dim(merged_spatial_filtered_endo_pass2)

saveRDS(merged_spatial_filtered_endo_pass2, "/scratch/aoill/projects/heart_transplant/new/merged_spatial_filtered_endo_pass2.rds")


clustered_obj_endo_pass2 <- read_rds("/scratch/aoill/projects/heart_transplant/new/clustered_obj_endo_pass2_NC100_NN20_PC20_2024_06_13.rds")


# Fix cluster types to character
clustered_obj_endo_pass2$subcluster <- as.character(clustered_obj_endo_pass2$subcluster)
clustered_obj_endo_pass2$leiden_res0.1 <- as.character(clustered_obj_endo_pass2$leiden_res0.1)
clustered_obj_endo_pass2$leiden_res0.2 <- as.character(clustered_obj_endo_pass2$leiden_res0.2)
clustered_obj_endo_pass2$leiden_res0.3 <- as.character(clustered_obj_endo_pass2$leiden_res0.3)
clustered_obj_endo_pass2$leiden_res0.5 <- as.character(clustered_obj_endo_pass2$leiden_res0.5)
clustered_obj_endo_pass2$leiden_res1.0 <- as.character(clustered_obj_endo_pass2$leiden_res1.0)
clustered_obj_endo_pass2$leiden_res1.5 <- as.character(clustered_obj_endo_pass2$leiden_res1.5)
clustered_obj_endo_pass2$leiden_res2.0 <- as.character(clustered_obj_endo_pass2$leiden_res2.0)


DimPlot(clustered_obj_endo_pass2, 
        reduction = "umap",
        raster = F,
        group.by = "subcluster",
        label = T
) + NoLegend()

DimPlot(clustered_obj_endo_pass2, 
        reduction = "umap",
        raster = F,
        group.by = "leiden_res0.5",
        label = T
) + NoLegend()
#LabelClusters(p, id = "ident",  fontface = "bold", color = "red")

DimPlot(clustered_obj_endo_pass2, 
        reduction = "umap",
        raster = F,
        group.by = "leiden_res0.5",
        split.by = "leiden_res0.5",
        ncol = 4)


markers <- presto::wilcoxauc(clustered_obj_endo_pass2, group_by = "leiden_res0.5", 
                             assay = "data", seurat_assay = "RNA")
topmarkers <- top_markers(markers, n = 10, auc_min = 0.6, pct_in_min = 50, 
                          pct_out_max = 100)


topmarkers


topmarkers <- top_markers(markers)
topmarkers

DefaultAssay(clustered_obj_endo_pass2) <- "RNA" 
clustered_obj_endo_pass2 <- NormalizeData(clustered_obj_endo_pass2)
VariableFeatures(clustered_obj_endo_pass2) <- rownames(clustered_obj_endo_pass2)
clustered_obj_endo_pass2 <- ScaleData(clustered_obj_endo_pass2)

DotPlot(clustered_obj_endo_pass2, 
        features = c(
          "EPCAM", "PECAM1",
          "VWF", "BMX", # Ednocardial cells
          "MYH11",  # Vascular smooth muscle cells (VSMCs) 
          "PDGFRB", "ACTA2", # Pericytes
          "DCN", "C7", "FBLN1", "LTBP2", 
          "OGN", "PDGFRA", # Fibroblasts
          "ADIPOQ", # Adipocytes
          "TTN", "DES", "FHL2", "S100A1", # cardiomyocytes
          "PTPRC",
          #"CD163", "MRC1", "MARCO", "FCGR1A", # Macrophage
          #"LYZ", "CD14", "FCGR3A",  "MS4A7", # Monocyte
          #"FCER1A", #mo-DCs, cDC1, or pDC
          #"CD1A", "CD1C", "CD209", "ITGAM", # monocyte derived DCs
          #"ITGAX", "XCR1", # cDC1
          #"GNLY", "NKG7", "KLRB1", #"KLRC1", #NK 
          #"KLRD1", #NK
          #"CD3E", "CD3D", "TRAC", "CD8A", "CD4", "FOXP3", # T
          #"MS4A2", # mast
          #"MS4A1", # B
          #"CD79A", # B and plasma
          #"TNFRSF17", "DERL3", # plasma
          #"LILRA4", "CCR7", # pDCs
          "MKI67"
        ), 
        group.by = "leiden_res0.5", 
        scale = TRUE) + coord_flip()



FeaturePlot(clustered_obj_endo_pass2, 
            features = c(
              "DCN", "C7", "FBLN1", "LTBP2", 
              "OGN", "PDGFRA", # Fibroblasts
              "CHST15"
            ),
            reduction = "umap",
            raster=FALSE,
            ncol = 3
)

FeaturePlot(clustered_obj_endo_pass2, 
            features = c(
              "PECAM1", "PTPRC"
            ),
            reduction = "umap",
            raster=FALSE,
            ncol = 1
)
VlnPlot(clustered_obj_endo_pass2, 
        features = c(
          "PECAM1"
        ),
        group.by = "leiden_res0.5",
        pt.size = 0,
        ncol = 1, raster=FALSE)

FeaturePlot(clustered_obj_endo_pass2, 
            features = c(
              "PECAM1", "VWF"
            ),
            reduction = "umap",
            raster=FALSE,
            ncol = 2
)

FeaturePlot(clustered_obj_endo_pass2, 
            features = c(
              "PECAM1", "PTPRC", "LYZ"
            ),
            reduction = "umap",
            raster=FALSE,
            ncol = 2
)
VlnPlot(clustered_obj_endo_pass2, 
        features = c(
          "PECAM1", "PTPRC", "LYZ"
        ),
        group.by = "leiden_res0.5",
        pt.size = 0,
        ncol = 3, raster=FALSE)


FeaturePlot(clustered_obj_endo_pass2, 
            features = c(
              #"DLL4", "NOTCH4", "EPHB4", "NRP1",
              #"CDH5", "EPHB2", "EFNB2", "HEY2",
              "DES", "TTN", "FHL2", "S100A1"
              #"MYBPC1" # cardiomyocytes
              #"ABCC11"
            ),
            reduction = "umap",
            raster=FALSE,
            ncol = 2
)

VlnPlot(clustered_obj_endo_pass2, 
        features = c(
          "PECAM1", "VWF",
          "DES", "TTN", "FHL2", "S100A1"
        ),
        group.by = "subcluster",
        pt.size = 0,
        ncol = 3, raster=FALSE)

FeaturePlot(clustered_obj_endo_pass2, 
            features = c(
              "DCN", "C7", "FBLN1", "LTBP2", 
              "OGN", "PDGFRA"# Fibroblasts
            ),
            reduction = "umap",
            raster=FALSE,
            ncol = 3
)

FeaturePlot(clustered_obj_endo_pass2, 
            features = c(
              "MKI67"
            ),
            reduction = "umap",
            raster=FALSE,
            ncol = 1
)

# fibroblasts
#POSTN
FeaturePlot(clustered_obj_endo_pass2, 
            features = c(
              "POSTN", "TFPI", "SELP", "TMEM100"
            ),
            reduction = "umap",
            raster=FALSE,
            ncol = 4
)

# cluster 2 top
FeaturePlot(clustered_obj_endo_pass2, 
            features = c(
              "CXCL11","CXCL9",  "FKBP1A", "UBD", 
              "CTSS",   "VCAM1",  "PLA1A",  "TAP1",   "ICAM1",  "CXCL10"
            ),
            reduction = "umap",
            raster=FALSE,
            ncol = 5
)
# Venous Endothelial Cells	
# COUP-TFII (NR2F2), VEGFR3 (Vascular endothelial growth factor receptor 3), CLDN5 (Claudin-5), NRP2 (Neuropilin-2)
FeaturePlot(clustered_obj_endo_pass2, 
            features = c(
              "NR2F2", "VEGFR3", "CLDN5", "NRP2"
            ),
            reduction = "umap",
            raster=FALSE,
            ncol = 4
)


FeaturePlot(clustered_obj_endo_pass2, 
            features = c(
              "PROX1", "LYVE1", "PDPN", "VEGFR3"
            ),
            reduction = "umap",
            raster=FALSE,
            ncol = 2
)

VlnPlot(clustered_obj_endo_pass2, 
        features = c(
          "PROX1", "LYVE1", "PDPN"
        ),
        group.by = "leiden_res0.5",
        pt.size = 0,
        ncol = 2, raster=FALSE)


FeaturePlot(clustered_obj_endo_pass2, 
            features = c(
              "CD34", "CDH5", "CLDN5", "VWF"
            ),
            reduction = "umap",
            raster=FALSE,
            ncol = 2
)


FeaturePlot(clustered_obj_endo_pass2, 
            features = c(
              "VCAM1", "CD31", "CDH5", "NRG1",
              "VWF", "PECAM1", "BMX", "NPR3", "TIE1", "FLK1", "ESM1", "THBD"
            ),
            reduction = "umap",
            raster=FALSE,
            ncol = 3
)


# Endocardial cells:
#VWF, PECAM1, BMX, NPR3
FeaturePlot(clustered_obj_endo_pass2, 
            features = c(
              "VWF", "PECAM1", "BMX", "NPR3",
              "POSTN", "SELP"
            ),
            reduction = "umap",
            raster=FALSE,
            ncol = 3
)
VlnPlot(clustered_obj_endo_pass2, 
        features = c(
          "VWF", "PECAM1", "BMX",
          "POSTN", "SELP" # SELP = activated
        ),
        group.by = "leiden_res0.5",
        pt.size = 0,
        ncol = 3, raster=FALSE)


### FINAL Endothelial object annotations ----
# subcluster
# Cardiomyocytes: 0,0

# Proliferating endothelial: 2
# Lymphatic endothelial: 5
# Activated endothelial: 4
# BMX+ Activated endothelial: 6
# Endothelial: 0,1, 0,2, 0,3, 1,0, 1,1, 3


# Add annotations to full object ----
## 1. Immune pass 2 object ----
# clustered_obj_imm_pass2
# subcluster2
# Endothelial: 13
# Fibroblasts: 7,0, 7,1
# Pericytes: 7,2
# These are final annotations for these clusters

cell_cluster_info1 <- clustered_obj_imm_pass2@meta.data %>%
  dplyr::select(X, subcluster2)

cell_cluster_info1$ct_first_pass <- "NA"

cell_cluster_info1$ct_first_pass[
  cell_cluster_info1$subcluster2 %in% c("13")] <- "Endothelial"

cell_cluster_info1$ct_first_pass[
  cell_cluster_info1$subcluster2 %in% c("7,0", "7,1")] <- "Fibroblasts"

cell_cluster_info1$ct_first_pass[
  cell_cluster_info1$subcluster2 %in% c("7,2")] <- "Pericytes"

table(cell_cluster_info1$ct_first_pass)

cell_cluster_info1_df <- cell_cluster_info1 %>% dplyr::select(X, ct_first_pass)
table(cell_cluster_info1_df$ct_first_pass)

cell_cluster_info1_df_final <- cell_cluster_info1_df %>%
  filter(ct_first_pass != "NA")
table(cell_cluster_info1_df_final$ct_first_pass)


## 2. Meyloid pass 2 object ----
# clustered_obj_meyloid
# subcluster
#Cardiomyocytes: 3
#Mast: 9
#pDC: 8
#Proliferating myeloid: 6,2
#Mature DC (mDCs): 6,1 
#cDC2: 5
#cDC1: 6,0
#Macrophages: 0, 1, 2, 4, 7
# NEW- Endothelial: 10

cell_cluster_info2 <- clustered_obj_meyloid@meta.data %>%
  dplyr::select(X, subcluster)

cell_cluster_info2$ct_first_pass <- "NA"

cell_cluster_info2$ct_first_pass[
  cell_cluster_info2$subcluster %in% c("3")] <- "Cardiomyocytes"

cell_cluster_info2$ct_first_pass[
  cell_cluster_info2$subcluster %in% c("9")] <- "Mast"

cell_cluster_info2$ct_first_pass[
  cell_cluster_info2$subcluster %in% c("8")] <- "pDC"

cell_cluster_info2$ct_first_pass[
  cell_cluster_info2$subcluster %in% c("6,2")] <- "Proliferating myeloid"

cell_cluster_info2$ct_first_pass[
  cell_cluster_info2$subcluster %in% c("6,1")] <- "mDC"

cell_cluster_info2$ct_first_pass[
  cell_cluster_info2$subcluster %in% c("5")] <- "cDC2"

cell_cluster_info2$ct_first_pass[
  cell_cluster_info2$subcluster %in% c("6,0")] <- "cDC1"

cell_cluster_info2$ct_first_pass[
  cell_cluster_info2$subcluster %in% c("0", "1", "2", "4", "7")] <- "Macrophages"

cell_cluster_info2$ct_first_pass[
  cell_cluster_info2$subcluster %in% c("10")] <- "Endothelial"

table(cell_cluster_info2$ct_first_pass)

cell_cluster_info2_df <- cell_cluster_info2 %>% dplyr::select(X, ct_first_pass)
table(cell_cluster_info2_df$ct_first_pass)

cell_cluster_info2_df_final <- cell_cluster_info2_df %>%
  filter(ct_first_pass != "NA")
table(cell_cluster_info2_df_final$ct_first_pass)

nrow(clustered_obj_meyloid@meta.data)
nrow(cell_cluster_info2_df_final)

## 3. Lymphoid pass 2 object ----
# clustered_obj_lymphoid
# subcluster
#B cells: 6
#Plasma: 0
#NK: 7,0, 7,2

cell_cluster_info3 <- clustered_obj_lymphoid@meta.data %>%
  dplyr::select(X, subcluster)

cell_cluster_info3$ct_first_pass <- "NA"

cell_cluster_info3$ct_first_pass[
  cell_cluster_info3$subcluster %in% c("6")] <- "B cells"

cell_cluster_info3$ct_first_pass[
  cell_cluster_info3$subcluster %in% c("0")] <- "Plasma"

cell_cluster_info3$ct_first_pass[
  cell_cluster_info3$subcluster %in% c("7,0", "7,2")] <- "NK"


table(cell_cluster_info3$ct_first_pass)

cell_cluster_info3_df <- cell_cluster_info3 %>% dplyr::select(X, ct_first_pass)
table(cell_cluster_info3_df$ct_first_pass)

cell_cluster_info3_df_final <- cell_cluster_info3_df %>%
  filter(ct_first_pass != "NA")
table(cell_cluster_info3_df_final$ct_first_pass)

nrow(clustered_obj_lymphoid@meta.data)
nrow(clustered_obj_T@meta.data)
nrow(cell_cluster_info3_df_final)
nrow(cell_cluster_info4_df_final)


## 4. T cell pass 2 object ----
# clustered_obj_T
# subcluster
#Treg: 5,2
#Proliferating T cells: 2
#CD4+ T cells: 0, 1,0, 1,1, 1,2, 5,0, 5,1 
#CD8+ T cells: 1,3, 1,4, 3, 4, 6

cell_cluster_info4 <- clustered_obj_T@meta.data %>%
  dplyr::select(X, subcluster)

cell_cluster_info4$ct_first_pass <- "NA"

cell_cluster_info4$ct_first_pass[
  cell_cluster_info4$subcluster %in% c("5,2")] <- "Treg"

cell_cluster_info4$ct_first_pass[
  cell_cluster_info4$subcluster %in% c("2")] <- "Proliferating T cells"

cell_cluster_info4$ct_first_pass[
  cell_cluster_info4$subcluster %in% c("0", "1,0", "1,1", "1,2", "5,0", "5,1")] <- "CD4+ T cells"

cell_cluster_info4$ct_first_pass[
  cell_cluster_info4$subcluster %in% c("1,3", "1,4", "3", "4", "6")] <- "CD8+ T cells"

table(cell_cluster_info4$ct_first_pass)

cell_cluster_info4_df <- cell_cluster_info4 %>% dplyr::select(X, ct_first_pass)
table(cell_cluster_info4_df$ct_first_pass)

cell_cluster_info4_df_final <- cell_cluster_info4_df %>%
  filter(ct_first_pass != "NA")
table(cell_cluster_info4_df_final$ct_first_pass)


## 5. Mesenchymal pass 2 object ----
# clustered_obj_mes_pass2
# subcluster

#Fibroblast associated CD8+ T cells: 0
#Cardiomyocytes: 4,0, 4,3
#Endothelial: 3,0
#Macrophages: 1 - label Macrophage (mes obj)

#vSMCs: 5, 3,1
#Fibroblasts: 2,0, 2,2, 2,3
#POSTN+ Fibroblasts: 2,1 
#Proliferating pericytes: 4,1  
#Myofibroblasts: 6
#Pericytes: 4,2, 4,4, 4,5, 4,6



cell_cluster_info5 <- clustered_obj_mes_pass2@meta.data %>%
  dplyr::select(X, subcluster)

cell_cluster_info5$ct_first_pass <- "NA"

cell_cluster_info5$ct_first_pass[
  cell_cluster_info5$subcluster %in% c("0")] <- "Fibroblast associated CD8+ T cells"

cell_cluster_info5$ct_first_pass[
  cell_cluster_info5$subcluster %in% c("4,0", "4,3")] <- "Cardiomyocytes"

cell_cluster_info5$ct_first_pass[
  cell_cluster_info5$subcluster %in% c("3,0")] <- "Endothelial"

cell_cluster_info5$ct_first_pass[
  cell_cluster_info5$subcluster %in% c("1")] <- "Macrophages (Mes obj)"

cell_cluster_info5$ct_first_pass[
  cell_cluster_info5$subcluster %in% c("5")] <- "vSMCs"

cell_cluster_info5$ct_first_pass[
  cell_cluster_info5$subcluster %in% c("3,1")] <- "vSMC associated endothelial"

cell_cluster_info5$ct_first_pass[
  cell_cluster_info5$subcluster %in% c("2,0", "2,2", "2,3")] <- "Fibroblasts"

cell_cluster_info5$ct_first_pass[
  cell_cluster_info5$subcluster %in% c("2,1")] <- "POSTN+ Fibroblasts"

cell_cluster_info5$ct_first_pass[
  cell_cluster_info5$subcluster %in% c("4,1")] <- "Proliferating pericytes"

cell_cluster_info5$ct_first_pass[
  cell_cluster_info5$subcluster %in% c("6")] <- "Myofibroblasts"

cell_cluster_info5$ct_first_pass[
  cell_cluster_info5$subcluster %in% c("4,2", "4,4", "4,5", "4,6")] <- "Pericytes"


table(cell_cluster_info5$ct_first_pass)

cell_cluster_info5_df <- cell_cluster_info5 %>% dplyr::select(X, ct_first_pass)
table(cell_cluster_info5_df$ct_first_pass)

cell_cluster_info5_df_final <- cell_cluster_info5_df %>%
  filter(ct_first_pass != "NA")
table(cell_cluster_info5_df_final$ct_first_pass)


## 6. Endothelial pass 2 object ----
# clustered_obj_endo_pass2
# subcluster
# Cardiomyocytes: 0,0

# Proliferating endothelial: 2
# Lymphatic endothelial: 5
# Activated endothelial: 4
# BMX+ Activated endothelial: 6
# Endothelial: 0,1, 0,2, 0,3, 1,0, 1,1, 3



cell_cluster_info6 <- clustered_obj_endo_pass2@meta.data %>%
  dplyr::select(X, subcluster)

cell_cluster_info6$ct_first_pass <- "NA"

cell_cluster_info6$ct_first_pass[
  cell_cluster_info6$subcluster %in% c("0,0")] <- "Cardiomyocytes"

cell_cluster_info6$ct_first_pass[
  cell_cluster_info6$subcluster %in% c("2")] <- "Proliferating endothelial"

cell_cluster_info6$ct_first_pass[
  cell_cluster_info6$subcluster %in% c("5")] <- "Lymphatic endothelial"

cell_cluster_info6$ct_first_pass[
  cell_cluster_info6$subcluster %in% c("4")] <- "Activated endothelial"

cell_cluster_info6$ct_first_pass[
  cell_cluster_info6$subcluster %in% c("6")] <- "BMX+ Activated endothelial"

cell_cluster_info6$ct_first_pass[
  cell_cluster_info6$subcluster %in% c("0,1", "0,2", "0,3", "1,0", "1,1", "3")] <- "Endothelial"


table(cell_cluster_info6$ct_first_pass)

cell_cluster_info6_df <- cell_cluster_info6 %>% dplyr::select(X, ct_first_pass)
table(cell_cluster_info6_df$ct_first_pass)

cell_cluster_info6_df_final <- cell_cluster_info6_df %>%
  filter(ct_first_pass != "NA")
table(cell_cluster_info6_df_final$ct_first_pass)


## 7. Cardiomyocytes object ----
# clustered_obj_card
# leiden_res0.4
# Cardiomyocytes: "0" "1" "3" "2"


cell_cluster_info7 <- clustered_obj_card@meta.data %>%
  dplyr::select(X, leiden_res0.4)

cell_cluster_info7$ct_first_pass <- "NA"

cell_cluster_info7$ct_first_pass[
  cell_cluster_info7$leiden_res0.4 %in% c("0", "1", "3", "2")] <- "Cardiomyocytes"



table(cell_cluster_info7$ct_first_pass)

cell_cluster_info7_df <- cell_cluster_info7 %>% dplyr::select(X, ct_first_pass)
table(cell_cluster_info7_df$ct_first_pass)

cell_cluster_info7_df_final <- cell_cluster_info7_df %>%
  filter(ct_first_pass != "NA")
table(cell_cluster_info7_df_final$ct_first_pass)


## 8. Cardiomyocytes from immune pass 1 object ----
# clustered_obj_imm
# leiden_res0.7
# Cardiomyocytes: "4"


cell_cluster_info8 <- clustered_obj_imm@meta.data %>%
  dplyr::select(X, leiden_res0.7)

cell_cluster_info8$ct_first_pass <- "NA"

cell_cluster_info8$ct_first_pass[
  cell_cluster_info8$leiden_res0.7 %in% c("4")] <- "Cardiomyocytes"



table(cell_cluster_info8$ct_first_pass)

cell_cluster_info8_df <- cell_cluster_info8 %>% dplyr::select(X, ct_first_pass)
table(cell_cluster_info8_df$ct_first_pass)

cell_cluster_info8_df_final <- cell_cluster_info8_df %>%
  filter(ct_first_pass != "NA")
table(cell_cluster_info8_df_final$ct_first_pass)


## 9. Cardiomyocytes from mesenchymal pass 1 object ----
# clustered_obj_mes
# leiden_res0.5
# Cardiomyocytes: "4"

cell_cluster_info9 <- clustered_obj_mes@meta.data %>%
  dplyr::select(X, leiden_res0.5)

cell_cluster_info9$ct_first_pass <- "NA"

cell_cluster_info9$ct_first_pass[
  cell_cluster_info9$leiden_res0.5 %in% c("4")] <- "Cardiomyocytes"



table(cell_cluster_info9$ct_first_pass)

cell_cluster_info9_df <- cell_cluster_info9 %>% dplyr::select(X, ct_first_pass)
table(cell_cluster_info9_df$ct_first_pass)

cell_cluster_info9_df_final <- cell_cluster_info9_df %>%
  filter(ct_first_pass != "NA")
table(cell_cluster_info9_df_final$ct_first_pass)


## 10. Cardiomyocytes from endothelial pass 1 object ----
# clustered_obj_endo
# subcluster2
# Cardiomyocytes: "0"

cell_cluster_info10 <- clustered_obj_endo@meta.data %>%
  dplyr::select(X, subcluster2)

cell_cluster_info10$ct_first_pass <- "NA"

cell_cluster_info10$ct_first_pass[
  cell_cluster_info10$subcluster2 %in% c("0")] <- "Cardiomyocytes"



table(cell_cluster_info10$ct_first_pass)

cell_cluster_info10_df <- cell_cluster_info10 %>% dplyr::select(X, ct_first_pass)
table(cell_cluster_info10_df$ct_first_pass)

cell_cluster_info10_df_final <- cell_cluster_info10_df %>%
  filter(ct_first_pass != "NA")
table(cell_cluster_info10_df_final$ct_first_pass)




## Merge all cell ids together ----
cell_cluster_info_all <- rbind(cell_cluster_info1_df_final, cell_cluster_info2_df_final, 
                               cell_cluster_info3_df_final, cell_cluster_info4_df_final,
                               cell_cluster_info5_df_final, cell_cluster_info6_df_final,
                               cell_cluster_info7_df_final, cell_cluster_info8_df_final,
                               cell_cluster_info9_df_final, cell_cluster_info10_df_final)

nrow(cell_cluster_info_all)

nrow(clustered_obj@meta.data)
# 164165


## join to full object ----
clustered_obj@meta.data$ct_first_pass <- NULL
clustered_obj_md <- clustered_obj@meta.data
nrow(clustered_obj_md)
#clustered_obj@meta.data <- clustered_obj_md
clustered_obj_md_all <- left_join(clustered_obj_md, cell_cluster_info_all)
nrow(clustered_obj_md_all)
#View(clustered_obj_md_all)
clustered_obj@meta.data$ct_first_pass <- clustered_obj_md_all$ct_first_pass


saveRDS(clustered_obj, "/scratch/aoill/projects/heart_transplant/new/clustered_obj_annotated_06_17_2024.rds")

# Plot results ----
# Define the number of colors you want
nb.cols <- length(unique(clustered_obj@meta.data$ct_first_pass))
#randomcoloR::distinctColorPalette(nb.cols)
mycolors <- randomcoloR::distinctColorPalette(nb.cols)
#mycolors_save2 <- mycolors
a <- DimPlot(clustered_obj, group.by = "ct_first_pass", 
             raster = F,
             reduction = "umap") +
  scale_color_manual(values = mycolors) + NoLegend() + coord_fixed()
pa <- LabelClusters(a, id = "ct_first_pass", box = TRUE, size = 3.5, label.size = 0.1, 
                    box.padding = 0.5, seed = 309,
                    #color = c(rep("black", 1), rep("white", 1), rep("white", 2), rep("black", 2), 
                    #          rep("white", 2), rep("black", 2), rep("white", 1), 
                    #          rep("black", 1)),
                    alpha = 0.8
)
pa

DimPlot(clustered_obj, group.by = "ct_first_pass", 
        split.by = "ct_first_pass", 
        raster = F,
        reduction = "umap",
        ncol = 7) +
  scale_color_manual(values = mycolors) + NoLegend()



DimPlot(clustered_obj, 
        group.by = "ct_first_pass", 
        split.by = "biopsy_timing", 
        raster = F,
        reduction = "umap",
        ncol = 2
) +
  scale_color_manual(values = mycolors) + NoLegend() 


# Dotplot heatmap ----
## Top markers ----
# Finding cluster markers for a heatmap/dotplot
markers <- presto::wilcoxauc(clustered_obj, group_by = "ct_first_pass")
markers <- top_markers(markers, n = 5, auc_min = 0.6, pct_in_min = 50, 
                       pct_out_max = 100)

markers 

all_markers <- markers %>%
  dplyr::select(-rank) %>% 
  unclass() %>% 
  stack() %>%
  pull(values) %>%
  unique() %>%
  .[!is.na(.)]

length(all_markers)


markers <- presto::wilcoxauc(clustered_obj, group_by = "ct_first_pass")
markers_2 <- top_markers(markers, n = 5)

#View(markers_2) 
all_markers <- markers_2 %>%
  dplyr::select(-rank) %>% 
  unclass() %>% 
  stack() %>%
  pull(values) %>%
  unique() %>%
  .[!is.na(.)]

length(all_markers)



# log-normalizing and scaling all features in the RNA assay. Scaling so that
# all features can be visualized using the same color scale
DefaultAssay(clustered_obj) <- "RNA"
clustered_obj <- NormalizeData(clustered_obj)
VariableFeatures(clustered_obj) <- rownames(clustered_obj)
clustered_obj <- ScaleData(clustered_obj)

p <- DotPlot(clustered_obj, features = all_markers, 
             group.by = "ct_first_pass", scale = TRUE) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
p

# Reproducing the dotplot using ComplexHeatmap
df <- p$data
head(df)

# Adding log-transformed values
#df$log.avg.exp <- log10(df$avg.exp)

# (Scaled) expression levels
exp_mat <- df %>% 
  dplyr::select(-pct.exp, -avg.exp) %>%  
  pivot_wider(names_from = id, values_from = avg.exp.scaled) %>% 
  as.data.frame()

row.names(exp_mat) <- exp_mat$features.plot
exp_mat <- exp_mat[,-1] %>% as.matrix()

# The percentage of cells expressing a feature
percent_mat <- df %>% 
  dplyr::select(-avg.exp, -avg.exp.scaled) %>%  
  pivot_wider(names_from = id, values_from = pct.exp) %>% 
  as.data.frame() 

row.names(percent_mat) <- percent_mat$features.plot  
percent_mat <- percent_mat[,-1] %>% as.matrix()

dim(exp_mat)
dim(percent_mat)

# Get an idea of the ranges of the matrix
quantile(exp_mat, c(0.1, 0.5, 0.9, 0.99))

# Any value that is greater than 2 will be mapped to yellow
col_fun = circlize::colorRamp2(c(-1, 0, 2), viridis(20)[c(1,10, 20)])

## Adding annotation colors
#mycolors_lst <- c("B cells" = "#AAE94E", "Cardiomyocytes" = "#D995DB", "CD8+ T cells" = "#E4DD3C",
#                  "cDC1" = "#67E8CA", "Endocardial cells" = "#8375D9", "Endothelial cells" = "#D0E7E3",
#                  "Fibroblast-like T cells" = "#B2EDAC", "Fibroblasts" = "#79D5E5", "Lymphatic ECs" = "#D652D1",
#                  "M2" = "#E6944E", "Macrophages" = "#D6618F", "Mast" = "#63EA7B", "mDC" = "#7BA2D8", 
#                  "Myofibroblasts" = "#913EE4", "NK" = "#E29F95", "pDC" = "#DCC5DF", 
#                  "Pericytes" = "#E85946", "Plasma"= "#B89E6C", "POSTN+ Fibroblasts" = "#83B4A7",
#                  "Proliferating CD8+ T cells" = "#95858D", "Proliferating pericytes" = "#DFDA7A",
#                  "Treg" = "#E3DCB5", "vSMCs" = "#64B76C")

cluster_anno <- colnames(exp_mat)
column_ha <- HeatmapAnnotation( # HeatmapAnnotation / rowAnnotation
  cluster = cluster_anno,
  #col = mycolors_lst,
  na_col = "lightgrey"
)

layer_fun = function(j, i, x, y, w, h, fill){
  grid.rect(x = x, y = y, width = w, height = h, 
            gp = gpar(col = NA, fill = NA))
  grid.points(x = x, y = y, 
              gp = gpar(col = col_fun(pindex((exp_mat), i, j))), # t(exp_mat)
              size = pindex((percent_mat), i, j)/100 * unit(4, "mm"), # t(percent_mat)
              pch = 16
  )
}

lgd_list = list(
  Legend( labels = c(0,0.25,0.5,0.75,1), title = "pt", type = "points", pch = 16, size = c(0,0.25,0.5,0.75,1) * unit(4,"mm"),
          legend_gp = gpar(col = "black")))

hp <- Heatmap((exp_mat), # t(exp_mat)
              heatmap_legend_param=list(title="Scaled expression"),
              column_title = "cluster_obj", 
              col = col_fun,
              rect_gp = gpar(type = "none"),
              layer_fun = layer_fun,
              row_names_gp = gpar(fontsize = 8),
              column_names_gp = gpar(fontsize = 8),
              #row_km = 4,
              column_km = 4,
              top_annotation = column_ha, # top_annotation / left_annotation
              border = "black")

draw(hp, annotation_legend_list = lgd_list)


# Extract macrophages and look at expression
clustered_obj_mac <- subset(clustered_obj, subset = ct_first_pass == "Macrophages")

DimPlot(clustered_obj_mac, group.by = "ct_first_pass", 
        raster = F,
        reduction = "umap")


FeaturePlot(clustered_obj_mac, 
            features = c(
              "PTPRC", "CD68", "CD163", "CD14", "LYZ"
            ),
            reduction = "umap",
            raster=FALSE,
            ncol = 5
)


FeaturePlot(clustered_obj_mac, 
            features = c(
              "TTN", "FHL2", "DES", "S100A1"
            ),
            reduction = "umap",
            raster=FALSE,
            ncol = 4
)


DimPlot(clustered_obj_mac, group.by = "leiden_res0.5", 
        raster = F,
        reduction = "umap")


# Re-cluster immune cells (immune pass 3) ----
imm_pass3 <- clustered_obj@meta.data %>% 
  filter(ct_first_pass %in% c("CD8+ T cells", "NK", "Macrophages (Mes obj)", 
                              "CD4+ T cells", "pDC", "Macrophages", 
                              "Fibroblast associated CD8+ T cells", 
                              "Proliferating T cells", "Treg", "B cells", "Mast", 
                              "cDC1", "Plasma", "cDC2", "mDC", 
                              "Proliferating myeloid")) %>% 
  dplyr::pull(cell_id)



FeaturePlot(clustered_obj, 
            features = c(
              "nCount_Xenium", "nFeature_Xenium"
            ),
            reduction = "umap",
            raster=FALSE,
            ncol = 2
)
VlnPlot(clustered_obj, 
        features = c(
          "nCount_Xenium", "nFeature_Xenium"
        ),
        group.by = "ct_first_pass",
        raster=FALSE,
        ncol = 1, pt.size = 0)



# Re-cluster all macrophages (pass 3 annotations) ----
cell_id_macro_pass3 <- clustered_obj@meta.data %>%  
  filter(ct_first_pass %in% c("Macrophages (Mes obj)", "Macrophages")) %>% 
  dplyr::pull(cell_id)

merged_spatial_filtered_macro_pass3 <- subset(merged_spatial_filtered_filtered, 
                                              subset = cell_id %in% cell_id_macro_pass3)
dim(merged_spatial_filtered_macro_pass3)

saveRDS(merged_spatial_filtered_macro_pass3, "/scratch/aoill/projects/heart_transplant/new/merged_spatial_filtered_macro_pass3.rds")


clustered_obj_macro_pass3 <- read_rds("/scratch/aoill/projects/heart_transplant/new/clustered_obj_macro_pass3_NC100_NN20_PC20_2024_06_18.rds")



clustered_obj_macro_pass3$subcluster <- as.character(clustered_obj_macro_pass3$subcluster)
clustered_obj_macro_pass3$leiden_res0.1 <- as.character(clustered_obj_macro_pass3$leiden_res0.1)
clustered_obj_macro_pass3$leiden_res0.2 <- as.character(clustered_obj_macro_pass3$leiden_res0.2)
clustered_obj_macro_pass3$leiden_res0.3 <- as.character(clustered_obj_macro_pass3$leiden_res0.3)
clustered_obj_macro_pass3$leiden_res0.5 <- as.character(clustered_obj_macro_pass3$leiden_res0.5)
clustered_obj_macro_pass3$leiden_res0.7 <- as.character(clustered_obj_macro_pass3$leiden_res0.7)
clustered_obj_macro_pass3$leiden_res0.75 <- as.character(clustered_obj_macro_pass3$leiden_res0.75)
clustered_obj_macro_pass3$leiden_res1.0 <- as.character(clustered_obj_macro_pass3$leiden_res1.0)
clustered_obj_macro_pass3$leiden_res1.5 <- as.character(clustered_obj_macro_pass3$leiden_res1.5)
clustered_obj_macro_pass3$leiden_res2.0 <- as.character(clustered_obj_macro_pass3$leiden_res2.0)


# Cell numbers
#0    1    2    3    4    5  6,0  6,1  6,2  6,3 
#7536  311 3326 7016  278 8511 1578 1279  255 1040 

DimPlot(clustered_obj_macro_pass3, 
        reduction = "umap",
        raster = F,
        group.by = "subcluster",
        label = T
) + NoLegend()

DimPlot(clustered_obj_macro_pass3, 
        reduction = "umap",
        raster = F,
        group.by = "subcluster",
        split.by = "subcluster",
        ncol = 4)

FeaturePlot(clustered_obj_macro_pass3, 
            features = c(
              "nCount_Xenium",
              "nFeature_Xenium"
            ),
            reduction = "umap",
            raster=FALSE,
            ncol = 2
)

VlnPlot(clustered_obj_macro_pass3, 
        features = c( "nCount_Xenium",
                      "nFeature_Xenium"
        ),
        group.by = "leiden_res0.7",
        pt.size = 0,
        ncol = 2, raster=FALSE)

FeaturePlot(clustered_obj_macro_pass3, 
            features = c(
              "PTPRC"#, # Immune
              #"EPCAM", # Epithelial
              #"PECAM1" # Endothelial
            ),
            reduction = "umap",
            raster=FALSE,
            ncol = 1
)
VlnPlot(clustered_obj_macro_pass3, 
        features = c("PTPRC"#, # Immune
                     #"EPCAM", # Epithelial
                     #"PECAM1" # Endothelial
        ),
        group.by = "subcluster2",
        pt.size = 0,
        ncol = 1, raster=FALSE)


FeaturePlot(clustered_obj_macro_pass3, 
            features = c(
              "PTPRC", "TTN", "DES", "FHL2", "S100A1"
            ),
            reduction = "umap",
            raster=FALSE,
            ncol = 2
)
VlnPlot(clustered_obj_macro_pass3, 
        features = c("PTPRC", "TTN", "DES", "FHL2", "S100A1"
        ),
        group.by = "subcluster",
        pt.size = 0,
        ncol = 2, raster=FALSE)

DefaultAssay(clustered_obj_macro_pass3) <- "RNA" 
clustered_obj_macro_pass3 <- NormalizeData(clustered_obj_macro_pass3)
VariableFeatures(clustered_obj_macro_pass3) <- rownames(clustered_obj_macro_pass3)
clustered_obj_macro_pass3 <- ScaleData(clustered_obj_macro_pass3)

DotPlot(clustered_obj_macro_pass3, 
        features = c(#"PECAM1",
          "PTPRC",
          "CD163", "MRC1", "MARCO", "FCGR1A", # Macrophage
          "LYZ", "CD14", "FCGR3A",  "MS4A7", # Monocyte
          "CD68",
          "FCER1A", #mo-DCs, cDC1, or pDC
          "CD1A", "CD1C", "CD209", "ITGAM", # monocyte derived DCs
          "ITGAX", "XCR1", # cDC1
          "GNLY", "NKG7", "KLRB1", #"KLRC1", #NK 
          "KLRD1", #NK
          "CD3E", "CD3D", "TRAC", "CD8A", "CD4", "FOXP3", # T
          "MS4A2", # mast
          "MS4A1", # B
          "CD79A", # B and plasma
          "TNFRSF17", "DERL3", # plasma
          "LILRA4", "CCR7", # pDCs
          "MKI67"
        ), 
        group.by = "subcluster", 
        scale = TRUE) + coord_flip()

DotPlot(clustered_obj_macro_pass3, 
        features = c(
          #"TSPY3", "PSG2", "RTKN2", "ELF5", "SOX2", "CELF4",
          "EPCAM", "PECAM1",
          "VWF", "BMX", # Ednocardial cells
          "MYH11",  # Vascular smooth muscle cells (VSMCs) 
          "PDGFRB", "ACTA2", # Pericytes
          "DCN", "C7", "FBLN1", "LTBP2", 
          "OGN", "PDGFRA", # Fibroblasts
          "ADIPOQ", # Adipocytes
          "TTN", "DES", "FHL2", "S100A1", # cardiomyocytes
          "PTPRC", "CD68",
          "CD163", "MRC1", "MARCO", "FCGR1A", # Macrophage
          "LYZ", "CD14", "FCGR3A",  "MS4A7", # Monocyte
          "FCER1A", #mo-DCs, cDC1, or pDC
          "CD1A", "CD1C", "CD209", "ITGAM", # monocyte derived DCs
          "ITGAX", "XCR1", # cDC1
          "GNLY", "NKG7", "KLRB1", #"KLRC1", #NK 
          "KLRD1", #NK
          "CD3E", "CD3D", "TRAC", "CD8A", "CD4", "FOXP3", # T
          "MS4A2", # mast
          "MS4A1", # B
          "CD79A", # B and plasma
          "TNFRSF17", "DERL3", # plasma
          "LILRA4", "CCR7", # pDCs
          "MKI67"
        ), 
        group.by = "subcluster", 
        scale = TRUE) + coord_flip() #+ 
#theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

# non-immune
FeaturePlot(clustered_obj_macro_pass3, 
            features = c(
              "PECAM1", "VWF", # Endothelial cells
              "MYH11",  # Vascular smooth muscle cells (VSMCs) 
              "PDGFRB", "ACTA2", # Pericytes
              "DCN", # Fibroblasts
              "TTN", "FHL2", "DES"
            ),
            reduction = "umap",
            raster=FALSE,
            ncol = 3
)

FeaturePlot(clustered_obj_macro_pass3, 
            features = c(
              "FHL2"
            ),
            reduction = "umap",
            split.by = "leiden_res0.7",
            raster=FALSE,
            ncol = 3
) + patchwork::plot_layout(ncol = 3, nrow = 3)

# T cells/NK
FeaturePlot(clustered_obj_macro_pass3, 
            features = c(
              #"GNLY", "NKG7", "KLRB1", #"KLRC1", #NK 
              #"KLRD1", #NK
              "CD3E", "CD3D", "TRAC", "CD8A", "CD4", "FOXP3" # T
              # "IL7R", "CCR7", # T
            ),
            reduction = "umap",
            raster=FALSE,
            ncol = 5
)

FeaturePlot(clustered_obj_macro_pass3, 
            features = c(
              "CD3E"
            ),
            reduction = "umap",
            split.by = "subcluster",
            raster=FALSE,
            ncol = 3
) + patchwork::plot_layout(ncol = 4, nrow = 3)

VlnPlot(clustered_obj_macro_pass3, 
        features = c(
          #"GNLY", "NKG7", "KLRB1", #"KLRC1", #NK 
          #"KLRD1", #NK
          "CD3E", "CD3D", "TRAC", "CD8A", "CD4", "FOXP3" # T
          # "IL7R", "CCR7", # T
          
        ),
        group.by = "subcluster",
        pt.size = 0,
        ncol = 4, raster=FALSE)



# Mono/mac
FeaturePlot(clustered_obj_macro_pass3, 
            features = c(
              #"CD80", "CD86", "CD206",
              "SPP1", "MARCO",
              "CD68", "CD163", "MRC1", "FCGR1A", # Macrophage
              "LYZ", "CD14", "FCGR3A",  "MS4A7" # Monocyte
            ),
            reduction = "umap",
            raster=FALSE,
            ncol = 5
)

VlnPlot(clustered_obj_macro_pass3, 
        features = c(
          "SPP1", "MARCO",
          "CD68", "CD163", "MRC1", "FCGR1A", # Macrophage
          "LYZ", "CD14", "FCGR3A",  "MS4A7" # Monocyte
          
        ),
        group.by = "subcluster",
        pt.size = 0,
        ncol = 5, raster=FALSE)



# moDCs
FeaturePlot(clustered_obj_macro_pass3, 
            features = c(
              "FCER1A", #mo-DCs, cDC1, or pDC
              "CD1A", "CD1C",	"MRC1", "CD209", "ITGAM" # monocyte derived DCs
              
            ),
            reduction = "umap",
            raster=FALSE,
            ncol = 4 
)
VlnPlot(clustered_obj_macro_pass3, 
        features = c(
          "FCER1A", #mo-DCs, cDC1, or pDC
          "CD1A", "CD1C",	"MRC1", "CD209", "ITGAM" # monocyte derived DCs
          
        ),
        group.by = "subcluster",
        pt.size = 0,
        ncol = 3, raster=FALSE)

# DC1 and 2
FeaturePlot(clustered_obj_macro_pass3, 
            features = c(
              "CD8A", "ITGAX", "XCR1", # cDC1
              "CD1C", "ITGAM" # cDC2 or mo-DC
            ),
            reduction = "umap",
            raster=FALSE,
            ncol = 3 
)
VlnPlot(clustered_obj_macro_pass3, 
        features = c(
          "CD8A", "ITGAX", "XCR1", # cDC1
          "CD1C", "ITGAM" # cDC2 or mo-DC
          
        ),
        group.by = "subcluster",
        pt.size = 0,
        ncol = 3, raster=FALSE)


# pDCs
FeaturePlot(clustered_obj_macro_pass3, 
            features = c(
              "LILRA4", "CCR7" # pDCs
            ),
            reduction = "umap",
            raster=FALSE,
            ncol = 2 
)
VlnPlot(clustered_obj_macro_pass3, 
        features = c(
          "LILRA4", "CCR7" # pDCs
        ),
        group.by = "subcluster",
        pt.size = 0,
        ncol = 2, raster=FALSE)


# B and plasma
FeaturePlot(clustered_obj_macro_pass3, 
            features = c(
              "MS4A1", # B
              "CD79A", # B and plasma
              "TNFRSF17", "DERL3" # plasma
            ),
            reduction = "umap",
            raster=FALSE,
            ncol = 4 
)
VlnPlot(clustered_obj_macro_pass3, 
        features = c(
          "MS4A1", # B
          "CD79A", # B and plasma
          "TNFRSF17", "DERL3" # plasma       
        ),
        group.by = "subcluster",
        pt.size = 0,
        ncol = 2, raster=FALSE)

FeaturePlot(clustered_obj_macro_pass3, 
            features = c(
              "BANK1"
            ),
            reduction = "umap",
            raster=FALSE,
            ncol = 2 
)

# Mast
FeaturePlot(clustered_obj_macro_pass3, 
            features = c(
              "MS4A2" # mast
            ),
            reduction = "umap",
            raster=FALSE 
)
VlnPlot(clustered_obj_macro_pass3, 
        features = c(
          "MS4A2" # mast      
        ),
        group.by = "subcluster",
        pt.size = 0, raster=FALSE)


FeaturePlot(clustered_obj_macro_pass3, 
            features = c(
              "MS4A2", # mast
              "SLC18A2"#, "VWA5A" # granulocyte
            ),
            reduction = "umap",
            raster=FALSE 
)

# Proliferating
FeaturePlot(clustered_obj_macro_pass3, 
            features = c(
              "MKI67"
            ),
            reduction = "umap",
            raster=FALSE
)

FeaturePlot(clustered_obj_macro_pass3, 
            features = c(
              "TTN", "FHL2"
            ),
            reduction = "umap",
            raster=FALSE
)

# Top markers
markers <- presto::wilcoxauc(clustered_obj_macro_pass3, group_by = "subcluster", 
                             assay = "data", seurat_assay = "RNA")
topmarkers <- top_markers(markers)
topmarkers <- top_markers(markers, n = 5, auc_min = 0.6, pct_in_min = 50, 
                          pct_out_max = 100)


topmarkers

## NEW annotations for cells that were originally labeled as macrophages ----
# Cardiomyocytes: 6,0
# Pericytes: 6,1
# Proliferating T cells: 6,2
# Fibroblasts: 6,3 
# Adipocytes: 4
# CD8+ T cells: 2
# SPP1+ Macrophages: 1
# Macrophages: 0, 3

# 5 not annotated here

## Re-cluster cluster 5 ----
cell_id_macro_c5_pass3 <- clustered_obj_macro_pass3@meta.data %>%  
  filter(leiden_res0.7 %in% c("5")) %>% 
  dplyr::pull(cell_id)

merged_spatial_filtered_macro_c5_pass3 <- subset(merged_spatial_filtered_filtered, 
                                                 subset = cell_id %in% cell_id_macro_c5_pass3)
dim(merged_spatial_filtered_macro_c5_pass3)

saveRDS(merged_spatial_filtered_macro_c5_pass3, "/scratch/aoill/projects/heart_transplant/new/merged_spatial_filtered_macro_c5_pass3.rds")


clustered_obj_macro_c5_pass3 <- read_rds("/scratch/aoill/projects/heart_transplant/new/clustered_obj_macro_c5_pass3_NC100_NN20_PC20_2024_06_18.rds")

clustered_obj_macro_c5_pass3$leiden_res0.1 <- as.character(clustered_obj_macro_c5_pass3$leiden_res0.1)
clustered_obj_macro_c5_pass3$leiden_res0.2 <- as.character(clustered_obj_macro_c5_pass3$leiden_res0.2)
clustered_obj_macro_c5_pass3$leiden_res0.3 <- as.character(clustered_obj_macro_c5_pass3$leiden_res0.3)
clustered_obj_macro_c5_pass3$leiden_res0.5 <- as.character(clustered_obj_macro_c5_pass3$leiden_res0.5)
clustered_obj_macro_c5_pass3$leiden_res0.7 <- as.character(clustered_obj_macro_c5_pass3$leiden_res0.7)
clustered_obj_macro_c5_pass3$leiden_res0.75 <- as.character(clustered_obj_macro_c5_pass3$leiden_res0.75)
clustered_obj_macro_c5_pass3$leiden_res1.0 <- as.character(clustered_obj_macro_c5_pass3$leiden_res1.0)
clustered_obj_macro_c5_pass3$leiden_res1.5 <- as.character(clustered_obj_macro_c5_pass3$leiden_res1.5)
clustered_obj_macro_c5_pass3$leiden_res2.0 <- as.character(clustered_obj_macro_c5_pass3$leiden_res2.0)


table(clustered_obj_macro_c5_pass3@meta.data$leiden_res0.7)

#0    1    2    3 
#2001 1397 3516 1597 

DimPlot(clustered_obj_macro_c5_pass3, 
        reduction = "umap",
        raster = F,
        group.by = "leiden_res0.7",
        label = T
) + NoLegend()

DimPlot(clustered_obj_macro_c5_pass3, 
        reduction = "umap",
        raster = F,
        group.by = "leiden_res0.7",
        split.by = "leiden_res0.7",
        ncol = 4)

FeaturePlot(clustered_obj_macro_c5_pass3, 
            features = c(
              "nCount_Xenium",
              "nFeature_Xenium"
            ),
            reduction = "umap",
            raster=FALSE,
            ncol = 2
)

VlnPlot(clustered_obj_macro_c5_pass3, 
        features = c( "nCount_Xenium",
                      "nFeature_Xenium"
        ),
        group.by = "leiden_res0.7",
        pt.size = 0,
        ncol = 2, raster=FALSE)

FeaturePlot(clustered_obj_macro_c5_pass3, 
            features = c(
              "PTPRC"#, # Immune
              #"EPCAM", # Epithelial
              #"PECAM1" # Endothelial
            ),
            reduction = "umap",
            raster=FALSE,
            ncol = 1
)
VlnPlot(clustered_obj_macro_c5_pass3, 
        features = c("PTPRC"#, # Immune
                     #"EPCAM", # Epithelial
                     #"PECAM1" # Endothelial
        ),
        group.by = "leiden_res0.7",
        pt.size = 0,
        ncol = 1, raster=FALSE)


FeaturePlot(clustered_obj_macro_c5_pass3, 
            features = c(
              "PTPRC", "TTN", "DES", "FHL2", "S100A1"
            ),
            reduction = "umap",
            raster=FALSE,
            ncol = 2
)
VlnPlot(clustered_obj_macro_c5_pass3, 
        features = c("PTPRC", "TTN", "DES", "FHL2", "S100A1"
        ),
        group.by = "leiden_res0.7",
        pt.size = 0,
        ncol = 2, raster=FALSE)

DefaultAssay(clustered_obj_macro_c5_pass3) <- "RNA" 
clustered_obj_macro_c5_pass3 <- NormalizeData(clustered_obj_macro_c5_pass3)
VariableFeatures(clustered_obj_macro_c5_pass3) <- rownames(clustered_obj_macro_c5_pass3)
clustered_obj_macro_c5_pass3 <- ScaleData(clustered_obj_macro_c5_pass3)

DotPlot(clustered_obj_macro_c5_pass3, 
        features = c(#"PECAM1",
          "PTPRC",
          "CD163", "MRC1", "MARCO", "FCGR1A", # Macrophage
          "LYZ", "CD14", "FCGR3A",  "MS4A7", # Monocyte
          "CD68",
          "FCER1A", #mo-DCs, cDC1, or pDC
          "CD1A", "CD1C", "CD209", "ITGAM", # monocyte derived DCs
          "ITGAX", "XCR1", # cDC1
          "GNLY", "NKG7", "KLRB1", #"KLRC1", #NK 
          "KLRD1", #NK
          "CD3E", "CD3D", "TRAC", "CD8A", "CD4", "FOXP3", # T
          "MS4A2", # mast
          "MS4A1", # B
          "CD79A", # B and plasma
          "TNFRSF17", "DERL3", # plasma
          "LILRA4", "CCR7", # pDCs
          "MKI67"
        ), 
        group.by = "leiden_res0.7", 
        scale = TRUE) + coord_flip()

DotPlot(clustered_obj_macro_c5_pass3, 
        features = c(
          #"TSPY3", "PSG2", "RTKN2", "ELF5", "SOX2", "CELF4",
          "EPCAM", "PECAM1",
          "VWF", "BMX", # Ednocardial cells
          "MYH11",  # Vascular smooth muscle cells (VSMCs) 
          "PDGFRB", "ACTA2", # Pericytes
          "DCN", "C7", "FBLN1", "LTBP2", 
          "OGN", "PDGFRA", # Fibroblasts
          "ADIPOQ", # Adipocytes
          "TTN", "DES", "FHL2", "S100A1", # cardiomyocytes
          "PTPRC", "CD68",
          "CD163", "MRC1", "MARCO", "FCGR1A", # Macrophage
          "LYZ", "CD14", "FCGR3A",  "MS4A7", # Monocyte
          "FCER1A", #mo-DCs, cDC1, or pDC
          "CD1A", "CD1C", "CD209", "ITGAM", # monocyte derived DCs
          "ITGAX", "XCR1", # cDC1
          "GNLY", "NKG7", "KLRB1", #"KLRC1", #NK 
          "KLRD1", #NK
          "CD3E", "CD3D", "TRAC", "CD8A", "CD4", "FOXP3", # T
          "MS4A2", # mast
          "MS4A1", # B
          "CD79A", # B and plasma
          "TNFRSF17", "DERL3", # plasma
          "LILRA4", "CCR7", # pDCs
          "MKI67"
        ), 
        group.by = "leiden_res0.7", 
        scale = TRUE) + coord_flip() #+ 
#theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

# non-immune
FeaturePlot(clustered_obj_macro_c5_pass3, 
            features = c(
              "PECAM1", "VWF", # Endothelial cells
              "MYH11",  # Vascular smooth muscle cells (VSMCs) 
              "PDGFRB", "ACTA2", # Pericytes
              "DCN", # Fibroblasts
              "TTN", "FHL2", "DES"
            ),
            reduction = "umap",
            raster=FALSE,
            ncol = 3
)

FeaturePlot(clustered_obj_macro_c5_pass3, 
            features = c(
              "FHL2"
            ),
            reduction = "umap",
            split.by = "leiden_res0.7",
            raster=FALSE,
            ncol = 3
) + patchwork::plot_layout(ncol = 3, nrow = 3)

# T cells/NK
FeaturePlot(clustered_obj_macro_c5_pass3, 
            features = c(
              #"GNLY", "NKG7", "KLRB1", #"KLRC1", #NK 
              #"KLRD1", #NK
              "CD3E", "CD3D", "TRAC", "CD8A", "CD4", "FOXP3" # T
              # "IL7R", "CCR7", # T
            ),
            reduction = "umap",
            raster=FALSE,
            ncol = 4
)

FeaturePlot(clustered_obj_macro_c5_pass3, 
            features = c(
              "CD3E"
            ),
            reduction = "umap",
            split.by = "leiden_res0.7",
            raster=FALSE,
            ncol = 3
) + patchwork::plot_layout(ncol = 4, nrow = 3)

VlnPlot(clustered_obj_macro_c5_pass3, 
        features = c(
          #"GNLY", "NKG7", "KLRB1", #"KLRC1", #NK 
          #"KLRD1", #NK
          "CD3E", "CD3D", "TRAC", "CD8A", "CD4", "FOXP3" # T
          # "IL7R", "CCR7", # T
          
        ),
        group.by = "leiden_res0.7",
        pt.size = 0,
        ncol = 4, raster=FALSE)



# Mono/mac
FeaturePlot(clustered_obj_macro_c5_pass3, 
            features = c(
              #"CD80", "CD86", "CD206",
              "SPP1", "MARCO",
              "CD68", "CD163", "MRC1", "FCGR1A", # Macrophage
              "LYZ", "CD14", "FCGR3A",  "MS4A7" # Monocyte
            ),
            reduction = "umap",
            raster=FALSE,
            ncol = 5
)

VlnPlot(clustered_obj_macro_c5_pass3, 
        features = c(
          "SPP1", "MARCO",
          "CD68", "CD163", "MRC1", "FCGR1A", # Macrophage
          "LYZ", "CD14", "FCGR3A",  "MS4A7" # Monocyte
          
        ),
        group.by = "leiden_res0.7",
        pt.size = 0,
        ncol = 5, raster=FALSE)



# moDCs
FeaturePlot(clustered_obj_macro_c5_pass3, 
            features = c(
              "FCER1A", #mo-DCs, cDC1, or pDC
              "CD1A", "CD1C",	"MRC1", "CD209", "ITGAM" # monocyte derived DCs
              
            ),
            reduction = "umap",
            raster=FALSE,
            ncol = 4 
)
VlnPlot(clustered_obj_macro_c5_pass3, 
        features = c(
          "FCER1A", #mo-DCs, cDC1, or pDC
          "CD1A", "CD1C",	"MRC1", "CD209", "ITGAM" # monocyte derived DCs
          
        ),
        group.by = "leiden_res0.7",
        pt.size = 0,
        ncol = 3, raster=FALSE)

# DC1 and 2
FeaturePlot(clustered_obj_macro_c5_pass3, 
            features = c(
              "CD8A", "ITGAX", "XCR1", # cDC1
              "CD1C", "ITGAM" # cDC2 or mo-DC
            ),
            reduction = "umap",
            raster=FALSE,
            ncol = 3 
)
VlnPlot(clustered_obj_macro_c5_pass3, 
        features = c(
          "CD8A", "ITGAX", "XCR1", # cDC1
          "CD1C", "ITGAM" # cDC2 or mo-DC
          
        ),
        group.by = "leiden_res0.7",
        pt.size = 0,
        ncol = 3, raster=FALSE)


# pDCs
FeaturePlot(clustered_obj_macro_c5_pass3, 
            features = c(
              "LILRA4", "CCR7" # pDCs
            ),
            reduction = "umap",
            raster=FALSE,
            ncol = 2 
)
VlnPlot(clustered_obj_macro_c5_pass3, 
        features = c(
          "LILRA4", "CCR7" # pDCs
        ),
        group.by = "leiden_res0.7",
        pt.size = 0,
        ncol = 2, raster=FALSE)


# B and plasma
FeaturePlot(clustered_obj_macro_c5_pass3, 
            features = c(
              "MS4A1", # B
              "CD79A", # B and plasma
              "TNFRSF17", "DERL3" # plasma
            ),
            reduction = "umap",
            raster=FALSE,
            ncol = 4 
)
VlnPlot(clustered_obj_macro_c5_pass3, 
        features = c(
          "MS4A1", "BANK1", # B
          "CD79A", # B and plasma
          "TNFRSF17", "DERL3" # plasma       
        ),
        group.by = "leiden_res0.7",
        pt.size = 0,
        ncol = 3, raster=FALSE)

FeaturePlot(clustered_obj_macro_c5_pass3, 
            features = c(
              "BANK1"
            ),
            reduction = "umap",
            raster=FALSE,
            ncol = 2 
)

# Mast
FeaturePlot(clustered_obj_macro_c5_pass3, 
            features = c(
              "MS4A2" # mast
            ),
            reduction = "umap",
            raster=FALSE 
)
VlnPlot(clustered_obj_macro_c5_pass3, 
        features = c(
          "MS4A2" # mast      
        ),
        group.by = "leiden_res0.7",
        pt.size = 0, raster=FALSE)


FeaturePlot(clustered_obj_macro_c5_pass3, 
            features = c(
              "MS4A2", # mast
              "SLC18A2"#, "VWA5A" # granulocyte
            ),
            reduction = "umap",
            raster=FALSE 
)

# Proliferating
FeaturePlot(clustered_obj_macro_c5_pass3, 
            features = c(
              "MKI67"
            ),
            reduction = "umap",
            raster=FALSE
)

FeaturePlot(clustered_obj_macro_c5_pass3, 
            features = c(
              "TTN", "FHL2", "MYBPC1",
              "COL1A2", "DCN", "COL1A1"
            ),
            reduction = "umap",
            raster=FALSE,
            ncol = 3
)

# Top markers
markers <- presto::wilcoxauc(clustered_obj_macro_c5_pass3, group_by = "leiden_res0.7", 
                             assay = "data", seurat_assay = "RNA")
topmarkers <- top_markers(markers)
topmarkers <- top_markers(markers, n = 5, auc_min = 0.6, pct_in_min = 50, 
                          pct_out_max = 100)


topmarkers

### NEW annotations for cells that were originally labeled as macrophages ----
# Macrophages: 1
# CD8+ T cells: 3
# Cardiomyocyte: 2
# Fibroblast: 0

# Re-do cell labeling on all cell UMAP ----
# Add annotations to full object ----
## 1. Immune pass 2 object ----
# clustered_obj_imm_pass2
# subcluster2
# Endothelial: 13
# Fibroblasts: 7,0, 7,1
# Pericytes: 7,2
# These are final annotations for these clusters

cell_cluster_info1 <- clustered_obj_imm_pass2@meta.data %>%
  dplyr::select(X, subcluster2)

cell_cluster_info1$ct_second_pass <- "NA"

cell_cluster_info1$ct_second_pass[
  cell_cluster_info1$subcluster2 %in% c("13")] <- "Endothelial"

cell_cluster_info1$ct_second_pass[
  cell_cluster_info1$subcluster2 %in% c("7,0", "7,1")] <- "Fibroblasts"

cell_cluster_info1$ct_second_pass[
  cell_cluster_info1$subcluster2 %in% c("7,2")] <- "Pericytes"

table(cell_cluster_info1$ct_second_pass)

cell_cluster_info1_df <- cell_cluster_info1 %>% dplyr::select(X, ct_second_pass)
table(cell_cluster_info1_df$ct_second_pass)

cell_cluster_info1_df_final <- cell_cluster_info1_df %>%
  filter(ct_second_pass != "NA")
table(cell_cluster_info1_df_final$ct_second_pass)


## 2. Meyloid pass 2 object ----
# clustered_obj_meyloid
# subcluster
#Cardiomyocytes: 3
#Mast: 9
#pDC: 8
#Proliferating myeloid: 6,2
#Mature DC (mDCs): 6,1 
#cDC2: 5
#cDC1: 6,0
#Macrophages: 0, 1, 2, 4, 7
# NEW- Endothelial: 10

cell_cluster_info2 <- clustered_obj_meyloid@meta.data %>%
  dplyr::select(X, subcluster)

cell_cluster_info2$ct_second_pass <- "NA"

cell_cluster_info2$ct_second_pass[
  cell_cluster_info2$subcluster %in% c("3")] <- "Cardiomyocytes"

cell_cluster_info2$ct_second_pass[
  cell_cluster_info2$subcluster %in% c("9")] <- "Mast"

cell_cluster_info2$ct_second_pass[
  cell_cluster_info2$subcluster %in% c("8")] <- "pDC"

cell_cluster_info2$ct_second_pass[
  cell_cluster_info2$subcluster %in% c("6,2")] <- "Proliferating DCs"

cell_cluster_info2$ct_second_pass[
  cell_cluster_info2$subcluster %in% c("6,1")] <- "mDC"

cell_cluster_info2$ct_second_pass[
  cell_cluster_info2$subcluster %in% c("5")] <- "cDC2"

cell_cluster_info2$ct_second_pass[
  cell_cluster_info2$subcluster %in% c("6,0")] <- "cDC1"

#cell_cluster_info2$ct_second_pass[
#  cell_cluster_info2$subcluster %in% c("0", "1", "2", "4", "7")] <- "Macrophages" # Remove and replace!!!

cell_cluster_info2$ct_second_pass[
  cell_cluster_info2$subcluster %in% c("10")] <- "Endothelial"

table(cell_cluster_info2$ct_second_pass)

cell_cluster_info2_df <- cell_cluster_info2 %>% dplyr::select(X, ct_second_pass)
table(cell_cluster_info2_df$ct_second_pass)

cell_cluster_info2_df_final <- cell_cluster_info2_df %>%
  filter(ct_second_pass != "NA")
table(cell_cluster_info2_df_final$ct_second_pass)

nrow(clustered_obj_meyloid@meta.data)
nrow(cell_cluster_info2_df_final)

## 3. Lymphoid pass 2 object ----
# clustered_obj_lymphoid
# subcluster
#B cells: 6
#Plasma: 0
#NK: 7,0, 7,2

cell_cluster_info3 <- clustered_obj_lymphoid@meta.data %>%
  dplyr::select(X, subcluster)

cell_cluster_info3$ct_second_pass <- "NA"

cell_cluster_info3$ct_second_pass[
  cell_cluster_info3$subcluster %in% c("6")] <- "B cells"

cell_cluster_info3$ct_second_pass[
  cell_cluster_info3$subcluster %in% c("0")] <- "Plasma"

cell_cluster_info3$ct_second_pass[
  cell_cluster_info3$subcluster %in% c("7,0", "7,2")] <- "NK"


table(cell_cluster_info3$ct_second_pass)

cell_cluster_info3_df <- cell_cluster_info3 %>% dplyr::select(X, ct_second_pass)
table(cell_cluster_info3_df$ct_second_pass)

cell_cluster_info3_df_final <- cell_cluster_info3_df %>%
  filter(ct_second_pass != "NA")
table(cell_cluster_info3_df_final$ct_second_pass)

nrow(clustered_obj_lymphoid@meta.data)
nrow(clustered_obj_T@meta.data)
nrow(cell_cluster_info3_df_final)
nrow(cell_cluster_info4_df_final)


## 4. T cell pass 2 object ----
# clustered_obj_T
# subcluster
#Treg: 5,2
#Proliferating T cells: 2
#CD4+ T cells: 1,0, 1,2,
#CD8+ T cells: 1,3, 1,4, 3, 4, 6
#   0, 1,1, 5,0, 5,1

cell_cluster_info4 <- clustered_obj_T@meta.data %>%
  dplyr::select(X, subcluster)

cell_cluster_info4$ct_second_pass <- "NA"

cell_cluster_info4$ct_second_pass[
  cell_cluster_info4$subcluster %in% c("5,2")] <- "Treg"

cell_cluster_info4$ct_second_pass[
  cell_cluster_info4$subcluster %in% c("2")] <- "Proliferating T cells"

cell_cluster_info4$ct_second_pass[
  cell_cluster_info4$subcluster %in% c("1,0", "1,2")] <- "CD4+ T cells"

cell_cluster_info4$ct_second_pass[
  cell_cluster_info4$subcluster %in% c("1,3", "1,4", "3", "4", "6",
                                       "0", "1,1", "5,0", "5,1")] <- "CD8+ T cells"

table(cell_cluster_info4$ct_second_pass)

cell_cluster_info4_df <- cell_cluster_info4 %>% dplyr::select(X, ct_second_pass)
table(cell_cluster_info4_df$ct_second_pass)

cell_cluster_info4_df_final <- cell_cluster_info4_df %>%
  filter(ct_second_pass != "NA")
table(cell_cluster_info4_df_final$ct_second_pass)


## 5. Mesenchymal pass 2 object ----
# clustered_obj_mes_pass2
# subcluster

#Fibroblast associated CD8+ T cells: 0
#Cardiomyocytes: 4,0, 4,3
#Endothelial: 3,0
#Macrophages: 1 - label Macrophage (mes obj)

#vSMCs: 5, 3,1
#Fibroblasts: 2,0, 2,2, 2,3
#POSTN+ Fibroblasts: 2,1 
#Proliferating pericytes: 4,1  
#Myofibroblasts: 6
#Pericytes: 4,2, 4,4, 4,5, 4,6



cell_cluster_info5 <- clustered_obj_mes_pass2@meta.data %>%
  dplyr::select(X, subcluster)

cell_cluster_info5$ct_second_pass <- "NA"

cell_cluster_info5$ct_second_pass[
  cell_cluster_info5$subcluster %in% c("0")] <- "CD8+ T cells" # Fibroblast associated CD8+ T cells

cell_cluster_info5$ct_second_pass[
  cell_cluster_info5$subcluster %in% c("4,0", "4,3")] <- "Cardiomyocytes"

cell_cluster_info5$ct_second_pass[
  cell_cluster_info5$subcluster %in% c("3,0")] <- "Endothelial"

#cell_cluster_info5$ct_second_pass[
#  cell_cluster_info5$subcluster %in% c("1")] <- "Macrophages (Mes obj)" # Remove and replace!!!

cell_cluster_info5$ct_second_pass[
  cell_cluster_info5$subcluster %in% c("5")] <- "vSMCs"

cell_cluster_info5$ct_second_pass[
  cell_cluster_info5$subcluster %in% c("3,1")] <- "Endothelial" # "vSMC associated endothelial

cell_cluster_info5$ct_second_pass[
  cell_cluster_info5$subcluster %in% c("2,0", "2,2", "2,3")] <- "Fibroblasts"

cell_cluster_info5$ct_second_pass[
  cell_cluster_info5$subcluster %in% c("2,1")] <- "POSTN+ Fibroblasts"

cell_cluster_info5$ct_second_pass[
  cell_cluster_info5$subcluster %in% c("4,1")] <- "Proliferating pericytes"

cell_cluster_info5$ct_second_pass[
  cell_cluster_info5$subcluster %in% c("6")] <- "Myofibroblasts"

cell_cluster_info5$ct_second_pass[
  cell_cluster_info5$subcluster %in% c("4,2", "4,4", "4,5", "4,6")] <- "Pericytes"


table(cell_cluster_info5$ct_second_pass)

cell_cluster_info5_df <- cell_cluster_info5 %>% dplyr::select(X, ct_second_pass)
table(cell_cluster_info5_df$ct_second_pass)

cell_cluster_info5_df_final <- cell_cluster_info5_df %>%
  filter(ct_second_pass != "NA")
table(cell_cluster_info5_df_final$ct_second_pass)


## 6. Endothelial pass 2 object ----
# clustered_obj_endo_pass2
# subcluster
# Cardiomyocytes: 0,0

# Proliferating endothelial: 2
# Lymphatic endothelial: 5
# Activated endothelial: 4
# BMX+ Activated endothelial: 6
# Endothelial: 0,1, 0,2, 0,3, 1,0, 1,1, 3



cell_cluster_info6 <- clustered_obj_endo_pass2@meta.data %>%
  dplyr::select(X, subcluster)

cell_cluster_info6$ct_second_pass <- "NA"

cell_cluster_info6$ct_second_pass[
  cell_cluster_info6$subcluster %in% c("0,0")] <- "Cardiomyocytes"

cell_cluster_info6$ct_second_pass[
  cell_cluster_info6$subcluster %in% c("2")] <- "Proliferating endothelial"

cell_cluster_info6$ct_second_pass[
  cell_cluster_info6$subcluster %in% c("5")] <- "Lymphatic endothelial"

cell_cluster_info6$ct_second_pass[
  cell_cluster_info6$subcluster %in% c("4")] <- "Activated endothelial"

cell_cluster_info6$ct_second_pass[
  cell_cluster_info6$subcluster %in% c("6")] <- "BMX+ Activated endothelial"

cell_cluster_info6$ct_second_pass[
  cell_cluster_info6$subcluster %in% c("0,1", "0,2", "0,3", "1,0", "1,1", "3")] <- "Endothelial"


table(cell_cluster_info6$ct_second_pass)

cell_cluster_info6_df <- cell_cluster_info6 %>% dplyr::select(X, ct_second_pass)
table(cell_cluster_info6_df$ct_second_pass)

cell_cluster_info6_df_final <- cell_cluster_info6_df %>%
  filter(ct_second_pass != "NA")
table(cell_cluster_info6_df_final$ct_second_pass)


## 7. Cardiomyocytes object ----
# clustered_obj_card
# leiden_res0.4
# Cardiomyocytes: "0" "1" "3" "2"


cell_cluster_info7 <- clustered_obj_card@meta.data %>%
  dplyr::select(X, leiden_res0.4)

cell_cluster_info7$ct_second_pass <- "NA"

cell_cluster_info7$ct_second_pass[
  cell_cluster_info7$leiden_res0.4 %in% c("0", "1", "3", "2")] <- "Cardiomyocytes"



table(cell_cluster_info7$ct_second_pass)

cell_cluster_info7_df <- cell_cluster_info7 %>% dplyr::select(X, ct_second_pass)
table(cell_cluster_info7_df$ct_second_pass)

cell_cluster_info7_df_final <- cell_cluster_info7_df %>%
  filter(ct_second_pass != "NA")
table(cell_cluster_info7_df_final$ct_second_pass)


## 8. Cardiomyocytes from immune pass 1 object ----
# clustered_obj_imm
# leiden_res0.7
# Cardiomyocytes: "4"


cell_cluster_info8 <- clustered_obj_imm@meta.data %>%
  dplyr::select(X, leiden_res0.7)

cell_cluster_info8$ct_second_pass <- "NA"

cell_cluster_info8$ct_second_pass[
  cell_cluster_info8$leiden_res0.7 %in% c("4")] <- "Cardiomyocytes"



table(cell_cluster_info8$ct_second_pass)

cell_cluster_info8_df <- cell_cluster_info8 %>% dplyr::select(X, ct_second_pass)
table(cell_cluster_info8_df$ct_second_pass)

cell_cluster_info8_df_final <- cell_cluster_info8_df %>%
  filter(ct_second_pass != "NA")
table(cell_cluster_info8_df_final$ct_second_pass)


## 9. Cardiomyocytes from mesenchymal pass 1 object ----
# clustered_obj_mes
# leiden_res0.5
# Cardiomyocytes: "4"

cell_cluster_info9 <- clustered_obj_mes@meta.data %>%
  dplyr::select(X, leiden_res0.5)

cell_cluster_info9$ct_second_pass <- "NA"

cell_cluster_info9$ct_second_pass[
  cell_cluster_info9$leiden_res0.5 %in% c("4")] <- "Cardiomyocytes"



table(cell_cluster_info9$ct_second_pass)

cell_cluster_info9_df <- cell_cluster_info9 %>% dplyr::select(X, ct_second_pass)
table(cell_cluster_info9_df$ct_second_pass)

cell_cluster_info9_df_final <- cell_cluster_info9_df %>%
  filter(ct_second_pass != "NA")
table(cell_cluster_info9_df_final$ct_second_pass)


## 10. Cardiomyocytes from endothelial pass 1 object ----
# clustered_obj_endo
# subcluster2
# Cardiomyocytes: "0"

cell_cluster_info10 <- clustered_obj_endo@meta.data %>%
  dplyr::select(X, subcluster2)

cell_cluster_info10$ct_second_pass <- "NA"

cell_cluster_info10$ct_second_pass[
  cell_cluster_info10$subcluster2 %in% c("0")] <- "Cardiomyocytes"



table(cell_cluster_info10$ct_second_pass)

cell_cluster_info10_df <- cell_cluster_info10 %>% dplyr::select(X, ct_second_pass)
table(cell_cluster_info10_df$ct_second_pass)

cell_cluster_info10_df_final <- cell_cluster_info10_df %>%
  filter(ct_second_pass != "NA")
table(cell_cluster_info10_df_final$ct_second_pass)


## 11. Macrophages pass 3 object ----
# clustered_obj_macro_pass3
# subcluster

# Cardiomyocytes: 6,0
# Pericytes: 6,1
# Proliferating T cells: 6,2
# Fibroblasts: 6,3 
# Adipocytes: 4
# CD8+ T cells: 2
# SPP1+ Macrophages: 1
# Macrophages: 0, 3

cell_cluster_info11 <- clustered_obj_macro_pass3@meta.data %>%
  dplyr::select(X, subcluster)

cell_cluster_info11$ct_second_pass <- "NA"

cell_cluster_info11$ct_second_pass[
  cell_cluster_info11$subcluster %in% c("6,0")] <- "Cardiomyocytes"

cell_cluster_info11$ct_second_pass[
  cell_cluster_info11$subcluster %in% c("6,1")] <- "Pericytes"

cell_cluster_info11$ct_second_pass[
  cell_cluster_info11$subcluster %in% c("6,2")] <- "Proliferating T cells"

cell_cluster_info11$ct_second_pass[
  cell_cluster_info11$subcluster %in% c("6,3")] <- "Fibroblasts"

cell_cluster_info11$ct_second_pass[
  cell_cluster_info11$subcluster %in% c("4")] <- "Adipocytes"

cell_cluster_info11$ct_second_pass[
  cell_cluster_info11$subcluster %in% c("2")] <- "CD8+ T cells"

cell_cluster_info11$ct_second_pass[
  cell_cluster_info11$subcluster %in% c("1")] <- "SPP1+ Macrophages"

cell_cluster_info11$ct_second_pass[
  cell_cluster_info11$subcluster %in% c("0", "3")] <- "Macrophages"

table(cell_cluster_info11$ct_second_pass)

cell_cluster_info11_df <- cell_cluster_info11 %>% dplyr::select(X, ct_second_pass)
table(cell_cluster_info11_df$ct_second_pass)

cell_cluster_info11_df_final <- cell_cluster_info11_df %>%
  filter(ct_second_pass != "NA")
table(cell_cluster_info11_df_final$ct_second_pass)


## 12. Macrophages pass 3 object, cluster 5 ----
# clustered_obj_macro_c5_pass3
# leiden_res0.7

# Macrophages: 1
# CD8+ T cells: 3
# Cardiomyocytes: 2
# Fibroblasts: 0

cell_cluster_info12 <- clustered_obj_macro_c5_pass3@meta.data %>%
  dplyr::select(X, leiden_res0.7)

cell_cluster_info12$ct_second_pass <- "NA"

cell_cluster_info12$ct_second_pass[
  cell_cluster_info12$leiden_res0.7 %in% c("1")] <- "Macrophages"

cell_cluster_info12$ct_second_pass[
  cell_cluster_info12$leiden_res0.7 %in% c("3")] <- "CD8+ T cells"

cell_cluster_info12$ct_second_pass[
  cell_cluster_info12$leiden_res0.7 %in% c("2")] <- "Cardiomyocytes"

cell_cluster_info12$ct_second_pass[
  cell_cluster_info12$leiden_res0.7 %in% c("0")] <- "Fibroblasts"

table(cell_cluster_info12$ct_second_pass)

cell_cluster_info12_df <- cell_cluster_info12 %>% dplyr::select(X, ct_second_pass)
table(cell_cluster_info12_df$ct_second_pass)

cell_cluster_info12_df_final <- cell_cluster_info12_df %>%
  filter(ct_second_pass != "NA")
table(cell_cluster_info12_df_final$ct_second_pass)



## Merge all cell ids together ----
cell_cluster_info_all <- rbind(cell_cluster_info1_df_final, cell_cluster_info2_df_final, 
                               cell_cluster_info3_df_final, cell_cluster_info4_df_final,
                               cell_cluster_info5_df_final, cell_cluster_info6_df_final,
                               cell_cluster_info7_df_final, cell_cluster_info8_df_final,
                               cell_cluster_info9_df_final, cell_cluster_info10_df_final,
                               cell_cluster_info11_df_final, cell_cluster_info12_df_final)

nrow(cell_cluster_info_all)

nrow(clustered_obj@meta.data)
# 164165


## Full join to full object ----
clustered_obj@meta.data$ct_second_pass <- NULL
clustered_obj_md <- clustered_obj@meta.data
nrow(clustered_obj_md)
#clustered_obj@meta.data <- clustered_obj_md
clustered_obj_md_all <- left_join(clustered_obj_md, cell_cluster_info_all)
nrow(clustered_obj_md_all)
#View(clustered_obj_md_all)
clustered_obj@meta.data$ct_second_pass <- clustered_obj_md_all$ct_second_pass



saveRDS(clustered_obj, "/scratch/aoill/projects/heart_transplant/new/clustered_obj_annotated_second_pass_06_19_2024.rds")

