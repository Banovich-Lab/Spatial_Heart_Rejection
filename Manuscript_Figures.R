# Load libraries ----
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(speckle)
library(Seurat)
library(ggalluvial)
library(ComplexHeatmap)
library(presto)
library(factoextra)
library(FactoMineR)
library(tidyr)
library(circlize)  


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


# Color pallet ----
# Define a custom palette for 28 cell types
mycolors <- c(
  # Immune
  "Macrophages" = "#6a0572", 
  "Mast" = "#E281CA",
  "cDC1" = "#D2A5E1",
  "cDC2" = "#CE6C7F",
  "SPP1+ Macrophages" = "#DC3F92",
  "mDC" = "#D035E9",
  "Proliferating DCs" = "#CA5BD6",
  "pDC" = "#E2B2BF",
  "NK" = "#C1C8E6",
  "CD8+ T cells" = "#6FA4E2", 
  "B cells" = "#76CADE", 
  "CD4+ T cells" = "#B6EBDB",
  "Proliferating T cells" = "#5EE8D7",
  "Treg" = "#0000ff", 
  "Plasma" = "#6699ff", 
  
  
  # Endothelial
  "Endothelial" = "#ae2012",
  "BMX+ Activated endothelial" = "#CD9764",
  "Proliferating endothelial" = "#E4DF8A",
  "Lymphatic endothelial" = "#DDBB44",
  "Activated endothelial" = "#e09074", 
  
  # Mesenchymal
  # "#74E8A0"
  "Fibroblasts" = "#C4E3AF",
  "Cardiomyocytes" = "#606d39", 
  "POSTN+ Fibroblasts" = "#76EA57",
  "Pericytes" = "#74E8A0",
  "vSMCs" = "#62B293",
  "Adipocytes" = "#52953e",
  "Myofibroblasts" = "#A4D367",
  "Proliferating pericytes" = "#9e9f3f" 
)


# Read in object ----
# This will be used for plotting
clustered_obj <- readRDS("/scratch/aoill/projects/heart_transplant/GEO/heart_spatial_obj.rds")

# Make cell groups file for explorer ----
cell_groups <- clustered_obj@meta.data %>% dplyr::select(X, ct_second_pass)
rownames(cell_groups) <- NULL
colnames(cell_groups) <- c("cell_id", "group")
head(cell_groups)
write.csv(cell_groups, "/scratch/aoill/projects/heart_transplant/00_final/cell_types.csv",
          quote = F, row.names = F)

# Main figures ----
# Figure 1A ----
# Generated using biorender

# Figure 1B ----
# Alluvial plot
spatial.alluvial <- read_csv("/home/aoill/projects/heart_transplant/spatial.demos_manual.update_including_nonbiopsies_ADDED_MISSING_POST_10.11.24.csv")

# Remove autopsy patient
spatial.alluvial <- spatial.alluvial %>%
  filter(biopsy_timing != "autopsy")

# Remove the one patient that was not in seurat object (should be Patient 16)
spatial.alluvial <- spatial.alluvial %>%
  filter(patient_id != "Patient_16")

# How many peds and adults
spatial.alluvial %>%
  filter(age_at_biopsy < 18) %>%
  summarise(N = unique(mrn))

spatial.alluvial %>%
  filter(age_at_biopsy >= 18) %>%
  summarise(N = unique(mrn))


spatial.alluvial <- spatial.alluvial %>%
  mutate(
    rejection = case_when(cellular_grading > 1 & antibody_grading == "pAMR0" ~ "ACR",
                          cellular_grading > 1 & antibody_grading == "pAMR1-h" ~ "ACR",
                          cellular_grading > 1 & antibody_grading == "pAMR1-H" ~ "ACR",
                          cellular_grading > 1 & is.na(antibody_grading) ~ "ACR",
                          cellular_grading > 1 & antibody_grading == "pAMR1-i" ~ "Mixed\nrejection",
                          cellular_grading > 1 & antibody_grading == "pAMR2" ~ "Mixed\nrejection",
                          cellular_grading < 2 & antibody_grading == "pAMR1-i" ~ "AMR",
                          cellular_grading < 2 & antibody_grading == "pAMR2" ~ "AMR",
                          cellular_grading < 2 & antibody_grading == "pAMR1-h" ~ "No\nrejection",
                          cellular_grading < 2 & antibody_grading == "pAMR0" ~ "No\nrejection")
  )

# Fix some spelling errors in the treatments
spatial.alluvial <- spatial.alluvial %>%
  mutate(
    treatment = case_when(treatment == "methylpre" ~ "methylpred",
                          treatment == "methypred" ~ "methylpred",
                          TRUE ~ treatment)
  )


# How many different treatments did patients receieve?
length(unique(spatial.alluvial$treatment)) # 14 unique treatments
table(spatial.alluvial$treatment)


spatial.alluvial <- spatial.alluvial %>%
  mutate(
    treatment_group = case_when(
      treatment == "none" ~ "No therapy",
      treatment == "oral_steroids" | treatment == "methylpred" ~ "Steroids",
      treatment == "methylpred, bortezomib" ~ "B cell\ntargeted\ntherapies",
      treatment == "methylpred, eculizumab, PLEX, IVIG, Rituximab" ~ "B cell\ntargeted\ntherapies",
      treatment == "methylpred, IVIG" ~ "B cell\ntargeted\ntherapies",
      treatment == "methylpred, IVIG, Rituximab" ~ "B cell\ntargeted\ntherapies",
      treatment == "methylpred, PLEX, IVIG" ~ "B cell\ntargeted\ntherapies",
      treatment == "methylpred, PLEX, IVIG, Rituximab" ~ "B cell\ntargeted\ntherapies",
      treatment == "methylpred, thymoglobulin" ~ "Cytolytic\ntherapy",
      treatment == "methylpred, thymoglobulin, IVIG, Rituximab" ~ "Combination\ntherapies",
      treatment == "methylpred, thymoglobulin, PLEX, IVIG" ~ "Combination\ntherapies",
      treatment == "methylpred, thymoglobulin, PLEX, IVIG, Rituximab" ~ "Combination\ntherapies",
      treatment == "methylpred, thymoglobulin, Rituximab" ~ "Combination\ntherapies"
    )
  )



# Only select columns of import
df <- spatial.alluvial %>%
  select(patient_id, biopsy_timing, treatment_group, rejection)

# Only keep unique rows
df <- df %>%
  group_by(patient_id, biopsy_timing) %>%
  distinct(.keep_all = TRUE)

# Wide format so that pre- and post-biopsy rejection grading is on same row
df_wide <- df %>%
  spread(biopsy_timing, rejection)

# Set different rejection grades and treatments as factors for ordering
df_wide <- df_wide %>%
  mutate(
    pre = factor(pre, levels = c("ACR",
                                 "AMR",
                                 "Mixed\nrejection",
                                 "No\nrejection")),
    post = factor(post, levels = c("ACR",
                                   "AMR",
                                   "No\nrejection")),
    treatment_group = factor(treatment_group, levels = c("Steroids",
                                                         "Cytolytic\ntherapy",
                                                         "B cell\ntargeted\ntherapies",
                                                         "Combination\ntherapies",
                                                         "No therapy"))
  )

# Change to long form - otherwise can't color both sets of rectangles
df_long <- df_wide %>%
  to_lodes_form(key = "biopsy_timing", axes = c("pre", "treatment_group", "post")) %>%
  drop_na()

# Plot the alluvial plot with treatment_group in the middle
p1 <- ggplot(df_long,
             aes(x = biopsy_timing, stratum = stratum, alluvium = patient_id, fill = stratum, label = stratum)) +
  geom_flow(stat = "alluvium", lode.guidance = "forward", curve_type = "quintic") +
  geom_stratum(width = 0.5) +
  geom_text(stat = "stratum", aes(label = after_stat(stratum)),
            size = 3.5) +
  scale_x_discrete(limits = c("pre", "treatment_group", "post"),
                   labels = c("Rejection type\n pre-treatment",
                              "Treatment",
                              "Rejection type\n post-treatment")) +
  scale_fill_brewer(palette = "Spectral", direction = 1) +
  cowplot::theme_cowplot() +
  labs(x = "",
       y = "Count") + Seurat::NoLegend()



pdf("/home/aoill/plots/00_heart_manuscript/Alluvial_plot.pdf",
    width = 5.7, height = 8.8)
p1
dev.off()



# Figure 1C ----
# UMAP of all cells 
a <- DimPlot(clustered_obj, group.by = "ct_second_pass", 
             raster = F,
             reduction = "umap") +
  scale_color_manual(values = mycolors) +
  NoLegend() + 
  coord_fixed() + 
  ggtitle(NULL) +  
  xlab("UMAP 1") + ylab("UMAP 2")

pa <- LabelClusters(a, id = "ct_second_pass", box = TRUE, size = 3.5, label.size = 0.1, 
                    box.padding = 0.5, seed = 309,
                    alpha = 0.8#, max.overlaps=Inf
) 

ggsave("/home/aoill/plots/00_heart_manuscript/revisions_01/UMAP.png", plot = pa, width = 180, height = 180, units = "mm", dpi = 800)

pdf("/home/aoill/plots/00_heart_manuscript/revisions_01/UMAP.pdf")
pa
dev.off()


# Figure 1D ----
# Bar plot with the number of cell per cell type faceted by lineage
# Extract metadata from Seurat object
endothelial_cells <- c("Activated endothelial", "BMX+ Activated endothelial", 
                       "Endothelial", "Lymphatic endothelial", "Proliferating endothelial")

stromal_cells <- c("Cardiomyocytes", "Adipocytes", "Fibroblasts", 
                       "Myofibroblasts", "Pericytes", "POSTN+ Fibroblasts", 
                       "Proliferating pericytes", "vSMCs")

immune_cells <- c("B cells", "CD4+ T cells", "CD8+ T cells", "cDC1", "cDC2", 
                  "Macrophages", "Mast", "mDC", "NK", "pDC", "Plasma", 
                  "Proliferating DCs", "Proliferating T cells", 
                  "SPP1+ Macrophages", "Treg")



# Create a new column in the metadatato classify each cell type into its lineage2
clustered_obj@meta.data <- clustered_obj@meta.data %>%
  mutate(lineage2 = case_when(
    ct_second_pass %in% endothelial_cells ~ "Endothelial",
    ct_second_pass %in% stromal_cells ~ "Stromal",
    ct_second_pass %in% immune_cells ~ "Immune",
    TRUE ~ "Other"
  ))


# Bar plot with the number of cell per cell type faceted by lineage2
# Extract metadata from Seurat object
metadata <- clustered_obj@meta.data

# Summarize the number of cells per cell type and lineage2
cell_counts <- metadata %>%
  dplyr::group_by(lineage2, ct_second_pass) %>%
  dplyr::summarise(cell_count = n(), .groups = "drop")

# Reorder cell types (ct_second_pass) by cell count within each lineage2
cell_counts <- cell_counts %>%
  group_by(lineage2) %>%
  mutate(ct_second_pass = factor(ct_second_pass, levels = ct_second_pass[order(-cell_count)])) %>%
  ungroup()


# Convert lineage2 to a factor for proper ordering
cell_counts$lineage2 <- factor(cell_counts$lineage2)

# Add formatted cell count labels
cell_counts <- cell_counts %>%
  mutate(label = ifelse(cell_count >= 1000, paste0(round(cell_count / 1000, 1), "k"), as.character(cell_count)))



pb <- ggplot(cell_counts, 
             aes(x = ct_second_pass, 
                 y = cell_count, 
                 fill = ct_second_pass)) +
  geom_bar(stat = "identity") +  # Set bar width
  geom_text(aes(label = label), vjust = -0.5, size = 3.5) +  # Add labels above bars
  facet_grid(~ lineage2, scales = "free_x", space = "free") +  # Proportional facet widths
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +  # Add space above the bars
  scale_fill_manual(values=mycolors) + 
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 11),  
    axis.text.y = element_text(size = 11), 
    axis.title.y = element_text(size = 11),
    strip.text = element_text(size = 13),  
    axis.title.x = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "none"
  ) +
  labs(
    y = "Number of Cells",
    title = "",
    fill = "lineage2"
  ) 


# Change facet colors
bar_plot_table <- ggplot_gtable(ggplot_build(pb))
striprt <- which(grepl("strip-r", bar_plot_table$layout$name) | 
                   grepl("strip-t", bar_plot_table$layout$name))
fills <- c("#ae2012", "#6FA4E2", "#606d39", NA, NA) # Lineages & Sample Types
colors <- c(rep("black", 3), rep(NA, 2))
font_colors <- c(rep("white", 3), rep("black", 2))
k <- 1
for (i in striprt) {
  j <- which(grepl("rect", bar_plot_table$grobs[[i]]$grobs[[1]]$childrenOrder))
  bar_plot_table$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- fills[k]
  bar_plot_table$grobs[[i]]$grobs[[1]]$children[[j]]$gp$col <- colors[k]
  bar_plot_table$grobs[[i]]$grobs[[1]]$children[[j+1]]$children[[1]]$gp$col <- font_colors[k]
  k <- k+1
}

gridExtra::grid.arrange(bar_plot_table)

ggsave("/home/aoill/plots/00_heart_manuscript/revisions_01/Figure_01D.png", plot = gridExtra::grid.arrange(bar_plot_table), width = 280, height = 100, units = "mm", dpi = 800)


pdf("/home/aoill/plots/00_heart_manuscript/revisions_01/Figure_01D.pdf",
  width = 13, 
  height = 4)
gridExtra::grid.arrange(bar_plot_table)
dev.off()



# Figure 1E ----
# Dot plot heatmap with select genes
baby_marker_genes <- rev(c(
  "FHL2", "PECAM1", "POSTN", "BMX", "LYVE1", "ACKR1", "PDGFRA", "SEMA3C", "PDGFRB",
  "MYH11", "ADIPOQ", "PTPRC", "CD3E", "CD8A", "CD4", "FOXP3", "GNLY", "MS4A1", "MZB1", "CD68", "CD163", "SPP1", 
  "MS4A2", "ITGAX", "XCR1", "CD1C", "CD83", "LILRA4", "MKI67"))

dotplot <- DotPlot(clustered_obj, group.by = "ct_second_pass", features = baby_marker_genes)

# Extract DotPlot data
dotplot_data <- dotplot$data

# Scaled expression levels
exp_mat <- dotplot_data %>% 
  select(-pct.exp, -avg.exp) %>%  
  pivot_wider(names_from = features.plot, values_from = avg.exp.scaled) %>% 
  column_to_rownames(var = "id") %>%
  as.matrix() 

# The percentage of cells expressing a feature
percent_mat <- dotplot_data %>% 
  select(-avg.exp, -avg.exp.scaled) %>%  
  pivot_wider(names_from = features.plot, values_from = pct.exp) %>% 
  column_to_rownames(var = "id") %>%
  as.matrix() 

# Ordering cell types
exp_mat_reordered <- exp_mat[c(
                               
                               "Endothelial", "Lymphatic endothelial", "Activated endothelial", 
                               "BMX+ Activated endothelial", "Proliferating endothelial", 
                               
                               "Cardiomyocytes", "Fibroblasts", "POSTN+ Fibroblasts", 
                               "Myofibroblasts", "Pericytes", "Proliferating pericytes", "vSMCs", "Adipocytes", 
                               
                               "CD8+ T cells", "CD4+ T cells", "Treg", 
                               "Proliferating T cells", "NK", "B cells", "Plasma", 
                               
                               
                               "Macrophages", "SPP1+ Macrophages", "Mast", "cDC1", "cDC2", 
                               "mDC", "pDC", "Proliferating DCs"), ]

# Ordering cell types
percent_mat_reordered <- percent_mat[c( 
                                       
                                       "Endothelial", "Lymphatic endothelial", "Activated endothelial", 
                                       "BMX+ Activated endothelial", "Proliferating endothelial", 
                                       
                                       "Cardiomyocytes", "Fibroblasts", "POSTN+ Fibroblasts", 
                                       "Myofibroblasts", "Pericytes", "Proliferating pericytes", 
                                       "vSMCs", "Adipocytes", 
                                       
                                       "CD8+ T cells", "CD4+ T cells", "Treg", 
                                       "Proliferating T cells", "NK", "B cells", "Plasma", 
                                       
                                       
                                       "Macrophages", "SPP1+ Macrophages", "Mast", "cDC1", "cDC2", 
                                       "mDC", "pDC", "Proliferating DCs"), ]

# Any value that is greater than 2 will be mapped to bright red
col_fun = circlize::colorRamp2(c(-2, 0, 2), colorspace::diverge_hsv(3))

# Updating lineage order
ct_col_order <- factor(c(
  rep("Endothelial", 5),
  rep("Stromal", 8),
  rep("Immune", 15)),
  levels = c("Endothelial", "Stromal", "Immune"))

# Creating a layer to add to plot
layer_fun = function(j, i, x, y, w, h, fill) {
  grid.points(x = x, y = y, 
              gp = gpar(col = col_fun(pindex(t(exp_mat_reordered), i, j))), 
              size = pindex(t(percent_mat_reordered), i, j)/100 * unit(1, "mm"),
              pch = 19)
}

# Create legend list for proportion and expression
lgd_list = list(
  Legend(labels = c(0,0.25,0.5,0.75,1), title = "Proportion", type = "points", 
         #pch = 19, size = c(0,0.25,0.5,0.75,1) * unit(1.4, "mm"),
         pch = 19, size = c(0,0.25,0.5,0.75,1) * unit(1, "mm"),
         legend_gp = gpar(col = "black"), direction = "horizontal", ncol = 1,
         title_position = "topcenter", background = NA, 
         grid_height = unit(2, "mm"),
         grid_width = unit(4, "mm"),
         title_gp = gpar(fontsize = 4), labels_gp = gpar(fontsize = 4),
         #legend_height = unit(5, "mm"),
         legend_height = unit(2, "mm"),
         legend_width = unit(10, "mm")
  )
)

# Heatmap parameters
hp <- Heatmap(t(exp_mat_reordered),
              column_title_gp = gpar(fontsize = 4),
              column_title_side = "bottom",
              column_split = ct_col_order,
              col = col_fun,
              rect_gp = gpar(type = "none"),
              layer_fun = layer_fun,
              row_names_gp = gpar(fontsize = 4, fontface = "italic"),
              row_names_side = "right",
              column_names_gp = gpar(fontsize = 4),
              column_dend_height = unit(2, "mm"),
              column_dend_side = "bottom",
              border = "black",
              column_names_rot = 90,
              heatmap_legend_param = list(title = "Scaled\nexpression",
                                          labels = c("-2", "-1", "0", "1", "2"),
                                          legend_direction = "horizontal",
                                          title_position = "topcenter",
                                          title_gp = gpar(fontsize = 4),
                                          labels_gp = gpar(fontsize = 4),
                                          background = NA,
                                          grid_width = unit(10, "mm"),
                                          grid_height = unit(2, "mm"),
                                          legend_height = unit(2, "mm"),
                                          legend_width = unit(10, "mm")),
              cluster_rows = FALSE, 
              cluster_columns = FALSE,
              cluster_row_slices = FALSE,
              show_heatmap_legend = TRUE)

ht <- draw(hp, 
           heatmap_legend_list = lgd_list,
           heatmap_legend_side = "bottom",
           annotation_legend_side = "bottom")


pdf("/home/aoill/plots/00_heart_manuscript/revisions_01/Figure_01E_new.pdf", 
    width = 2.75, height = 3.25)
ht
dev.off()


# Figure 1F ----
# Cell type proportion bar plot stratified by biopsy level rejection category
# cell_prop_by_rej_type IS SUPPLEMENTARY TABLE 1. SEE Propeller.R
# FOR HOW THIS FILE WAS GENERATED
# Read in propeller results
cell_prop_by_rej_type <- read_csv("/home/aoill/projects/heart_transplant/propeller_pre-treatment_by_biopsy_rejection_type_correct.csv")


# Transform the data into long format
desired_cell_type_order <- c(
  
  # Endothelial
  "Endothelial", "Activated endothelial", "BMX+ Activated endothelial", 
  "Proliferating endothelial", "Lymphatic endothelial", 
  
  # Immune
  "CD8+ T cells", "Macrophages", "NK", "Plasma", 
  "cDC2", "CD4+ T cells", "B cells", "Proliferating T cells", 
  "cDC1", "Treg", "SPP1+ Macrophages", "Mast", "pDC", "mDC", 
  "Proliferating DCs", 
  
  # Stromal
  "Cardiomyocytes", "Fibroblasts", "Pericytes", "POSTN+ Fibroblasts", "vSMCs", 
  "Myofibroblasts", "Adipocytes", "Proliferating pericytes"
)


cell_prop_long <- cell_prop_by_rej_type %>%
  pivot_longer(
    cols = starts_with("PropMean"),
    names_to = "rejection_type",
    values_to = "proportion"
  ) %>%
  mutate(
    rejection_type = str_replace(rejection_type, "PropMean\\.", ""),
    significance = ifelse(FDR < 0.10, "***", "")
  ) %>%
  dplyr::rename(
    cell_type = `...1` # Use dplyr::rename if masking is an issue
  ) %>%
  mutate(
    rejection_type = case_when(
      rejection_type == "no_rejection" ~ "No rejection",
      rejection_type == "mixed_rejection" ~ "Mixed rejection",
      rejection_type == "cellular_rejection" ~ "ACR",
      rejection_type == "antibody_rejection." ~ "AMR"
    ),
    cell_type = factor(cell_type, levels = desired_cell_type_order)
  ) %>%
  dplyr::select(-BaselineProp, -Fstatistic, -P.Value)


p1 <- cell_prop_long %>%
  ggplot(aes(x = rejection_type, y = proportion, fill = cell_type)) +
  geom_bar(stat = "identity", position = "fill") +
  scale_fill_manual(values = mycolors) + # This will use your defined colors
  cowplot::theme_cowplot() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 16),
    axis.text.y = element_text(size = 16),
    axis.title.x = element_text(size = 16),
    axis.title.y = element_text(size = 16),
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 12),
    panel.grid = element_blank()
  ) +
  scale_y_continuous(labels = scales::percent) +
  labs(
    title = "",
    x = "",
    y = "Proportion of Cells",
    fill = "Cell Type"
  ) +
  guides(fill = guide_legend(ncol = 1))

p1


pdf("/home/aoill/plots/00_heart_manuscript/revisions_01/Figure_01F.pdf", 
    width = 7.25, height = 8.75)
p1
dev.off()


# Figure 2 ----
# SEE calc_cellproximity_heart_spatial.R AND Fisher_exact_ct_wPlot_heart_spatial.R
# FOR HOW INPUT FILES WERE GENERATED.
# Figure 2A ----
enrich_score <- readRDS("/home/aoill/projects/heart_transplant/00_final/cell_prox_CR_biopsies_pre.rds")

pval <- reshape2::dcast(prox_to ~ prox_ct, value.var = 'fdr', data = enrich_score) 
or <- reshape2::dcast(prox_to ~ prox_ct, value.var = 'logOR', data = enrich_score)

rownames(pval) <- pval$prox_to
pval$prox_to <- NULL

rownames(or) <- or$prox_to
or$prox_to <- NULL

ct_pval<- pval
ct_pval[is.na(ct_pval)] = 1

ct_or<- or
ct_or[is.na(ct_or)] = 0
ct_or[is.infinite(as.matrix(ct_or))] <- 0

mean(rownames(ct_or) == rownames(ct_pval))
mean(colnames(ct_or) == colnames(ct_pval))

ct_or[ct_pval > 0.05] = 0 # replace non-sig results with 0


calc_heatmap_w <- function(df){width <- ncol(df) * 0.24; return(width)}
calc_heatmap_h <- function(df){height <- nrow(df) * 0.18; return(height)}

# plot odds ratio
min01 <- quantile(as.matrix(ct_or), 0.01)
max99 <- quantile(as.matrix(ct_or), 0.99)
col_fun = colorRamp2(c(min01, 0, max99), c("#377EB8", "white", "#E41A1C"))

hp <- Heatmap(t(ct_or), name = 'log(OR)', 
              row_title = 'Starting Cell', 
              row_title_side = 'right', 
              row_names_side = 'left', 
              row_dend_side = 'right',  
              row_title_gp = gpar(fontsize = 15), 
              
              column_title = 'Nearest Neighbor Cell', 
              column_title_side = "bottom",
              column_names_side = 'top', 
              column_names_rot = 45,
              
              column_dend_side = 'bottom',  
              
              column_title_gp = gpar(fontsize = 15),
              
              rect_gp = gpar(col = "#b8b7b6", lwd = 0.5),
              width = unit(calc_heatmap_w(ct_or), "in"), 
              height = unit(calc_heatmap_h(ct_or), "in"),
              na_col = 'white', 
              col = col_fun, 
              #row_dend_width = unit(0.5, 'cm'), 
              #column_dend_height = unit(0.5, 'cm')
) 


pdf("/home/aoill/plots/Figure_02A.pdf",
    width = 12,
    height = 10)
draw(hp)
dev.off()


# Figure 2C ----
enrich_score <- readRDS("/home/aoill/projects/heart_transplant/00_final/cell_prox_AMR_biopsies_pre.rds")


# Plotting
pval <- reshape2::dcast(prox_to ~ prox_ct, value.var = 'fdr', data = enrich_score) 
or <- reshape2::dcast(prox_to ~ prox_ct, value.var = 'logOR', data = enrich_score)

rownames(pval) <- pval$prox_to
pval$prox_to <- NULL

rownames(or) <- or$prox_to
or$prox_to <- NULL

ct_pval<- pval
ct_pval[is.na(ct_pval)] = 1

ct_or<- or
ct_or[is.na(ct_or)] = 0
ct_or[is.infinite(as.matrix(ct_or))] <- 0


mean(rownames(ct_or) == rownames(ct_pval))
mean(colnames(ct_or) == colnames(ct_pval))

ct_or[ct_pval > 0.05] = 0


calc_heatmap_w <- function(df){width <- ncol(df) * 0.24; return(width)}
calc_heatmap_h <- function(df){height <- nrow(df) * 0.18; return(height)}

# plot odds ratio
min01 <- quantile(as.matrix(ct_or), 0.01)
max99 <- quantile(as.matrix(ct_or), 0.99)
col_fun = colorRamp2(c(min01, 0, max99), c("#377EB8", "white", "#E41A1C"))

hp <- Heatmap(t(ct_or), name = 'log(OR)', 
              row_title = 'Starting Cell', 
              row_title_side = 'right', 
              row_names_side = 'left', 
              row_dend_side = 'right',  
              row_title_gp = gpar(fontsize = 15), 
              
              column_title = 'Nearest Neighbor Cell', 
              column_title_side = "bottom",
              column_names_side = 'top', 
              column_names_rot = 45,
              
              column_dend_side = 'bottom',  
              
              column_title_gp = gpar(fontsize = 15),
              
              rect_gp = gpar(col = "#b8b7b6", lwd = 0.5),
              width = unit(calc_heatmap_w(ct_or), "in"), 
              height = unit(calc_heatmap_h(ct_or), "in"),
              na_col = 'white', 
              col = col_fun, 
              #row_dend_width = unit(0.5, 'cm'), 
              #column_dend_height = unit(0.5, 'cm')
) 


pdf("/home/aoill/plots/Figure_02C.pdf",
    width = 12,
    height = 10)
draw(hp)
dev.off()


# Figure 3A ----
# Generated in biorender


# Figure 3B ----
# Transcriptional similarity analysis heatmap
# SEE Transcriptional_similarity_analysis.R FOR HOW INPUT FILE WAS GENERATED.
fe_results_all <- readRDS("/scratch/aoill/projects/heart_transplant/transcriptional_similarity_results.rds")

pval <- reshape2::dcast(cell_cell_type ~ neighbor, value.var = 'fdr', data = fe_results_all) 
or <- reshape2::dcast(cell_cell_type ~ neighbor, value.var = 'logOR', data = fe_results_all)
# some ORs are -INF so might need to replace -INF with 0

rownames(pval) <- pval$cell_cell_type
pval$cell_cell_type <- NULL

rownames(or) <- or$cell_cell_type
or$cell_cell_type <- NULL

ct_pval<- pval

ct_or<- or


# Replace Inf and -Inf with placeholders for plotting
max_val <- max(ct_or[is.finite(as.matrix(ct_or))], na.rm = TRUE)
min_val <- min(ct_or[is.finite(as.matrix(ct_or))], na.rm = TRUE)
ct_or[is.infinite(as.matrix(ct_or)) & ct_or > 0] <- max(ct_or[is.finite(as.matrix(ct_or))], na.rm = TRUE) + .1  # Replace Inf
ct_or[is.infinite(as.matrix(ct_or)) & ct_or < 0] <- min(ct_or[is.finite(as.matrix(ct_or))], na.rm = TRUE) - .1  # Replace -Inf

#ct_or

# Define significance threshold
sig_threshold <- 0.05

# Create a matrix for annotations (TRUE = significant, FALSE = not significant)
significance <- ct_pval < sig_threshold


# reorder matricies
# Reorder columns: 3, 2, 1, 0
ct_or <- ct_or[, c("3", "2", "1", "0")]
significance <- significance[, c("3", "2", "1", "0")]



## Subset by source cell grade ----

### Source biopsy grade 0 ----
# Select rows where rownames start with "0_"
ct_or_0 <- ct_or[grep("^0_", rownames(ct_or)), ]
significance_0 <- significance[grep("^0_", rownames(significance)), ]

# drop neighbor results that are 0
ct_or_0 <- ct_or_0[, colnames(ct_or_0) != "0"]
significance_0 <- significance_0[, colnames(significance_0) != "0"]



ct_or_0 <- ct_or_0[, c("1", "2", "3")]
significance_0 <- significance_0[, c("1", "2", "3")]

ht1 <-
  Heatmap(ct_or_0, name = 'log(OR)', 
          row_title = '', 
          row_title_side = 'left', 
          row_title_gp = gpar(fontsize = 8, fontface = "bold"), 
          column_names_gp = gpar(fontsize = 8),
          row_names_gp = gpar(fontsize = 8),
          column_title = "Source biopsy grade 0", 
          column_names_side = 'bottom', 
          #column_dend_side = 'bottom',  
          column_title_gp = gpar(fontsize = 8, fontface = "bold"),
          rect_gp = gpar(col = "black", lwd = 0.5),
          #width = unit(calc_heatmap_w(ct_or), "in"), 
          #height = unit(calc_heatmap_h(ct_or), "in"),
          na_col = 'white', 
          col = col_fun, 
          cluster_rows = FALSE,
          cluster_columns = FALSE,
          #column_names_rot = 45, 
          row_dend_width = unit(0.5, 'cm'), 
          column_dend_height = unit(0.5, 'cm'),
          heatmap_legend_param = list(
            title = "log(OR)",
            at = c(min_val - 1, min_val, 0, max_val, max_val + 1),  # Define breaks
            labels = c("-Inf", round(min_val, digits = 2), "0", round(max_val, digits = 2), "Inf"),
            title_gp = gpar(fontsize = 8),  # Legend title font size
            labels_gp = gpar(fontsize = 8)  # Legend labels font size
          ),
          cell_fun = function(j, i, x, y, width, height, fill) {
            if (significance_0[i, j]) {
              grid.text("*", x, y, gp = gpar(fontsize = 16, col = "black"))
            }
          }
  )

### Source biopsy grade 1 ----
ct_or_1 <- ct_or[grep("^1_", rownames(ct_or)), ]
significance_1 <- significance[grep("^1_", rownames(significance)), ]

# drop neighbor results that are 0
ct_or_1 <- ct_or_1[, colnames(ct_or_1) != "0"]
significance_1 <- significance_1[, colnames(significance_1) != "0"]


ct_or_1 <- ct_or_1[, c("1", "2", "3")]
significance_1 <- significance_1[, c("1", "2", "3")]

ht2 <- Heatmap(ct_or_1, name = 'log(OR)', 
               row_title = '', 
               row_title_side = 'right', 
               row_title_gp = gpar(fontsize = 8, fontface = "bold"), 
               column_names_gp = gpar(fontsize = 8),
               row_names_gp = gpar(fontsize = 8),
               column_title = "Source biopsy grade 1", 
               column_names_side = 'bottom', 
               #column_dend_side = 'bottom',  
               column_title_gp = gpar(fontsize = 8, fontface = "bold"),
               rect_gp = gpar(col = "black", lwd = 0.5),
               #width = unit(calc_heatmap_w(ct_or), "in"), 
               #height = unit(calc_heatmap_h(ct_or), "in"),
               na_col = 'white', 
               col = col_fun, 
               cluster_rows = FALSE,
               cluster_columns = FALSE,
               #column_names_rot = 45, 
               row_dend_width = unit(0.5, 'cm'), 
               column_dend_height = unit(0.5, 'cm'),
               heatmap_legend_param = list(
                 title = "log(OR)",
                 at = c(min_val - 1, min_val, 0, max_val, max_val + 1),  # Define breaks
                 labels = c("-Inf", round(min_val, digits = 2), "0", round(max_val, digits = 2), "Inf"),
                 title_gp = gpar(fontsize = 8),  # Legend title font size
                 labels_gp = gpar(fontsize = 8)  # Legend labels font size
               ),
               cell_fun = function(j, i, x, y, width, height, fill) {
                 if (significance_1[i, j]) {
                   grid.text("*", x, y, gp = gpar(fontsize = 16, col = "black"))
                 }
               }
)



### Source biopsy grade 2 ----
ct_or_2 <- ct_or[grep("^2_", rownames(ct_or)), ]
significance_2 <- significance[grep("^2_", rownames(significance)), ]
# drop neighbor results that are 0
ct_or_2 <- ct_or_2[, colnames(ct_or_2) != "0"]
significance_2 <- significance_2[, colnames(significance_2) != "0"]

ct_or_2 <- ct_or_2[, c("1", "2", "3")]
significance_2 <- significance_2[, c("1", "2", "3")]

ht3 <- Heatmap(ct_or_2, name = 'log(OR)', 
               row_title = '', 
               row_title_side = 'right', 
               row_title_gp = gpar(fontsize = 8, fontface = "bold"), 
               column_names_gp = gpar(fontsize = 8),
               row_names_gp = gpar(fontsize = 8),
               column_title = "Source biopsy grade 2", 
               column_names_side = 'bottom', 
               #column_dend_side = 'bottom',  
               column_title_gp = gpar(fontsize = 8, fontface = "bold"),
               rect_gp = gpar(col = "black", lwd = 0.5),
               #width = unit(calc_heatmap_w(ct_or), "in"), 
               #height = unit(calc_heatmap_h(ct_or), "in"),
               na_col = 'white', 
               col = col_fun, 
               cluster_rows = FALSE,
               cluster_columns = FALSE,
               #column_names_rot = 45, 
               row_dend_width = unit(0.5, 'cm'), 
               column_dend_height = unit(0.5, 'cm'),
               heatmap_legend_param = list(
                 title = "log(OR)",
                 at = c(min_val - 1, min_val, 0, max_val, max_val + 1),  # Define breaks
                 labels = c("-Inf", round(min_val, digits = 2), "0", round(max_val, digits = 2), "Inf"),
                 title_gp = gpar(fontsize = 8),  # Legend title font size
                 labels_gp = gpar(fontsize = 8)  # Legend labels font size
               ),
               cell_fun = function(j, i, x, y, width, height, fill) {
                 if (significance_2[i, j]) {
                   grid.text("*", x, y, gp = gpar(fontsize = 16, col = "black"))
                 }
               }
)


### Source biopsy grade 3 ----
ct_or_3 <- ct_or[grep("^3_", rownames(ct_or)), ]
significance_3 <- significance[grep("^3_", rownames(significance)), ]

# drop neighbor results that are 0
ct_or_3 <- ct_or_3[, colnames(ct_or_3) != "0"]
significance_3 <- significance_3[, colnames(significance_3) != "0"]


ct_or_3_orig <- ct_or_3
significance_3_orig <- significance_3

rownames(ct_or_3) <- gsub("3_", "", rownames(ct_or_3))
rownames(significance_3) <- gsub("3_", "", rownames(significance_3))


ct_or_3 <- ct_or_3[, c("1", "2", "3")]
significance_3 <- significance_3[, c("1", "2", "3")]

ht4 <- Heatmap(ct_or_3, name = 'log(OR)', 
               row_title = 'source', 
               row_title_side = 'right', 
               row_title_gp = gpar(fontsize = 8, fontface = "bold"), 
               column_names_gp = gpar(fontsize = 8),
               row_names_gp = gpar(fontsize = 8),
               column_title = "Source biopsy grade 3", 
               column_names_side = 'bottom', 
               #column_dend_side = 'bottom',  
               column_title_gp = gpar(fontsize = 8, fontface = "bold"),
               rect_gp = gpar(col = "black", lwd = 0.5),
               #width = unit(calc_heatmap_w(ct_or), "in"), 
               #height = unit(calc_heatmap_h(ct_or), "in"),
               na_col = 'white', 
               col = col_fun, 
               cluster_rows = FALSE,
               cluster_columns = FALSE,
               #column_names_rot = 45, 
               row_dend_width = unit(0.5, 'cm'), 
               column_dend_height = unit(0.5, 'cm'),
               heatmap_legend_param = list(
                 title = "log(OR)",
                 at = c(min_val - 1, min_val, 0, max_val, max_val + 1),  # Define breaks
                 labels = c("-Inf", round(min_val, digits = 2), "0", round(max_val, digits = 2), "Inf"),
                 title_gp = gpar(fontsize = 5),  # Legend title font size
                 labels_gp = gpar(fontsize = 5)  # Legend labels font size
               ),
               cell_fun = function(j, i, x, y, width, height, fill) {
                 if (significance_3[i, j]) {
                   grid.text("*", x, y, gp = gpar(fontsize = 16, col = "black"))
                 }
               }
)


# Combine
ht_list = ht1 + ht2 + ht3 + ht4


draw(ht_list, 
     row_title = "Source cell type", 
     row_title_side = "right",
     row_title_gp = gpar(fontsize = 5, fontface = "bold"),
     column_title = "Nearest neighbor biopsy grade", 
     column_title_side = "bottom",
     column_title_gp = gpar(fontsize = 5),
     
)

pdf("/home/aoill/plots/00_heart_manuscript/Figure_03.pdf",
    width = 7.5,
    height = 7)
draw(ht_list, 
     row_title = "Source cell type", 
     row_title_side = "right",
     row_title_gp = gpar(fontsize = 8, fontface = "bold"),
     column_title = "Nearest neighbor biopsy grade", 
     column_title_side = "bottom",
     column_title_gp = gpar(fontsize = 8),
     
)
dev.off()


# Figure 3C ----
# Pseudo-bulked expression PCA 
# SEE Expression_PCA.R FOR HOW INPUT FILE WAS GENERATED
smpl.pca <- readRDS("/home/aoill/projects/heart_transplant/00_final/smpl.pca_CR_pre_pseudo_all_together.rds")

p1 <- fviz_pca_ind(smpl.pca,
                   label = "none", # hide individual labels
                   habillage = as.factor(mean_expression_per_sample_info$patient_rejection_type), # color by groups
                   addEllipses = F, # Concentration ellipses
                   title = "",
                   pointsize = 3,
                   mean.point = FALSE
) +
  labs(color = "Rejection type", shape = "Rejection type")  +
  coord_fixed() + theme_classic() +
  theme(
    axis.text = element_text(size = 14),
    axis.title = element_text(size = 16),
    legend.text = element_text(size = 14),  
    legend.title = element_text(size = 16)
  ) 

p1



pdf("/home/aoill/plots/00_heart_manuscript/PCA_expression_pseudo_all_cts_patient_level_rejection_type.pdf",
    width = 7,
    height = 7)
p1
dev.off()


# Figure 3D ----
# Expression correlation heat map
clustered_obj_pre <- subset(clustered_obj, subset = biopsy_timing == "pre")

# Pseudo-bulk
# Extract the normalized expression data (data is log-transformed counts)
expr_matrix <- as.data.frame(as.matrix(clustered_obj_pre@assays$RNA@data))
expr_matrix_t <- t(expr_matrix)
expr_matrix_t <- as.data.frame(expr_matrix_t)

# Add sample information to the data
expr_matrix_t$patient_id <- clustered_obj_pre@meta.data$patient_id

# Summarize by Sample to calculate the mean expression for each gene and sample
mean_expression_per_sample <- expr_matrix_t %>%
  group_by(patient_id) %>%
  summarise(across(everything(), mean, na.rm = TRUE))

# Convert back to a data frame 
mean_expression_per_sample <- as.data.frame(mean_expression_per_sample)
rownames(mean_expression_per_sample) <- mean_expression_per_sample$patient_id
mean_expression_per_sample$patient_id <- NULL
#head(t(mean_expression_per_sample))
avg_expression_df <- t(mean_expression_per_sample)

# Add in sample info for plotting
# make a col called age_cat for pediatric or adult
clustered_obj_pre@meta.data$age_cat <- ifelse(clustered_obj_pre@meta.data$transplant_age < 18, "Pediatric", "Adult")

sample_data <- clustered_obj_pre@meta.data %>% 
  dplyr::select(patient_id, patient_rejection_type, patient_cellular_grading, patient_antibody_grading,
                age_cat) %>%
  unique()
#outcome_data$Sample <- gsub("_", "-", outcome_data$Sample)
head(sample_data)
rownames(sample_data) <- NULL

# add to expression matrix
mean_expression_per_sample$patient_id <- rownames(mean_expression_per_sample)
head(mean_expression_per_sample)
mean_expression_per_sample_info <- left_join(mean_expression_per_sample, sample_data,
                                             by = "patient_id")

# add rownames back in
rownames(mean_expression_per_sample_info) <- mean_expression_per_sample_info$patient_id
head(mean_expression_per_sample_info)

# run correlation
cor_matrix <- cor(avg_expression_df, method = "pearson")

# Set sample IDs as row names if not already
sample_data_tmp <- sample_data %>%
  column_to_rownames(var = "patient_id")  # Ensure patient_id colnames in cor_matrix

# Make sure the order of cor_matrix is the same as annotation_df
annotation_df <- annotation_df[colnames(cor_matrix), , drop = FALSE]

colnames(annotation_df) <- c("Rejection Type", "ACR Grade", "AMR Grade", "Age Category")

# Define colors for annotations
anno_colors <- list(
  `Age Category` = c("Pediatric" = "blue", "Adult" = "red"),
  `ACR Grade` = c("0" = "#28b463", "1" = "#ca6f1e", "2" = "#cc0000", 
                  "3" = "#0b5394", "NA" = "grey"),
  #`AMR Grade` = c("pAMR1-h" = "#DAF7A6", "pAMR1-i" = "#FFC300", 
  #                    "pAMR0" = "#FF5733", "pAMR2" = "#C70039", "NA" = "grey"),
  `AMR Grade` = c("pAMR0" = "#28b463", "pAMR1-h" = "#ca6f1e", 
                  "pAMR1-i" = "#cc0000", "pAMR2" = "#0b5394", "NA" = "grey"),
  
  `Rejection Type` = c("cellular_rejection" = "green", 
                       "mixed_rejection" = "purple", 
                       "antibody_rejection" = "orange")
)

# Create a HeatmapAnnotation object
col_anno <- HeatmapAnnotation(df = annotation_df, 
                              col = anno_colors)

ht <- Heatmap(cor_matrix, 
              name = "Correlation",
              top_annotation = col_anno,  # Add sample metadata annotations
              #col = colorRamp2(c(-1, 0, 1), c("blue", "white", "red")),  # Correlation color scale
              col = colorRamp2(c(min(cor_matrix), max(cor_matrix)), c("white", "blue")),  # Correlation color scale
              clustering_distance_rows = "euclidean",
              clustering_distance_columns = "euclidean",
              clustering_method_rows = "complete",
              clustering_method_columns = "complete",
              show_column_names = F,
              show_row_names = F)


pdf("/home/aoill/plots/00_heart_manuscript/pre_treatment_expression_correlation_heatmap_small_20250225.pdf",
    width = 8.25, height = 6)
draw(ht)
dev.off()


# Figure 4A ----
# DGE analysis pre-treatment samples AMR vs ACR
# THIS IS SUPPLEMENTARY TABLE 3. SEE T_test_pre_AMR_CR.R
# FOR HOW THIS FILE WAS GENERATED
# Volcano plot
acr_v_amr <- read_csv("/home/aoill/projects/heart_transplant/00_final/revisions/t_test_AMR_vs_CR_results.csv")


vp <- ggplot(acr_v_amr, aes(x = Log2FC, y = -log10(padj), color = cell_type, alpha = padj < 0.05)) +
  geom_point() +
  cowplot::theme_cowplot() +
  scale_alpha_manual(values = c(0.1, 1)) +
  scale_color_manual(values = mycolors, name = "Cell Type") +
  labs(x = "LogFC", y = "-log10(Adj. p-value)") +
  scale_x_continuous(breaks = seq(-4, 4, by = 2)) +
  geom_segment(data = NULL, inherit.aes = FALSE,
               aes(x = 0.25, y = 6, xend = 4, yend = 6),
               arrow = arrow(type = "closed", length = unit(0.1, "inches")),
               color = "black") +
  annotate("text", x = 1.8, y = 6.4, label = "Upregulated in ACR", size = 4.5, hjust = 0.5) +
  geom_segment(data = NULL, inherit.aes = FALSE,
               aes(x = -0.25, y = 6, xend = -4, yend = 6),
               arrow = arrow(type = "closed", length = unit(0.1, "inches")),
               color = "black") +
  annotate("text", x = -1.8, y = 6.4, label = "Upregulated in AMR", size = 4.5, hjust = 0.5) + 
  geom_hline(yintercept=-log10(0.05), linetype="dashed")


pdf("/home/aoill/plots/00_heart_manuscript/revisions_01/t_test_volcano_AMRvsACR.pdf",
    width = 9, height = 6)
vp
dev.off()


# Figure 4B ----
# DGE analysis pre-treatment samples AMR vs ACR
# Bar plot of the number of significant DEGs among cell types
fdr.acr_v_amr <- acr_v_amr %>%
  filter(padj < 0.05)

pdf("/home/aoill/plots/00_heart_manuscript/revisions_01/t_test_num_genes_per_ct_AMRvsACR_new.pdf",
    width = 8, height = 5
)
ggplot(fdr.acr_v_amr %>%
         count(cell_type) %>%
         arrange(desc(n)),
       aes(x = reorder(cell_type, -n), y = n, fill = cell_type)) + 
  geom_bar(stat = "identity") + 
  scale_fill_manual(values = mycolors) +
  labs(x = "Cell Type", y = "Number of DEGs") +
  cowplot::theme_cowplot() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  Seurat::NoLegend()
dev.off()


# Figure 4C ----
# DGE analysis pre-treatment samples AMR vs ACR
# Gene set enrichment analysis (GSEA) depicting results of hallmark pathways 
# that showed significance at an FDR < 0.1.
# gsea_results.acr_v_amr IS SUPPLEMENTARY TABLE 4. SEE GSEA_from_DEG_pre_AMR_CR.R
# FOR HOW THIS FILE WAS GENERATED 
gsea_results.acr_v_amr <- read.csv("/home/aoill/projects/heart_transplant/00_final/20250805_gsea_ACR_v_AMR.csv")

# Filter pathways with padj < 0.1
filtered_data <- gsea_results.acr_v_amr %>%
  filter(padj < 0.1) %>%
  mutate(pathway = fct_reorder(pathway, NES))

# Create lollipop plot
pdf("/home/aoill/plots/00_heart_manuscript/revisions_01/GSEA_AMR_ACR_lollipop_plot_revised_01.pdf",
    width = 8, height = 4.5)
ggplot(filtered_data, aes(x = NES, y = pathway, color = cell_type)) +
  geom_segment(aes(x = 0, xend = NES, y = pathway, yend = pathway), color = "grey") +
  geom_point(size = 5) +
  scale_color_manual(values = mycolors, name = "Cell Type") +
  cowplot::theme_cowplot() +
  labs(
    x = "Normalized Enrichment Score (NES)",
    y = "Pathway"
  ) +
  theme(
    legend.position = "right",
    axis.text.y = element_text(size = 12),
    axis.text.x = element_text(size = 12),
    panel.grid.major.y = element_blank(), # Remove y grid lines for cleaner look
    panel.grid.minor = element_blank()
  )
dev.off()


# Figure 4D ----
# DGE ACR pre-treatment resolution vs persistent (responders vs nonresponders)
# THIS IS SUPPLEMENTARY TABLE 7. SEE lm_pre_ACR_responders_nonresponders.R
# FOR HOW THIS FILE WAS GENERATED
acr_resolved_v_persistent <- read_csv("/home/aoill/projects/heart_transplant/00_final/revisions/lm_ACR_resp_vs_non_resp_results.csv")

# Visualize with limited x-axis
vp2 <- ggplot(acr_resolved_v_persistent, aes(x = estimate, y = -log10(padj), 
                                      color = cell_type, alpha = padj < 0.05)) +
  geom_point() +
  cowplot::theme_cowplot() +
  scale_alpha_manual(values = c(0.1, 1)) +
  scale_color_manual(values = mycolors, name = "Cell Type") +
  labs(x = "β", y = "-log10(Adj. p-value)") +
  scale_x_continuous(limits = c(-1, 1)) +
  geom_segment(data = NULL, inherit.aes = FALSE,
               aes(x = 0.25, y = 3.5, xend = 0.9, yend = 3.5),
               arrow = arrow(type = "closed", length = unit(0.1, "inches")),
               color = "black") +
  annotate("text", x = 0.5, y = 3.8, label = "Upregulated in Non-Responders", size = 4.5, hjust = 0.5) +
  geom_segment(data = NULL, inherit.aes = FALSE,
               aes(x = -0.25, y = 3.5, xend = -0.9, yend = 3.5),
               arrow = arrow(type = "closed", length = unit(0.1, "inches")),
               color = "black") +
  annotate("text", x = -0.5, y = 3.8, label = "Upregulated in Responders", size = 4.5, hjust = 0.5) + NoLegend() +
  geom_hline(yintercept=-log10(0.05), linetype="dashed")



vp2

pdf("/home/aoill/plots/00_heart_manuscript/revisions_01/lm_volcano_ACR_resp_nonresp_revised_01.pdf",
    width = 7, height = 9)
print(vp2)
dev.off()

# Figure 4E ----
# DGE ACR pre-treatment resolution vs persistent (responders vs nonresponders)
# Bar plot of the number of significant DEGs among cell types
acr_resolved_v_persistent_fdr_05 <- acr_resolved_v_persistent %>%
  filter(adj.P.Val < 0.05)

pdf("/home/aoill/plots/00_heart_manuscript/revisions_01/DGE_num_genes_per_ct_ACR_resp_nonresp_revised_01.pdf",
    #width = 8, 
    width = 5, height = 4
)
ggplot(acr_resolved_v_persistent_fdr_05 %>%
         count(cell_type) %>%
         arrange(desc(n)),
       aes(x = reorder(cell_type, -n), y = n, fill = cell_type)) + 
  geom_bar(stat = "identity") + 
  scale_fill_manual(values = mycolors) +
  labs(x = "Cell Type", y = "Number of DEGs") +
  cowplot::theme_cowplot() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  Seurat::NoLegend()
dev.off() 


# Figure 5 ----
# Figure 5A ----
# Not generated with code

# Figure 5B ----
# SEE GLMM_CAV_*.R FOR HOW FILES WERE GENERATED 
# THESE ARE SIPPLEMENTAL TABLES 9-13.
glmm_acr_pre_cav <- read_csv("/home/aoill/projects/heart_transplant/00_final/revisions/GLMM_CAV_grade_CR_pre_treatment_biopsies_revisions_fixed.csv")
glmm_amr_post_cav <- read_csv("/home/aoill/projects/heart_transplant/00_final/revisions/GLMM_CAV_grade_AMR_post_treatment_biopsies_revisions_fixed.csv")
glmm_amr_pre_cav <- read_csv("/home/aoill/projects/heart_transplant/00_final/revisions/GLMM_CAV_grade_AMR_pre_treatment_biopsies_revisions_fixed.csv")
glmm_acr_post_cav <- read_csv("/home/aoill/projects/heart_transplant/00_final/revisions/GLMM_CAV_grade_CR_post_treatment_biopsies_revisions_fixed.csv")
glmm_mixed_pre_cav <- read_csv("/home/aoill/projects/heart_transplant/00_final/revisions/GLMM_CAV_grade_Mixed_pre_treatment_biopsies_revisions_fixed.csv")

# Create the list of all genes from the GLMM models
gene.list.ALL <- bind_rows(glmm_acr_pre_cav, glmm_acr_post_cav,
                           glmm_amr_pre_cav, glmm_amr_post_cav,
                           glmm_mixed_pre_cav)

# Highlight the top 20 genes (by absolute value of beta)
top_genes <- gene.list.ALL %>%
  filter(padj < 0.05) %>%
  mutate(composite_score = abs(beta) * -log10(padj)) %>% 
  arrange(desc(composite_score)) %>%
  slice_head(n = 20)

# Visualize
pdf("/home/aoill/plots/00_heart_manuscript/revisions_01/CAV_GLMM_all_volcano_revised_01.pdf",
    width = 10, height = 6
)
ggplot(gene.list.ALL, aes(x = beta, 
                          y = -log10(padj), 
                          color = ct, 
                          alpha = padj < 0.05)) +
  geom_point() +
  cowplot::theme_cowplot() +
  ggrepel::geom_text_repel(
    data = top_genes,
    aes(label = genei),
    size = 4,
    max.overlaps = 40,
    box.padding = 0.5,
    point.padding = 0.5,
    segment.color = "black",
    segment.size = 0.4
  ) +
  scale_alpha_manual(values = c(0.05, 1)) +
  scale_color_manual(values = mycolors, name = "Cell Type") +
  labs(x = "β", y = "-log10(Adj. p-value)")
dev.off()


# Figure 5C ----
# Choosing the top 50 genes that are all FDR < 0.05
top_genes.lollipop <- gene.list.ALL %>%
  filter(padj < 0.05) %>%
  arrange(desc(abs(beta))) %>%
  distinct(genei, .keep_all = TRUE) %>%
  slice_head(n = 50)
length(unique(top_genes.lollipop$genei)) # 50 unique genes

# Visualize
pdf("/home/aoill/plots/00_heart_manuscript/revisions_01/betas_CAV_top50_lollipop_revised_01.pdf",
    width = 7, height = 9
)
ggplot(top_genes.lollipop, aes(x = reorder(genei, desc(beta)), y = beta, fill = ct)) +
  geom_segment(aes(xend = genei, y = 0, yend = beta), color = "black") +
  geom_point(size = 5, shape = 21, color = "black") +
  coord_flip() +
  scale_fill_manual(values = mycolors) +
  cowplot::theme_cowplot() +
  labs(x = "Gene", y = "β", fill = "Cell Type")
dev.off()



# Extended Data Figures ----
# Extended Data Figure 1 ----
# Top markers + canonical markers
# Order info for heat map
new_ct_order2 <- c("Endothelial", "Lymphatic endothelial", "Activated endothelial", 
                   "BMX+ Activated endothelial", "Proliferating endothelial", 
                   
                   "Cardiomyocytes", "Fibroblasts", "POSTN+ Fibroblasts", 
                   "Myofibroblasts", "Pericytes", "Proliferating pericytes", "vSMCs", "Adipocytes", 
                   
                   "CD8+ T cells", "CD4+ T cells", "Treg", 
                   "Proliferating T cells", "NK", "B cells", "Plasma", 
                   
                   
                   "Macrophages", "SPP1+ Macrophages", "Mast", "cDC1", "cDC2", 
                   "mDC", "pDC", "Proliferating DCs"
)

# Finding cluster markers for a heatmap/dotplot
markers <- presto::wilcoxauc(clustered_obj, group_by = "ct_second_pass")
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

all_markers <- unique(c(all_markers, 
                        "CD8A", "CD4","BMX", "CD68", "SPP1", "CD1C", "PTPRC", 
                        "PECAM1", "EPCAM"))


# log-normalizing and scaling all features in the RNA assay. Scaling so that
# all features can be visualized using the same color scale
DefaultAssay(clustered_obj) <- "RNA"
clustered_obj <- NormalizeData(clustered_obj)
VariableFeatures(clustered_obj) <- rownames(clustered_obj)
clustered_obj <- ScaleData(clustered_obj)

p <- DotPlot(clustered_obj, features = all_markers, 
             group.by = "ct_second_pass", scale = TRUE) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
#p

# Reproducing the dotplot using ComplexHeatmap
df <- p$data
head(df)

# (Scaled) expression levels
exp_mat <- df %>% 
  dplyr::select(-pct.exp, -avg.exp) %>%  
  tidyr::pivot_wider(names_from = id, values_from = avg.exp.scaled) %>% 
  as.data.frame() %>%
  dplyr::select(features.plot, new_ct_order2)

row.names(exp_mat) <- exp_mat$features.plot
exp_mat <- exp_mat[,-1] %>% as.matrix()

# The percentage of cells expressing a feature
percent_mat <- df %>% 
  dplyr::select(-avg.exp, -avg.exp.scaled) %>%  
  tidyr::pivot_wider(names_from = id, values_from = pct.exp) %>% 
  as.data.frame() %>%
  dplyr::select(features.plot, new_ct_order2)

row.names(percent_mat) <- percent_mat$features.plot  
percent_mat <- percent_mat[,-1] %>% as.matrix()

dim(exp_mat)
dim(percent_mat)

# Get an idea of the ranges of the matrix
quantile(exp_mat, c(0.1, 0.5, 0.9, 0.99), na.rm=T)

# Replace NA with 0 in matrix
exp_mat[!is.finite(exp_mat)] <- 0
quantile(exp_mat, c(0.1, 0.5, 0.9, 0.99))

# Any value that is greater than 2 will be mapped to yellow
col_fun = circlize::colorRamp2(c(-1, 0, 2), viridis(20)[c(1,10, 20)])



mycolors_lst <- c(list(cluster = mycolors))
cluster_anno <- colnames(exp_mat)
column_ha <- HeatmapAnnotation( # HeatmapAnnotation / rowAnnotation
  cluster = cluster_anno,
  col = mycolors_lst,
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
  Legend( labels = c(0,0.25,0.5,0.75,1), 
          title = "Proportion", 
          type = "points", 
          pch = 16, 
          size = c(0,0.25,0.5,0.75,1) * unit(4,"mm"),
          legend_gp = gpar(col = "black"),
          background = NA,
          title_gp = gpar(fontsize = 8),
          labels_gp = gpar(fontsize = 8)
  )
)


ct_col_order <- factor(c(
  rep("Endothelial", 5),
  rep("Stromal", 8),
  rep("Immune", 15)),
  levels = c("Endothelial", "Stromal", "Immune"))

hp <- Heatmap((exp_mat), # t(exp_mat)
              heatmap_legend_param=list(ncol = 1, 
                                        title="Scaled Exp.",
                                        legend_direction = "vertical",
                                        title_position = "topcenter",
                                        title_gp = gpar(fontsize = 8),
                                        labels_gp = gpar(fontsize = 8),
                                        background = NA#,
                                        #grid_width = unit(2, "mm"),
                                        #grid_height = unit(6, "mm"),
                                        #legend_height = unit(10, "mm"),
                                        #legend_width = unit(2, "mm")
              ),
              
              #column_title = " ", 
              col = col_fun,
              rect_gp = gpar(type = "none"),
              layer_fun = layer_fun,
              
              row_names_gp = gpar(fontsize = 8),
              row_title_side = "right",
              row_names_side = "left",
              
              column_title_gp = gpar(fontsize = 9),
              column_names_gp = gpar(fontsize = 8),
              column_names_rot = 35,
              border = "black",
              cluster_columns = F,
              column_split = ct_col_order
)

pdf("/home/aoill/plots/00_heart_manuscript/revisions_01/Extended_data_01_new.pdf",
    width = 7, height = 10)
draw(hp,
     heatmap_legend_list = lgd_list,
     heatmap_legend_side = "right",
     annotation_legend_side = "right",
     align_heatmap_legend = "global_center"
)
dev.off()


# Supplementary Figures ----

# Supplementary Figure 2 ----
# Box plot comparing our data to Heart cell atlas
# SEE scRNA_heart_props.ipynb FOR HOW INPUT DATA WAS GENERATED
# Read in proportion data (gathered in python because data was in H5AD format)
heart_atlas_props <- read.csv("/scratch/aoill/projects/heart_transplant/scrna_atlas/donor_lineage_proportions.csv")

# Add a column called data_set
heart_atlas_props$data_set <- "heart_cell_atlas"

head(heart_atlas_props)
# change Donor.ID to Sample to match spatial data set (for plotting)
colnames(heart_atlas_props) <- c("Sample", "lineage", "Proportion", "data_set")
head(heart_atlas_props)

# Get cell type proportions for each biopsy and then in the data_set col, label
# each biopsy spatial_rejection, spatial_no_rejection
clustered_obj <- readRDS("/scratch/aoill/projects/heart_transplant/GEO/heart_spatial_obj.rds")

# pull out meta 
conf_clustered_obj_meta <- clustered_obj@meta.data

# Remove rows with NA in the 'biopsy_rejection_type' column
conf_clustered_obj_meta <- conf_clustered_obj_meta %>%
  filter(!is.na(biopsy_rejection_type))

unique(conf_clustered_obj_meta$biopsy_rejection_type)

# Add data_set column for plotting
conf_clustered_obj_meta <- conf_clustered_obj_meta %>%
  mutate(data_set = ifelse(biopsy_rejection_type == "no_rejection", 
                           "spatial_no_rejection", 
                           "spatial_rejection"))

# check result
conf_clustered_obj_meta %>%
  dplyr::select(biopsy_rejection_type, data_set) %>%
  tail()


# shorten the data to only the necessary columns. This just makes it easier for
# me when I am checking the files
conf_clustered_obj_meta_short <- conf_clustered_obj_meta %>%
  dplyr::select(Sample, patient_id, patient_biopsy, ct_second_pass, biopsy_rejection_type, data_set, lineage)

head(conf_clustered_obj_meta_short)
# replace mesenchymal with stromal in lineage colunm
conf_clustered_obj_meta_short$lineage <- str_replace_all(conf_clustered_obj_meta_short$lineage, "Mesenchymal", "Stromal")


# Get cell type proportions for each Sample and save to a data frame along with
# the lineage and data_set columns

# Calculate lineage proportions for each Sample
lineage_proportions <- conf_clustered_obj_meta_short %>%
  group_by(Sample, lineage, data_set) %>%  # Group by Sample, lineage, and data_set
  summarise(count = n(), .groups = "drop") %>%  # Count the number of cells
  group_by(Sample, data_set) %>%  # Group by Sample and data_set to calculate proportions
  mutate(total_count = sum(count),
         Proportion = count / total_count) %>%  # Calculate proportions
  select(Sample, lineage, Proportion, data_set)  # Select relevant columns

lineage_proportions <- as.data.frame(lineage_proportions)
# check result
head(lineage_proportions)
head(heart_atlas_props)


# combine data
all_props <- rbind(heart_atlas_props, lineage_proportions)


all_props$data_set <- factor(
  all_props$data_set, 
  levels = c("spatial_rejection", "spatial_no_rejection", "heart_cell_atlas")
)

# Remove neural because we do not have a compariable cell type
all_props_filter <- all_props %>% 
  filter(lineage != "Neural")

# Plot all results
pdf("/home/aoill/plots/00_heart_manuscript/revisions_01/atlas_comparison_boxplot_new.pdf",
    width = 7, height = 2.75)
ggplot(all_props_filter, aes(x = data_set, y = Proportion, fill = data_set)) +
  geom_boxplot() +
  facet_wrap(~ lineage, scales = "free_y", nrow = 1) +  # Facet by lineage
  theme_classic() +
  ylim(0, 1) +
  labs(
    x = "Data Set",
    y = "Proportion"
  ) +
  theme(
    #axis.text.x = element_text(angle = 45, hjust = 1), 
    strip.text = element_text(size = 10, face = "bold"),
    axis.text=element_text(size=8),
    axis.title=element_text(size=10)
  ) +
  scale_fill_brewer(palette = "Set2") +
  scale_x_discrete(labels= c("Spatial\n(rejection)", "Spatial\n(no rejection)", "HCA")) +
  NoLegend()
dev.off()



# Supplementary Figure 3 ----
# Scree plot 
smpl.pca <- readRDS("/home/aoill/projects/heart_transplant/00_final/smpl.pca_CR_pre_pseudo_all_together.rds")

# Extract eigenvalues/variances
get_eig(smpl.pca)
p2 <- fviz_screeplot(smpl.pca, addlabels = TRUE, ylim = c(0, 50))
ggsave("/home/aoill/plots/00_heart_manuscript/PCA_scree_expression_pseudo_all_cts_patient_level_rejection_type.png", plot = p2, width = 160, height = 100, units = "mm", dpi = 800)


# Supplementary Figure 4 ----
# PC loadings of the top 20 genes by their absolute loadings in each of the first 3 PCs.
smpl.pca <- readRDS("/home/aoill/projects/heart_transplant/00_final/smpl.pca_CR_pre_pseudo_all_together.rds")
# Get PCA loadings (coordinates of variables)
loadings <- as.data.frame(smpl.pca$var$coord) 

# Add gene names
loadings$gene <- rownames(loadings)

# Reshape for plotting (long format)
library(tidyr)
loadings_long <- pivot_longer(loadings, cols = starts_with("Dim"), 
                              names_to = "PC", values_to = "loading")

# Clean PC names
loadings_long$PC <- gsub("Dim.", "PC", loadings_long$PC)

# Select top 20 genes for each of the first 3 PCs
top_genes <- loadings_long %>%
  filter(PC %in% c("PC1", "PC2", "PC3")) %>%  # Keep only PC1, PC2, PC3
  group_by(PC) %>%
  arrange(desc(abs(loading))) %>%  # Sort by absolute loading
  slice(1:20) %>%  # Take top 20 genes per PC
  ungroup()

# Plot as bar plots
p4 <- ggplot(top_genes, aes(x = reorder(gene, loading), y = loading, fill = loading > 0)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c("blue", "red")) +  # Negative = blue, Positive = red
  coord_flip() +  # Flip to make it horizontal
  facet_wrap(~PC, scales = "free_y") +  # One panel per PC
  theme_classic() +
  labs(x = "Gene", y = "Loading") +
  theme(legend.position = "none", strip.background = element_rect(fill = "grey"))  # Remove legend

ggsave("/home/aoill/plots/00_heart_manuscript/expression_PCA_loadings_barplot.png", 
       plot = p4, width = 200, height = 150, units = "mm", dpi = 800)


# Save loading
#write.csv(loadings, "/home/aoill/projects/heart_transplant/00_final/expression_PCA_loadings.csv", 
#          quote = F, row.names = F)


# Supplementary Figure 5 ----
# Pseudo-bulked expression for each cell type separately 
# By cell type
clustered_obj_pre <- subset(clustered_obj, subset = biopsy_timing == "pre")

cts_to_analyze <- levels(as.factor(clustered_obj_pre@meta.data$ct_second_pass))
pca_plots <- list()  # Empty list to store plots


custom_colors <- c("antibody_rejection" = "#F8766D", 
                   "cellular_rejection" = "#00BFC4",  
                   "mixed_rejection" = "#7CAE00")
# Define custom shapes 
custom_shapes <- c("antibody_rejection" = 16,
                   "cellular_rejection" = 17,
                   "mixed_rejection" = 15)

for (ctn in 1:length(cts_to_analyze)) {
  ct <- cts_to_analyze[ctn]
  print(ct)
  
  # subset object to cell type
  obj_CT <- subset(clustered_obj_pre, subset = ct_second_pass == ct)
  
  # Extract the normalized expression data (data is log-transformed counts)
  expr_matrix <- as.data.frame(as.matrix(obj_CT@assays$RNA@data))
  expr_matrix_t <- t(expr_matrix)
  expr_matrix_t <- as.data.frame(expr_matrix_t)
  
  #all_genes <- colnames(expr_matrix_t)
  
  # Add sample information to the data
  expr_matrix_t$patient_id <- obj_CT@meta.data$patient_id
  
  # Summarize by Sample to calculate the mean expression for each gene and sample
  mean_expression_per_sample <- expr_matrix_t %>%
    group_by(patient_id) %>%
    summarise(across(everything(), mean, na.rm = TRUE))
  
  # Convert back to a data frame 
  mean_expression_per_sample <- as.data.frame(mean_expression_per_sample)
  rownames(mean_expression_per_sample) <- mean_expression_per_sample$patient_id
  mean_expression_per_sample$patient_id <- NULL
  #head(t(mean_expression_per_sample))
  avg_expression_df <- t(mean_expression_per_sample)
  
  # mean_expression_per_sample will be used for PCA but need to add rejection type
  # info for plotting
  mean_expression_per_sample$patient_id <- rownames(mean_expression_per_sample)
  
  tmp <- obj_CT@meta.data %>%
    dplyr::select(patient_id, patient_biopsy, biopsy_timing, #rejection_type, 
                  sex, donor_sex, #cellular_grading, antibody_grading, 
                  treatment,
                  patient_cellular_grading, patient_antibody_grading, patient_rejection_type) %>%
    unique()
  mean_expression_per_sample_info <- left_join(mean_expression_per_sample, tmp)
  rownames(mean_expression_per_sample_info) <- rownames(mean_expression_per_sample)
  
  
  ## All samples, pesudo-bulked by patient and biopsy timing across all cell types together
  
  smpl.pca <- PCA(mean_expression_per_sample_info[,1:477], graph = F)
  
  
  # Generate PCA plot
  pca_plots[[ct]] <- fviz_pca_ind(smpl.pca,
                                  label = "none", # hide individual labels
                                  habillage = as.factor(mean_expression_per_sample_info$patient_rejection_type), # color by groups
                                  addEllipses = F, # Concentration ellipses
                                  title = ct
  ) +
    labs(color = "Rejection type", shape = "Rejection type") +
    scale_color_manual(values = custom_colors) +
    scale_shape_manual(values = custom_shapes) +
    theme(plot.title = element_text(hjust = 0.5))
  
  
}

# Arrange all plots into a grid
pdf("/home/aoill/plots/00_heart_manuscript/pca_exp_cell_type.pdf",
    width = 17,
    height = 19)
ggarrange(plotlist = pca_plots, ncol = 5, nrow = ceiling(length(pca_plots) / 5),
          common.legend = TRUE, legend = "bottom")
dev.off()


# Supplementary Figure 7 ----
# Distribution of TMA grade in ACR patients whose rejection resolved 
# (responders) and those whose rejection persisted (nonresponders) 
# Subset data
conf_clustered_obj_pre <- subset(clustered_obj, subset = biopsy_timing == "pre" )
conf_clustered_obj_pre_CR <- subset(conf_clustered_obj_pre, subset = patient_rejection_type == "cellular_rejection")

# Pull out columns of interest into a data frame
patient_info_df <- conf_clustered_obj_pre_CR@meta.data %>%
  dplyr::select(patient_id, Sample, biopsy_cellular_grading, patient_cellular_grading, 
                resolution_broad) %>% 
  unique() %>%
  na.omit() 
rownames(patient_info_df) <- NULL
head(patient_info_df)

# Biopsy level grading
p1 <- ggplot(patient_info_df, aes(x=biopsy_cellular_grading, fill=resolution_broad)) +
  geom_histogram(alpha=0.5, position="identity", binwidth = 1) +
  xlab("TMA cellular grade") +
  theme_classic()

# Plot
pdf("/home/aoill/plots/00_heart_manuscript/revisions_01/resolution_dist_TMA_grade_ACR.pdf",
    width = 5, height = 4)
print(p1)
dev.off()

# also save as a png
ggsave("/home/aoill/plots/00_heart_manuscript/revisions_01/resolution_dist_TMA_grade_ACR.png", 
       plot = p1, width = 110, height = 72, units = "mm", dpi = 800)



