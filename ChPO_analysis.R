############################################################
## ChP RNA-seq Analysis
##
## Groups:
##   Br  -> CC
##   ChP -> ChP
##   M   -> ChPO
##
## Colors:
##   CC = Green
##   ChP   = Pink
##   ChPO  = Blue
############################################################


############################################################
## 1. Load libraries
############################################################

library(Rsubread)
library(limma)
library(edgeR)
library(DESeq2)
library(ggplot2)
library(ggforce)
library(pheatmap)
library(grid)


############################################################
## 2. Set working directory
############################################################

directory <- "/Volumes/AtsadangSSD/Pin_re/"
setwd(directory)


############################################################
## 3. Create output folder
############################################################

plot_dir <- file.path(getwd(), "Plots_check")

if (!dir.exists(plot_dir)) {
  dir.create(
    plot_dir,
    recursive = TRUE
  )
}

cat(
  "\nPlots will be saved to:\n",
  plot_dir,
  "\n"
)


############################################################
## 4. Read raw count data
############################################################

df <- read.table(
  "rawcount.csv",
  header = TRUE,
  sep = "\t",
  check.names = FALSE
)

head(df)


############################################################
## 5. Clean original sample names
############################################################

old_names <- colnames(df)[3:ncol(df)]

## Remove ".sorted.bam"
old_names <- gsub(
  "\\.sorted\\.bam$",
  "",
  old_names
)

## Remove SRR information
old_names <- sub(
  "_SRR.*",
  "",
  old_names
)

old_names


############################################################
## 6. Assign sample groups
##
## Br  -> CC
## ChP -> ChP
## M   -> ChPO
############################################################

sample_group <- ifelse(
  grepl("^Br", old_names, ignore.case = TRUE),
  "CC",
  ifelse(
    grepl("^ChP", old_names, ignore.case = TRUE),
    "ChP",
    ifelse(
      grepl("^M", old_names, ignore.case = TRUE),
      "ChPO",
      NA
    )
  )
)


############################################################
## Check group assignment
############################################################

sample_check <- data.frame(
  Original_name = old_names,
  Group = sample_group
)

print(sample_check)

if (any(is.na(sample_group))) {
  stop(
    "Some sample names could not be assigned to CC, ChP, or ChPO."
  )
}


############################################################
## 7. Create continuous sample numbers
##
## CC1, CC2, ...
## ChP1, ChP2, ...
## ChPO1, ChPO2, ...
############################################################

running_number <- ave(
  seq_along(sample_group),
  sample_group,
  FUN = seq_along
)

new_names <- paste0(
  sample_group,
  running_number
)


############################################################
## Check old and new names
############################################################

sample_name_table <- data.frame(
  Original_name = old_names,
  New_name = new_names,
  Group = sample_group
)

print(sample_name_table)


############################################################
## Save sample name conversion
############################################################

write.csv(
  sample_name_table,
  file = file.path(
    plot_dir,
    "Sample_name_conversion.csv"
  ),
  row.names = FALSE
)


############################################################
## Replace names
############################################################

colnames(df)[3:ncol(df)] <- new_names

colnames(df)


############################################################
## 8. Aggregate counts by GeneID
############################################################

df1 <- aggregate(
  df[, 3:ncol(df)],
  by = list(
    GeneID = df$GeneID
  ),
  FUN = sum
)

head(df1)


############################################################
## 9. Create count matrix
############################################################

rownames(df1) <- df1$GeneID

countData <- df1[, -1]

countData <- as.matrix(
  countData
)

mode(countData) <- "integer"

dim(countData)

head(countData)

colnames(countData)


############################################################
## 10. Create DESeq2 sample information
############################################################

colData <- data.frame(
  row.names = colnames(countData),
  
  group = ifelse(
    grepl("^CC", colnames(countData)),
    "CC",
    
    ifelse(
      grepl("^ChPO", colnames(countData)),
      "ChPO",
      "ChP"
    )
  )
)


############################################################
## Set group order
############################################################

colData$group <- factor(
  colData$group,
  levels = c(
    "CC",
    "ChP",
    "ChPO"
  )
)

colData

table(
  colData$group
)


############################################################
## 11. Create DESeq2 dataset
############################################################

dds <- DESeqDataSetFromMatrix(
  countData = countData,
  colData = colData,
  design = ~ group
)


############################################################
## 12. Filter low-count genes
############################################################

dds <- dds[
  rowSums(counts(dds)) >= 10,
]


############################################################
## 13. Run DESeq2
############################################################

dds <- DESeq(
  dds
)


############################################################
## 14. Variance Stabilizing Transformation
############################################################

vsd <- vst(
  dds,
  blind = TRUE
)

vsd_mat <- assay(
  vsd
)


############################################################
## 15. MASTER group colors
############################################################

group_colors <- c(
  "CC" = "#7CAE00",   # Green
  "ChP"   = "#F8766D",   # Pink
  "ChPO"  = "#00BFC4"    # Blue
)


############################################################
## PCA shapes
############################################################

group_legend_labels <- c(
  "CC" = "CC, cerebral cortex",
  "ChP" = "ChP, choroid plexus",
  "ChPO" = "ChPO, choroid plexus organoids"
)

group_shapes <- c(
  "CC" = 16,
  "ChP"   = 15,
  "ChPO"  = 17
)


############################################################
## Heatmap annotation
############################################################

annotation <- data.frame(
  Group = factor(
    group_legend_labels[as.character(colData$group)],
    levels = unname(group_legend_labels)
  )
)

rownames(annotation) <- rownames(
  colData
)


############################################################
## Heatmap annotation colors
############################################################

annotation_colors <- list(
  Group = setNames(unname(group_colors), unname(group_legend_labels))
)


############################################################
## 16. Helper function for large viewing window
############################################################

# Output-first rendering: avoid mixing RStudio graphics and grid viewports.
# Heatmaps are constructed silently and written to dedicated PDF/PNG devices.
# Open the files in Plots/ to inspect them.
open_plot_window <- function(width = 12, height = 10) {
  invisible(NULL)  # retained for compatibility; deliberately opens no device
}

safe_pheatmap <- function(..., filename = NA, width = 13, height = 10) {
  hm <- pheatmap::pheatmap(..., silent = TRUE)
  if (length(filename) == 1L && !is.na(filename)) {
    extension <- tolower(tools::file_ext(filename))
    if (extension == "pdf") {
      grDevices::pdf(filename, width = width, height = height,
                     family = "Helvetica", onefile = FALSE)
    } else if (extension == "png") {
      grDevices::png(filename, width = width, height = height,
                     units = "in", res = 600, bg = "white")
    } else {
      stop("Heatmap output must be .pdf or .png: ", filename)
    }
    tryCatch({
      grid::grid.newpage()
      grid::grid.draw(hm$gtable)
    }, finally = grDevices::dev.off())
    message("Saved heatmap: ", filename)
  }
  invisible(hm)
}


############################################################
## 17. PCA
############################################################

pcaData <- plotPCA(
  vsd,
  intgroup = "group",
  returnData = TRUE
)

percentVar <- round(
  100 *
    attr(
      pcaData,
      "percentVar"
    ),
  1
)

pcaData$name <- rownames(
  pcaData
)


############################################################
## PCA plot - larger, clearer typography
## Use polygon ellipses rather than geom_mark_ellipse()
## to avoid ggforce/grid depth(NULL) rendering errors.
############################################################

make_pca_ellipse <- function(data, radius = 3.2, padding = 8,
                             n_points = 120) {
  pieces <- lapply(split(data, data$group), function(sub) {
    xy <- as.matrix(sub[, c("PC1", "PC2")])
    if (nrow(xy) < 3L || any(!is.finite(xy))) return(NULL)
    center <- colMeans(xy)
    covariance <- stats::cov(xy)
    eig <- eigen(covariance, symmetric = TRUE)
    # Floor small eigenvalues for stable contours with few samples.
    spread <- sqrt(pmax(eig$values, 1)) * radius + padding
    angles <- seq(0, 2 * pi, length.out = n_points)
    unit_circle <- rbind(cos(angles), sin(angles))
    coords <- eig$vectors %*% (diag(spread, 2L) %*% unit_circle)
    data.frame(
      PC1 = coords[1, ] + center[1],
      PC2 = coords[2, ] + center[2],
      group = as.character(sub$group[1])
    )
  })
  do.call(rbind, pieces)
}

pca_ellipse_data <- make_pca_ellipse(pcaData)

p <- ggplot(
  pcaData,
  aes(x = PC1, y = PC2, color = group, fill = group, shape = group)
) +
  geom_polygon(
    data = pca_ellipse_data,
    aes(x = PC1, y = PC2, group = group, fill = group),
    inherit.aes = FALSE,
    alpha = 0.15,
    color = NA,
    show.legend = FALSE
  ) +
  geom_point(size = 7, stroke = 1.2) +
  scale_color_manual(values = group_colors, labels = group_legend_labels) +
  scale_fill_manual(values = group_colors, labels = group_legend_labels) +
  scale_shape_manual(values = group_shapes, labels = group_legend_labels) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  ggtitle("PCA of VST-normalized counts") +
  theme_bw(base_size = 23) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5, size = 24),
    axis.title.x = element_text(
      size = 25, face = "bold", margin = margin(t = 12)
    ),
    axis.title.y = element_text(
      size = 25, face = "bold", margin = margin(r = 12)
    ),
    axis.text.x = element_text(size = 20, color = "black"),
    axis.text.y = element_text(size = 20, color = "black"),
    axis.line = element_line(linewidth = 0.6, color = "black"),
    panel.grid.minor = element_blank(),
    legend.title = element_blank(),
    legend.key = element_blank(),
    legend.text = element_text(size = 23, color = "black")
  ) +
  coord_fixed(xlim = c(-150, 150), ylim = c(-150, 150))

############################################################
## Show PCA
############################################################

## Preview PCA by opening Plots/PCA_VST.pdf after saving.
## No RStudio device is opened (prevents mixed grid/graphics errors).

############################################################
## Save PCA as PDF (vector)
############################################################

ggsave(
  filename = file.path(plot_dir, "PCA_VST2.pdf"),
  plot = p,
  width = 12,
  height = 7
)

############################################################
## Save PCA as PNG (600 dpi)
############################################################

ggsave(
  filename = file.path(plot_dir, "PCA_VST2.png"),
  plot = p,
  width = 12,
  height = 7,
  dpi = 600,
  bg = "white"
)


############################################################
## 17B. CENTROID CORRELATION AND SAMPLE CORRELATION
## Uses the existing CC / ChP / ChPO labels and colors
############################################################

## Match sample identifiers to columns of the VST matrix
brain_samples <- rownames(colData)[colData$group == "CC"]
chp_samples   <- rownames(colData)[colData$group == "ChP"]
org_samples   <- rownames(colData)[colData$group == "ChPO"]

stopifnot(
  all(c(brain_samples, chp_samples, org_samples) %in% colnames(vsd_mat)),
  all(c("CC", "ChP", "ChPO") %in% as.character(colData$group))
)

## Group centroids: mean VST expression for each gene
brain_centroid <- rowMeans(vsd_mat[, brain_samples, drop = FALSE])
chp_centroid   <- rowMeans(vsd_mat[, chp_samples, drop = FALSE])
org_centroid   <- rowMeans(vsd_mat[, org_samples, drop = FALSE])

centroid_mat <- cbind(
  CC = brain_centroid,
  ChP   = chp_centroid,
  ChPO  = org_centroid
)

centroid_pearson  <- cor(centroid_mat, method = "pearson")
centroid_spearman <- cor(centroid_mat, method = "spearman")

cat("\nPearson correlation between group centroids:\n")
print(round(centroid_pearson, 4))
cat("\nSpearman correlation between group centroids:\n")
print(round(centroid_spearman, 4))

write.csv(centroid_pearson,
          file.path(plot_dir, "Centroid_correlation_Pearson.csv"))
write.csv(centroid_spearman,
          file.path(plot_dir, "Centroid_correlation_Spearman.csv"))

## Select up to 2,000 most variable genes across all samples
## (plotPCA uses 500 by default, for comparison)
gene_var <- apply(vsd_mat, 1, var)
top_n <- min(2000L, sum(is.finite(gene_var)))
if (top_n < 2L) stop("Not enough variable genes for correlation heatmap")
top_genes <- names(sort(gene_var[is.finite(gene_var)],
                        decreasing = TRUE))[seq_len(top_n)]

## Keep sample names in the same order as existing figures
sample_order <- c(brain_samples, chp_samples, org_samples)
corMat_hvg <- cor(vsd_mat[top_genes, sample_order, drop = FALSE],
                  method = "pearson")

## Reuse existing sample labels, group colors and heatmap palette
cor_annotation <- annotation[sample_order, , drop = FALSE]
cor_colors <- colorRampPalette(c(
  "#4575B4", "#D9EDF7", "#F5F5B5", "#F4A259", "#D73027"
))(100)

## Correlations generally high: use actual range so differences are visible
cor_range <- range(corMat_hvg, na.rm = TRUE)
cor_breaks <- seq(min(cor_range[1], 0.999999),
                  max(cor_range[2], 1.000001), length.out = 101)

## One plotting function for screen / PDF / PNG (identical formatting)
plot_correlation_heatmap <- function(filename = NA) {
  safe_pheatmap(
    corMat_hvg,
    color = cor_colors,
    breaks = cor_breaks,
    annotation_col = cor_annotation,
    annotation_row = cor_annotation,
    annotation_colors = annotation_colors,
    cluster_rows = TRUE,
    cluster_cols = TRUE,
    clustering_method = "complete",
    border_color = "#999999",
    show_rownames = TRUE,
    show_colnames = TRUE,
    angle_col = 90,
    fontsize = 16,
    fontsize_row = 14,
    fontsize_col = 14,
    treeheight_row = 90,
    treeheight_col = 90,
    cellwidth = 32,
    cellheight = 32,
    display_numbers = TRUE,
    number_format = "%.2f",
    fontsize_number = 8,
    legend = TRUE,
    annotation_legend = TRUE,
    annotation_names_col = TRUE,
    annotation_names_row = TRUE,
    main = "Pearson correlation (top 2000 variable genes)",
    filename = filename,
    width = 15,
    height = 10
  )
}

# Construct silently; save both files using a fresh device for each.
plot_correlation_heatmap(file.path(plot_dir, "Pearson_correlation_top2000_2.pdf"))
plot_correlation_heatmap(file.path(plot_dir, "Pearson_correlation_top2000_2.png"))
write.csv(corMat_hvg,
          file.path(plot_dir, "Pearson_correlation_top2000_matrix.csv"))

## Permutation test: ChP-ChPO similarity minus CC-ChPO similarity
set.seed(123)
n_perm <- 1000L
obs_diff <- mean(corMat_hvg[chp_samples, org_samples, drop = FALSE]) -
  mean(corMat_hvg[brain_samples, org_samples, drop = FALSE])

## Permute sample group labels while preserving each group size
perm_diffs <- replicate(n_perm, {
  shuffled <- sample(as.character(colData$group))
  names(shuffled) <- rownames(colData)
  chp_s <- names(shuffled)[shuffled == "ChP"]
  brain_s <- names(shuffled)[shuffled == "CC"]
  org_s <- names(shuffled)[shuffled == "ChPO"]
  mean(corMat_hvg[chp_s, org_s, drop = FALSE]) -
    mean(corMat_hvg[brain_s, org_s, drop = FALSE])
})

## Two-sided permutation p-value (+1 correction)
p_value <- (sum(abs(perm_diffs) >= abs(obs_diff)) + 1) /
  (n_perm + 1)

cat("\nObserved correlation difference (ChP-ChPO minus CC-ChPO): ",
    obs_diff, "\n", sep = "")
cat("Two-sided permutation p-value: ", p_value, "\n", sep = "")

write.csv(data.frame(
  observed_difference = obs_diff,
  permutation_p_value = p_value,
  number_of_permutations = n_perm,
  selected_genes = length(top_genes)
), file.path(plot_dir, "Correlation_permutation_result.csv"),
row.names = FALSE)

############################################################
## 18. SAMPLE-TO-SAMPLE DISTANCE
############################################################

sampleDists <- dist(
  t(vsd_mat)
)

sampleDistMatrix <- as.matrix(
  sampleDists
)

rownames(sampleDistMatrix) <- colnames(
  vsd_mat
)

colnames(sampleDistMatrix) <- colnames(
  vsd_mat
)


############################################################
## Distance heatmap colors
############################################################

dist_colors <- colorRampPalette(
  c(
    "#4575B4",
    "#D9EDF7",
    "#F2E394",
    "#F4A259",
    "#D73027"
  )
)(100)


############################################################
## Display sample-to-sample heatmap
############################################################

open_plot_window(
  12,
  10
)

# Preview omitted: open the corresponding saved PDF or PNG in Plots/.


############################################################
## Save sample-to-sample heatmap PDF
############################################################

safe_pheatmap(
  sampleDistMatrix,
  
  color = dist_colors,
  
  annotation_col = annotation,
  annotation_row = annotation,
  
  annotation_colors = annotation_colors,
  
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  
  border_color = "#999999",
  
  show_rownames = TRUE,
  show_colnames = TRUE,
  
  angle_col = 90,
  
  fontsize = 16,
  fontsize_row = 14,
  fontsize_col = 14,
  
  treeheight_row = 90,
  treeheight_col = 90,
  
  cellwidth = 32,
  cellheight = 32,
  
  legend = TRUE,
  annotation_legend = TRUE,
  
  annotation_names_col = TRUE,
  annotation_names_row = TRUE,
  
  main = "Sample-to-sample distance",
  
  filename = file.path(
    plot_dir,
    "Sample_to_sample_distance.pdf"
  ),
  
  width = 15,
  height = 10
)


############################################################
## Save sample-to-sample heatmap PNG
############################################################

safe_pheatmap(
  sampleDistMatrix,
  
  color = dist_colors,
  
  annotation_col = annotation,
  annotation_row = annotation,
  
  annotation_colors = annotation_colors,
  
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  
  border_color = "#999999",
  
  show_rownames = TRUE,
  show_colnames = TRUE,
  
  angle_col = 90,
  
  fontsize = 16,
  fontsize_row = 14,
  fontsize_col = 14,
  
  treeheight_row = 90,
  treeheight_col = 90,
  
  cellwidth = 32,
  cellheight = 32,
  
  legend = TRUE,
  annotation_legend = TRUE,
  
  annotation_names_col = TRUE,
  annotation_names_row = TRUE,
  
  main = "Sample-to-sample distance",
  
  filename = file.path(
    plot_dir,
    "Sample_to_sample_distance.png"
  ),
  
  width = 15,
  height = 10
)


############################################################
## 19. HEATMAP COLOR SCALE
##
## Same style for both marker heatmaps
############################################################

marker_colors <- colorRampPalette(
  c(
    "#4575B4",
    "#91BFDB",
    "#E0F3F8",
    "#F5F5B5",
    "#F1D38A",
    "#F39C5A",
    "#D73027"
  )
)(100)

marker_breaks <- seq(
  -2.3,
  2.3,
  length.out = 101
)


############################################################
## 20. ORIGINAL CHOROID PLEXUS MARKER HEATMAP
############################################################

markers_original <- c(
  "TTR",
  "AQP1",
  "CLIC6",
  "KCNJ13",
  "OTX2",
  "IGFBP7",
  "KRT18",
  "KRT19",
  "EPCAM",
  "CLDN1",
  "CLDN2",
  "CLDN3",
  "TJP1",
  "MSX1",
  "LMX1A",
  "FOXJ1",
  "RSPH1",
  "DNAH5"
)


############################################################
## Genes found / missing
############################################################

markers_original_present <- intersect(
  markers_original,
  rownames(vsd)
)

markers_original_missing <- setdiff(
  markers_original,
  rownames(vsd)
)

cat(
  "\nOriginal ChP markers found:\n"
)

print(
  markers_original_present
)

cat(
  "\nOriginal ChP markers missing:\n"
)

print(
  markers_original_missing
)


############################################################
## Extract expression
############################################################

marker_original_mat <- assay(vsd)[
  markers_original_present,
  ,
  drop = FALSE
]


############################################################
## Row Z-score
############################################################

marker_original_scaled <- t(
  scale(
    t(marker_original_mat)
  )
)


############################################################
## Remove NA genes
############################################################

marker_original_scaled <- marker_original_scaled[
  complete.cases(marker_original_scaled),
  ,
  drop = FALSE
]


############################################################
## DISPLAY ORIGINAL MARKER HEATMAP
############################################################

open_plot_window(
  12,
  10
)

# Preview omitted: open the corresponding saved PDF or PNG in Plots/.


############################################################
## Save ORIGINAL marker heatmap PDF
############################################################

safe_pheatmap(
  marker_original_scaled,
  
  color = marker_colors,
  breaks = marker_breaks,
  
  scale = "none",
  
  annotation_col = annotation,
  
  annotation_colors = annotation_colors,
  
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  
  clustering_method = "complete",
  
  border_color = "#999999",
  
  show_rownames = TRUE,
  show_colnames = TRUE,
  
  angle_col = 90,
  
  fontsize = 16,
  fontsize_row = 16,
  fontsize_col = 14,
  
  treeheight_row = 90,
  treeheight_col = 90,
  
  cellwidth = 32,
  cellheight = 30,
  
  legend = TRUE,
  annotation_legend = TRUE,
  
  annotation_names_col = TRUE,
  
  main = "Choroid Plexus Marker Genes",
  
  filename = file.path(
    plot_dir,
    "Choroid_Plexus_marker_genes_original.pdf"
  ),
  
  width = 13,
  height = 10
)


############################################################
## Save ORIGINAL marker heatmap PNG
############################################################

safe_pheatmap(
  marker_original_scaled,
  
  color = marker_colors,
  breaks = marker_breaks,
  
  scale = "none",
  
  annotation_col = annotation,
  
  annotation_colors = annotation_colors,
  
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  
  clustering_method = "complete",
  
  border_color = "#999999",
  
  show_rownames = TRUE,
  show_colnames = TRUE,
  
  angle_col = 90,
  
  fontsize = 16,
  fontsize_row = 16,
  fontsize_col = 14,
  
  treeheight_row = 90,
  treeheight_col = 90,
  
  cellwidth = 32,
  cellheight = 30,
  
  legend = TRUE,
  annotation_legend = TRUE,
  
  annotation_names_col = TRUE,
  
  main = "Choroid Plexus Marker Genes",
  
  filename = file.path(
    plot_dir,
    "Choroid_Plexus_marker_genes_original2.png"
  ),
  
  width = 15,
  height = 10
)


############################################################
## 21. ONE HEATMAP WITH TWO GENE SECTIONS
## 1. ChP-associated genes
## 2. CC-associated genes
############################################################

##########################################################
## Define gene sets
##########################################################

ChP_markers <- c(
  "TTR",
  "CLIC6",
  "FOLR1",
  "AQP1",
  "KCNE2",
  "CLDN2",
  "OTX2",
  "KRT18",
  "NME5",
  "LUM",
  "DCN",
  "RBP1",
  "PLTP",
  "CA2",
  "SLC13A4",
  "TRPM3",
  "IFI6",
  "SLC39A12",
  "WFIKKN2",
  "C1QTNF5",
  "PRODH",
  "FHIT",
  "MFRP",
  "ETV3L"
)

length(ChP_markers)
# 24

Brain_markers <- c(
  "PMCH",
  "GRM4",
  "GABRA6",
  "SLC6A3",
  "MOG",
  "HAPLN2",
  "BARHL1",
  "OXT",
  "NEUROD6",
  "AVP",
  "FGF3",
  "HCRT"
)

##########################################################
## Combine genes in desired order
## ChP first, Brain second
##########################################################

markers_combined <- c(
  ChP_markers,
  Brain_markers
)

##########################################################
## Keep only genes found in dataset, preserving order
##########################################################

markers_combined_present <- markers_combined[
  markers_combined %in% rownames(vsd)
]

markers_combined_missing <- setdiff(
  markers_combined,
  rownames(vsd)
)

cat("\n========================================\n")
cat("COMBINED ChP + BRAIN HEATMAP\n")
cat("========================================\n")

cat("\nMarkers found:\n")
print(markers_combined_present)

cat("\nMarkers NOT found:\n")
print(markers_combined_missing)

##########################################################
## Extract expression
##########################################################

marker_combined_mat <- assay(vsd)[
  markers_combined_present,
  ,
  drop = FALSE
]

##########################################################
## Row Z-score
##########################################################

marker_combined_scaled <- t(
  scale(
    t(marker_combined_mat)
  )
)

##########################################################
## Remove rows with NA after scaling
##########################################################

marker_combined_scaled <- marker_combined_scaled[
  complete.cases(marker_combined_scaled),
  ,
  drop = FALSE
]

##########################################################
## Optional: cluster rows WITHIN each section only
##########################################################

cluster_within_group <- function(mat) {
  if (nrow(mat) <= 1) {
    return(rownames(mat))
  }
  rownames(mat)[
    hclust(dist(mat), method = "complete")$order
  ]
}

chp_mat <- marker_combined_scaled[
  rownames(marker_combined_scaled) %in% ChP_markers,
  ,
  drop = FALSE
]

brain_mat <- marker_combined_scaled[
  rownames(marker_combined_scaled) %in% Brain_markers,
  ,
  drop = FALSE
]

chp_order <- cluster_within_group(chp_mat)
brain_order <- cluster_within_group(brain_mat)

final_row_order <- c(chp_order, brain_order)

marker_combined_scaled <- marker_combined_scaled[
  final_row_order,
  ,
  drop = FALSE
]

##########################################################
## Row annotation for gene section
##########################################################

row_annotation <- data.frame(
  GeneSet = c(
    rep("ChP-associated", length(chp_order)),
    rep("CC-associated", length(brain_order))
  )
)

rownames(row_annotation) <- rownames(marker_combined_scaled)

##########################################################
## Row gap position
##########################################################

gap_position <- length(chp_order)

##########################################################
## Row annotation colors
##########################################################

annotation_colors2 <- list(
  Group = setNames(unname(group_colors), unname(group_legend_labels)),
  GeneSet = c(
    "ChP-associated" = "#F4A6B7",
    "CC-associated" = "#BDBDBD"
  )
)

##########################################################
## Display final genes used
##########################################################

cat("\nGenes included in combined heatmap:\n")
print(rownames(marker_combined_scaled))

##########################################################
## DISPLAY HEATMAP
##########################################################

open_plot_window(
  12,
  12
)

# Preview omitted: open the corresponding saved PDF or PNG in Plots/.

##########################################################
## Save PDF
##########################################################

safe_pheatmap(
  marker_combined_scaled,
  
  color = marker_colors,
  breaks = marker_breaks,
  
  scale = "none",
  
  annotation_col = annotation,
  annotation_row = row_annotation,
  
  annotation_colors = annotation_colors2,
  
  cluster_rows = FALSE,
  cluster_cols = TRUE,
  
  clustering_method = "complete",
  
  gaps_row = gap_position,
  
  border_color = "#999999",
  
  show_rownames = TRUE,
  show_colnames = TRUE,
  
  angle_col = 90,
  
  fontsize = 16,
  fontsize_row = 15,
  fontsize_col = 14,
  
  treeheight_row = 0,
  treeheight_col = 90,
  
  cellwidth = 32,
  cellheight = 26,
  
  legend = TRUE,
  annotation_legend = TRUE,
  annotation_names_col = TRUE,
  annotation_names_row = TRUE,
  
  main = "ChP-associated marker genes",
  
  filename = file.path(
    plot_dir,
    "ChP_Brain_marker_genes_two_sections.pdf"
  ),
  
  width = 13,
  height = 12
)


##########################################################
## Save PNG - 600 dpi
##########################################################

# Create heatmap without opening a graphics device
hm <- pheatmap::pheatmap(
  marker_combined_scaled,
  
  color = marker_colors,
  breaks = marker_breaks,
  
  scale = "none",
  
  annotation_col = annotation,
  annotation_row = row_annotation,
  
  annotation_colors = annotation_colors2,
  
  cluster_rows = FALSE,
  cluster_cols = TRUE,
  
  clustering_method = "complete",
  
  gaps_row = gap_position,
  
  border_color = "#999999",
  
  show_rownames = TRUE,
  show_colnames = TRUE,
  
  angle_col = 90,
  
  fontsize = 20,
  fontsize_row = 20,
  fontsize_col = 18,
  
  treeheight_row = 0,
  treeheight_col = 90,
  
  cellwidth = 38,
  cellheight = 32,
  
  legend = TRUE,
  annotation_legend = TRUE,
  annotation_names_col = TRUE,
  annotation_names_row = TRUE,
  
  main = "ChP associated marker genes",
  
  silent = TRUE
)

# Save PNG at 600 dpi
png(
  filename = file.path(
    plot_dir,
    "ChP_Brain_marker_genes_two_sections2.png"
  ),
  
  width = 17,
  height = 20,
  units = "in",
  res = 600,
  bg = "white"
)

grid::grid.newpage()
grid::grid.draw(hm$gtable)

dev.off()
