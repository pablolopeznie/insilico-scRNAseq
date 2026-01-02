if (!requireNamespace("pacman", quietly = TRUE)) {
  install.packages("pacman")
}

# Loading required packages to the session
pacman::p_load(Seurat, SeuratObject, Biobase, ggplot2, patchwork, glmGamPoi,
               celldex, SingleR, SingleCellExperiment, dplyr)

# Load healthy donor and patient datasets
load('./healthy_CITE_seurat_preprocessed.rda')
load('./patients_CITE_seu_filt.rda')

full_list <- c(filtered_seurat_controls, seurat_list)

# Normalize each individual dataset with SCTransform (vst:glmGamPoi)
hd.pat.integrated <- lapply(full_list, Seurat::SCTransform) # method = "glmGamPoi") -> not used in server

# Find common features to perform integration
integ.feat <- SelectIntegrationFeatures(object.list = hd.pat.integrated, nfeatures = 2000) # reduced from 3000

# Prepare the object SCT-normalized for integration
hd.pat.integrated <- PrepSCTIntegration(object.list = hd.pat.integrated, anchor.features = integ.feat)

# Run PCA to summarize the biology varibility in linear components
hd.pat.integrated <- lapply(hd.pat.integrated, Seurat::RunPCA, features = integ.feat)

# Find most variable features to perform integration
hd.pat.anchors <- FindIntegrationAnchors(
  object.list = hd.pat.integrated, 
  normalization.method = 'SCT',
  anchor.features = integ.feat, 
  reduction = 'rpca',
  k.anchor = 20) # reduced from 50

# Perform integration of all the datasets into a single one
hd.pat.integrated <- IntegrateData(anchorset = hd.pat.anchors, normalization.method = 'SCT')

# Visualization of integrated datasets
hd.pat.integrated <- RunPCA(hd.pat.integrated, assay = 'integrated')

# Run UMAP to visualize data according to non-linear relationship
hd.pat.integrated <- RunUMAP(hd.pat.integrated, reduction = "pca", dims = 1:15, assay = 'integrated')

save(hd.pat.integrated, file = './hd_pat_integrated.rda')
