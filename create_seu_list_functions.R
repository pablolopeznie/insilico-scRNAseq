library(Seurat)
library(Biobase)

#1 Function to create Seurat objects from Healthy Donors
create_seurat_hd <- function(path_gex, path_adt, patient_id, ds) {
  # Import GEX and ADT data from source directories
  gex <- Read10X(data.dir = path_gex)
  
  # Adjust ADT since it is not in the standard format
  adt <- ReadMtx(
    mtx = paste0(path_adt, sep = '/', 'matrix.mtx.gz'),
    features = paste0(path_adt, sep = '/', 'features.tsv.gz'),
    cells = paste0(path_adt, sep = '/', 'barcodes.tsv.gz'),
    feature.column = 1
  )
  
  # Homogenize cell barcodes between GEX and ADT
  adt_ids <- paste0(colnames(adt), sep = '-', '1')
  colnames(adt) <- adt_ids
  
  # Extract common cell barcodes between GEX and ADT
  common_barcodes <- intersect(colnames(adt), colnames(gex))
  
  # Filter initial objects
  gex_filtered <- gex[, common_barcodes]
  adt_filtered <- adt[, common_barcodes]
  
  # Create Seurat object with Gene Expression Data
  seurat_obj <- CreateSeuratObject(
    counts = gex_filtered,
    project = 'AML_Pre-Tr'
  )
  
  # Add counts from Antibody-derived tags
  seurat_obj[['ADT']] <- CreateAssayObject(counts = adt_filtered)
  
  # Specify Patient ID and Disease State
  seurat_obj$`patient.id` <- patient_id
  seurat_obj$`disease.state` <- ds
  seurat_obj$`response` <- NA
  seurat_obj$`AML.type` <- NA
  seurat_obj$`age.range` <- NA
  seurat_obj$`sex` <- NA
  
  # Notify user about completion
  message('Seurat object for ', patient_id, ' created successfully.')
  
  # Return Seurat object
  return(seurat_obj)
}

#2 Function to create Seurat objects from Patient Data
create_seurat <- function(
    path_cite, patient_id, ds, response, AML_type, age_range, sex) {
  # Import CITE-seq joint data (ADT+GEX)
  cite_list <- Read10X(data.dir = path_cite)
  
  # Create Seurat object with imported GEX Data
  seurat_obj <- CreateSeuratObject(
    counts = cite_list$`Gene Expression`,
    project = 'AML_Pre-Tr'
  )
  
  # Add ADT as an additional layer
  seurat_obj[['ADT']] <- CreateAssayObject(counts = cite_list$`Antibody Capture`)
  
  # Add study covariables to the object
  seurat_obj$`patient.id` <- patient_id
  seurat_obj$`disease.state` <- ds
  seurat_obj$`response` <- response
  seurat_obj$`AML.type` <- AML_type
  seurat_obj$`age.range` <- age_range
  seurat_obj$`sex` <- sex
  
  # Notify user about completion
  message('Seurat object for ', patient_id, ' created successfully.')
  
  # Return object
  return(seurat_obj)
}