rm(list = ls())
pacman::p_load(Biobase, GEOquery, Seurat)
source('./create_seu_list_functions.R')

## Trabajaremos con datos de scRNA-seq de un estudio de AML
## Emplearemos el paquete Seurat para ello
## Se trata de un estudio de tipo CITE-seq, por lo que tenemos:
  ## datos de expresión (scRNA-seq)
  ## datos de antígenos de superficie (ADT)

#1 Download supplementary files from GEO accession
GEOacc <- 'GSE271406'
GEOquery::getGEOSuppFiles(GEOacc)

#2 Erase files which do not correspond to scRNA-seq experiment 
file.remove('./GSE271406/GSE271406_count.csv.gz')
file.remove('./GSE271406/GSE271406_HD_feature_reference.csv.gz')
file.remove('./GSE271406/GSE271406_PDAML_feature_reference.csv.gz')
file.remove('./GSE271406/GSE271406_RAW.tar')

#3 Create individual directories for each sample
dir.files <- list.files('./GSE271406', full.names = TRUE)
file.names <- basename(dir.files)

sample.ids <- sub(
  pattern = "^GSE271406_(.+?)_(barcodes|features|matrix).*",
  replacement = "\\1",
  x = file.names
)

for (sample in unique(sample.ids)) {
  # Create dir for each sample
  sample.dir <- file.path('./GSE271406', sample)
  dir.create(sample.dir)
  
  sample.files <- dir.files[sample.ids == sample]
  
  # Move files to their corresponding directory
  file.rename(from = sample.files, to = file.path(sample.dir, basename(sample.files)))
}

# Rename files to Seurat standard format
sample.dirs <- list.dirs('./GSE271406', full.names = TRUE, recursive = FALSE)

for (d in sample.dirs) {
  files <- list.files(d, full.names = TRUE)
  
  for (f in files) {
    if (grepl("barcodes", f)) {
      file.rename(f, file.path(d, "barcodes.tsv.gz"))}
    if (grepl("features", f)) {
      file.rename(f, file.path(d, "features.tsv.gz"))}
    if (grepl("matrix", f)) {
      file.rename(f, file.path(d, "matrix.mtx.gz"))}
  }
}

#4 Create a list of Seurat objects for patients
# Read metadata from samples.csv
samples <- read.csv(
  file = './samples.csv',
  header = TRUE,
  quote = "\'"
  
)

# Iterate over the sample sheet and generate a list with Seurat objects
seurat_list <- lapply(1:nrow(samples), function(i) {
  create_seurat(
    path_cite = samples$path_cite[i],
    patient_id = samples$patient_id[i],
    ds = samples$ds[i],
    response = samples$response[i],
    AML_type = samples$AML_type[i],
    age_range = samples$age_range[i],
    sex = samples$sex[i]
    
  )
})

# Save the list for further analysis
save(seurat_list, file = './patients_CITE_seurat.rda')

#5 Get healthy donor samples
list.hd.files <- c('GSM8376300', 
                   'GSM8376301',
                   'GSM8376302',
                   'GSM8376303',
                   'GSM8376304',
                   'GSM8376305')

for (gsm in list.hd.files) {
  GEOquery::getGEOSuppFiles(gsm)
}
    
for (d in paste0('./',list.hd.files)) {
  
  files <- list.files(d, full.names = TRUE)
  
  for (f in files) {
    
    if (grepl("barcodes", f)) {
      file.rename(f, file.path(d, "barcodes.tsv.gz"))
    }
    
    if (grepl("features", f)) {
      file.rename(f, file.path(d, "features.tsv.gz"))
    }
    
    if (grepl("matrix", f)) {
      file.rename(f, file.path(d, "matrix.mtx.gz"))
    }
  }
}

# Create Seurat objects for healthy donors
hd1_seurat <- create_seurat_hd(
  path_gex = './GSM8376300/',
  path_adt = './GSM8376303/',
  patient_id = 'HD1',
  ds = 'healthy_donor'
)

hd2_seurat <- create_seurat_hd(
  path_gex = './GSM8376301/',
  path_adt = './GSM8376304/',
  patient_id = 'HD2',
  ds = 'healthy_donor'
)

hd3_seurat <- create_seurat_hd(
  path_gex = './GSM8376302/',
  path_adt = './GSM8376305/',
  patient_id = 'HD3',
  ds = 'healthy_donor'
)

# Create list for healthy donor datasets
seurat_controls <- list(hd1_seurat, hd2_seurat, hd3_seurat)
save(seurat_controls, file = './healthy_CITE_seurat.rda')

