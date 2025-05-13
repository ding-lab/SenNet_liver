## Alla Karpova
### merge Xenium objects from sennet project
# v1.2 05.23.2024  - added ability to merge all panels together

suppressMessages(library(Signac))
suppressMessages(library(Seurat))
suppressMessages(library(GenomeInfoDb))
suppressMessages(library(ggplot2))
suppressMessages(library(RColorBrewer))

suppressMessages(library(tidyverse))
set.seed(1234)

suppressMessages(library(dplyr))
suppressMessages(library(data.table))

suppressMessages(library(EnsDb.Hsapiens.v86))
suppressMessages(library(GenomicRanges))
suppressMessages(library(future))
suppressMessages(library(optparse))
suppressMessages(library(BSgenome.Hsapiens.UCSC.hg38))

suppressMessages(library(googlesheets4))
suppressMessages(library(stringr))
suppressMessages(library(harmony))
library(doParallel)
options(future.globals.maxSize= 1610612736) #1.5Gb

############## FUNCTIONS #####################

option_list = list(
  make_option(c("-o", "--output"),
              type="character",
              default="./", 
              help="output folder path",
              metavar="character"),
  make_option(c("-e", "--extra"),
              type="character",
              default="./", 
              help="add unique string identifier for your data",
              metavar="character"),
  make_option(c("-m","--metadata.file"),
              type="character",
              default=NULL,
              help = "path to metadats file with cell types, make cell barcodes in the 1st column",
              metavar="character"),
  make_option(c("-d","--data.type"),
              type="character",
              default='combo',
              help = "data.type to use for merging, use 5sc for 5'scRNA-seq, 3sc for 3'scRNA seq, combo for multiome",
              metavar="character"),
  make_option(c("-t","--tissue"),
              type="character",
              default='liver',
              help = "liver, skin, bone",
              metavar="character"),
  make_option(c("-s","--species"),
              type="character",
              default='Human',
              help = "Human or Mouse",
              metavar="character"),
  make_option(c("-v","--panel.version"),
              type="character",
              default=NULL,
              help = "8UT8B9_hMulti or 8P9U3X_hMulti",
              metavar="character")
  
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

# read in initial arguments

out_path <- opt$output
add_filename <- opt$extra
meta.path <- opt$metadata.file
data.type <- opt$data.type
tissue <- opt$tissue
species <- opt$species
ver <- opt$panel.version

if (data.type == "Xenium" | grepl('enium', data.type)) {
  dt.tofilter = "Xenium"
}
cat(paste('Data type to use', dt.tofilter))

if (tissue == 'bone') {
  tis.tofilter = "bone marrow"
} else {
  tis.tofilter = tissue
}
cat(paste('Tissue to use', tis.tofilter))

select <- dplyr::select

dir.create(out_path, showWarnings = F, recursive = T)
setwd(out_path)


###### load google sheet and extract samples from there ########
gs4_deauth()
samples <- read_sheet("https://docs.google.com/spreadsheets/d/1VeWme__vvVHAhHaQB3wCvAGuq-w3WrhZ5cT-7Mh7Sr0/edit#gid=0", 
                      sheet = "mCRC data", trim_ws = T)

samples <- samples %>% dplyr::filter(`Include in downstream` == 'Yes')
samples <- samples %>% dplyr::filter(`Data type` %in% dt.tofilter)
samples <- samples %>% dplyr::filter(Tissue == tis.tofilter)
samples <- samples %>% dplyr::filter(Species == species)
if(!is.null(ver)) {
  samples <- samples %>% dplyr::filter(`Xenium panel version` == ver)
}

samples <- samples %>% dplyr::select(`Patient ID`, Piece_ID, `Sample name`, `Data type`, Age, Treatment, Tissue, `Rds objects`, `Xenium panel version`)

samples.id <- samples$`Sample name` %>% as.character()

cat (paste("Samples found:" ,length(samples.id), '\n'))


# my.metadata <- fread(meta.path, data.table = F) %>% 
#   data.frame(row.names = 1, check.rows = F, check.names = F) %>%
#   dplyr::select(-seurat_clusters)

if (!file.exists(paste0(length(samples.id),"_Merged_not_normalized_",add_filename,".rds"))) {
  cat('creating object \n')
  paths <- samples$`Rds objects`
  paths
  
  # make the list of objects
  registerDoParallel(cores=10)
  cat ('Reading in objects\n')
  
  obj <- foreach (s=samples.id, 
                  p = paths, 
                  pid = samples$`Patient ID`, 
                  piece = samples$Piece_ID, 
                  dt = samples$`Data type`, 
                  a = samples$Age, 
                  tr = samples$Treatment,
                  tis = samples$Tissue, 
                  v=samples$`Xenium panel version`, .combine=c) %dopar% {
                    print(s)
                    obj=readRDS(p) 
                    print(paste('opened', s))
                    if (dt.tofilter == 'Combo RNA' & length(dt.tofilter) == 1) {
                      DefaultAssay(obj) <- assay.towork
                      if (!file.exists(Fragments(obj)[[1]]@path)) stop("Urgh, this sample object can't locate fragments file")
                      obj <- DietSeurat(obj, assays = c(assay.towork, 'RNA'))
                    } else if (dt.tofilter=='Xenium') {
                      DefaultAssay(obj) <- 'Xenium'
                      obj <- DietSeurat(obj, assays = 'Xenium')
                    } else {
                      DefaultAssay(obj) <- 'RNA'
                      obj <- DietSeurat(obj, assays = 'RNA')
                    }
                    
                    obj$Sample = s
                    obj$Patient_ID = pid
                    obj$Piece_ID = piece
                    obj$Data_type = dt
                    obj$Age = as.numeric(a)
                    obj$Treatment = tr
                    obj$Tissue = tis
                    obj$Panel_version = v
                    return(obj)
                  }
  stopImplicitCluster()
  
  combined <- merge(x = obj[[1]], y = obj[-1], add.cell.ids = samples.id)  
  saveRDS(combined,  paste0(length(samples.id),"_Merged_not_normalized_",add_filename,".rds"))
  
} else {
  combined <- readRDS(paste0(length(samples.id),"_Merged_not_normalized_",add_filename,".rds"))
}

cat('normalizing Xenium\n')
DefaultAssay(combined) <- 'Xenium'


combined <- combined %>%
  SCTransform(
    assay = 'Xenium',
    return.only.var.genes = TRUE, 
    verbose = T) %>%
  RunPCA(assay = 'SCT', do.print = FALSE, verbose = T) %>%
  RunUMAP(dims = 1:30, verbose = T) %>%
  RunHarmony('Patient_ID', reduction = 'pca', assay.use = 'SCT') %>%
  FindNeighbors(reduction = "harmony", dims = 1:30) %>%
  FindClusters(resolution = 1.5, verbose = FALSE) %>%
  RunUMAP(reduction = "harmony",reduction.name = 'umap.harmony', reduction.key = 'harmonyUMAP_',  dims = 1:30)


cat('saving the object...\n')
saveRDS(combined,  paste0(length(samples.id),"_Merged_normalized_",add_filename,".rds"))

#combined <- AddMetaData(combined, my.metadata)
combined$Age.group <- case_when(combined$Age >= 60 ~ 'Old', 
                                combined$Age >= 40 ~ 'Middle age', 
                                TRUE ~ 'Young')

DimPlot(combined,  group.by = "seurat_clusters", reduction = 'umap.harmony', label = TRUE, label.size = 2.5, repel = TRUE)
ggsave(paste0(length(samples.id),"_Merged_clusters_", add_filename, ".pdf"),height=10,width=11)
DimPlot(combined,  group.by = "Sample", reduction = 'umap.harmony', label = TRUE, label.size = 2.5, repel = TRUE)
ggsave(paste0(length(samples.id),"_Merged_Sample_", add_filename, ".pdf"),height=10,width=12)
DimPlot(combined,  group.by = "Age.group", reduction = 'umap.harmony', label = TRUE, label.size = 2.5, repel = TRUE)
ggsave(paste0(length(samples.id),"_Merged_Age.group_", add_filename, ".pdf"),height=10,width=11)


#DimPlot(combined, reduction = "umap.rna", group.by = "cell_type", label = TRUE, label.size = 2.5, repel = TRUE)
#ggsave(paste0(length(samples.id),"_Merged_cell_type_", add_filename, ".pdf"),height=10,width=11)

fwrite(combined@meta.data, 
       paste0(length(samples.id),"_Merged_normalized_",add_filename,".metadata.tsv"), sep = '\t', row.names = T)

cat('saving the object...\n')
saveRDS(combined,  paste0(length(samples.id),"_Merged_normalized_",add_filename,".rds"))



