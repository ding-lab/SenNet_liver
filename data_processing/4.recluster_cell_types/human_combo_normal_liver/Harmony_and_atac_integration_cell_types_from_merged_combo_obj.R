# Recluster cell types 
###libraries
##################
library(future)

#plan("multicore", workers = 20)
#options(future.globals.maxSize = 100 * 1024 ^ 3)

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
suppressMessages(library(doParallel))
suppressMessages(library(harmony))


normalize_rna_harmony <- function(obj, dims=30, column = 'Patient_ID') {
  DefaultAssay(obj) <- "SCT"
  obj <- obj %>%
    RunPCA(assay = 'SCT', do.print = FALSE) %>% 
    RunHarmony(column, reduction = 'pca', assay.use = 'SCT') %>%
    RunUMAP(reduction = "harmony",reduction.name = 'umap.harmony', reduction.key = 'harmonyUMAP_',  dims = 1:dims)
  
  return(obj)
  
}


integrate_atac <- function (int.sub.f, column = 'Patient_ID') {

  DefaultAssay(int.sub.f) <- 'ATAC_merged'
  
  int.sub.f@meta.data$Batches <- int.sub.f@meta.data[[column]]
  batch.counts <- table(int.sub.f@meta.data$Batches)
  
  bad.count <- names(batch.counts)[batch.counts < 100]
  best.count <- names(batch.counts)[which.max(batch.counts)]
  
  int.sub.f@meta.data$Batches <- case_when( int.sub.f@meta.data$Batches %in% bad.count ~ best.count,
                                              TRUE ~ int.sub.f@meta.data$Batches)
  
  atac.split <- SplitObject(int.sub.f, split.by = 'Batches')
  
  atac.split <- map(atac.split, function(obj) {
    obj <- FindTopFeatures(obj, min.cutoff = 50) %>%
      RunTFIDF() %>%
      RunSVD(reduction.key = 'LSI_',
             reduction.name = 'lsi',
             irlba.work = 400)
    return(obj)
  })
  
  #######integration############
  plan("multiprocess", workers = 10)
  options(future.globals.maxSize = 100 * 1024^3)
  
  cat('FindIntegrationAnchors\n')
  integration.anchors <- FindIntegrationAnchors(
    object.list = atac.split,
    anchor.features = rownames(int.sub.f),
    reduction = "rlsi",
    k.filter=200,
    dims = 2:50
  )
  
  # integrate LSI embeddings
  cat('IntegrateEmbeddings\n')
  integrated <- IntegrateEmbeddings(
    anchorset = integration.anchors,
    reductions = int.sub.f[["lsi"]],
    new.reduction.name = "integrated_lsi",
    dims.to.integrate = 1:50, 
    k.weight = 30
  )
  
  # create a new UMAP using the integrated embeddings
  integrated <- RunUMAP(integrated, 
                        reduction = "integrated_lsi", 
                        dims = 2:50, 
                        reduction.name = "atac.umap", reduction.key = "atacUMAP_")
  
  
  #Annotation(integrated) <- annotations.f
  return(integrated)
}

normalize_multiome_with_integration_harmony <- function(obj,dims = 50, harm.column='Patient_ID', atac.int.column = 'Patient_ID') {
  if (!file.exists('Intermediate.rds')) {
    obj <- normalize_rna_harmony(obj, column = harm.column)
    obj <- integrate_atac(int.sub.f = obj, column = atac.int.column)
    saveRDS(obj, 'Intermediate.rds')
  } else {
    obj <- readRDS('Intermediate.rds')
  }
  DefaultAssay(obj) <- 'SCT'
  obj <- FindMultiModalNeighbors(obj, 
                                 reduction.list = list("harmony", "integrated_lsi"), 
                                 dims.list = list(1:dims, 2:dims))
  obj <- RunUMAP(obj, nn.name = "weighted.nn", 
                 reduction.name = "wnn.umap", 
                 reduction.key = "wnnUMAP_")
  obj <- FindClusters(obj, graph.name = "wsnn", algorithm = 4,  resolution=1.2, verbose = T)
  
  return(obj)
}


###options###
######################
option_list = list(
  make_option(c("-i", "--input.object"),
              type="character",
              default=NULL, 
              help="path to merged object",
              metavar="character"),
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
              metavar="character")
  
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser)
###################################


# read in initial arguments
input.path <- opt$input.object
out_path <- opt$output
add_filename <- opt$extra
meta.path <- opt$metadata.file


dir.create(out_path, showWarnings = F)
setwd(out_path)

filter <- dplyr::filter
select <- dplyr::select

panc.my <- readRDS(input.path)


if(!is.null(meta.path)) {
  my.metadata <- fread(meta.path, data.table = F, header = TRUE) %>% 
    data.frame(row.names = 1, check.rows = F, check.names = F)
  panc.my <- AddMetaData(panc.my, my.metadata)
}


int <- normalize_multiome_with_integration_harmony(panc.my, dims = 30)

cat('saving the object...\n')
saveRDS(int,  paste0( add_filename,"_harmony_atac_integrated.rds"))




