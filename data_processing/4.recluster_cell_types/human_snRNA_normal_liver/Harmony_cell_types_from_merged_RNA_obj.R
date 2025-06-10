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

runHarmonyNormalization <- function(obj, dims=30, column = 'Patient_ID') {
  
  # obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^MT-")
  # s.genes <- cc.genes$s.genes
  # g2m.genes <- cc.genes$g2m.genes
  # 
  # obj <- CellCycleScoring(obj, s.features = s.genes, g2m.features = g2m.genes, set.ident = F)
  # 
  # obj <- obj %>%
  #   SCTransform(
  #     assay = 'RNA',
  #     vars.to.regress =  c("nCount_RNA", "percent.mt", "S.Score", "G2M.Score"),
  #     conserve.memory = T,
  #     return.only.var.genes = T,
  #     verbose = FALSE
  #   )
  obj <- obj %>%
    #RunPCA(assay = 'SCT', do.print = FALSE) %>%
    RunHarmony(column, reduction = 'pca', assay.use = 'SCT') %>%
    FindNeighbors(reduction = "harmony", dims = 1:dims) %>%
    FindClusters(verbose = FALSE, resolution = 1.2, algorithm = 4,
                 method='igraph') %>%
    RunUMAP(reduction = "harmony",reduction.name = 'umap.harmony', reduction.key = 'harmonyUMAP_',  dims = 1:dims)
  
  #obj <- NormalizeData(obj, assay = 'RNA')
  
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


int <- runHarmonyNormalization(panc.my)

cat('saving the object...\n')
saveRDS(int,  paste0( add_filename,"_harmony.rds"))

DimPlot(int, reduction='umap.harmony', group.by = 'Patient_ID')
ggsave(glue::glue("Dimplot_{add_filename}_Patient.id.pdf"), width = 8, height = 5)

DimPlot(int, reduction='umap.harmony', group.by = 'seurat_clusters')
ggsave(glue::glue("Dimplot_{add_filename}_seurat_clusters.pdf"), width = 6.5, height = 5)

