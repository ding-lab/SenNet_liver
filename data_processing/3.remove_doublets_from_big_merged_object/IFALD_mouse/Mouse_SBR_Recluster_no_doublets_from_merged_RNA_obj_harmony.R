# Remove doublets in mouse samples
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



################################

#####################################
####### FUNCTIONS ##################
####################################



############################################

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
              metavar="character"),
  make_option(c("-c","--cell_type_column"),
              type="character",
              default='cell_type',
              help = "column in the metadata with most recent cell types",
              metavar="character")
  
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser)
###################################

NormalizeRNA <- function(obj){
  ######## Normalize RNA
  DefaultAssay(obj) <- 'RNA'
  cat('normalizing RNA\n')
  obj <- obj %>%
    NormalizeData(assay = 'RNA') %>%
    SCTransform(
      assay = 'RNA',
      vars.to.regress = c("nCount_RNA", 
                          "percent.mt"
      ),
      return.only.var.genes = TRUE, verbose = T) %>%
    RunPCA(assay = 'SCT', do.print = FALSE, verbose = T) %>%
    RunUMAP(dims = 1:50,reduction.name = 'umap.rna', reduction.key = 'rnaUMAP_', verbose = T) 
  
  return(obj)
}

runHarmonyNormalization <- function(obj, dims=30, column = 'Mouse_ID') {
  DefaultAssay(obj) <- 'SCT'
  obj <- obj %>%
    RunHarmony(column, reduction = 'pca', assay.use = 'SCT') %>%
    FindNeighbors(reduction = "harmony", dims = 1:dims) %>%
    FindClusters(verbose = FALSE, resolution = 1, 
                 algorithm = 4,
                 method='igraph') %>%
    RunUMAP(reduction = "harmony",reduction.name = 'umap.harmony', reduction.key = 'harmonyUMAP_',  dims = 1:dims)
  
  return(obj)
  
}


# read in initial arguments
input.path <- opt$input.object
out_path <- opt$output
add_filename <- opt$extra
meta.path <- opt$metadata.file
cell_column <- opt$cell_type_column

dir.create(out_path, showWarnings = F)
setwd(out_path)

select <- dplyr::select
filter <- dplyr::filter
my.metadata <- fread(meta.path, data.table = F, header = TRUE) %>% 
  data.frame(row.names = 1, check.rows = F, check.names = F) 

panc.my <- readRDS(input.path)
panc.my <- AddMetaData(panc.my, my.metadata)


panc.my <- subset(x = panc.my, cells = rownames(dplyr::filter(panc.my@meta.data, !grepl('Doubl', .data[[cell_column]]))))

panc.my <- NormalizeRNA(panc.my)
panc.my@meta.data$Mouse_ID <- as.character(panc.my@meta.data$Mouse_ID)
panc.my <- runHarmonyNormalization(panc.my)

saveRDS(panc.my,  paste0(add_filename,".rds"))

DimPlot(panc.my, group.by = 'Mouse_ID', reduction='umap.harmony', label = TRUE)
ggsave(glue::glue("Dimplot_Mouse_ID_{add_filename}_umap.harmony.pdf"),
       height=5,width=8,useDingbats=FALSE,limitsize = FALSE)

DimPlot(panc.my, group.by = 'seurat_clusters', reduction='umap.harmony', label = TRUE)
ggsave(glue::glue("Dimplot_seurat_clusters_{add_filename}_umap.harmony.pdf"),
       height=5,width=8,useDingbats=FALSE,limitsize = FALSE)






