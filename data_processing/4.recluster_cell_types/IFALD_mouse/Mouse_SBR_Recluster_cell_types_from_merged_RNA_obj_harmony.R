# Recluster cell types from a combo object and run chromvar mouse data
###libraries
##################
#library(future)

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

suppressMessages(library(EnsDb.Mmusculus.v79))
suppressMessages(library(BSgenome.Mmusculus.UCSC.mm10))
suppressMessages(library(GenomicRanges))
#suppressMessages(library(future))
suppressMessages(library(optparse))


suppressMessages(library(googlesheets4))
suppressMessages(library(stringr))
suppressMessages(library(harmony))

require(magrittr)
library(ggplot2)
library(RColorBrewer)
library(JASPAR2020)
library(TFBSTools)
library(patchwork)
library(chromVAR)

################################

#####################################
####### FUNCTIONS ##################
####################################
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

DefaultAssay(panc.my) <- 'RNA'
panc.my <- DietSeurat(panc.my, assays = c('RNA'))


cell.types.oi <- c( 'Hepatocytes', 'Cholangiocytes', 'HSCs|fibro',  'LSECs'
                    )

cell.types.in.object <- unique(as.character(unlist(panc.my[[cell_column]])))
cell.types.touse <- intersect(cell.types.oi, cell.types.in.object)
print(cell.types.touse)

print(dim(panc.my))


cell.types.oi %>% walk (function(ct) {
  
    if(!file.exists(paste0(add_filename,"_",make.names(ct), ".rds"))) {
      print(ct)
      int.sub <- subset(x = panc.my, 
                        cells = rownames(dplyr::filter(panc.my@meta.data, 
                                                       grepl(ct, .data[[cell_column]])
                        )
                        )
      )
      print(dim(int.sub))
      
      ######## Normalize RNA
      DefaultAssay(int.sub) <- 'RNA'
      cat('normalizing RNA\n')
      int.sub[["percent.mt"]] <- PercentageFeatureSet(int.sub, pattern = "^mt-") 
      int.sub <- int.sub %>%
        NormalizeData(assay = 'RNA') 
      
      int.sub <- int.sub %>%
        SCTransform(
          assay = 'RNA',
          vars.to.regress = c("nCount_RNA", "percent.mt"),
          return.only.var.genes = TRUE, verbose = F) %>%
        RunPCA(assay = 'SCT', do.print = FALSE, verbose = F) %>%
        RunUMAP(dims = 1:50,reduction.name = 'umap.rna', reduction.key = 'rnaUMAP_', verbose = F) 
    
      int.sub@meta.data$Mouse_ID <- as.character(int.sub@meta.data$Mouse_ID)
      int.sub <- runHarmonyNormalization(int.sub)
      
      cat('saving the object...\n')
      saveRDS(int.sub,  paste0(add_filename,"_",make.names(ct), ".rds"))
      
      
      DimPlot(int.sub, group.by = 'Mouse_ID', reduction='umap.rna', label = TRUE)
      ggsave(glue::glue("Dimplot_Mouse_ID_{ct}.pdf"),
             height=5,width=8,useDingbats=FALSE,limitsize = FALSE)
      
      DimPlot(int.sub, group.by = 'Mouse_ID', reduction='umap.harmony', label = TRUE)
      ggsave(glue::glue("Dimplot_Mouse_ID_{ct}_umap.harmony.pdf"),
             height=5,width=8,useDingbats=FALSE,limitsize = FALSE)
      
      DimPlot(int.sub, group.by = 'seurat_clusters', reduction='umap.harmony', label = TRUE)
      ggsave(glue::glue("Dimplot_seurat_clusters_{ct}_umap.harmony.pdf"),
             height=5,width=8,useDingbats=FALSE,limitsize = FALSE)
      
    } 
    
}
)








