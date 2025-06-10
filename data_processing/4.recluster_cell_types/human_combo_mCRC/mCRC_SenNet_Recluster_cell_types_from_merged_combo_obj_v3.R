# Recluster cell types from a combo object and run chromvar
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

suppressMessages(library(EnsDb.Hsapiens.v100))
suppressMessages(library(GenomicRanges))
#suppressMessages(library(future))
suppressMessages(library(optparse))
suppressMessages(library(BSgenome.Hsapiens.UCSC.hg38))

suppressMessages(library(googlesheets4))
suppressMessages(library(stringr))
#suppressMessages(library(doParallel))

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
runAllChromvar <- function(obj, assay = 'ATAC_merged') {
  DefaultAssay(obj) <- assay
  # Get a list of motif position frequency matrices from the JASPAR database
  pfm <- getMatrixSet(
    x = JASPAR2020,
    opts = list(species = 9606, all_versions = FALSE)
  )
  
  # Scan the DNA sequence of each peak for the presence of each motif
  motif.matrix <- CreateMotifMatrix(
    features = granges(obj),
    pwm = pfm,
    genome = 'BSgenome.Hsapiens.UCSC.hg38',
    use.counts = FALSE
  )
  
  # Create a new Motif object to store the results
  motif <- CreateMotifObject(
    data = motif.matrix,
    pwm = pfm
  )
  
  # Add the Motif object to the assay
  obj <- SetAssayData(
    object = obj,
    assay = assay,
    slot = 'motifs',
    new.data = motif
  )
  
  cat('doing chromvar\n')
  obj <- RegionStats(object = obj, genome = BSgenome.Hsapiens.UCSC.hg38)
  
  obj <- RunChromVAR(
    object = obj,
    genome = BSgenome.Hsapiens.UCSC.hg38
  )
  
  DefaultAssay(obj) <- 'chromvar'
  obj@assays$chromvar@scale.data <- obj@assays$chromvar@data
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
              metavar="character"),
  make_option(c("-a", "--assay"),
              type="character",
              default="ATAC_merged", 
              help="which assay should be used to merge objects? ATAC_merged, peaks",
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
assay.towork <- opt$assay

dir.create(out_path, showWarnings = F)
setwd(out_path)

select <- dplyr::select
filter <- dplyr::filter

my.metadata <- fread(meta.path, data.table = F, header = TRUE) %>% 
  data.frame(row.names = 1, check.rows = F, check.names = F)

panc.my <- readRDS(input.path)
panc.my <- AddMetaData(panc.my, my.metadata)

DefaultAssay(panc.my) <- 'RNA'
panc.my <- DietSeurat(panc.my, assays = c('RNA', assay.towork))



annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v100,  standard.chromosomes = TRUE)
genome(annotations) <- "NA"
seqlevelsStyle(annotations) <- 'UCSC' # instead of USCS because it keeps return error https://github.com/stuart-lab/signac/issues/826
genome(annotations) <- "hg38"

cell.types.oi <- c('Portal endo|Vascular endo|LSECs', 'Hepatocytes', 'Cholangiocytes', "fibroblasts|stellate", 
                   'macs')

cell.types.in.object <- unique(as.character(unlist(panc.my[[cell_column]])))
cell.types.touse <- intersect(cell.types.oi, cell.types.in.object)
print(cell.types.touse)

print(dim(panc.my))

DefaultAssay(panc.my)=assay.towork
Annotation(panc.my) <- annotations


cell.types.oi %>% walk (function(ct) {
  if(!file.exists(paste0(add_filename,"_",make.names(ct), ".chromvar.rds"))) {
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
      int.sub[["percent.mt"]] <- PercentageFeatureSet(int.sub, pattern = "^MT-") 
      int.sub <- int.sub %>%
        NormalizeData(assay = 'RNA') %>%
        CellCycleScoring(s.features = cc.genes.updated.2019$s.genes, g2m.features = cc.genes.updated.2019$g2m.genes, set.ident = F)
      
      int.sub <- int.sub %>%
        SCTransform(
          assay = 'RNA',
          vars.to.regress = c("nCount_RNA", "percent.mt", "S.Score", "G2M.Score"),
          return.only.var.genes = TRUE, verbose = F) %>%
        RunPCA(assay = 'SCT', do.print = FALSE, verbose = F) %>%
        RunUMAP(dims = 1:50,reduction.name = 'umap.rna', reduction.key = 'rnaUMAP_', verbose = F) 
      
      ######## Normalize ATAC
      DefaultAssay(int.sub) <- "ATAC_merged"
      
      Annotation(int.sub) <- annotations
      
      cat('normalizing ATAC\n')
      int.sub <- int.sub %>% 
        RunTFIDF() %>%
        FindTopFeatures(min.cutoff = 'q10') %>%
        RunSVD(verbose = F) %>%
        RunUMAP(reduction = 'lsi', 
                dims = 2:50, reduction.name = "umap.atac", reduction.key = "atacUMAP_", verbose = F)
      
      # do wnn analysis
      cat('doing WNN\n')
      int.sub <- FindMultiModalNeighbors(int.sub, 
                                         reduction.list = list("pca", "lsi"), 
                                         dims.list = list(1:50, 2:50)) %>%
        RunUMAP( nn.name = "weighted.nn", 
                 reduction.name = "wnn.umap", 
                 reduction.key = "wnnUMAP_")  %>%
        FindClusters(graph.name = "wsnn", algorithm = 1, resolution = 2, verbose = T)
      
      cat('saving the object...\n')
      saveRDS(int.sub,  paste0(add_filename,"_",make.names(ct), ".rds"))
      
      DimPlot(int.sub, group.by = cell_column, reduction='wnn.umap', label = TRUE)
      ggsave(glue::glue("Dimplot_{cell_column}_{ct}.pdf"),
             height=5,width=7,useDingbats=FALSE,limitsize = FALSE)
      
      DimPlot(int.sub, group.by = 'Patient_ID', reduction='wnn.umap', label = TRUE)
      ggsave(glue::glue("Dimplot_Patient_ID_{ct}.pdf"),
             height=5,width=8,useDingbats=FALSE,limitsize = FALSE)
      
      DimPlot(int.sub, group.by = 'seurat_clusters', reduction='wnn.umap', label = TRUE)
      ggsave(glue::glue("Dimplot_seurat_clusters_{ct}.pdf"),
             height=5,width=8,useDingbats=FALSE,limitsize = FALSE)
      
    } else {
      int.sub <- readRDS(paste0(add_filename,"_",make.names(ct), ".rds"))
      
    }
    
    #run chromvar
    int.sub <- runAllChromvar(int.sub, assay = assay.towork)
    
    cat('saving the object...\n')
    saveRDS(int.sub,  paste0(add_filename,"_",make.names(ct), ".chromvar.rds"))
    
  }
}
)








