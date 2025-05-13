## Alla Karpova
## 02.15.2024 Alla added species filter
### merge objects of one data type from sennet project
# v2.1 05.30.2023 

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
suppressMessages(library(future))
suppressMessages(library(optparse))

suppressMessages(library(googlesheets4))
suppressMessages(library(stringr))
suppressMessages(library(doParallel))
library(harmony)

############## FUNCTIONS #####################
runHarmonyNormalization <- function(obj, dims=30, column = 'Mouse_ID') {
  
  obj <- obj %>%
    RunHarmony(column, reduction = 'pca', assay.use = 'SCT') %>%
    FindNeighbors(reduction = "harmony", dims = 1:dims) %>%
    FindClusters(verbose = FALSE, resolution = 1, algorithm = 4,
                 method='igraph') %>%
    RunUMAP(reduction = "harmony",reduction.name = 'umap.harmony', reduction.key = 'harmonyUMAP_',  dims = 1:dims)
  
  #obj <- NormalizeData(obj, assay = 'RNA')
  
  return(obj)
  
}



option_list = list(
  make_option(c("-i", "--input.folder"),
              type="character",
              default=NULL, 
              help="path to folder with cancer level merged rds objects",
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
              metavar="character")
  
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

# read in initial arguments
input.path <- opt$input.folder
out_path <- opt$output
add_filename <- opt$extra
meta.path <- opt$metadata.file
data.type <- opt$data.type
tissue <- opt$tissue
species <- opt$species

if (data.type == '5sc') {
  dt.tofilter = "5' scRNA"
} else if (data.type == '3sc') {
  dt.tofilter = "3' scRNA"
} else if (data.type == 'combo' | grepl('multi|Multi', data.type)) {
  dt.tofilter = "Combo RNA"
} else if (data.type == 'snRNA' ) {
  dt.tofilter = "snRNA-seq"
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
                      sheet = "mouse_IFALD_data", trim_ws = T)

samples <- samples %>% dplyr::filter(`Include in downstream` == 'Yes')
samples <- samples %>% dplyr::filter(`Data type` == dt.tofilter)

samples <- samples %>% dplyr::select(`Correct.Mouse ID`,`Correct.Sample name`, `Sequencing batch`,Cohort, `Timepoint (weeks)`, Tissue,Sex, `Rds objects`)

samples.id <- samples$`Correct.Sample name` %>% as.character()

cat (paste("Samples found:" ,length(samples.id), '\n'))

cat('creating object \n')
paths <- samples$`Rds objects`
paths

# my.metadata <- fread(meta.path, data.table = F) %>% 
#   data.frame(row.names = 1, check.rows = F, check.names = F) %>%
#   dplyr::select(-seurat_clusters)

if (!file.exists(paste0(length(samples.id),"_Mouse_Merged_not_normalized_",add_filename,".rds"))) {
  
  # make the list of objects
  registerDoParallel(cores=10)
  cat ('Reading in objects\n')
  
  obj <- foreach (s=samples.id, 
                  p = paths, 
                  pid = samples$`Correct.Mouse ID`, 
                  sq = samples$`Sequencing batch`, 
                  coh = samples$Cohort,
                  a = samples$`Timepoint (weeks)`, 
                  tis = samples$Tissue, 
                  sex = samples$Sex,
                  .combine=c) %dopar% {
    print(s)
    obj=readRDS(p) 
    print(paste('opened', s))

    DefaultAssay(obj) <- 'RNA'
    obj <- DietSeurat(obj, assays = 'RNA')

    
    obj@meta.data$Sample_ID = s
    obj@meta.data$Mouse_ID = pid
    obj@meta.data$Sequencing_batch = sq
    obj@meta.data$Data_type = 'snRNA-seq'
    obj@meta.data$Cohort = coh
    obj@meta.data$Timepoint.weeks = a
    obj@meta.data$Tissue = tis
    obj@meta.data$Sex = sex
    return(obj)
  }
  stopImplicitCluster()
}


  cat('doing normalization of reg RNA object\n')
  
  combined <- merge(x = obj[[1]], y = obj[-1], add.cell.ids = samples.id)  
  
  
  cat('normalizing RNA\n')
  DefaultAssay(combined) <- 'RNA'
  
  combined <- combined %>%
    NormalizeData(assay = 'RNA') 
  
  combined[["percent.mt"]] <- PercentageFeatureSet(combined, pattern = "^mt-") 
  print(head(combined@meta.data))
  #print(head(combined$G2M.Score))
  
  combined <- combined %>%
    SCTransform(
      assay = 'RNA',
      vars.to.regress = c("nCount_RNA", 
                          "percent.mt"
                          
      ),
      return.only.var.genes = FALSE, verbose = T) %>%
    RunPCA(assay = 'SCT', do.print = FALSE, verbose = T) %>%
    RunUMAP(dims = 1:30,reduction.name = 'umap.rna', reduction.key = 'rnaUMAP_', verbose = T) %>%
    FindNeighbors( dims = 1:30) %>%
    FindClusters(resolution = 1, verbose = FALSE)
  
  cat('saving the object...\n')
  saveRDS(combined,  paste0(length(samples.id),"_Mouse_Merged_normalized_",add_filename,".rds"))
  
  #combined <- AddMetaData(combined, my.metadata)
  
  DimPlot(combined, reduction = "umap.rna", group.by = "seurat_clusters", label = TRUE, label.size = 2.5, repel = TRUE)
  ggsave(paste0(length(samples.id),"_Mouse_Merged_clusters_", add_filename, ".pdf"),height=10,width=11)
  DimPlot(combined, reduction = "umap.rna", group.by = "Sample_ID", label = TRUE, label.size = 2.5, repel = TRUE)
  ggsave(paste0(length(samples.id),"_Mouse_Merged_Sample_ID_", add_filename, ".pdf"),height=10,width=12)
   
  DimPlot(combined, reduction = "umap.rna", group.by = "Timepoint.weeks", label = TRUE, label.size = 2.5, repel = TRUE)
  ggsave(paste0(length(samples.id),"_Mouse_Merged_Timepoint.weeks_", add_filename, ".pdf"),height=10,width=12)
  
  
  combined <- combined %>% runHarmonyNormalization()
  
  saveRDS(combined,  paste0(length(samples.id),"_Mouse_Merged_normalized_",add_filename,".rds"))
  
  DimPlot(combined, reduction = "umap.harmony", group.by = "Sample_ID", label = TRUE, label.size = 2.5, repel = TRUE)
  ggsave(paste0(length(samples.id),"_Mouse_Merged_harmony_Sample_ID_", add_filename, ".pdf"),height=10,width=12)
  
  DimPlot(combined, reduction = "umap.harmony", group.by = "Timepoint.weeks", label = TRUE, label.size = 2.5, repel = TRUE)
  ggsave(paste0(length(samples.id),"_Mouse_Merged_harmony_Timepoint.weeks_", add_filename, ".pdf"),height=10,width=12)
  
  
  fwrite(cbind(Embeddings(combined, reduction = "umap.rna"),
               combined@meta.data), 
         paste0(length(samples.id),"_Mouse_Merged_normalized_",add_filename,".metadata.tsv"), sep = '\t', row.names = T)
  


