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
    CellCycleScoring(s.features = cc.genes.updated.2019$s.genes, g2m.features = cc.genes.updated.2019$g2m.genes, set.ident = F) %>%
    SCTransform(
      assay = 'RNA',
      vars.to.regress = c("nCount_RNA", 
                          "percent.mt", 
                          "S.Score", 
                          "G2M.Score"
      ),
      return.only.var.genes = TRUE, verbose = T) %>%
    RunPCA(assay = 'SCT', do.print = FALSE, verbose = T) %>%
    RunUMAP(dims = 1:50,reduction.name = 'umap.rna', reduction.key = 'rnaUMAP_', verbose = T) %>%
    FindNeighbors( dims = 1:50) %>%
    FindClusters(resolution = 2, verbose = FALSE)
  
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


#add clinical info
clinical <- fread('/diskmnt/Projects/SenNet_analysis/Main.analysis/clinical_data/Cohort_full_clinical_v4.csv', data.table = F) 
clinical$Patient_ID <- paste0('SN', sprintf("%03d", as.numeric(clinical$`Participant ID`)), 'H1')
clinical$Gender <- factor(clinical$Gender, levels = c('Woman', 'Man'))

clinical <- clinical %>%
  mutate(
    BMI.cat = case_when(`Body Mass Index` < 25 ~ 'Healthy weight',
                        `Body Mass Index` < 30 ~ 'Overweight',
                        `Body Mass Index` < 35 ~ 'Obese',
                        `Body Mass Index` >= 35 ~ 'Extremely obese'),
    BMI.cat = factor(BMI.cat, levels =c('Healthy weight','Overweight','Obese','Extremely obese')),
    Age.group = case_when(Age < 40 ~ 'Young',
                          Age >= 60 ~ 'Old',
                          TRUE ~ 'Middle age'),
    Age.group =factor(Age.group, levels = c('Young', 'Middle age', 'Old')),
    Alcohol = str_split_fixed(`Alcohol History`, ' ', 2)[,1],
    Smoking = str_split_fixed(`Smoking History`, ' ', 2)[,1],
    Race = str_split_fixed(Race, ' ', 2)[,1],
    History.of.Cancer = `History of Cancer`)

head(clinical)

panc.my@meta.data <- panc.my@meta.data %>% 
  rownames_to_column(var='B') %>% 
  select(-Age, -Age.group) %>%
  left_join(clinical, by='Patient_ID') %>%
  column_to_rownames(var = 'B')


panc.my <- subset(x = panc.my, cells = rownames(dplyr::filter(panc.my@meta.data, !grepl('Doubl', .data[[cell_column]]))))

panc.my <- NormalizeRNA(panc.my)

saveRDS(panc.my,  paste0(add_filename,".rds"))








