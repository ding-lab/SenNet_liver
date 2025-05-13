#!/usr/bin/env Rscript --vanilla
#v1.3 Alla removed scrublet cutoff threshold option
# v1.2 Alla added adding scrublet annotation and not removing doublets
# v1.1 Alla added regression of cell cycle scores
# load required libraries
library(optparse)
library(Seurat)
library(tidyverse)
library(Matrix)
library(RColorBrewer)
library(ggplot2)
library(data.table)

# create user options
option_list = list(
  make_option(c("-i", "--input"),
              type="character",
              default=NULL,
              help="path to data folder (e.g. cellranger output's raw matrices folder)",
              metavar="character"),
  make_option(c("--pre_filter"),
              type="integer",
              default=300,
              help="min number of reads per cell to prefilter",
              metavar="integer"),
  make_option(c("--nfeature_min"),
              type="integer",
              default=200,
              help="nFeature_RNA min value for filtering",
              metavar="integer"),
  make_option(c("--nfeature_max"),
              type="integer",
              default=10000,
              help="nFeature_RNA max value for filtering",
              metavar="integer"),
  make_option(c("--ncount_min"),
              type="integer",
              default=1000,
              help="nCount_RNA min value for filtering",
              metavar="integer"),
  make_option(c("--ncount_max"),
              type="integer",
              default=80000,
              help="nCount_RNA max value for filtering",
              metavar="integer"),
  make_option(c("--mito_max"),
              type="double",
              default=0.1,
              help="maximum allowed mitochondrial fraction",
              metavar="double"),
  make_option(c("-o", "--output"),
              type="character",
              default="./",
              help="output folder path",
              metavar="character"),
  make_option(c("-s","--sample_id"),
              type="character",
              default="single_cell_study",
              help="Name of your sample",
              metavar="character"),
  make_option(c("--pc_num"),
              type="integer",
              default=30,
              help="number of principal components to use",
              metavar="integer"),
  make_option(c("--scrublet"),
              type="character",
              default=NULL,
              help="path to scrublet folder output",
              metavar="character") #,
 # make_option(c("--scrub.cutoff"),
 #             type="double",
 #             default=0.34,
 #             help="cutoff for doublet calls",
 #             metavar="double")
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

# complain if there's no data
if (is.null(opt$input)){
  print_help(opt_parser)
  stop("Path to data is required (--input).n", call.=FALSE)
}

# read in initial arguments
sample_id <- opt$sample_id
out_path <- opt$output
matrix_dir = opt$input


# make output dir if it doesn't exist
out_path=paste(out_path, '/', sample_id,"/",sep="")
dir.create(out_path, showWarnings = F, recursive = T)
setwd(out_path)

# get direct paths to data
barcode.path <- paste0(matrix_dir, "/barcodes.tsv.gz")
features.path <- paste0(matrix_dir, "/features.tsv.gz")
matrix.path <- paste0(matrix_dir, "/matrix.mtx.gz")

# read in matrix
cat('read in matrix\n')
input <- readMM(file = matrix.path)
feature.names = read.delim(features.path, header = FALSE,stringsAsFactors = FALSE)
barcode.names = read.delim(barcode.path, header = FALSE,stringsAsFactors = FALSE)
colnames(input) = barcode.names$V1
rownames(input) = feature.names$V2

# pre-filter and create initial seurat object
cat('prefilter\n')
bc_sums <- colSums(input)
bg_id <- names(bc_sums[bc_sums >= opt$pre_filter])
panc = CreateSeuratObject(counts = input[,bg_id],project=opt$sample_id,min.cells = 0)

### QC
# get percent mitochondrial content
cat('get percent mitochondrial content\n')
panc[["percent.mt"]] <- PercentageFeatureSet(panc, pattern = "^MT-") / 100


if (!is.null(opt$scrublet)){
  cat('add scrublet\n')
  scrub.path <- list.files(path = opt$scrublet, pattern = 'output_table.csv', full.names = T)
  print(scrub.path)
  scrub.table <- fread(scrub.path, header = T) %>% data.frame(row.names = 1)
 # scrub.table$predicted_doublet <- ifelse(scrub.table$doublet_score>= opt$scrub.cutoff, 'True', 'False')
  panc <- AddMetaData(panc, scrub.table)
}

# plot pre-filter metadata
#panc$percent.mito<-percent.mito
pdf(paste("QC_in_sample_",sample_id, ".pdf", sep=""), width=15, height=9)
VlnPlot(object = panc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
dev.off()

# plot metadata associations
pdf(paste0("FeatureScatter_in_sample_",sample_id,".pdf",sep=""),width=12,height=7)
plot1 <- FeatureScatter(object = panc, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(object = panc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))
dev.off()

# filter step
panc<-subset(x = panc, subset = nFeature_RNA > opt$nfeature_min & 
               nFeature_RNA < opt$nfeature_max & 
               nCount_RNA > opt$ncount_min & 
               nCount_RNA < opt$ncount_max & 
               percent.mt<opt$mito_max)
panc %>% dim

# plot post-filter metadata
pdf(paste("After_QC_in_sample_",sample_id, ".pdf", sep=""), width=15, height=9)
VlnPlot(object = panc, features = c("nCount_RNA", "nFeature_RNA", "percent.mt"), ncol = 3)
dev.off()

# Run the standard workflow for visualization and clustering

s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
panc <- NormalizeData(panc, assay = 'RNA')
panc <- CellCycleScoring(panc, s.features = s.genes, g2m.features = g2m.genes, set.ident = F)
print(head(panc@meta.data))
panc <- SCTransform(panc, 
                    vars.to.regress = c("nCount_RNA","percent.mt","S.Score", "G2M.Score"),return.only.var.genes = F)
panc <- RunPCA(panc, npcs = opt$pc_num, verbose = FALSE)

# t-SNE and Clustering
panc <- RunUMAP(panc, reduction = "pca", dims = 1:opt$pc_num)
panc <- FindNeighbors(panc, reduction = "pca", dims = 1:opt$pc_num)
panc <- FindClusters(panc, resolution = 0.5)

# plot the clusters
pdf(paste0("DimPlot_",sample_id,".pdf"),useDingbats=FALSE)
DimPlot(object = panc, reduction = "umap",label=TRUE,label.size=6)
dev.off()

# plot scrublet
if (!is.null(opt$scrublet)){
  pdf(paste0("DimPlot_predicted_doublet_",sample_id,".pdf"),useDingbats=FALSE)
  DimPlot(object = panc, group.by = 'predicted_doublet',cols = c('black', 'yellow'), 
          reduction = "umap",label=TRUE,label.size=6)
  dev.off()
}

# save object so far
saveRDS(panc,file = paste(sample_id, "_processed.rds", sep=""))
