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

suppressMessages(library(EnsDb.Hsapiens.v86))
suppressMessages(library(GenomicRanges))
suppressMessages(library(future))
suppressMessages(library(optparse))
suppressMessages(library(BSgenome.Hsapiens.UCSC.hg38))

suppressMessages(library(googlesheets4))
suppressMessages(library(stringr))
suppressMessages(library(doParallel))


############## FUNCTIONS #####################
iterative_removal <- function(all_peaks.f) {
  #print(cancer.type)
  #all_peaks.f <- all_peaks.f[order(score.norm, decreasing = T), ]
  #just load existing peaks if any
  if (file.exists(paste0('peaks/',length(samples.id),'_recentered_final.',add_filename,'.tsv'))) {
    recentered_final.f <- fread(paste0('peaks/',length(samples.id),'_recentered_final.',add_filename,'.tsv'))
  } else {
    
    recentered_p=StringToGRanges(regions = all_peaks.f$new_peak, sep = c("-", "-"))
    
    cat(paste0('finding overlapping peaks \n'))
    overlapping=as.data.table(x = findOverlaps(query = recentered_p, 
                                               subject = recentered_p)) # find which peaks overlap
    overlapping=overlapping[queryHits!=subjectHits,]
    overlapping.peak.number <- unique(x = overlapping$queryHits) #these are numbers of overlapping peaks that denote their position in all_peaks.f table
    recentered_non_overlapping=all_peaks.f[-overlapping.peak.number,] # select peaks that are not overlapping as non-overlapping peaks
    fwrite(recentered_non_overlapping,paste0('peaks/',length(samples.id),'_recentered_nonOverlapping.',add_filename,'.tsv'),
           sep='\t',row.names=FALSE)
    if (length(overlapping.peak.number)>0) {
      tmp <- data.table(chr = all_peaks.f$seqnames[overlapping.peak.number], 
                        num = overlapping.peak.number)
      overlapping.peak.number.split <- split(tmp, by = 'chr', keep.by = T) #split peaks by chromosome 
      registerDoParallel(cores=25)
      #this is where iterative removal of peaks is done
      best_in_overlapping_num <- foreach(peak.numbers=overlapping.peak.number.split) %dopar% {
        cat('removing overlapping peaks in each chromosome\n')
        iterative_removal_core (peak.numbers = peak.numbers, overlapping.f = overlapping)
      }
      stopImplicitCluster()
      best_in_overlapping_num <- do.call('c', best_in_overlapping_num) #combine best peak numbers from all chromosomes
      best_in_overlapping_cancer <- all_peaks.f[best_in_overlapping_num,] #extract peaks themselves
      fwrite(best_in_overlapping_cancer,paste0('peaks/',length(samples.id),'_recentered_Overlapping.',add_filename,'.tsv'),
             sep='\t',row.names=FALSE)
      recentered_final.f=rbindlist(list(recentered_non_overlapping,best_in_overlapping_cancer))
    } else {
      recentered_final.f=recentered_non_overlapping
    }
    final.overlaps <-  recentered_final.f$new_peak %>% 
      unique %>% 
      StringToGRanges %>% 
      countOverlaps
    if (sum(final.overlaps>1)>0) {
      stop("Execution stopped. Overlapping peaks remained")
    }
 
  }
  return(recentered_final.f)
}


# this works like a charm
iterative_removal_core <- function(peak.numbers, overlapping.f) {
  chr = peak.numbers$chr[1]
  running.vector <- peak.numbers$num
  peaks.to.trash <- NULL
  peaks.to.keep <- NULL
  while (length(running.vector) != 0) {
    n <- running.vector[1] # this is the first and the best peak since peaks are sorted by scores
    neighbor.peaks.num.discard <- overlapping.f[queryHits==n, subjectHits] #find positions of other peaks overlapping with the first one 
    running.vector <- setdiff(running.vector, neighbor.peaks.num.discard) # remove them from the list of peaks
    running.vector <- setdiff(running.vector, n)
    peaks.to.keep <- c(peaks.to.keep, n) # add this peak to the keeping list
    peaks.to.trash <- unique(c(peaks.to.trash, neighbor.peaks.num.discard)) # add neighbors to the list of peaks to discard
  }
  cat('done\n')
  return(peaks.to.keep)
}

load_peaks <- function(sample, input.path.f){
  peaks=fread(paste0(input.path.f,'/', sample,"/recentered_final.filtered",sample,".tsv"))
  peaks$Sample=sample
  
  peaks$new_peak = paste(peaks$seqnames, peaks$recentered_start, peaks$recentered_end, sep = '-') #make sure all peaks are written with -
  total.score.per.mil <- sum(peaks$neg_log10qvalue_summit)/1000000 # this is scaling factor for MACS2 score
  peaks$score.norm <- peaks$neg_log10qvalue_summit / total.score.per.mil # normalize peak score in each sample aka score per million
  return(peaks)
}

getFeatureMatrix <- function (obj, peaks, pro_n, assay.towork.f) {
  frag <- Fragments(obj[[assay.towork.f]])
  cat('Making a large count matrix...\n')
  matrix.counts <- FeatureMatrix(
    fragments = frag,
    features = peaks,
    process_n = pro_n,
    sep = c("-","-"),
    cells = colnames(obj)
  )
  return(matrix.counts)
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
              metavar="character"),
  make_option(c("-a", "--assay"),
              type="character",
              default="ATAC_MACS2", 
              help="which assay should be used to merge objects? ex: X500peaksMACS2, peaks, ATAC_MACS2",
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
assay.towork <- opt$assay
species <- opt$species

if (data.type == '5sc') {
  dt.tofilter = "5' scRNA"
} else if (data.type == '3sc') {
  dt.tofilter = "3' scRNA"
} else if (data.type == 'combo' | grepl('multi|Multi', data.type)) {
  dt.tofilter = "Combo RNA"
} else if (data.type == "3' snRNA" | grepl('snRNA', data.type)) {
  dt.tofilter = c("Combo RNA", "3' snRNA")
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
                      sheet = "Patient single cells data", trim_ws = T)

samples <- samples %>% dplyr::filter(`Include in downstream` == 'Yes')
samples <- samples %>% dplyr::filter(`Data type` %in% dt.tofilter)
samples <- samples %>% dplyr::filter(Tissue == tis.tofilter)
samples <- samples %>% dplyr::filter(Species == species)
samples <- samples %>% dplyr::select(`Patient ID`,`Sample name`, `Data type`, Age, Tissue, `Rds objects`)

samples.id <- samples$`Sample name` %>% as.character()

cat (paste("Samples found:" ,length(samples.id), '\n'))

cat('creating object \n')
paths <- samples$`Rds objects`
paths

my.metadata <- fread(meta.path, data.table = F) %>% 
  data.frame(row.names = 1, check.rows = F, check.names = F) %>%
  dplyr::select(-seurat_clusters)

if (!file.exists(paste0(length(samples.id),"_Merged_not_normalized_",add_filename,".rds"))) {
    
  # make the list of objects
  registerDoParallel(cores=10)
  cat ('Reading in objects\n')
  
  obj <- foreach (s=samples.id, p = paths, pid = samples$`Patient ID`, dt = samples$`Data type`, a = samples$Age, tis = samples$Tissue, .combine=c) %dopar% {
    print(s)
    obj=readRDS(p) 
    print(paste('opened', s))
    if (dt.tofilter == 'Combo RNA' & length(dt.tofilter) == 1) {
      DefaultAssay(obj) <- assay.towork
      if (!file.exists(Fragments(obj)[[1]]@path)) stop("Urgh, this sample object can't locate fragments file")
      obj <- DietSeurat(obj, assays = c(assay.towork, 'RNA'))
    } else {
      DefaultAssay(obj) <- 'RNA'
      obj <- DietSeurat(obj, assays = 'RNA')
    }
    
    obj$Sample = s
    obj$Patient_ID = pid
    obj$Data_type = dt
    obj$Age = as.numeric(a)
    obj$Tissue = tis
    return(obj)
  }
  stopImplicitCluster()
}

if (dt.tofilter == 'Combo RNA' & length(dt.tofilter) == 1) {
  if (!file.exists(paste0(length(samples.id),"_Merged_not_normalized_",add_filename,".rds"))) {
    

    cat('doing normalization of combo object\n')
    dir.create('peaks')
  
    ############# create a common set peak
    ###########################################
    # load peaks
    registerDoParallel(cores=10)
    all_peaks <- foreach (s=samples.id,  p = paths) %dopar% {
      load_peaks (sample = s, input.path.f = input.path)
    }
    stopImplicitCluster()
    
    all_peaks <- rbindlist(all_peaks) # peaks from all samples
    all_peaks <- all_peaks[order(score.norm, decreasing = T), ] # order peaks by normalized scores, this is essential for filtering out overlapping peaks
    fwrite(all_peaks, paste0('peaks/',length(samples.id),"_sample_MACS2_peaks_",add_filename,".tsv"),
           sep='\t',row.names=FALSE)
    
    recentered_final <- iterative_removal(all_peaks)
    fwrite(recentered_final, paste0('peaks/',length(samples.id),'_recentered_final.',add_filename,'.tsv'),sep='\t',
           row.names=FALSE)
    recentered_p=StringToGRanges(unique(recentered_final$new_peak), sep = c("-", "-"))
    
    ###
    plan("multicore", workers = 20)
    options(future.globals.maxSize = 250 * 1024^3) # for 250 Gb RAM
    
    peak.number <- length(unique(recentered_final$new_peak))
    n.peaks <- round(peak.number/20)
    
    matrix.counts <- map (obj, ~getFeatureMatrix(.x, recentered_p, n.peaks, assay.towork))
    
    registerDoParallel(cores=10)
    cat ('creating assays on common peaks\n')
    obj <- foreach (obj = obj, co = matrix.counts, .combine=c) %dopar% {
      DefaultAssay(obj) <- 'RNA'
      frag = Fragments(obj[[assay.towork]])
      obj[[assay.towork]] <- NULL
      obj[['ATAC_merged']] <- CreateChromatinAssay(counts = co,
                                                   fragments=frag, 
                                                   min.cells = -1, 
                                                   min.features = -1)
      
      return(obj)
    }
    stopImplicitCluster()
    
    
    combined <- merge(x = obj[[1]], y = obj[-1], add.cell.ids = samples.id)  
    combined <- AddMetaData(combined, my.metadata)
    combined <- subset(combined, TSS.enrichment >= 4) # should remove crappy atac cells that form a cloud on the UMAP
    saveRDS(combined,  paste0(length(samples.id),"_Merged_not_normalized_",add_filename,".rds"))
  } else {
    combined <- readRDS(paste0(length(samples.id),"_Merged_not_normalized_",add_filename,".rds"))
  }
  
  # add metadata
  ifelse('passed_filters' %in% colnames(combined@meta.data),
         total_fragments_cell <- combined$passed_filters, 
         total_fragments_cell <- combined$atac_fragments)
  
  peak.counts <- colSums(x = GetAssayData(combined, slot = 'counts'))
  frip <- peak.counts *100 / total_fragments_cell
  combined <- AddMetaData(object = combined, metadata = frip, col.name = 'pct_read_in_peaks_ATAC_merged')
  combined <- AddMetaData(object = combined, metadata = peak.counts, col.name = 'peak_region_fragments_ATAC_merged')
  
  ######## Normalize RNA
  DefaultAssay(combined) <- 'RNA'
  cat('normalizing RNA\n')
  combined[["percent.mt"]] <- PercentageFeatureSet(combined, pattern = "^MT-") 
  combined <- combined %>%
    NormalizeData(assay = 'RNA') %>%
    CellCycleScoring(s.features = cc.genes$s.genes, g2m.features = cc.genes$g2m.genes, set.ident = F)
  combined <- combined %>%
    SCTransform(
      assay = 'RNA',
      vars.to.regress = c("nCount_RNA", "percent.mt", "S.Score", "G2M.Score"),
      return.only.var.genes = TRUE, verbose = F) %>%
    RunPCA(assay = 'SCT', do.print = FALSE, verbose = F) %>%
    RunUMAP(dims = 1:30,reduction.name = 'umap.rna', reduction.key = 'rnaUMAP_', verbose = F) 
  
  ######## Normalize ATAC
  DefaultAssay(combined) <- "ATAC_merged"
  
  cat('add gene annotation\n')
  annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
  genome(annotations) <- "NA"
  seqlevelsStyle(annotations) <- 'UCSC' # instead of USCS because it keeps return error https://github.com/stuart-lab/signac/issues/826
  genome(annotations) <- "hg38"
  #saveRDS(annotations, 'Annotations.EnsDb.Hsapiens.v86.rds')
  
  Annotation(combined) <- annotations
  
  cat('normalizing ATAC\n')
  combined <- combined %>% 
    RunTFIDF() %>%
    FindTopFeatures(min.cutoff = 'q5') %>%
    RunSVD(verbose = F) %>%
    RunUMAP(reduction = 'lsi', 
            dims = 2:30, reduction.name = "umap.atac", reduction.key = "atacUMAP_", verbose = F)
  
  # do wnn analysis
  cat('doing WNN\n')
  combined <- FindMultiModalNeighbors(combined, 
                                      reduction.list = list("pca", "lsi"), 
                                      dims.list = list(1:30, 2:30)) %>%
    RunUMAP( nn.name = "weighted.nn", 
             reduction.name = "wnn.umap", 
             reduction.key = "wnnUMAP_")  %>%
    FindClusters(graph.name = "wsnn", algorithm = 3,resolution = 2, verbose = T)
  
  cat('saving the object...\n')
  saveRDS(combined,  paste0(length(samples.id),"_Merged_normalized_",add_filename,".rds"))
  
  combined$Age.group <- case_when(combined$Age >= 50 ~ 'Old', 
                                  TRUE ~ 'Young')
  
  p1 <- DimPlot(combined, reduction = "umap.rna", group.by = "seurat_clusters", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("RNA")
  p2 <- DimPlot(combined, reduction = "umap.atac", group.by = "seurat_clusters", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("ATAC")
  p3 <- DimPlot(combined, reduction = "wnn.umap", group.by = "seurat_clusters", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("WNN")
  p1 + p2 + p3 & NoLegend() & theme(plot.title = element_text(hjust = 0.5))
  ggsave(paste0(length(samples.id),"_Merged_clusters_", add_filename, ".pdf"),height=10,width=30)
  
  
  p1 <- DimPlot(combined, reduction = "umap.rna", group.by = "Sample", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("RNA")
  p2 <- DimPlot(combined, reduction = "umap.atac", group.by = "Sample", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("ATAC")
  p3 <- DimPlot(combined, reduction = "wnn.umap", group.by = "Sample", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("WNN")
  p1 + p2 + p3 & NoLegend() & theme(plot.title = element_text(hjust = 0.5))
  ggsave(paste0(length(samples.id),"_Merged_Sample_", add_filename, ".pdf"),height=10,width=30)
  
  
  p1 <- DimPlot(combined, reduction = "umap.rna", group.by = "Age.group", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("RNA")
  p2 <- DimPlot(combined, reduction = "umap.atac", group.by = "Age.group", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("ATAC")
  p3 <- DimPlot(combined, reduction = "wnn.umap", group.by = "Age.group", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("WNN")
  p1 + p2 + p3 & NoLegend() & theme(plot.title = element_text(hjust = 0.5))
  ggsave(paste0(length(samples.id),"_Merged_Age.group_", add_filename, ".pdf"),height=10,width=30)
  
  
  
  p1 <- DimPlot(combined, reduction = "umap.rna", group.by = "cell_type", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("RNA")
  p2 <- DimPlot(combined, reduction = "umap.atac", group.by = "cell_type", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("ATAC")
  p3 <- DimPlot(combined, reduction = "wnn.umap", group.by = "cell_type", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("WNN")
  p1 + p2 + p3 & NoLegend() & theme(plot.title = element_text(hjust = 0.5))
  ggsave(paste0(length(samples.id),"_Merged_cell_type_", add_filename, ".pdf"),height=10,width=30)
  
  
  fwrite(cbind(Embeddings(combined, reduction = "wnn.umap"),
               Embeddings(combined, reduction = "umap.rna"),
               Embeddings(combined, reduction = "umap.atac"),
               combined@meta.data), 
         paste0(length(samples.id),"_Merged_normalized_",add_filename,".metadata.tsv"), sep = '\t', row.names = T)

  
  } else {
  cat('doing normalization of reg RNA object\n')
    
  combined <- merge(x = obj[[1]], y = obj[-1], add.cell.ids = samples.id)  
  

  cat('normalizing RNA\n')
  DefaultAssay(combined) <- 'RNA'
  
  combined <- combined %>%
    NormalizeData(assay = 'RNA') %>%
    CellCycleScoring(s.features = cc.genes$s.genes, g2m.features = cc.genes$g2m.genes, set.ident = F)
  
  combined[["percent.mt"]] <- PercentageFeatureSet(combined, pattern = "^MT-") 
  print(head(combined@meta.data))
  #print(head(combined$G2M.Score))
  
  combined <- combined %>%
    SCTransform(
      assay = 'RNA',
      vars.to.regress = c("nCount_RNA", 
                          "percent.mt", 
                          "S.Score", 
                          "G2M.Score"
                          ),
      return.only.var.genes = TRUE, verbose = T) %>%
    RunPCA(assay = 'SCT', do.print = FALSE, verbose = T) %>%
    RunUMAP(dims = 1:30,reduction.name = 'umap.rna', reduction.key = 'rnaUMAP_', verbose = T) %>%
    FindNeighbors( dims = 1:30) %>%
    FindClusters(resolution = 2, verbose = FALSE)
  
  cat('saving the object...\n')
  saveRDS(combined,  paste0(length(samples.id),"_Merged_normalized_",add_filename,".rds"))
  
  combined <- AddMetaData(combined, my.metadata)
  combined$Age.group <- case_when(combined$Age >= 50 ~ 'Old', 
                                  TRUE ~ 'Young')
  
  DimPlot(combined, reduction = "umap.rna", group.by = "seurat_clusters", label = TRUE, label.size = 2.5, repel = TRUE)
  ggsave(paste0(length(samples.id),"_Merged_clusters_", add_filename, ".pdf"),height=10,width=11)
  DimPlot(combined, reduction = "umap.rna", group.by = "Sample", label = TRUE, label.size = 2.5, repel = TRUE)
  ggsave(paste0(length(samples.id),"_Merged_Sample_", add_filename, ".pdf"),height=10,width=12)
  DimPlot(combined, reduction = "umap.rna", group.by = "Age.group", label = TRUE, label.size = 2.5, repel = TRUE)
  ggsave(paste0(length(samples.id),"_Merged_Age.group_", add_filename, ".pdf"),height=10,width=11)
  
 
  DimPlot(combined, reduction = "umap.rna", group.by = "cell_type", label = TRUE, label.size = 2.5, repel = TRUE)
  ggsave(paste0(length(samples.id),"_Merged_cell_type_", add_filename, ".pdf"),height=10,width=11)
  
  fwrite(cbind(Embeddings(combined, reduction = "umap.rna"),
               combined@meta.data), 
         paste0(length(samples.id),"_Merged_normalized_",add_filename,".metadata.tsv"), sep = '\t', row.names = T)
  
  cat('saving the object...\n')
  saveRDS(combined,  paste0(length(samples.id),"_Merged_normalized_",add_filename,".rds"))
  
}

