#ATAC-seq rn7 analysis
set.seed(1234)
#Load libraries
library(Seurat)
library(SeuratObject)
library(SeuratData)
library(Signac)
library(ggplot2)
library(devtools)
library(dplyr)
library(GenomeInfoDb)
library(rtracklayer)

#Filtering datasets
#Read h5 files
FS_counts <- Read10X_h5(filename = '/Users/ethanwan/Desktop/rn7_ATAC/FS/filtered_peak_bc_matrix.h5')
FC_counts <- Read10X_h5(filename = '/Users/ethanwan/Desktop/rn7_ATAC/FC/filtered_peak_bc_matrix.h5')
MS_counts <- Read10X_h5(filename = '/Users/ethanwan/Desktop/rn7_ATAC/MS/filtered_peak_bc_matrix.h5')
MC_counts <- Read10X_h5(filename = '/Users/ethanwan/Desktop/rn7_ATAC/mC/filtered_peak_bc_matrix.h5')

#Load rn7 GTF file for annotation
rn7_gtf <- rtracklayer::import("/Users/ethanwan/Desktop/Rattus_norvegicus.mRatBN7.2.105.gtf")
rn7_gtf <- rn7_gtf[rn7_gtf$type == "gene"]
rn7_gtf <- keepStandardChromosomes(rn7_gtf, pruning.mode = "coarse")

#Metadata creation for all groups
FS_metadata <- read.csv(
  file = '/Users/ethanwan/Desktop/rn7_ATAC/FS/singlecell.csv',
  header = TRUE,
  row.names = 1
)
FC_metadata <- read.csv(
  file = '/Users/ethanwan/Desktop/rn7_ATAC/FC/singlecell.csv',
  header = TRUE,
  row.names = 1
)
MS_metadata <- read.csv(
  file = '/Users/ethanwan/Desktop/rn7_ATAC/MS/singlecell.csv',
  header = TRUE,
  row.names = 1
)
MC_metadata <- read.csv(
  file = '/Users/ethanwan/Desktop/rn7_ATAC/MC/singlecell.csv',
  header = TRUE,
  row.names = 1
)

#Chromatin Assay creation for all groups
FS_assay <- CreateChromatinAssay(
  counts = FS_counts,
  sep = c(":","-"),
  genome = "rn7",
  fragments = "/Users/ethanwan/Desktop/rn7_ATAC/FS/fragments.tsv.gz",
  min.cells = 10
)
FC_assay <- CreateChromatinAssay(
  counts = FC_counts,
  sep = c(":","-"),
  genome = "rn7",
  fragments = "/Users/ethanwan/Desktop/rn7_ATAC/FC/fragments.tsv.gz",
  min.cells = 10
)
MS_assay <- CreateChromatinAssay(
  counts = MS_counts,
  sep = c(":","-"),
  genome = "rn7",
  fragments = "/Users/ethanwan/Desktop/rn7_ATAC/MS/fragments.tsv.gz",
  min.cells = 10
)
MC_assay <- CreateChromatinAssay(
  counts = MC_counts,
  sep = c(":","-"),
  genome = "rn7",
  fragments = "/Users/ethanwan/Desktop/rn7_ATAC/MC/fragments.tsv.gz",
  min.cells = 10
)

#Creation of Seurat objects for all groups
Fem_Sal <- CreateSeuratObject(
  counts = FS_assay,
  assay = 'peaks',
  meta.data = FS_metadata
)
Fem_Coc <- CreateSeuratObject(
  counts = FC_assay,
  assay = 'peaks',
  meta.data = FC_metadata
)
Male_Sal <- CreateSeuratObject(
  counts = MS_assay,
  assay = 'peaks',
  meta.data = MS_metadata
)
Male_Coc <- CreateSeuratObject(
  counts = MC_assay,
  assay = 'peaks',
  meta.data = MC_metadata
)

#Adding gene annotations to all groups
Annotation(Fem_Sal) <- rn7_gtf
Annotation(Fem_Coc) <- rn7_gtf
Annotation(Male_Sal) <- rn7_gtf
Annotation(Male_Coc) <- rn7_gtf

#QC metrics
Fem_Sal <- NucleosomeSignal(object = Fem_Sal)
Fem_Coc <- NucleosomeSignal(object = Fem_Coc)
Male_Sal <- NucleosomeSignal(object = Male_Sal)
Male_Coc <- NucleosomeSignal(object = Male_Coc)

Fem_Sal <- TSSEnrichment(object = Fem_Sal, fast = FALSE, verbose = TRUE)
Fem_Coc <- TSSEnrichment(object = Fem_Coc, fast = FALSE, verbose = TRUE)
Male_Sal <- TSSEnrichment(object = Male_Sal, fast = FALSE, verbose = TRUE)
Male_Coc <- TSSEnrichment(object = Male_Coc, fast = FALSE, verbose = TRUE)

Fem_Sal$pct_reads_in_peaks <- Fem_Sal$peak_region_fragments / Fem_Sal$passed_filters * 100
Fem_Sal$blacklist_ratio <- Fem_Sal$blacklist_region_fragments/Fem_Sal$peak_region_fragments

Fem_Coc$pct_reads_in_peaks <- Fem_Coc$peak_region_fragments / Fem_Coc$passed_filters * 100
Fem_Coc$blacklist_ratio <- Fem_Coc$blacklist_region_fragments / Fem_Coc$peak_region_fragments

Male_Sal$pct_reads_in_peaks <- Male_Sal$peak_region_fragments / Male_Sal$passed_filters * 100
Male_Sal$blacklist_ratio <- Male_Sal$blacklist_region_fragments / Male_Sal$peak_region_fragments

Male_Coc$pct_reads_in_peaks <- Male_Coc$peak_region_fragments / Male_Coc$passed_filters * 100
Male_Coc$blacklist_ratio <- Male_Coc$blacklist_region_fragments / Male_Coc$peak_region_fragments

#TSS Enrichment plots
Fem_Sal$high.tss <- ifelse(Fem_Sal$TSS.enrichment > 2, 'High', 'Low')
Fem_Sal_TSS <- TSSPlot(Fem_Sal, group.by = 'high.tss') + NoLegend() 
ggsave(Fem_Sal_TSS, file = '/Users/ethanwan/Desktop/rn7_ATAC/FS/Fem_Sal_TSS.pdf', width = 16, height = 12)

Fem_Coc$high.tss <- ifelse(Fem_Coc$TSS.enrichment > 2, 'High', 'Low')
Fem_Coc_TSS <- TSSPlot(Fem_Coc, group.by = 'high.tss') + NoLegend() 
ggsave(Fem_Coc_TSS, file = '/Users/ethanwan/Desktop/rn7_ATAC/FC/Fem_Coc_TSS.pdf', width = 16, height = 12)

Male_Sal$high.tss <- ifelse(Male_Sal$TSS.enrichment > 2, 'High', 'Low')
Male_Sal_TSS <- TSSPlot(Male_Sal, group.by = 'high.tss') + NoLegend() 
ggsave(Male_Sal_TSS, file = '/Users/ethanwan/Desktop/rn7_ATAC/MS/Male_Sal_TSS.pdf', width = 16, height = 12)

Male_Coc$high.tss <- ifelse(Male_Coc$TSS.enrichment > 2, 'High', 'Low')
Male_Coc_TSS <- TSSPlot(Male_Coc, group.by = 'high.tss') + NoLegend()
ggsave(Male_Coc_TSS, file = '/Users/ethanwan/Desktop/rn7_ATAC/MC/Male_Coc_TSS.pdf', width = 16, height = 12)

#Vln Plots
Fem_Sal_Vln <- VlnPlot(
  object = Fem_Sal,
  features = c('pct_reads_in_peaks', 'peak_region_fragments', 'TSS.enrichment', 'blacklist_ratio', 'nucleosome_signal'),
  pt.size = 0.1,
  ncol = 5
)
ggsave(Fem_Sal_Vln, file = '/Users/ethanwan/Desktop/rn7_ATAC/FS/VlnPlot.pdf', width = 16, height = 12)

Fem_Coc_Vln <- VlnPlot(
  object = Fem_Coc,
  features = c('pct_reads_in_peaks', 'peak_region_fragments', 'TSS.enrichment', 'blacklist_ratio', 'nucleosome_signal'),
  pt.size = 0.1,
  ncol = 5
)
ggsave(Fem_Coc_Vln, file = '/Users/ethanwan/Desktop/rn7_ATAC/FC/VlnPlot.pdf', width = 16, height = 12)

Male_Sal_Vln <- VlnPlot(
  object = Male_Sal,
  features = c('pct_reads_in_peaks', 'peak_region_fragments', 'TSS.enrichment', 'blacklist_ratio', 'nucleosome_signal'),
  pt.size = 0.1,
  ncol = 5
)
ggsave(Male_Sal_Vln, file = '/Users/ethanwan/Desktop/rn7_ATAC/MS/VlnPlot.pdf', width = 16, height = 12)

Male_Coc_Vln <- VlnPlot(
  object = Male_Coc,
  features = c('pct_reads_in_peaks', 'peak_region_fragments', 'TSS.enrichment', 'blacklist_ratio', 'nucleosome_signal'),
  pt.size = 0.1,
  ncol = 5
)
ggsave(Male_Coc_Vln, file = '/Users/ethanwan/Desktop/rn7_ATAC/MC/VlnPlot.pdf', width = 16, height = 12)

#Subsetting after QC metrics are evaluated
FS <- subset(Fem_Sal, 
             subset = peak_region_fragments > 1000 &
               peak_region_fragments < 50000 &
               pct_reads_in_peaks > 10 &
               nucleosome_signal < 4 &
               TSS.enrichment > 2)

FC <- subset(Fem_Coc, 
             subset = peak_region_fragments > 1000 &
               peak_region_fragments < 50000 &
               pct_reads_in_peaks > 10 &
               nucleosome_signal < 4 &
               TSS.enrichment > 2)

MS <- subset(Male_Sal, 
             subset = peak_region_fragments > 1000 &
               peak_region_fragments < 50000 &
               pct_reads_in_peaks > 10 &
               nucleosome_signal < 4 &
               TSS.enrichment > 2)

MC <- subset(Male_Coc, 
             subset = peak_region_fragments > 1000 &
               peak_region_fragments < 50000 &
               pct_reads_in_peaks > 10 &
               nucleosome_signal < 4 &
               TSS.enrichment > 2)

#Normalization and linear dimension reduction
FS <- RunTFIDF(FS)
FC <- RunTFIDF(FC)
MS <- RunTFIDF(MS)
MC <- RunTFIDF(MC)

FS <- FindTopFeatures(FS, min.cutoff = 'q0')
FC <- FindTopFeatures(FC, min.cutoff = 'q0')
MS <- FindTopFeatures(MS, min.cutoff = 'q0')
MC <- FindTopFeatures(MC, min.cutoff = 'q0')

FS <- RunSVD(object = FS)
FC <- RunSVD(object = FC)
MS <- RunSVD(object = MS)
MC <- RunSVD(object = MC)

FS$Stim <- 'Saline'
FC$Stim <- 'Cocaine'
MS$Stim <- 'Saline'
MC$Stim <- 'Cocaine'

FS$Sex <- 'Female'
FC$Sex <- 'Female'
MS$Sex <- 'Male'
MC$Sex <- 'Male'

FS$dataset <- 'FS'
FC$dataset <- 'FC'
FS$dataset <- 'MS'
FS$dataset <- 'MC'

#Clustering and dimensionality reduction for individual objects
FS <- RunUMAP(object = FS, reduction = 'lsi', dims = 2:30)
FS <- FindNeighbors(object = FS, reduction = 'lsi', dims = 2:30)
FS <- FindClusters(object = FS, verbose = FALSE, algorithm = 1)
DimPlot(object = FS, label = TRUE) + NoLegend()

FC <- RunUMAP(object = FC, reduction = 'lsi', dims = 2:30)
FC <- FindNeighbors(object = FC, reduction = 'lsi', dims = 2:30)
FC <- FindClusters(object = FC, verbose = FALSE, algorithm = 1)
DimPlot(object = FC, label = TRUE) + NoLegend()

MS <- RunUMAP(object = MS, reduction = 'lsi', dims = 2:30)
MS <- FindNeighbors(object = MS, reduction = 'lsi', dims = 2:30)
MS <- FindClusters(object = MS, verbose = FALSE, algorithm = 1)
DimPlot(object = MS, label = TRUE) + NoLegend()

MC <- RunUMAP(object = MC, reduction = 'lsi', dims = 2:30)
MC <- FindNeighbors(object = MC, reduction = 'lsi', dims = 2:30)
MC <- FindClusters(object = MC, verbose = FALSE, algorithm = 1)
DimPlot(object = MC, label = TRUE) + NoLegend()

#Add RNA assay and create gene activities
gene.activities_FS <- GeneActivity(FS)
gene.activities_FC <- GeneActivity(FC)
gene.activities_MS <- GeneActivity(MS)
gene.activities_MC <- GeneActivity(MC)

FS[['RNA']] <- CreateAssayObject(counts = gene.activities_FS)
FC[['RNA']] <- CreateAssayObject(counts = gene.activities_FC)
MS[['RNA']] <- CreateAssayObject(counts = gene.activities_MS)
MC[['RNA']] <- CreateAssayObject(counts = gene.activities_MC)

saveRDS(FS, file = '/Users/ethanwan/Desktop/rn7_ATAC/FS/FS_filtered.rds')
saveRDS(FC, file = '/Users/ethanwan/Desktop/rn7_ATAC/FC/FC_filtered.rds')
saveRDS(MS, file = '/Users/ethanwan/Desktop/rn7_ATAC/MS/MS_filtered.rds')
saveRDS(MC, file = '/Users/ethanwan/Desktop/rn7_ATAC/MC/MC_filtered.rds')

#Reading cellranger peaks to combine the samples into one object
FS_peaks <- read.table(
  file = "/Users/ethanwan/Desktop/rn7_ATAC/FS/peaks.bed",
  col.names = c("chr", "start", "end")
)
FC_peaks <- read.table(
  file = "/Users/ethanwan/Desktop/rn7_ATAC/FC/peaks.bed",
  col.names = c("chr", "start", "end")
)
MS_peaks <- read.table(
  file = "/Users/ethanwan/Desktop/rn7_ATAC/MS/peaks.bed",
  col.names = c("chr", "start", "end")
)
MC_peaks <- read.table(
  file = "/Users/ethanwan/Desktop/rn7_ATAC/MC/peaks.bed",
  col.names = c("chr", "start", "end")
)

#Converting peaks to GRanges
FS_range <- makeGRangesFromDataFrame(FS_peaks)
FC_range <- makeGRangesFromDataFrame(FC_peaks)
MS_range <- makeGRangesFromDataFrame(MS_peaks)
MC_range <- makeGRangesFromDataFrame(MC_peaks)

#Combine all peaks and filter out bad peaks, these will be our anchor features
combo_peaks <- reduce(c(FS_range, FC_range, MS_range, MC_range))
peakwidths <- width(combo_peaks)
combo_peaks <- combo_peaks[peakwidths < 10000 & peakwidths > 20]
combo_peaks

#Find cells with more than 1000 counts and also filter previous objects in the same way
FS_pass <- names(which(FS$passed_filters > 1000))
FC_pass <- names(which(FC$passed_filters > 1000))
MS_pass <- names(which(MS$passed_filters > 1000))
MC_pass <- names(which(MC$passed_filters > 1000))

#Quantifying peaks in individual datasets
FS_fragpath <- '/Users/ethanwan/Desktop/rn7_ATAC/FS/fragments.tsv.gz'
FC_fragpath <- '/Users/ethanwan/Desktop/rn7_ATAC/FC/fragments.tsv.gz'
MS_fragpath <- '/Users/ethanwan/Desktop/rn7_ATAC/MS/fragments.tsv.gz'
MC_fragpath <- '/Users/ethanwan/Desktop/rn7_ATAC/MC/fragments.tsv.gz'

#Creating fragment objects
FS_frags <- CreateFragmentObject(
  path = FS_fragpath,
  cells = FS_pass
)
FC_frags <- CreateFragmentObject(
  path = FC_fragpath,
  cells = FC_pass
)
MS_frags <- CreateFragmentObject(
  path = MS_fragpath,
  cells = MS_pass
)
MC_frags <- CreateFragmentObject(
  path = MC_fragpath,
  cells = MC_pass
)

#Quantifying peaks in each set
FS.counts <- FeatureMatrix(
  fragments = FS_frags, 
  features = combo_peaks,
  cells = FS_pass
)
FC.counts <- FeatureMatrix(
  fragments = FC_frags, 
  features = combo_peaks,
  cells = FC_pass
)
MS.counts <- FeatureMatrix(
  fragments = MS_frags, 
  features = combo_peaks,
  cells = MS_pass
)
MC.counts <- FeatureMatrix(
  fragments = MC_frags, 
  features = combo_peaks,
  cells = MC_pass
)

#Creating New Seurat objects to be integrated together 
FS.assay <- CreateChromatinAssay(FS.counts, fragments = FS_frags)
FS_reduced <- CreateSeuratObject(FS.assay, assay = 'ATAC')

FC.assay <- CreateChromatinAssay(FC.counts, fragments = FC_frags)
FC_reduced <- CreateSeuratObject(FC.assay, assay = 'ATAC')

MS.assay <- CreateChromatinAssay(MS.counts, fragments = MS_frags)
MS_reduced <- CreateSeuratObject(MS.assay, assay = 'ATAC')

MC.assay <- CreateChromatinAssay(MC.counts, fragments = MC_frags)
MC_reduced <- CreateSeuratObject(MC.assay, assay = 'ATAC')

#Adding some metadata information
FS_reduced$dataset <- 'FS'
FC_reduced$dataset <- 'FC'
MS_reduced$dataset <- 'MS'
MC_reduced$dataset <- 'MC'

FS_reduced$Sex <- 'Female'
FC_reduced$Sex <- 'Female'
MS_reduced$Sex <- 'Male'
MC_reduced$Sex <- 'Male'

FS_reduced$Stim <- 'Saline'
FC_reduced$Stim <- 'Cocaine'
MS_reduced$Stim <- 'Saline'
MC_reduced$Stim <- 'Cocaine'

#Saving objects
saveRDS(FS_reduced, file = '/Users/ethanwan/Desktop/rn7_ATAC/FS/FS_reduced.rds')
saveRDS(FC_reduced, file = '/Users/ethanwan/Desktop/rn7_ATAC/FC/FC_reduced.rds')
saveRDS(MS_reduced, file = '/Users/ethanwan/Desktop/rn7_ATAC/MS/MS_reduced.rds')
saveRDS(MC_reduced, file = '/Users/ethanwan/Desktop/rn7_ATAC/MC/MC_reduced.rds')
#Clear environment
rm(list = ls())

FS_reduced <- readRDS('/Users/ethanwan/Desktop/rn7_ATAC/FS/FS_reduced.rds')
FC_reduced <- readRDS('/Users/ethanwan/Desktop/rn7_ATAC/FC/FC_reduced.rds')
MS_reduced <- readRDS('/Users/ethanwan/Desktop/rn7_ATAC/MS/MS_reduced.rds')
MC_reduced <- readRDS('/Users/ethanwan/Desktop/rn7_ATAC/MC/MC_reduced.rds')

#Merging new datasets and running normalization
combined_reduced <- merge(x = FS_reduced, y = list(FC_reduced, MS_reduced, MC_reduced), add.cell.ids = c('FS', 'FC', 'MS', 'MC'))

#Loading original objects
FS <- readRDS(file = '/Users/ethanwan/Desktop/rn7_ATAC/FS/FS_filtered.rds')
FC <- readRDS(file = '/Users/ethanwan/Desktop/rn7_ATAC/FC/FC_filtered.rds')
MS <- readRDS(file = '/Users/ethanwan/Desktop/rn7_ATAC/MS/MS_filtered.rds')
MC <- readRDS(file = '/Users/ethanwan/Desktop/rn7_ATAC/MC/MC_filtered.rds')

DefaultAssay(FS) <- 'peaks'
DefaultAssay(FC) <- 'peaks'
DefaultAssay(MS) <- 'peaks'
DefaultAssay(MC) <- 'peaks'

FS <- RenameCells(object = FS,add.cell.id = "FS")
FC <- RenameCells(object = FC,add.cell.id = "FC")
MS <- RenameCells(object = MS,add.cell.id = "MS")
MC <- RenameCells(object = MC,add.cell.id = "MC")

#Creating merged counts for RNA assay creation
Merged_counts <- cbind(as.matrix(GetAssayData(object = FS,assay = "RNA")),
                       as.matrix(GetAssayData(object = FC,assay = "RNA")),
                       as.matrix(GetAssayData(object = MS,assay = "RNA")),
                       as.matrix(GetAssayData(object = MC,assay = "RNA")))

#Create new RNA assay in merged
combined_reduced[['RNA']] <- CreateAssayObject(counts = Merged_counts)

#Now normalize the RNA counts 
combined_reduced <- NormalizeData(object = combined_reduced,
                                  assay = "RNA",
                                  normalization.method = "LogNormalize",
                                  scale.factor = median(combined_reduced$nCount_RNA))

#Find variable features in the same manner that we found them in the Repeated scRNA analysis 
DefaultAssay(combined_reduced) <- "RNA"
combined_reduced <- FindVariableFeatures(combined_reduced, selection.method = "vst", nfeatures = 2000)

#Now scale the data 
combined_reduced <- ScaleData(combined_reduced, verbose = TRUE)  

#Change default assay to ATAC
DefaultAssay(combined_reduced) <- "ATAC"

#Run dimensionality reduction and build the umap
combined_reduced <- RunTFIDF(combined_reduced)
combined_reduced <- FindTopFeatures(combined_reduced, min.cutoff = 'q0')
combined_reduced <- RunSVD(combined_reduced)
combined_reduced <- RunUMAP(combined_reduced, dims = 2:30, reduction = 'lsi')
combined_reduced <- FindNeighbors(combined_reduced, dims = 2:30, reduction = 'lsi')
combined_reduced <- FindClusters(combined_reduced, resolution = 0.20)
DimPlot(combined_reduced, reduction = "umap")
saveRDS(combined_reduced, file = '/Users/ethanwan/Desktop/rn7_ATAC/combined_reduced.rds')
rm(list = ls())

#Integrating with repeated cocaine snRNAseq object
combined_reduced <- readRDS('/Users/ethanwan/Desktop/rn7_ATAC/combined_reduced.rds')
All_Groups_log_Doublets <- readRDS('/Users/ethanwan/Desktop/rn7_ATAC/All_Groups_log_Doublets.rds')
Annotation(combined_reduced) <- rn7_gtf

#Find the transfer anchors 
reduced_anchors <- FindTransferAnchors(reference       = All_Groups_log_Doublets,
                                       query           = combined_reduced,
                                       features        = VariableFeatures(All_Groups_log_Doublets),
                                       reference.assay = "RNA",
                                       query.assay     = "RNA",
                                       reduction       = "cca")
reduced_celltype_predictions <- TransferData(anchorset = reduced_anchors, 
                                             refdata = All_Groups_log_Doublets$CellType, 
                                             weight.reduction = combined_reduced[["lsi"]], 
                                             dims = 2:30)
combined_reduced <- AddMetaData(combined_reduced, metadata = reduced_celltype_predictions)

hist(combined_reduced$prediction.score.max,main = "Histogram of Predicted Cell Type Prediction Scores for Reduced CellRanger Peaks")
abline(v = 0.5, col = "red")

table(combined_reduced$prediction.score.max > 0.5) #10085 kept

#Keep only those cells with a prediction score greater than 50% 
combined_filtered <- subset(combined_reduced, subset = prediction.score.max > 0.5)
#Make the colors match
combined_filtered$predicted.id <- factor(combined_filtered$predicted.id, levels = levels(All_Groups_log_Doublets))

#Make the plots 
p1 <-  DimPlot(All_Groups_log_Doublets, group.by = "CellType", label = TRUE, repel = TRUE) +
  ggtitle("snRNA-seq cells") + 
  NoLegend()

p2 <- DimPlot(combined_filtered, group.by = "predicted.id", label = TRUE, repel = TRUE) + 
  ggtitle("snATAC-seq cells") + 
  NoLegend() +
  scale_colour_hue(drop = FALSE)

p1+p2

rn7_gtf_all <- rtracklayer::import("/Users/ethanwan/Desktop/Rattus_norvegicus.mRatBN7.2.105.gtf")
rn7_gtf_all <- keepStandardChromosomes(rn7_gtf_all, pruning.mode = "coarse")
rn7_gtf_all$tx_id <- rn7_gtf_all$transcript_id

Annotation(combined_filtered) <- rn7_gtf_all

Idents(combined_filtered) <- combined_filtered$predicted.id
CoveragePlot(
  object = combined_filtered,
  region = "3-116907322-116908266",
  annotation = TRUE,
  extend.upstream = 10000,
  extend.downstream = 60000
)

