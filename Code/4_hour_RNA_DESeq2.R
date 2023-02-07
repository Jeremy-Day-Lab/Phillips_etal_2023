#Script is based on NBI DESeq2 workflow by Lara Ianov which is based off DESeq2 vignette
#Set seed and working directory
set.seed(1234)
setwd(dir = "/Volumes/Day-Lab$/RobertPhillips/Bioinformatics/Phillips_etal_2023/4_hour_RNA/")

#load libraries
library(DESeq2)
library(ggplot2)
library(gridExtra)
library(Seurat)
library(biomaRt)
library(pheatmap) 


##################################
######4 hour RNA-seq dataset######
##################################

# Read in fixed_raw_counts.txt file, removes chromosome, start, and end columns
rawCounts <- read.table("fixed_raw_counts.txt",
                        header   = TRUE, 
                        row.names= 1)
rawCounts <- rawCounts[-c(1:3)]


# create coldata dataframe
colData <- data.frame(Sample    = colnames(rawCounts),
                      treatment = c(rep("KCl",4),rep("Veh",4)))

# Perfrom DDS
dds <- DESeqDataSetFromMatrix(countData = rawCounts,
                              colData   = colData,
                              design    = ~treatment)

dds


############## Exploring the data by visualization - GOOD QA-QC ################

# PCA plot
rld <- rlog(dds, blind=FALSE)
pca_data <- plotPCA(rld, intgroup=c("treatment"),returnData = TRUE)
pca_1 <- ggplot(data = pca_data,aes(x = PC1, y = PC2,color = treatment)) +
  geom_point(size = 3) +
  scale_color_manual(values = c("KCl" = "orange",
                                "Veh" = "gray56")) +
  theme_bw() +
  NoGrid() +
  labs(x     = "PC1: 97% variance",
       y     = "PC2: 1% variance",
       title = "4 hour 10mM KCl\nRNA-Seq") +
  theme(plot.title = element_text(hjust = 0.5))
ggsave(filename = "plots/PCA_4hours.pdf",plot = pca_1)


# Heatmap and hierarchical clustering with Euclidean distances
sampleDists <- dist(t(assay(rld))) 
library(RColorBrewer)
library(pheatmap)
sampleDistMatrix <- as.matrix(sampleDists)
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)


############## Differential expression ################
# For differential expression with outlier detection as default, 
# where second argument in contrast is the 'control', as this is the denominator.
# Thus all results of the fold change direction will be relative to numerator (treatment).

# DESeq function
dds <- DESeq(dds)

# Results function. Prints summary of differentially expressed genes in Group_1 vs Group_2 comparison.
res <- results(dds, contrast=c("treatment",
                               "KCl",
                               "Veh"), 
               alpha = 0.05,
               pAdjustMethod = "BH")

summary(res)
# out of 22408 with nonzero total read count
# adjusted p-value < 0.05
# LFC > 0 (up)       : 3648, 16%
# LFC < 0 (down)     : 3383, 15%
# outliers [1]       : 1, 0.0045%
# low counts [2]     : 6392, 29%
# (mean count < 4)
# [1] see 'cooksCutoff' argument of ?results
# [2] see 'independentFiltering' argument of ?results

sum(res$padj < 0.05, na.rm=TRUE)
# [1] 7031

# Order results by padj
resOrdered <- res[order(res$padj),]

#Remove NAs from adjusted p values 
resFiltered <- as.data.frame(resOrdered[!is.na(resOrdered$padj), ])

#Make the rownames a column that will be useful for later merging 
resFiltered$ensembl_gene_id <- row.names(resFiltered)

#Subset for genes with baseMean >= 50
resFiltered <- subset(resFiltered,subset=(baseMean >= 50))

#Filter for genes with a adj.p.val < 0.05 and |log2FC| > 0.5
resFiltered_sig <- subset(resFiltered, padj<0.05 & log2FoldChange > 0.5 | padj < 0.05 & log2FoldChange < (-0.5))

#Make the rownames a column that will be useful for later merging. 
resFiltered_sig$ensembl_gene_id <- row.names(resFiltered_sig)

############### Assign Gene IDs and coordinates to resOrdered file using biomaRt ##############
#For all genes
# Specify Ensembl anotation version. Use version which matches gtf file version used in primary pipeline.
ensembl <- useMart("ensembl",
                   dataset = "rnorvegicus_gene_ensembl") #Should pull the latest ensembl (This analysis was completed on 02/07/2023)

#Pull attributes
genes <- getBM(attributes = c("ensembl_gene_id", 
                              "external_gene_name",
                              "description",
                              "chromosome_name",
                              "start_position",
                              "end_position",
                              "strand"),
               filters = "ensembl_gene_id", 
               values  = resFiltered$ensembl,
               mart    = ensembl)

resFiltered <- merge(x  = as.data.frame(resFiltered),
                     y  = genes,
                     by = "ensembl_gene_id")


#Write out the tables
write.table(x         = resFiltered,
            file      = "Gene_Lists/FourHour_AllGenes.txt",
            sep       = "\t",
            col.names = TRUE,
            row.names = FALSE,
            quote     = FALSE)

write.table(x         = resFiltered_sig,
            file      = "Gene_Lists/FourHour_DEGs.txt",
            sep       = "\t",
            col.names = TRUE,
            row.names = FALSE,
            quote     = FALSE)

write.table(x         = resFiltered_sig$ensembl_gene_id,
            file      = "Gene_Lists/FourHour_DEGs_ensemblID.txt",
            sep       = "\t",
            col.names = TRUE,
            row.names = FALSE,
            quote     = FALSE)


# Volcano plot
four_hour_vp <- ggplot(data = as.data.frame(resFiltered),aes(x = log2FoldChange,y=abs(stat))) + 
  geom_point(colour = "grey",size = 3) +
  geom_point(data = subset(as.data.frame(resFiltered), padj<.05 & log2FoldChange>0.5),colour = "red",size = 3) + 
  geom_point(data = subset(as.data.frame(resFiltered), padj<.05 & log2FoldChange<(-0.5)),colour = "dodger blue",size = 3) + 
  geom_point(data = subset(resFiltered,subset=(external_gene_name == "Fos" | 
                                                 external_gene_name == "Fosb" |
                                                 external_gene_name == "Junb" |
                                                 external_gene_name == "Pdyn" |
                                                 external_gene_name == "Scg2" |
                                                 external_gene_name == "Slc7a1" |
                                                 external_gene_name == "Kcnf1" |
                                                 external_gene_name == "LRRTM1")),
             color = "black",
             shape= 1,
             size = 3) +
  xlim(c(-5.5,5.5)) +
  geom_vline(xintercept = 0.5,lty = 2,colour = "black") + 
  geom_vline(xintercept = -0.5,lty = 2,colour = "black") +
  geom_text(data = subset(resFiltered,subset=(external_gene_name == "Fos" | 
                                                external_gene_name == "Fosb" |
                                                external_gene_name == "Junb" |
                                                external_gene_name == "Pdyn" |
                                                external_gene_name == "Scg2" |
                                                external_gene_name == "Slc7a1"|
                                                external_gene_name == "Kcnf1" |
                                                external_gene_name == "LRRTM1"),size = 3),
            aes(label = external_gene_name),
            hjust=0, 
            vjust=0) +
  ggtitle("4 Hour 10mM KCl \n1,680 DEGs") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5)) +
  NoGrid() +
  annotate(geom = "text",x = -5.5,y = 72,label = nrow(subset(resFiltered_sig,subset=(log2FoldChange > 0))),color = "red") +
  annotate(geom = "text",x = -5.5,y = 70,label = nrow(subset(resFiltered_sig,subset=(log2FoldChange < 0))),color = "dodger blue") 

ggsave(plot = four_hour_vp,
       filename = "plots/4Hour_DEGs_VolcanoPlot.pdf",
       height   = 8,
       width    = 8)


sessionInfo()
# R version 4.0.3 (2020-10-10)
# Platform: x86_64-apple-darwin17.0 (64-bit)
# Running under: macOS Big Sur 10.16
# 
# Matrix products: default
# LAPACK: /Library/Frameworks/R.framework/Versions/4.0/Resources/lib/libRlapack.dylib
# 
# locale:
#   [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
# 
# attached base packages:
#   [1] parallel  stats4    stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
#   [1] RColorBrewer_1.1-2          pheatmap_1.0.12             biomaRt_2.46.3              SeuratObject_4.0.0          Seurat_4.0.1               
# [6] gridExtra_2.3               ggplot2_3.3.3               DESeq2_1.30.1               SummarizedExperiment_1.20.0 Biobase_2.50.0             
# [11] MatrixGenerics_1.2.1        matrixStats_0.58.0          GenomicRanges_1.42.0        GenomeInfoDb_1.26.7         IRanges_2.24.1             
# [16] S4Vectors_0.28.1            BiocGenerics_0.36.1        
# 
# loaded via a namespace (and not attached):
#   [1] BiocFileCache_1.14.0   plyr_1.8.6             igraph_1.2.6           lazyeval_0.2.2         splines_4.0.3          BiocParallel_1.24.1   
# [7] listenv_0.8.0          scattermore_0.7        digest_0.6.27          htmltools_0.5.4        fansi_0.4.2            magrittr_2.0.1        
# [13] memoise_2.0.1          tensor_1.5             cluster_2.1.0          ROCR_1.0-11            globals_0.14.0         annotate_1.68.0       
# [19] askpass_1.1            spatstat.sparse_2.0-0  prettyunits_1.1.1      colorspace_2.0-0       blob_1.2.1             rappdirs_0.3.3        
# [25] ggrepel_0.9.1          dplyr_1.0.5            crayon_1.4.1           RCurl_1.98-1.3         jsonlite_1.7.2         genefilter_1.72.1     
# [31] spatstat.data_2.1-0    survival_3.2-7         zoo_1.8-9              glue_1.4.2             polyclip_1.10-0        gtable_0.3.0          
# [37] zlibbioc_1.36.0        XVector_0.30.0         leiden_0.3.7           DelayedArray_0.16.3    future.apply_1.7.0     abind_1.4-5           
# [43] scales_1.1.1           DBI_1.1.1              miniUI_0.1.1.1         Rcpp_1.0.7             viridisLite_0.4.0      xtable_1.8-4          
# [49] progress_1.2.2         reticulate_1.18        spatstat.core_2.0-0    bit_4.0.4              htmlwidgets_1.5.3      httr_1.4.2            
# [55] ellipsis_0.3.2         ica_1.0-2              farver_2.1.0           pkgconfig_2.0.3        XML_3.99-0.6           uwot_0.1.10           
# [61] dbplyr_2.1.1           deldir_1.0-6           locfit_1.5-9.4         utf8_1.2.1             labeling_0.4.2         tidyselect_1.1.0      
# [67] rlang_0.4.10           reshape2_1.4.4         later_1.1.0.1          AnnotationDbi_1.52.0   munsell_0.5.0          tools_4.0.3           
# [73] cachem_1.0.4           generics_0.1.0         RSQLite_2.2.7          ggridges_0.5.3         stringr_1.4.0          fastmap_1.1.0         
# [79] goftest_1.2-2          bit64_4.0.5            fitdistrplus_1.1-3     purrr_0.3.4            RANN_2.6.1             pbapply_1.4-3         
# [85] future_1.21.0          nlme_3.1-149           mime_0.10              xml2_1.3.2             compiler_4.0.3         plotly_4.9.3          
# [91] curl_4.3               png_0.1-7              spatstat.utils_2.3-1   tibble_3.1.0           geneplotter_1.68.0     stringi_1.5.3         
# [97] lattice_0.20-41        Matrix_1.2-18          vctrs_0.3.8            pillar_1.5.1           lifecycle_1.0.0        spatstat.geom_2.4-0   
# [103] lmtest_0.9-38          RcppAnnoy_0.0.19       data.table_1.14.0      cowplot_1.1.1          bitops_1.0-6           irlba_2.3.3           
# [109] httpuv_1.5.5           patchwork_1.1.1        R6_2.5.0               promises_1.2.0.1       KernSmooth_2.23-17     parallelly_1.24.0     
# [115] codetools_0.2-16       MASS_7.3-53            assertthat_0.2.1       openssl_1.4.3          withr_2.4.1            sctransform_0.3.2     
# [121] GenomeInfoDbData_1.2.4 mgcv_1.8-33            hms_1.1.0              grid_4.0.3             rpart_4.1-15           tidyr_1.1.3           
# [127] Rtsne_0.15             shiny_1.6.0           
