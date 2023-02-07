#Script is based on NBI DESeq2 workflow by Lara Ianov which is based off DESeq2 vignette
#Set seed and working directory
set.seed(1234)
setwd(dir = "/Volumes/Day-Lab$/RobertPhillips/Bioinformatics/Phillips_etal_2023/1_hour_RNA/")
#load libraries
library(DESeq2)
library(ggplot2)
library(gridExtra)
library(Seurat)
library(biomaRt)
library(pheatmap) 

##################################
######1 hour RNA-seq dataset######
##################################

#Read in the counts dataframe 
raw_counts <- read.delim(file = "fixed_raw_counts.txt")
#The Ensembl gene IDs are currently in a column. To make this count table useful for DESeq2 we are going to make the ensembl gene IDs the 
# rownames 
rownames(raw_counts) <- raw_counts$GeneID
#Table contains three columns regarding positional information for the genes. We are going to remove these three columns and the geneID column 
raw_counts <- raw_counts[-c(1:4)]

#Now we need to remove the _S## suffix that is on each of the column names for the samples 
for(i in 1:ncol(raw_counts)){
  names(raw_counts)[i] <- strsplit(x = names(raw_counts)[i],split = "_")[[1]][1]
}

#Now we need to split the samples up by brain region. 
Str <- raw_counts[,grep(pattern = "S",names(raw_counts))]

# make a Coldata table
colData <- data.frame(Sample    = names(Str),
                      Treatment = c(rep("KCl",4),rep("Veh",4)))

# Perfrom DDS
dds <- DESeqDataSetFromMatrix(countData = Str,
                              colData   =  colData,
                              design    = ~Treatment)

############## Exploring the data by visualization - GOOD QA-QC ################

# PCA plot
rld <- rlog(dds, blind=FALSE)
pca_1 <- plotPCA(rld, intgroup=c("Treatment")) #to plot a simple PCA
#Outlier present in the PCA. save for before and after. 
ggsave(filename = "Plots/PCA_With_Outlier.pdf",plot = pca_1)


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

#The QC measures show that K3S is an obvious outlier. Thus, we will remove it
Str     <- Str[,-2]
colData <- colData[-2,]

#Now rerun the DESeqDataSetFromMatrix command 
dds <- DESeqDataSetFromMatrix(countData = Str,
                              colData   =  colData,
                              design    = ~Treatment)


# PCA plot
rld <- rlog(dds, blind=FALSE)
pca_data <- plotPCA(rld, intgroup=c("Treatment"),returnData = TRUE) #to plot a simple PCA
pca_2 <- ggplot(data = pca_data,aes(x = PC1, y = PC2,color = Treatment)) +
  geom_point(size = 3) +
  scale_color_manual(values = c("KCl" = "orange",
                                "Veh" = "gray56")) +
  theme_bw() +
  NoGrid() +
  labs(x     = "PC1: 92% variance",
       y     = "PC2: 5% variance",
       title = "1 hour 10mM KCl\nRNA-Seq") +
  theme(plot.title = element_text(hjust = 0.5))
ggsave(filename = "Plots/PCA_OutlierRemoved.pdf",plot = pca_2)


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

#Perform differential expression analysis via DESeq command
# DESeq function
dds <- DESeq(dds)

# Results function. Prints summary of differentially expressed genes in Group_1 vs Group_2 comparison.
res <- results(dds, contrast=c("Treatment", 
                               "KCl", 
                               "Veh"), 
               alpha         = 0.05, 
               pAdjustMethod = "BH")

summary(res)
# out of 21330 with nonzero total read count
# adjusted p-value < 0.05
# LFC > 0 (up)       : 536, 2.5%
# LFC < 0 (down)     : 210, 0.98%
# outliers [1]       : 2, 0.0094%
# low counts [2]     : 8082, 38%
# (mean count < 14)
# [1] see 'cooksCutoff' argument of ?results
# [2] see 'independentFiltering' argument of ?results

sum(res$padj < 0.05, na.rm=TRUE)
# [1] 746

# Order results by padj
resOrdered <- res[order(res$padj),]

#Remove NAs from adjusted p values 
resFiltered <- as.data.frame(resOrdered[!is.na(resOrdered$padj), ])

#Make the rownames a column that will be useful for later merging 
resFiltered$ensembl_gene_id <- row.names(resFiltered)

#Subset for genes with baseMean >= 50
resFiltered <- subset(resFiltered,subset=(baseMean >= 50))

#Filter for genes with a adj.p.val < 0.05 and |log2FC| > 0.5
resFiltered_sig <- subset(resFiltered, padj<.05 & log2FoldChange > 0.5| padj < 0.05 & log2FoldChange < (-0.5))

#Make the rownames a column that will be useful for later merging. 
resFiltered_sig$ensembl_gene_id <- row.names(resFiltered_sig)

#For all genes
# Specify Ensembl anotation version. Use version which matches gtf file version used in primary pipeline.
ensembl <- useMart("ensembl",
                   dataset = "rnorvegicus_gene_ensembl") #Should pull the latest ensembl (This analysis was completed on 02/07/2023)

#Pull wanted attributes
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

#Merge with filtered dataframe
resFiltered <- merge(x  = as.data.frame(resFiltered),
                     y  = genes,
                     by = "ensembl_gene_id")

#Write out the tables
write.table(x         = resFiltered,
            file      = "Gene_Lists/OneHour_AllGenes.txt",
            sep       = "\t",
            col.names = TRUE,
            row.names = FALSE,
            quote     = FALSE)

write.table(x         = resFiltered_sig,
            file      = "Gene_Lists/OneHour_DEGs.txt",
            sep       = "\t",
            col.names = TRUE,
            row.names = FALSE,
            quote     = FALSE)

write.table(x         = resFiltered_sig$ensembl_gene_id,
            file      = "Gene_Lists/OneHour_DEGs_ensemblID.txt",
            sep       = "\t",
            col.names = TRUE,
            row.names = FALSE,
            quote     = FALSE)


# Volcano plot
One_hour_vp <- ggplot(data = as.data.frame(resFiltered),aes(x = log2FoldChange,y=abs(stat))) + 
  geom_point(colour = "grey",size = 3) +
  geom_point(data = subset(as.data.frame(resFiltered), padj<.05 & log2FoldChange>0.5),colour = "red",size = 3) + 
  geom_point(data = subset(as.data.frame(resFiltered), padj<.05 & log2FoldChange<(-0.5)),colour = "dodger blue",size = 3) + 
  geom_point(data = subset(resFiltered,subset=(external_gene_name == "Fos" | 
                                                 external_gene_name == "Fosb" |
                                                 external_gene_name == "Fosl1" |
                                                 external_gene_name == "Fosl2" |
                                                 external_gene_name == "Junb" |
                                                 external_gene_name == "Npas4" |
                                                 external_gene_name == "Ccn1" |
                                                 external_gene_name == "Gadd45g" |
                                                 external_gene_name == "Nr4a1" |
                                                 external_gene_name == "Nr4a2" |
                                                 external_gene_name == "Jun")),
             color = "black",
             shape= 1,
             size = 3) +
  xlim(c(-4,4)) +
  geom_vline(xintercept = 0.5,lty = 2,colour = "black") + 
  geom_vline(xintercept = -0.5,lty = 2,colour = "black") +
  geom_text(data = subset(resFiltered,subset=(external_gene_name == "Fos" | 
                                                external_gene_name == "Fosb" |
                                                external_gene_name == "Fosl1" |
                                                external_gene_name == "Fosl2" |
                                                external_gene_name == "Junb" |
                                                external_gene_name == "Npas4" |
                                                external_gene_name == "Ccn1" |
                                                external_gene_name == "Gadd45g" |
                                                external_gene_name == "Nr4a1" |
                                                external_gene_name == "Nr4a2" |
                                                external_gene_name == "Jun")),
            aes(label = external_gene_name),
            hjust=0,
            vjust=0) +
  ggtitle("1 Hour 10mM KCl \n207 DEGs") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5)) +
  NoGrid() +
  annotate(geom = "text",x = -3.5,y = 30,label = nrow(subset(resFiltered_sig,subset=(log2FoldChange > 0))),color = "red") +
  annotate(geom = "text",x = -3.5,y = 29,label = nrow(subset(resFiltered_sig,subset=(log2FoldChange < 0))),color = "dodger blue") 

ggsave(plot     = One_hour_vp,
       filename = "Plots/1Hour_DEGs_VolcanoPlot.pdf",
       height   = 8,
       width    = 8)

sessionInfo()
# R version 4.2.2 (2022-10-31)
# Platform: aarch64-apple-darwin20 (64-bit)
# Running under: macOS Ventura 13.2
# 
# Matrix products: default
# LAPACK: /Library/Frameworks/R.framework/Versions/4.2-arm64/Resources/lib/libRlapack.dylib
# 
# locale:
#   [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
# 
# attached base packages:
#   [1] stats4    stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
#   [1] RColorBrewer_1.1-3          pheatmap_1.0.12             biomaRt_2.54.0              SeuratObject_4.1.3          Seurat_4.3.0               
# [6] gridExtra_2.3               ggplot2_3.4.0               DESeq2_1.38.3               SummarizedExperiment_1.28.0 Biobase_2.58.0             
# [11] MatrixGenerics_1.10.0       matrixStats_0.63.0          GenomicRanges_1.50.2        GenomeInfoDb_1.34.9         IRanges_2.32.0             
# [16] S4Vectors_0.36.1            BiocGenerics_0.44.0        
# 
# loaded via a namespace (and not attached):
#   [1] BiocFileCache_2.6.0    plyr_1.8.8             igraph_1.3.5           lazyeval_0.2.2         sp_1.6-0               splines_4.2.2         
# [7] BiocParallel_1.32.5    listenv_0.9.0          scattermore_0.8        digest_0.6.31          htmltools_0.5.4        fansi_1.0.4           
# [13] magrittr_2.0.3         memoise_2.0.1          tensor_1.5             cluster_2.1.4          ROCR_1.0-11            globals_0.16.2        
# [19] Biostrings_2.66.0      annotate_1.76.0        spatstat.sparse_3.0-0  prettyunits_1.1.1      colorspace_2.1-0       rappdirs_0.3.3        
# [25] blob_1.2.3             ggrepel_0.9.3          dplyr_1.1.0            crayon_1.5.2           RCurl_1.98-1.10        jsonlite_1.8.4        
# [31] progressr_0.13.0       spatstat.data_3.0-0    survival_3.4-0         zoo_1.8-11             glue_1.6.2             polyclip_1.10-4       
# [37] gtable_0.3.1           zlibbioc_1.44.0        XVector_0.38.0         leiden_0.4.3           DelayedArray_0.24.0    future.apply_1.10.0   
# [43] abind_1.4-5            scales_1.2.1           DBI_1.1.3              spatstat.random_3.1-3  miniUI_0.1.1.1         Rcpp_1.0.10           
# [49] viridisLite_0.4.1      xtable_1.8-4           progress_1.2.2         reticulate_1.28        bit_4.0.5              htmlwidgets_1.6.1     
# [55] httr_1.4.4             ellipsis_0.3.2         ica_1.0-3              farver_2.1.1           pkgconfig_2.0.3        XML_3.99-0.13         
# [61] uwot_0.1.14            dbplyr_2.3.0           deldir_1.0-6           locfit_1.5-9.7         utf8_1.2.3             labeling_0.4.2        
# [67] tidyselect_1.2.0       rlang_1.0.6            reshape2_1.4.4         later_1.3.0            AnnotationDbi_1.60.0   munsell_0.5.0         
# [73] tools_4.2.2            cachem_1.0.6           cli_3.6.0              generics_0.1.3         RSQLite_2.2.20         ggridges_0.5.4        
# [79] stringr_1.5.0          fastmap_1.1.0          goftest_1.2-3          bit64_4.0.5            fitdistrplus_1.1-8     purrr_1.0.1           
# [85] RANN_2.6.1             KEGGREST_1.38.0        pbapply_1.7-0          future_1.31.0          nlme_3.1-160           mime_0.12             
# [91] xml2_1.3.3             compiler_4.2.2         filelock_1.0.2         curl_5.0.0             plotly_4.10.1          png_0.1-8             
# [97] spatstat.utils_3.0-1   tibble_3.1.8           geneplotter_1.76.0     stringi_1.7.12         lattice_0.20-45        Matrix_1.5-1          
# [103] vctrs_0.5.2            pillar_1.8.1           lifecycle_1.0.3        spatstat.geom_3.0-6    lmtest_0.9-40          RcppAnnoy_0.0.20      
# [109] data.table_1.14.6      cowplot_1.1.1          bitops_1.0-7           irlba_2.3.5.1          httpuv_1.6.8           patchwork_1.1.2       
# [115] R6_2.5.1               promises_1.2.0.1       KernSmooth_2.23-20     parallelly_1.34.0      codetools_0.2-18       assertthat_0.2.1      
# [121] MASS_7.3-58.1          withr_2.5.0            sctransform_0.3.5      GenomeInfoDbData_1.2.9 parallel_4.2.2         hms_1.1.2             
# [127] grid_4.2.2             tidyr_1.3.0            Rtsne_0.16             spatstat.explore_3.0-6 shiny_1.7.4      



