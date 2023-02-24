#Set seed and working directory
set.seed(1234)
setwd(dir = "/Users/rphill3/4_hour_ATAC_local/")

#load libraries
library(DiffBind)
library(ggplot2)
library(Seurat)
library(edgeR)

#Create the dba object. Here the excel sheet contains narrow peak files from MACS2 which are in the BED6+4 format.
#Thus, with the narrowpeak file we can just use bed as our peakcaller
KCl_Veh_Str_4 <- dba(sampleSheet = "4_hour_ATAC.csv",
                     minOverlap = 0) #NArrow

#create a heatmap
plot(KCl_Veh_Str_4)

#Now we count reads to create a binding matrix
#The summits option recenters the peak by adding 250 bp up and downstream
KCl_Veh_Str_4 <- dba.count(KCl_Veh_Str_4,
                           bParallel = TRUE,
                           summits = 250,
                           bUseSummarizeOverlaps = TRUE,
                           score = DBA_SCORE_TMM_READS_FULL_CPM)

#Write out the per sample counts 
Counts <- as.data.frame(dba.peakset(KCl_Veh_Str_4, bRetrieve=TRUE))
write.table(x         = Counts,
            file      = "4_Hour_ATAC_Counts.txt",
            quote     = FALSE,
            col.names = TRUE,
            row.names = FALSE,
            sep       = "\t")


#Make the contrast. KCl first so that all of regions increased by KCl are positive. 
KCl_Veh_Str_4 <- dba.contrast(KCl_Veh_Str_4,contrast = c("Treatment","KCl","Veh"))
# Computing results names...

#Now run the DESeq2 analysis 
KCl_Veh_Str_4 <- dba.analyze(KCl_Veh_Str_4,method = DBA_DESEQ2)
# Applying Blacklist/Greylists...
# No genome detected.
# Normalize DESeq2 with defaults...
# Analyzing...
# gene-wise dispersion estimates
# mean-dispersion relationship
# final dispersion estimates

#Write out the DARs
FourHour_DARs <- as.data.frame(dba.report(KCl_Veh_Str_4, th = 1,bNormalized = FALSE))

#Write out the dataframe with stats for all DARs
write.table(x         = FourHour_DARs,
            file      = "Four_Hour_DARs_Rn7_all.txt",
            sep       = "\t",
            quote     = FALSE,
            row.names = FALSE,
            col.names = TRUE)

#Subset for significant
sig <- subset(FourHour_DARs,subset=(FDR <= 0.05))

#Write out the files as a bed file. 
write.table(x         = data.frame(chr = sig$seqnames,
                                   start = sig$start,
                                   end.  = sig$end),
            file      = "Four_Hour_DARs_Rn7.bed",
            sep       = "\t",
            quote     = FALSE,
            row.names = FALSE,
            col.names = FALSE)


#Ggplot Volcano plot 
p1 <- ggplot(data = FourHour_DARs,aes(x = Fold,y = -log10(FDR))) + 
  geom_point(colour = "grey") +
  geom_point(data = subset(FourHour_DARs,subset=(Fold > 0 & FDR <= 0.05)),colour = "red") +
  geom_point(data = subset(FourHour_DARs,subset=(Fold < 0 & FDR <= 0.05)),colour = "dodgerblue") +
  xlim(c(-3.5,3.5)) +
  theme_bw() +
  NoGrid() +
  ylim(c(0,8)) +
  geom_hline(yintercept = -log10(0.05),lty = 2) +
  labs(x = "log2(KCl)-log2(Veh)",
       title = paste("Four Hour 10mM KCl\n",nrow(subset(FourHour_DARs,subset=(FDR <= 0.05))),"DARs")) +
  annotate(geom = "text",x = -2,y = 8,label = nrow(subset(FourHour_DARs,subset=(FDR <= 0.05 & Fold > 0))),color = "red") +
  annotate(geom = "text",x = -2,y = 7.7,label = nrow(subset(FourHour_DARs,subset=(FDR <= 0.05 & Fold < 0))),color = "dodgerblue") +
  theme(plot.title = element_text(hjust = 0.5))

ggsave(filename = "Four_hour_KCl_ATAC_Volcanoplot_bNormalizedFALSE.pdf",plot = p1,height = 8,width = 8)


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
#   [1] edgeR_3.40.2                limma_3.54.1                SeuratObject_4.1.3          Seurat_4.3.0                ggplot2_3.4.0              
# [6] DiffBind_3.8.4              SummarizedExperiment_1.28.0 Biobase_2.58.0              MatrixGenerics_1.10.0       matrixStats_0.63.0         
# [11] GenomicRanges_1.50.2        GenomeInfoDb_1.34.9         IRanges_2.32.0              S4Vectors_0.36.1            BiocGenerics_0.44.0        
# 
# loaded via a namespace (and not attached):
#   [1] utf8_1.2.3               spatstat.explore_3.0-6   reticulate_1.28          tidyselect_1.2.0         RSQLite_2.2.20           AnnotationDbi_1.60.0    
# [7] htmlwidgets_1.6.1        grid_4.2.2               BiocParallel_1.32.5      Rtsne_0.16               munsell_0.5.0            codetools_0.2-18        
# [13] ragg_1.2.5               ica_1.0-3                interp_1.1-3             systemPipeR_2.4.0        future_1.31.0            miniUI_0.1.1.1          
# [19] withr_2.5.0              spatstat.random_3.1-3    colorspace_2.1-0         progressr_0.13.0         rstudioapi_0.14          ROCR_1.0-11             
# [25] tensor_1.5               listenv_0.9.0            labeling_0.4.2           bbmle_1.0.25             GenomeInfoDbData_1.2.9   mixsqp_0.3-48           
# [31] hwriter_1.3.2.1          polyclip_1.10-4          bit64_4.0.5              farver_2.1.1             coda_0.19-4              parallelly_1.34.0       
# [37] vctrs_0.5.2              generics_0.1.3           R6_2.5.1                 apeglm_1.20.0            invgamma_1.1             locfit_1.5-9.7          
# [43] bitops_1.0-7             spatstat.utils_3.0-1     cachem_1.0.6             DelayedArray_0.24.0      promises_1.2.0.1         BiocIO_1.8.0            
# [49] scales_1.2.1             gtable_0.3.1             globals_0.16.2           goftest_1.2-3            rlang_1.0.6              systemfonts_1.0.4       
# [55] splines_4.2.2            rtracklayer_1.58.0       lazyeval_0.2.2           spatstat.geom_3.0-6      yaml_2.3.7               reshape2_1.4.4          
# [61] abind_1.4-5              httpuv_1.6.8             tools_4.2.2              ellipsis_0.3.2           gplots_3.1.3             RColorBrewer_1.1-3      
# [67] ggridges_0.5.4           Rcpp_1.0.10              plyr_1.8.8               zlibbioc_1.44.0          purrr_1.0.1              RCurl_1.98-1.10         
# [73] deldir_1.0-6             pbapply_1.7-0            ashr_2.2-54              cowplot_1.1.1            zoo_1.8-11               ggrepel_0.9.3           
# [79] cluster_2.1.4            magrittr_2.0.3           data.table_1.14.6        scattermore_0.8          lmtest_0.9-40            RANN_2.6.1              
# [85] truncnorm_1.0-8          mvtnorm_1.1-3            SQUAREM_2021.1           amap_0.8-19              fitdistrplus_1.1-8       patchwork_1.1.2         
# [91] mime_0.12                xtable_1.8-4             XML_3.99-0.13            emdbook_1.3.12           jpeg_0.1-10              gridExtra_2.3           
# [97] compiler_4.2.2           bdsmatrix_1.3-6          tibble_3.1.8             KernSmooth_2.23-20       crayon_1.5.2             htmltools_0.5.4         
# [103] later_1.3.0              tidyr_1.3.0              geneplotter_1.76.0       DBI_1.1.3                MASS_7.3-58.1            ShortRead_1.56.1        
# [109] Matrix_1.5-1             cli_3.6.0                parallel_4.2.2           igraph_1.3.5             pkgconfig_2.0.3          GenomicAlignments_1.34.0
# [115] numDeriv_2016.8-1.1      sp_1.6-0                 plotly_4.10.1            spatstat.sparse_3.0-0    annotate_1.76.0          XVector_0.38.0          
# [121] stringr_1.5.0            digest_0.6.31            sctransform_0.3.5        RcppAnnoy_0.0.20         spatstat.data_3.0-0      Biostrings_2.66.0       
# [127] leiden_0.4.3             uwot_0.1.14              restfulr_0.0.15          GreyListChIP_1.30.0      shiny_1.7.4              Rsamtools_2.14.0        
# [133] gtools_3.9.4             rjson_0.2.21             lifecycle_1.0.3          nlme_3.1-160             jsonlite_1.8.4           viridisLite_0.4.1       
# [139] BSgenome_1.66.2          fansi_1.0.4              pillar_1.8.1             lattice_0.20-45          KEGGREST_1.38.0          fastmap_1.1.0           
# [145] httr_1.4.4               survival_3.4-0           glue_1.6.2               png_0.1-8                bit_4.0.5                stringi_1.7.12          
# [151] blob_1.2.3               textshaping_0.3.6        DESeq2_1.38.3            latticeExtra_0.6-30      caTools_1.18.2           memoise_2.0.1           
# [157] dplyr_1.1.0              irlba_2.3.5.1            future.apply_1.10.0     
