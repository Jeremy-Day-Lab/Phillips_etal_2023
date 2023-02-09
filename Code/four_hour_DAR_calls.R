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

#Now run the DESeq2 analysis 
KCl_Veh_Str_4 <- dba.analyze(KCl_Veh_Str_4)
#Here is the DESeq2 output in the console. 
# Applying Blacklist/Greylists...
# No genome detected.
# Normalize DESeq2 with defaults...
# Analyzing...
# gene-wise dispersion estimates
# mean-dispersion relationship
# final dispersion estimates

#Write out the DARs
FourHour_DARs <- as.data.frame(dba.report(KCl_Veh_Str_4, th = 1))

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
  xlim(c(-2.5,2.5)) +
  theme_bw() +
  NoGrid() +
  ylim(c(0,8)) +
  geom_hline(yintercept = -log10(0.05),lty = 2) +
  labs(x = "log2(KCl)-log2(Veh)",
       title = paste("Four Hour 10mM KCl\n",nrow(subset(FourHour_DARs,subset=(FDR <= 0.05))),"DARs")) +
  annotate(geom = "text",x = -2,y = 8,label = nrow(subset(FourHour_DARs,subset=(FDR <= 0.05 & Fold > 0))),color = "red") +
  annotate(geom = "text",x = -2,y = 7.7,label = nrow(subset(FourHour_DARs,subset=(FDR <= 0.05 & Fold < 0))),color = "dodgerblue") +
  theme(plot.title = element_text(hjust = 0.5))

ggsave(filename = "Four_hour_KCl_ATAC_Volcanoplot.pdf",plot = p1,height = 8,width = 8)


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
#   [1] edgeR_3.40.2                limma_3.54.1                SeuratObject_4.1.3          Seurat_4.3.0               
# [5] ggplot2_3.4.0               DiffBind_3.8.4              SummarizedExperiment_1.28.0 Biobase_2.58.0             
# [9] MatrixGenerics_1.10.0       matrixStats_0.63.0          GenomicRanges_1.50.2        GenomeInfoDb_1.34.9        
# [13] IRanges_2.32.0              S4Vectors_0.36.1            BiocGenerics_0.44.0        
# 
# loaded via a namespace (and not attached):
#   [1] plyr_1.8.8               igraph_1.3.5             lazyeval_0.2.2           sp_1.6-0                
# [5] splines_4.2.2            BiocParallel_1.32.5      listenv_0.9.0            scattermore_0.8         
# [9] amap_0.8-19              digest_0.6.31            invgamma_1.1             htmltools_0.5.4         
# [13] SQUAREM_2021.1           fansi_1.0.4              memoise_2.0.1            magrittr_2.0.3          
# [17] BSgenome_1.66.2          tensor_1.5               cluster_2.1.4            ROCR_1.0-11             
# [21] annotate_1.76.0          globals_0.16.2           Biostrings_2.66.0        systemPipeR_2.4.0       
# [25] spatstat.sparse_3.0-0    bdsmatrix_1.3-6          jpeg_0.1-10              colorspace_2.1-0        
# [29] blob_1.2.3               apeglm_1.20.0            ggrepel_0.9.3            dplyr_1.1.0             
# [33] crayon_1.5.2             RCurl_1.98-1.10          jsonlite_1.8.4           progressr_0.13.0        
# [37] spatstat.data_3.0-0      survival_3.4-0           zoo_1.8-11               glue_1.6.2              
# [41] polyclip_1.10-4          gtable_0.3.1             zlibbioc_1.44.0          XVector_0.38.0          
# [45] leiden_0.4.3             DelayedArray_0.24.0      future.apply_1.10.0      abind_1.4-5             
# [49] scales_1.2.1             mvtnorm_1.1-3            DBI_1.1.3                spatstat.random_3.1-3   
# [53] miniUI_0.1.1.1           Rcpp_1.0.10              viridisLite_0.4.1        xtable_1.8-4            
# [57] emdbook_1.3.12           reticulate_1.28          bit_4.0.5                truncnorm_1.0-8         
# [61] htmlwidgets_1.6.1        httr_1.4.4               gplots_3.1.3             RColorBrewer_1.1-3      
# [65] ellipsis_0.3.2           ica_1.0-3                farver_2.1.1             pkgconfig_2.0.3         
# [69] XML_3.99-0.13            uwot_0.1.14              deldir_1.0-6             locfit_1.5-9.7          
# [73] utf8_1.2.3               labeling_0.4.2           AnnotationDbi_1.60.0     tidyselect_1.2.0        
# [77] rlang_1.0.6              reshape2_1.4.4           later_1.3.0              cachem_1.0.6            
# [81] munsell_0.5.0            tools_4.2.2              cli_3.6.0                RSQLite_2.2.20          
# [85] generics_0.1.3           ggridges_0.5.4           stringr_1.5.0            fastmap_1.1.0           
# [89] goftest_1.2-3            yaml_2.3.7               bit64_4.0.5              fitdistrplus_1.1-8      
# [93] caTools_1.18.2           purrr_1.0.1              RANN_2.6.1               KEGGREST_1.38.0         
# [97] nlme_3.1-160             pbapply_1.7-0            future_1.31.0            mime_0.12               
# [101] compiler_4.2.2           plotly_4.10.1            png_0.1-8                spatstat.utils_3.0-1    
# [105] geneplotter_1.76.0       tibble_3.1.8             stringi_1.7.12           lattice_0.20-45         
# [109] Matrix_1.5-1             vctrs_0.5.2              pillar_1.8.1             lifecycle_1.0.3         
# [113] spatstat.geom_3.0-6      lmtest_0.9-40            RcppAnnoy_0.0.20         data.table_1.14.6       
# [117] cowplot_1.1.1            bitops_1.0-7             irlba_2.3.5.1            httpuv_1.6.8            
# [121] patchwork_1.1.2          rtracklayer_1.58.0       R6_2.5.1                 BiocIO_1.8.0            
# [125] latticeExtra_0.6-30      hwriter_1.3.2.1          promises_1.2.0.1         ShortRead_1.56.1        
# [129] KernSmooth_2.23-20       gridExtra_2.3            parallelly_1.34.0        codetools_0.2-18        
# [133] MASS_7.3-58.1            gtools_3.9.4             DESeq2_1.38.3            rjson_0.2.21            
# [137] withr_2.5.0              GenomicAlignments_1.34.0 sctransform_0.3.5        Rsamtools_2.14.0        
# [141] GenomeInfoDbData_1.2.9   parallel_4.2.2           grid_4.2.2               tidyr_1.3.0             
# [145] coda_0.19-4              GreyListChIP_1.30.0      ashr_2.2-54              Rtsne_0.16              
# [149] mixsqp_0.3-48            spatstat.explore_3.0-6   bbmle_1.0.25             numDeriv_2016.8-1.1     
# [153] shiny_1.7.4              interp_1.1-3             restfulr_0.0.15    
