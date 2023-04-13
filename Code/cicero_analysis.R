## Co-Acccesibility Analysis ##
set.seed(1234)
#Load libraries
library(Signac)
library(Seurat)
library(SeuratWrappers)
library(ggplot2)
library(patchwork)
library(monocle3)
library(cicero)

#Load ATAC object
ATAC <- readRDS(file='/Users/ethanwan/Desktop/coaccessibility_rn7/ATAC_filtered.rds')

#Create cicero object
atac.cds <- as.cell_data_set(x = ATAC)
atac.cicero <- make_cicero_cds(atac.cds, reduced_coordinates = reducedDims(atac.cds)$UMAP)

#Perform Cicero analysis
#Load in chromosome lengths
genome.df <- readRDS('/Users/ethanwan/Desktop/coaccessibility_rn7/genomedf.rds')

#Run cicero - takes a long time
connections <- run_cicero(atac.cicero, genomic_coords = genome.df, sample_num = 200)
#Find cis-co-accessible networks
cisco <- generate_ccans(connections)

#Add Seurat object links
links <- ConnectionsToLinks(conns = connections, ccans = cisco)
Links(ATAC) <- links
save(connections, cisco, links, file='/Users/ethanwan/Desktop/cicero_gear.RData')

#Pdyn coverage plot
Coverage <- CoveragePlot(
  object = ATAC,
  region = "3-116907322-116908266",
  annotation = TRUE,
  extend.upstream = 10000,
  extend.downstream = 60000,
  links = T
) 

## D1 Cicero analysis ## 
#Genome dataframe only contains chromosome lengths, so we can use the same one as earlier
#Same procedure as before
D1 <- subset(ATAC, predicted.id=='Drd1-MSN')
D1.cds <- as.cell_data_set(x = D1)
D1.cicero <- make_cicero_cds(D1.cds, reduced_coordinates = reducedDims(D1.cds)$UMAP)
connections <- run_cicero(D1.cicero, genomic_coords = genome.df, sample_num = 200)
cisco <- generate_ccans(connections)
links <- ConnectionsToLinks(conns = connections, ccans = cisco)
Links(D1) <- links
save(connections, cisco, links, file = '/Users/ethanwan/Desktop/D1_cicero_gear.RData')

#Plots
Coverage <- CoveragePlot(
  object = D1,
  region = "3-116907322-116908266",
  annotation = TRUE,
  extend.upstream = 10000,
  extend.downstream = 60000,
  links = F
) 
Links <- LinkPlot(object = D1,
                  region = "3-116907322-116908266", 
                  extend.upstream = 10000,
                  extend.downstream = 60000) + scale_color_gradientn(colors=c('lightgrey', 'red'))
CombineTracks(plotlist = list(Coverage, Links),
              heights = c(10,3))

#sessionInfo()
# R version 4.2.2 (2022-10-31)
# Platform: x86_64-apple-darwin17.0 (64-bit)
# Running under: macOS Ventura 13.0.1
# 
# Matrix products: default
# LAPACK: /Library/Frameworks/R.framework/Versions/4.2/Resources/lib/libRlapack.dylib
# 
# locale:
#   [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
# 
# attached base packages:
#   [1] grid      stats4    stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
#   [1] cicero_1.3.9                Gviz_1.42.1                 monocle3_1.3.1             
# [4] SingleCellExperiment_1.20.0 SummarizedExperiment_1.28.0 GenomicRanges_1.50.2       
# [7] GenomeInfoDb_1.34.9         IRanges_2.32.0              S4Vectors_0.36.1           
# [10] MatrixGenerics_1.10.0       matrixStats_0.63.0          Biobase_2.58.0             
# [13] BiocGenerics_0.44.0         patchwork_1.1.2             ggplot2_3.4.0              
# [16] SeuratWrappers_0.3.1        SeuratObject_4.1.3          Seurat_4.3.0               
# [19] Signac_1.9.0               
# 
# loaded via a namespace (and not attached):
#   [1] rappdirs_0.3.3           rtracklayer_1.58.0       scattermore_0.8          R.methodsS3_1.8.2       
# [5] tidyr_1.3.0              bit64_4.0.5              knitr_1.41               irlba_2.3.5.1           
# [9] DelayedArray_0.24.0      R.utils_2.12.2           data.table_1.14.6        rpart_4.1.19            
# [13] KEGGREST_1.38.0          RCurl_1.98-1.9           AnnotationFilter_1.22.0  generics_0.1.3          
# [17] GenomicFeatures_1.50.4   terra_1.6-47             cowplot_1.1.1            RSQLite_2.3.0           
# [21] RANN_2.6.1               VGAM_1.1-7               future_1.30.0            bit_4.0.5               
# [25] spatstat.data_3.0-0      xml2_1.3.3               httpuv_1.6.8             assertthat_0.2.1        
# [29] xfun_0.36                hms_1.1.2                evaluate_0.19            promises_1.2.0.1        
# [33] fansi_1.0.3              restfulr_0.0.15          progress_1.2.2           dbplyr_2.3.1            
# [37] igraph_1.3.5             DBI_1.1.3                htmlwidgets_1.6.1        spatstat.geom_3.0-3     
# [41] purrr_1.0.1              ellipsis_0.3.2           dplyr_1.1.0              backports_1.4.1         
# [45] biomaRt_2.54.0           deldir_1.0-6             vctrs_0.5.2              remotes_2.4.2           
# [49] ensembldb_2.22.0         ROCR_1.0-11              abind_1.4-5              cachem_1.0.6            
# [53] withr_2.5.0              ggforce_0.4.1            BSgenome_1.66.3          progressr_0.13.0        
# [57] checkmate_2.1.0          sctransform_0.3.5        GenomicAlignments_1.34.0 prettyunits_1.1.1       
# [61] goftest_1.2-3            cluster_2.1.4            lazyeval_0.2.2           crayon_1.5.2            
# [65] spatstat.explore_3.0-5   pkgconfig_2.0.3          labeling_0.4.2           tweenr_2.0.2            
# [69] nlme_3.1-160             ProtGenerics_1.30.0      nnet_7.3-18              rlang_1.0.6             
# [73] globals_0.16.2           lifecycle_1.0.3          miniUI_0.1.1.1           filelock_1.0.2          
# [77] BiocFileCache_2.6.1      rsvd_1.0.5               dichromat_2.0-0.1        polyclip_1.10-4         
# [81] lmtest_0.9-40            Matrix_1.5-1             boot_1.3-28              zoo_1.8-11              
# [85] base64enc_0.1-3          ggridges_0.5.4           png_0.1-8                viridisLite_0.4.1       
# [89] rjson_0.2.21             bitops_1.0-7             R.oo_1.25.0              KernSmooth_2.23-20      
# [93] Biostrings_2.66.0        blob_1.2.3               stringr_1.5.0            parallelly_1.33.0       
# [97] spatstat.random_3.0-1    jpeg_0.1-10              scales_1.2.1             memoise_2.0.1           
# [101] magrittr_2.0.3           plyr_1.8.8               ica_1.0-3                zlibbioc_1.44.0         
# [105] compiler_4.2.2           BiocIO_1.8.0             RColorBrewer_1.1-3       lme4_1.1-31             
# [109] fitdistrplus_1.1-8       Rsamtools_2.14.0         cli_3.6.0                XVector_0.38.0          
# [113] listenv_0.9.0            pbapply_1.7-0            htmlTable_2.4.1          Formula_1.2-5           
# [117] MASS_7.3-58.1            tidyselect_1.2.0         stringi_1.7.12           yaml_2.3.6              
# [121] latticeExtra_0.6-30      ggrepel_0.9.2            VariantAnnotation_1.44.1 fastmatch_1.1-3         
# [125] tools_4.2.2              future.apply_1.10.0      parallel_4.2.2           rstudioapi_0.14         
# [129] foreign_0.8-83           gridExtra_2.3            farver_2.1.1             Rtsne_0.16              
# [133] digest_0.6.31            BiocManager_1.30.19      FNN_1.1.3.1              shiny_1.7.4             
# [137] Rcpp_1.0.9               later_1.3.0              RcppAnnoy_0.0.20         httr_1.4.4              
# [141] AnnotationDbi_1.60.0     biovizBase_1.46.0        colorspace_2.0-3         XML_3.99-0.13           
# [145] tensor_1.5               reticulate_1.27          splines_4.2.2            uwot_0.1.14             
# [149] RcppRoll_0.3.0           spatstat.utils_3.0-1     sp_1.5-1                 plotly_4.10.1           
# [153] xtable_1.8-4             jsonlite_1.8.4           nloptr_2.0.3             glasso_1.11             
# [157] R6_2.5.1                 Hmisc_4.8-0              pillar_1.8.1             htmltools_0.5.4         
# [161] mime_0.12                glue_1.6.2               fastmap_1.1.0            minqa_1.2.5             
# [165] BiocParallel_1.32.5      codetools_0.2-18         utf8_1.2.2               lattice_0.20-45         
# [169] spatstat.sparse_3.0-0    tibble_3.1.8             curl_5.0.0               leiden_0.4.3            
# [173] interp_1.1-3             survival_3.4-0           rmarkdown_2.19           munsell_0.5.0           
# [177] GenomeInfoDbData_1.2.9   reshape2_1.4.4           gtable_0.3.1         
