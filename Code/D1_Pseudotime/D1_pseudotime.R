#Calculating pseudotime for Drd1 subcluster
set.seed(1234)
#Load libraries
library(Seurat)
library(SeuratData)
library(monocle3)
library(ggplot2)
library(patchwork)
library(magrittr)

#Load Acute cocaine rn7 object
All_Groups_log <- readRDS('/Users/ethanwan/Desktop/Acute_Decontaminated.rds')
#Validate UMAP
DimPlot(object = All_Groups_log,reduction = "umap")

#Subcluster Drd1 neurons
#Subset for only D1s
D1 <- subset(All_Groups_log, subset = CellType == 'Drd1-MSN')

#Run standard dimensionality reduction and clustering workflow from Seurat4
D1 <- ScaleData(D1,verbose = FALSE)
D1 <- RunPCA(D1,npcs = 17,verbose = FALSE)
D1 <- RunUMAP(D1, reduction = "pca", dims = 1:17)
D1 <- FindNeighbors(D1, reduction = "pca", dims = 1:17)
D1 <- FindClusters(D1, resolution = 0.2)

#Visualize D1s on UMAP
DimPlot(D1,reduction = 'umap',label = TRUE) + NoLegend()

#Find genes of interest 
FeaturePlot(D1, features = c('Ebf1','Htr4', 'Pdyn', 'Fos'))

#To run pseudotime, the Seurat4 object must be converted into a Monocle3 cds object
#To make cds object manually, we need to extract annotations, counts metadata, and barcodes
#SeuratWrappers contains documentation on a shorter way to do this, but it does not give as biologically believable results
#For example, when using SeuratWrappers, some Htr4 positive cells are assigned high pseudotimes whereas when creating the Monocle3 object manually, this does not occur

#Extract gene annotations
gene_annotation <- as.data.frame(rownames(D1@reductions[["pca"]]@feature.loadings),
                                 row.names = rownames(D1@reductions[["pca"]]@feature.loadings))
#Rename column
colnames(gene_annotation) <- "gene_short_name"

#Extract barcodes
cell_metadata <- as.data.frame(D1@assays[["RNA"]]@counts@Dimnames[[2]],
                               row.names = D1@assays[["RNA"]]@counts@Dimnames[[2]])
#Rename column
colnames(cell_metadata) <- "barcode"

#Extract counts data
New_matrix <- D1@assays[["RNA"]]@counts
New_matrix <- New_matrix[rownames(D1@reductions[["pca"]]@feature.loadings), ]
expression_matrix <- New_matrix

#Create a cds object from Seurat using annotations, counts, and barcodes
cds_from_seurat <- new_cell_data_set(expression_matrix,
                                     cell_metadata = cell_metadata,
                                     gene_metadata = gene_annotation)

#Monocle3 requires that partitions are calculated in conjunction with the clustering information
#However, we will not be using the partition data in this pseudotime analysis, so we are creating a filler partition element in the cds object to satisfy this requirement

#Create a dataframe that has the same number of rows as cds object - all values are "1" to stand in as partition data
recreate.partition <- c(rep(1, length(cds_from_seurat@colData@rownames)))
#Add barcodes as names 
names(recreate.partition) <- cds_from_seurat@colData@rownames
recreate.partition <- as.factor(recreate.partition)
#Add dummy dataframe to cds object 
cds_from_seurat@clusters@listData[["UMAP"]][["partitions"]] <- recreate.partition

#Add clustering data and UMAP dimensions from Seurat object to cds object so the UMAPs match
list_cluster <- D1@active.ident
names(list_cluster) <- D1@assays[["RNA"]]@data@Dimnames[[2]]
cds_from_seurat@clusters@listData[["UMAP"]][["clusters"]] <- list_cluster
cds_from_seurat@clusters@listData[["UMAP"]][["louvain_res"]] <- "NA"
cds_from_seurat@int_colData@listData$reducedDims@listData[["UMAP"]] <-D1@reductions[["umap"]]@cell.embeddings
cds_from_seurat@reduce_dim_aux$gene_loadings <- D1@reductions[["pca"]]@feature.loadings

#Perform pseudotime analysis
cds_from_seurat <- learn_graph(cds_from_seurat, use_partition = F)

#Choose root cell/origin
cds_from_seurat <- order_cells(cds_from_seurat, reduction_method = 'UMAP')

#Visualize on UMAP using Monocle3, should be the same as Seurat UMAP
plot_cells(cds_from_seurat, 
           color_cells_by = 'cluster',
           label_groups_by_cluster=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           graph_label_size=4)

#Visualize pseudotime on UMAP
plot_cells(cds_from_seurat,
           color_cells_by = "pseudotime",
           graph_label_size=5,
           show_trajectory_graph = TRUE)

#Write pseudotime data out as .csv
pseudo <- pseudotime(cds_from_seurat)
write.csv(pseudo, file = 'D1_rn7_pseudotime.csv', quote = FALSE)

#Save D1 object
saveRDS(D1, file = 'D1_rn7.rds')


#SessionInfo()
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
#   [1] stats4    stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
#   [1] magrittr_2.0.3              patchwork_1.1.2             ggplot2_3.4.0              
# [4] monocle3_1.3.1              SingleCellExperiment_1.20.0 SummarizedExperiment_1.28.0
# [7] GenomicRanges_1.50.2        GenomeInfoDb_1.34.6         IRanges_2.32.0             
# [10] S4Vectors_0.36.1            MatrixGenerics_1.10.0       matrixStats_0.63.0         
# [13] Biobase_2.58.0              BiocGenerics_0.44.0         SeuratData_0.2.2           
# [16] SeuratObject_4.1.3          Seurat_4.3.0               
# 
# loaded via a namespace (and not attached):
#   [1] minqa_1.2.5            Rtsne_0.16             colorspace_2.0-3       deldir_1.0-6          
# [5] ellipsis_0.3.2         ggridges_0.5.4         XVector_0.38.0         rstudioapi_0.14       
# [9] spatstat.data_3.0-0    leiden_0.4.3           listenv_0.9.0          ggrepel_0.9.2         
# [13] fansi_1.0.3            codetools_0.2-18       splines_4.2.2          polyclip_1.10-4       
# [17] jsonlite_1.8.4         nloptr_2.0.3           ica_1.0-3              cluster_2.1.4         
# [21] png_0.1-8              uwot_0.1.14            shiny_1.7.4            sctransform_0.3.5     
# [25] spatstat.sparse_3.0-0  compiler_4.2.2         httr_1.4.4             Matrix_1.5-1          
# [29] fastmap_1.1.0          lazyeval_0.2.2         cli_3.6.0              later_1.3.0           
# [33] htmltools_0.5.4        tools_4.2.2            igraph_1.3.5           GenomeInfoDbData_1.2.9
# [37] gtable_0.3.1           glue_1.6.2             RANN_2.6.1             reshape2_1.4.4        
# [41] dplyr_1.0.10           rappdirs_0.3.3         Rcpp_1.0.9             scattermore_0.8       
# [45] vctrs_0.5.1            spatstat.explore_3.0-5 nlme_3.1-160           progressr_0.13.0      
# [49] lmtest_0.9-40          spatstat.random_3.0-1  stringr_1.5.0          globals_0.16.2        
# [53] lme4_1.1-31            mime_0.12              miniUI_0.1.1.1         lifecycle_1.0.3       
# [57] irlba_2.3.5.1          goftest_1.2-3          terra_1.6-47           future_1.30.0         
# [61] zlibbioc_1.44.0        MASS_7.3-58.1          zoo_1.8-11             scales_1.2.1          
# [65] promises_1.2.0.1       spatstat.utils_3.0-1   parallel_4.2.2         RColorBrewer_1.1-3    
# [69] reticulate_1.27        pbapply_1.7-0          gridExtra_2.3          stringi_1.7.12        
# [73] boot_1.3-28            bitops_1.0-7           rlang_1.0.6            pkgconfig_2.0.3       
# [77] lattice_0.20-45        ROCR_1.0-11            purrr_1.0.1            tensor_1.5            
# [81] htmlwidgets_1.6.1      cowplot_1.1.1          tidyselect_1.2.0       parallelly_1.33.0     
# [85] RcppAnnoy_0.0.20       plyr_1.8.8             R6_2.5.1               generics_0.1.3        
# [89] DelayedArray_0.24.0    DBI_1.1.3              withr_2.5.0            pillar_1.8.1          
# [93] fitdistrplus_1.1-8     RCurl_1.98-1.9         survival_3.4-0         abind_1.4-5           
# [97] sp_1.5-1               tibble_3.1.8           future.apply_1.10.0    crayon_1.5.2          
# [101] KernSmooth_2.23-20     utf8_1.2.2             spatstat.geom_3.0-3    plotly_4.10.1         
# [105] grid_4.2.2             data.table_1.14.6      digest_0.6.31          xtable_1.8-4          
# [109] tidyr_1.2.1            httpuv_1.6.8           munsell_0.5.0          viridisLite_0.4.1    





