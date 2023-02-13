#Reconstructing Drd1 static gene regulatory network
set.seed(1234)
#Load libraries
library(Seurat)
library(SeuratObject)
library(SeuratData)
library(ggplot2)
library(epoch)

#Load in D1 subcluster and pseudotime data
D1 <- readRDS('D1_rn7.rds')
pseudo <- read.csv('D1_rn7_pseudotime.csv')

#Add pseudotime data to Seurat object
D1$pseudotime <- pseudo$x

#Checking to see if genes of interest are present
hg <- rownames(D1)
grep("Drd1", hg)
grep("Ebf1", hg) 
grep("Htr4", hg) 
grep("Fos", hg)  
grep("Pdyn", hg)

#Epoch requires RNA assay data and metadata to reconstruct a static network
#Extract needed data
expDat <- as.matrix(GetAssayData(D1, assay = "RNA", slot = "data"))
sampTab <- D1@meta.data

#This procedure is based on the outline provided by the Cahan Lab on their github page for epoch
#https://github.com/pcahan1/epoch

#Reconstructing static network
#Start by finding dynamically expressed genes
xdyn <- findDynGenes(expDat, sampTab, group_column='seurat_clusters', pseudotime_column="pseudotime")
#Filter by a significance threshold of p < 0.05
pThresh<-0.05
dgenes<-names(xdyn$genes)[xdyn$genes<pThresh]
#Remove all NA values
dyn.genes <- dgenes[!(is.na(dgenes))]

#Read in rattus norvegicus transcription factors from AnimalTFDB 4.0 database
TF <- read.table("/Users/ethanwan/Downloads/Rattus_norvegicus_TF.txt", sep="\t", header = T)
#Filter transcription factors that are found in the database and in our object
TF.1 <- intersect(rownames(expDat), unique(TF[,2]))

#Reconstruct network using Context of Likelihood Relatedness (CLR) algorithm
grn.1 <- reconstructGRN(expDat, TF.1, dyn.genes, method = "pearson", zThresh = 0)
#Perform crossweighting to refine
grn.1 <- crossweight(grn.1, expDat, xdyn, filter_thresh = 0)

#Adding interaction type, e.g. activation or repression
grn.1 <- grn.1[order(grn.1$weighted_score, decreasing = T), ]
grn.1.t <- add_interaction_type(grn.1)

#Save GRN data
write.table(grn.1.t, file = "D1.static.network.GRNs.TF.txt", quote = F, sep ="\t", row.names = FALSE)
save(expDat, TF, xdyn, sampTab, file="/Users/ethanwan/Downloads/D1_rn7_grn.RData")
load('D1_grn.RData')

#Dynamic network extraction
#Calculating epochs
hh <- sampTab[rownames(xdyn$cells), ]
table(rownames(hh) == rownames(xdyn$cells))
#Separating the cells based on pseudotime
xdyn <- define_epochs(xdyn, expDat, method = "pseudotime", num_epochs = 6)
dyn.m <- xdyn$cells
dyn.m <- dyn.m[rownames(D1@meta.data), ]
rownames(dyn.m) == rownames(D1@meta.data)
D1@meta.data$epoch <- dyn.m$epoch

#Generating dynamic network
#Change NA values to 1, NA values mess up the epoch assignments
xdyn$genes[which(is.na(xdyn$genes))] <- 1
#Assign epochs
epoch_assignments <- assign_epochs(expDat = expDat,dynRes =  xdyn)

#Compute GRN 
dynamic_grn <- epochGRN(grn.1.t, epoch_assignments)
gene_rank <- compute_pagerank(dynamic_grn, weight_column = 'weighted_score')
dynamic_grn <- add_interaction_type(dynamic_grn)

#Plotting dynamic networks fof Pdyn and Scg2, which is a known Fos target
genes_of_interest <- plot_targets_with_top_regulators_detail(dynamic_grn, c("Pdyn","Scg2"), epoch_assignments, weight_column = 'weighted_score', numTopRegulators = 20, declutter = TRUE, order = 'epoch6..epoch6')
ggsave(genes_of_interest, file = 'epoch6_top20.pdf', width = 8, height = 6)

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
#   [1] splines   stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
#   [1] epoch_0.0.0.9000   viridis_0.6.2      viridisLite_0.4.1  gridExtra_2.3      ggnetwork_0.5.10  
# [6] DescTools_0.99.47  GENIE3_1.20.0      RColorBrewer_1.1-3 pheatmap_1.0.12    minet_3.56.0      
# [11] gam_1.22           foreach_1.5.2      loomR_0.2.0        itertools_0.1-3    iterators_1.0.14  
# [16] hdf5r_1.3.7        R6_2.5.1           igraph_1.3.5       reshape2_1.4.4     ggplot2_3.4.0     
# [21] SeuratData_0.2.2   SeuratObject_4.1.3 Seurat_4.3.0      
# 
# loaded via a namespace (and not attached):
#   [1] Rtsne_0.16             colorspace_2.0-3       deldir_1.0-6           class_7.3-20          
# [5] ellipsis_0.3.2         ggridges_0.5.4         gld_2.6.6              proxy_0.4-27          
# [9] rstudioapi_0.14        spatstat.data_3.0-0    leiden_0.4.3           listenv_0.9.0         
# [13] ggrepel_0.9.2          bit64_4.0.5            mvtnorm_1.1-3          fansi_1.0.3           
# [17] codetools_0.2-18       rootSolve_1.8.2.3      polyclip_1.10-4        jsonlite_1.8.4        
# [21] ica_1.0-3              cluster_2.1.4          png_0.1-8              uwot_0.1.14           
# [25] shiny_1.7.4            sctransform_0.3.5      spatstat.sparse_3.0-0  compiler_4.2.2        
# [29] httr_1.4.4             Matrix_1.5-1           fastmap_1.1.0          lazyeval_0.2.2        
# [33] cli_3.6.0              later_1.3.0            htmltools_0.5.4        tools_4.2.2           
# [37] lmom_2.9               gtable_0.3.1           glue_1.6.2             RANN_2.6.1            
# [41] dplyr_1.0.10           rappdirs_0.3.3         Rcpp_1.0.9             scattermore_0.8       
# [45] cellranger_1.1.0       vctrs_0.5.1            spatstat.explore_3.0-5 nlme_3.1-160          
# [49] progressr_0.13.0       lmtest_0.9-40          spatstat.random_3.0-1  stringr_1.5.0         
# [53] globals_0.16.2         mime_0.12              miniUI_0.1.1.1         lifecycle_1.0.3       
# [57] irlba_2.3.5.1          goftest_1.2-3          future_1.30.0          MASS_7.3-58.1         
# [61] zoo_1.8-11             scales_1.2.1           promises_1.2.0.1       spatstat.utils_3.0-1  
# [65] parallel_4.2.2         expm_0.999-7           Exact_3.2              reticulate_1.27       
# [69] pbapply_1.7-0          stringi_1.7.12         e1071_1.7-12           boot_1.3-28           
# [73] rlang_1.0.6            pkgconfig_2.0.3        matrixStats_0.63.0     lattice_0.20-45       
# [77] ROCR_1.0-11            purrr_1.0.1            tensor_1.5             patchwork_1.1.2       
# [81] htmlwidgets_1.6.1      cowplot_1.1.1          bit_4.0.5              tidyselect_1.2.0      
# [85] parallelly_1.33.0      RcppAnnoy_0.0.20       plyr_1.8.8             magrittr_2.0.3        
# [89] generics_0.1.3         DBI_1.1.3              pillar_1.8.1           withr_2.5.0           
# [93] fitdistrplus_1.1-8     survival_3.4-0         abind_1.4-5            sp_1.5-1              
# [97] tibble_3.1.8           future.apply_1.10.0    crayon_1.5.2           KernSmooth_2.23-20    
# [101] utf8_1.2.2             spatstat.geom_3.0-3    plotly_4.10.1          readxl_1.4.1          
# [105] grid_4.2.2             data.table_1.14.6      digest_0.6.31          xtable_1.8-4          
# [109] tidyr_1.2.1            httpuv_1.6.8           munsell_0.5.0         


