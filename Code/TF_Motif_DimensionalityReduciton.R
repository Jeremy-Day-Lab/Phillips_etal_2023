#setseed and working directory
set.seed(1234)
setwd(dir = "/Volumes/Day-Lab$/RobertPhillips/Bioinformatics/Phillips_etal_2023/4_hour_ATAC/")

#load libraries needed for the analysis
library(umap)
library(ggplot2)
library(Seurat) #For NoGrid

#Read in the nmotifs matrix
nMotifs_DARs <- read.delim(file = "4hourDARs_nMotifs_1kb.txt",header = TRUE)

#Pull just the motif data and get rid of everything else
nMotifs_DARs_short <- nMotifs_DARs[,c(1,22:458)]

#Change the rownames to the first column which is the peak_id
rownames(nMotifs_DARs_short) <- nMotifs_DARs_short$PeakID..cmd.annotatePeaks.pl.4_hour_DARs_HOMER.bed.Rn7_use_this..size..1000.1000..nmotifs..m..data.user.rphill3.homer.motifs.original.all_motifs.motif.

#Get rid of the first column
nMotifs_DARs_short <- nMotifs_DARs_short[,-1]

# #alter the colnames to get rid of extraneous information
colnames(nMotifs_DARs_short) <- as.character(lapply(strsplit(x = colnames(nMotifs_DARs_short),split = "Homer"),"[",1))

#Normalize the dataframe by dividing each row by the total number of motifs found within that dataframe. 
nMotifs_DAR_normalized <- nMotifs_DARs_short/rowSums(nMotifs_DARs_short)

#Build a umap and visualize. 
#Here I will use default parameters, which are:
#n_neighbors = 15
#n_components = 2
#min_dist = 0.1
umap_vis <- umap(nMotifs_DAR_normalized,
                     preserve.seed = TRUE)

#Make a dataframe for visualization
umap_vis <- as.data.frame(umap_vis$layout)
colnames(umap_vis) <- c("UMAP_1","UMAP_2")

#Now plot. 
ggplot(data = umap_vis,aes(x = UMAP_1, y = UMAP_2)) +
  geom_point()

#Now generate clusterabble embeddings. 
#change configurations for clustering. 
#This section is based on https://umap-learn.readthedocs.io/en/latest/clustering.html
custom.config              <- umap.defaults
custom.config$min_dist     <- 1e-250 #Super small min_dist value to pack clusters tightly. 
custom.config$n_neighbors  <- 30 #Double n_neighbors to focus more on global structure, 
custom.config$n_components <- 10 #Calculate 10 umap components
umap_clustering <- umap(nMotifs_DAR_normalized,
                        config = custom.config,
                        preserve.seed = TRUE,
                        verbose = TRUE)

#Extract dataframe
umap_clustering_df <- as.data.frame(umap_clustering$layout)
colnames(umap_clustering_df) <- paste0("UMAP_",1:10)
umap_clustering_df$ID <- row.names(umap_clustering_df)

#Load the dbscan package for density based clustering. 
library(dbscan)

#Run the hdbscan algorithm
clustering <- hdbscan(x = as.matrix(umap_clustering_df[,c(1:10)]),minPts = 150)

#Add the cluster to the umap_vis dataframe
umap_vis$Cluster <- as.character(clustering$cluster)

#Now run the umap with added clusters. 
umap_with_legend <- ggplot(umap_vis,aes(x=UMAP_1,y=UMAP_2,color = Cluster)) +
  geom_point() +
  theme_bw() +
  NoGrid()
ggsave(filename = "TF_Motif_UMAP_withLegend.pdf",
       plot     = umap_with_legend,
       height   = 8,
       width    = 8)

#UMAP without legend
umap_without_legend <- ggplot(umap_vis,aes(x=UMAP_1,y=UMAP_2,color = Cluster)) +
  geom_point() +
  theme_bw() +
  NoGrid() +
  theme(legend.position = "none")
ggsave(filename = "TF_Motif_UMAP_withoutLegend.pdf",
       plot     = umap_without_legend,
       height   = 8,
       width    = 8)

#Now make plots in which the points are colored by presence or absence of a motif. 
#Make a binary matrix
nMotifs_binary <- ifelse(nMotifs_DARs_short > 0,
                         "Present",
                         "Absent")
nMotifs_binary <- as.data.frame(nMotifs_binary)

#Create an ID column for both the binary matrix and 
nMotifs_binary$ID <- row.names(nMotifs_binary)
umap_vis$ID <- row.names(umap_vis)

#Merge with the umap_vis dataframe
umap_vis_binary <- merge(x = umap_vis,
                         y = nMotifs_binary,
                         by = "ID")

for(i in colnames(umap_vis_binary)[5:441]){
  print(i)
  p <- ggplot(data = umap_vis_binary, aes(x = UMAP_1,y = UMAP_2)) +
    geom_point(aes_string(color = i)) +
    scale_color_manual(values = c("Absent" = "grey", "Present" = "red")) +
    theme_bw() +
    NoGrid() +
    labs(title = i) +
    theme(plot.title = element_text(hjust = 0.5),
          legend.position = "none") 
  ggsave(filename = paste0("/Users/rphill3/4_hour_ATAC_local/TF_UMAPs/",i,".pdf"),
         plot     = p,
         height   = 8, 
         width    = 8)
}

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
#   [1] stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
#   [1] dbscan_1.1-11      SeuratObject_4.1.3 Seurat_4.3.0       ggplot2_3.4.0      umap_0.2.10.0     
# 
# loaded via a namespace (and not attached):
#   [1] Rtsne_0.16             colorspace_2.1-0       deldir_1.0-6           ellipsis_0.3.2         ggridges_0.5.4        
# [6] rstudioapi_0.14        spatstat.data_3.0-0    leiden_0.4.3           listenv_0.9.0          farver_2.1.1          
# [11] ggrepel_0.9.3          RSpectra_0.16-1        fansi_1.0.4            codetools_0.2-18       splines_4.2.2         
# [16] polyclip_1.10-4        jsonlite_1.8.4         ica_1.0-3              cluster_2.1.4          png_0.1-8             
# [21] uwot_0.1.14            shiny_1.7.4            sctransform_0.3.5      spatstat.sparse_3.0-0  compiler_4.2.2        
# [26] httr_1.4.4             Matrix_1.5-1           fastmap_1.1.0          lazyeval_0.2.2         cli_3.6.0             
# [31] later_1.3.0            htmltools_0.5.4        tools_4.2.2            igraph_1.4.1           gtable_0.3.1          
# [36] glue_1.6.2             RANN_2.6.1             reshape2_1.4.4         dplyr_1.1.0            Rcpp_1.0.10           
# [41] scattermore_0.8        vctrs_0.5.2            spatstat.explore_3.0-6 nlme_3.1-160           progressr_0.13.0      
# [46] lmtest_0.9-40          spatstat.random_3.1-3  stringr_1.5.0          globals_0.16.2         mime_0.12             
# [51] miniUI_0.1.1.1         lifecycle_1.0.3        irlba_2.3.5.1          goftest_1.2-3          future_1.31.0         
# [56] MASS_7.3-58.1          zoo_1.8-11             scales_1.2.1           ragg_1.2.5             promises_1.2.0.1      
# [61] spatstat.utils_3.0-1   parallel_4.2.2         RColorBrewer_1.1-3     reticulate_1.28        pbapply_1.7-0         
# [66] gridExtra_2.3          stringi_1.7.12         rlang_1.0.6            pkgconfig_2.0.3        systemfonts_1.0.4     
# [71] matrixStats_0.63.0     lattice_0.20-45        ROCR_1.0-11            purrr_1.0.1            tensor_1.5            
# [76] patchwork_1.1.2        htmlwidgets_1.6.1      labeling_0.4.2         cowplot_1.1.1          tidyselect_1.2.0      
# [81] parallelly_1.34.0      RcppAnnoy_0.0.20       plyr_1.8.8             magrittr_2.0.3         R6_2.5.1              
# [86] generics_0.1.3         DBI_1.1.3              pillar_1.8.1           withr_2.5.0            fitdistrplus_1.1-8    
# [91] survival_3.4-0         abind_1.4-5            sp_1.6-0               tibble_3.1.8           future.apply_1.10.0   
# [96] crayon_1.5.2           KernSmooth_2.23-20     utf8_1.2.3             spatstat.geom_3.0-6    plotly_4.10.1         
# [101] grid_4.2.2             data.table_1.14.6      digest_0.6.31          xtable_1.8-4           tidyr_1.3.0           
# [106] httpuv_1.6.8           textshaping_0.3.6      openssl_2.0.5          munsell_0.5.0          viridisLite_0.4.1     
# [111] askpass_1.1   
