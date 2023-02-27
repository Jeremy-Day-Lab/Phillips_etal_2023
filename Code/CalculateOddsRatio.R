#set working directory and seed. 
setwd(dir = "/Volumes/Day-Lab$/RobertPhillips/Bioinformatics/Phillips_etal_2023/4_hour_ATAC/")
set.seed(1234)

#Load ggplot2 and seurat
library(ggplot2)
library(Seurat)

#Read in the vehicle peaks
Veh_Peaks <- read.delim(file = "Veh_Peaks_Merged_anno_genomiclocation.txt",sep = "\t",header = TRUE)

#Make a less detailed annotation column. 
Veh_Peaks$Short_Annotation <- as.character(lapply(strsplit(as.character(lapply(strsplit(Veh_Peaks$Annotation,split = "[(]"),"[",1))," "),"[",1))

#Read in the DAR peaks
DAR_Peaks <- read.delim(file = "Four_Hour_DARs_anno_genomiclocation.txt",sep = "\t",header = TRUE)

#Make a less detailed annotation column. 
DAR_Peaks$Short_Annotation <- as.character(lapply(strsplit(as.character(lapply(strsplit(DAR_Peaks$Annotation,split = "[(]"),"[",1))," "),"[",1))

#Now make two dataframes that are the percentage of total peaks that are in eahc genomic annotation. 
Veh_Peaks_df <- as.data.frame(table(Veh_Peaks$Short_Annotation))
Veh_Peaks_df$Peakset <- "Vehicle Peaks"

DAR_Peaks_df <- as.data.frame(table(DAR_Peaks$Short_Annotation))
DAR_Peaks_df$Peakset <- "DARs"

#merge the two dataframes
Peaks_df <- rbind(Veh_Peaks_df,DAR_Peaks_df)

#EMo
Peaks_OR <- as.data.frame(matrix(ncol = 2,nrow = length(unique(Peaks_df$Var1))))
colnames(Peaks_OR) <- unique(Peaks_df$Peakset)
rownames(Peaks_OR) <- unique(Peaks_df$Var1)
Peaks_OR$Annotation <- rownames(Peaks_OR)

#Calculate odds ratio for intergenic
#2x2
#Columns will be DARs and Vehicle peaks
#Rows will be Intergenic vs other
for(i in unique(Peaks_OR$Annotation)){
  print(i)
  Peaks_mat <- matrix(nrow = 2,ncol = 2)
  colnames(Peaks_mat) <- c(i,"Other")
  rownames(Peaks_mat) <- c("DARs","Veh_Peaks")
  Peaks_mat["DARs",i] <- subset(Peaks_df,subset=(Var1 == i & Peakset == "DARs"))$Freq
  Peaks_mat["Veh_Peaks",i] <- subset(Peaks_df,subset=(Var1 == i & Peakset == "Vehicle Peaks"))$Freq
  Peaks_mat["DARs","Other"] <- sum(subset(Peaks_df,subset=(Peakset == "DARs"))$Freq) - Peaks_mat["DARs",i]
  Peaks_mat["Veh_Peaks","Other"] <- sum(subset(Peaks_df,subset=(Peakset == "Vehicle Peaks"))$Freq) - Peaks_mat["Veh_Peaks",i]
  Peaks_OR[i,"DARs"] <- (Peaks_mat["DARs",i]*Peaks_mat["Veh_Peaks","Other"])/(Peaks_mat["Veh_Peaks",i]*Peaks_mat["DARs","Other"])
  Peaks_OR[i,"Vehicle Peaks"] <- (Peaks_mat["Veh_Peaks",i]*Peaks_mat["DARs","Other"])/(Peaks_mat["DARs",i]*Peaks_mat["Veh_Peaks","Other"])
}

#Write out the table
write.table(x         = Peaks_OR,
            file      = "Peakset_Oddsratios.txt",
            sep       = "\t",
            col.names = TRUE,
            row.names = FALSE,
            quote     = FALSE)


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
#   [1] SeuratObject_4.1.3 Seurat_4.3.0       ggplot2_3.4.0     
# 
# loaded via a namespace (and not attached):
#   [1] nlme_3.1-160           matrixStats_0.63.0     spatstat.sparse_3.0-0  RcppAnnoy_0.0.20       RColorBrewer_1.1-3    
# [6] httr_1.4.4             sctransform_0.3.5      tools_4.2.2            utf8_1.2.3             R6_2.5.1              
# [11] irlba_2.3.5.1          KernSmooth_2.23-20     DBI_1.1.3              uwot_0.1.14            lazyeval_0.2.2        
# [16] colorspace_2.1-0       withr_2.5.0            sp_1.6-0               tidyselect_1.2.0       gridExtra_2.3         
# [21] compiler_4.2.2         progressr_0.13.0       cli_3.6.0              spatstat.explore_3.0-6 plotly_4.10.1         
# [26] scales_1.2.1           lmtest_0.9-40          spatstat.data_3.0-0    ggridges_0.5.4         pbapply_1.7-0         
# [31] goftest_1.2-3          stringr_1.5.0          digest_0.6.31          spatstat.utils_3.0-1   pkgconfig_2.0.3       
# [36] htmltools_0.5.4        parallelly_1.34.0      fastmap_1.1.0          htmlwidgets_1.6.1      rlang_1.0.6           
# [41] rstudioapi_0.14        shiny_1.7.4            generics_0.1.3         zoo_1.8-11             jsonlite_1.8.4        
# [46] spatstat.random_3.1-3  ica_1.0-3              dplyr_1.1.0            magrittr_2.0.3         patchwork_1.1.2       
# [51] Matrix_1.5-1           Rcpp_1.0.10            munsell_0.5.0          fansi_1.0.4            abind_1.4-5           
# [56] reticulate_1.28        lifecycle_1.0.3        stringi_1.7.12         MASS_7.3-58.1          Rtsne_0.16            
# [61] plyr_1.8.8             grid_4.2.2             parallel_4.2.2         listenv_0.9.0          promises_1.2.0.1      
# [66] ggrepel_0.9.3          deldir_1.0-6           miniUI_0.1.1.1         lattice_0.20-45        cowplot_1.1.1         
# [71] splines_4.2.2          tensor_1.5             pillar_1.8.1           igraph_1.3.5           spatstat.geom_3.0-6   
# [76] future.apply_1.10.0    reshape2_1.4.4         codetools_0.2-18       leiden_0.4.3           glue_1.6.2            
# [81] data.table_1.14.6      png_0.1-8              vctrs_0.5.2            httpuv_1.6.8           polyclip_1.10-4       
# [86] gtable_0.3.1           RANN_2.6.1             purrr_1.0.1            tidyr_1.3.0            scattermore_0.8       
# [91] future_1.31.0          mime_0.12              xtable_1.8-4           later_1.3.0            survival_3.4-0        
# [96] viridisLite_0.4.1      tibble_3.1.8           cluster_2.1.4          globals_0.16.2         fitdistrplus_1.1-8    
# [101] ellipsis_0.3.2         ROCR_1.0-11           
