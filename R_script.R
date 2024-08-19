###My initial lod


### Set Working Directory
{
  setwd('/data/Greten_Lab/Rajiv/TAS/Zhang 2022 Immune phenotypic linkage')
  #setwd('/Users/trehanrs/Library/CloudStorage/OneDrive-NationalInstitutesofHealth/Treg/Local Work')
  #setwd("/Volumes/Greten_Lab/Rajiv/Treg")
  #setwd("/data/trehanrs/ZhangCRC2022/Zhang 2022 Immune phenotypic linkage/Uncompressed files")
  #setwd("~/Library/CloudStorage/OneDrive-NationalInstitutesofHealth/Obtained/Firouzeh/Treg")
}
### common r packages to load
{
  library(pacman)
  library(devtools)
  library(Seurat)
  library(dplyr)
  library(Matrix)
  library(base)
  library(methods)
  library(utils)
  library(stats)
  library(gdata)
  library(graphics)
  library(grDevices)
  library(ggplot2)
  library(cowplot)
  library(hdf5r)
  library(patchwork)
  library(future)
  library(future.batchtools)
  library(future.apply)
  library(doFuture)
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(AnnotationDbi)
  library(DESeq2)
  library(data.table)
  library(SingleCellExperiment)
  library(BiocParallel)
  library(BiocManager)
  library(reticulate) #use py install for python packages
  library(xlsx)
  library(scRepertoire)
  library(monocle3) #'cole-trapnell-lab/monocle3'
  library(plotly)
  #library(writexl)
}
{
  library(CellChat)
  library(patchwork)
  library(ComplexHeatmap)
  options(stringsAsFactors = FALSE)
  library(circlize)
  library(NMF)
  library(ggalluvial)
  library(devtools)
  library(ggalluvial)
  library(parallel)
  library(enrichR)
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(AnnotationDbi)
  library(AUCell)
  p_load(Seurat)
  p_load(anndata)
  p_load(SeuratDisk)
  library(monocle3)
  #p_load_gh('satijalab/seurat-wrappers')
  library(future)
  library(future.batchtools)
  library(future.apply)
  library(harmony)
  p_load('ClusterMap') 
  p_load("dittoSeq")
  library('SeuratWrappers') #remotes::install_github("satijalab/seurat-wrappers", "seurat5") #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
  #library(volcano3Ddata)
  library(volcano3D)
  library(DESeq2)
  library(limma)
  library(monocle3)
}
###Their packages 
{
  #!/usr/bin/env Rscript
  
  #library("ssVis")
  #library("Startrac")
  library("tictoc")
  library("ggpubr")
  library("ggplot2")
  library("ComplexHeatmap")
  library("RColorBrewer")
  library("circlize")
  library("data.table")
  library("plyr")
  library("ggpubr")
  library(EnhancedVolcano)
  library('scCustomize')
  library(ReactomePA)
  library(clusterProfiler)
  library(enrichplot)
  library(tidyverse)
  #source("./func.R")
  
  #RhpcBLASctl::omp_set_num_threads(1)
  #doParallel::registerDoParallel(10)
}
{
  calculatingforIndividualPatientsPerCluster <- function(x,gene,inputThresh){
    #gene <- deparse(substitute(gene))
    uniquePatients <- unique(x@meta.data$patient)
    meanCalc = list()
    #x@meta.data$TASBinary <- FALSE
    
    for (i in 1:length(uniquePatients)){
      #This works: perPatient = subset(ln_cd8,subset = patient=='patient09')
      if (sum(x@meta.data$patient!=uniquePatients[i])){
        #x@meta.data$TAS <- scale(x@meta.data$TAS)
        perPatient = subset(x,subset = patient==uniquePatients[i])
        Idents(perPatient) = perPatient@meta.data[["celltype_major"]]
        percent <- Percent_Expressing(seurat_object = perPatient,threshold = inputThresh ,features = gene, entire_object = FALSE)
        print(as.data.frame(percent))
        print(uniquePatients[i])
        meanCalc = c(meanCalc,as.data.frame(percent))
        
      }
    }
    
    return(meanCalc)
  }
  calculatingforIndividualPatientsSum <- function(x,gene,inputThresh){
    #gene <- deparse(substitute(gene))
    uniquePatients <- unique(x@meta.data$patient)
    meanCalc = list()
    #x@meta.data$TASBinary <- FALSE
    
    for (i in 1:length(uniquePatients)){
      #This works: perPatient = subset(ln_cd8,subset = patient=='patient09')
      if (sum(x@meta.data$patient!=uniquePatients[i])){
        #x@meta.data$TAS <- scale(x@meta.data$TAS)
        perPatient = subset(x,subset = patient==uniquePatients[i])
        Idents(perPatient) = perPatient@meta.data[["celltype_global"]]
        percent <- Percent_Expressing(seurat_object = perPatient,threshold = inputThresh ,features = gene, entire_object = TRUE)
        percent<-mean(percent$All_Cells)
        print(as.data.frame(percent))
        print(uniquePatients[i])
        meanCalc = c(meanCalc,as.data.frame(percent))
        
      }
    }
    
    return(meanCalc)
  }
  
  test2 <- calculatingforIndividualPatients(mt,'CD44',1)
  VlnPlot(pbmc,features = 'PDGFB')
  test2 <- calculatingforIndividualPatients(pt,'PDGFB',0)
  test2 <- calculatingforIndividualPatients(pt,'CXCL9',0)
  test2 <- calculatingforIndividualPatients(pt,'IL17A',0)
  test2 <- calculatingforIndividualPatients(pt,'CXCL10',0)
  test2 <- calculatingforIndividualPatients(pt,'TLR4',0)
  test2 <- calculatingforIndividualPatients(pt,'TLR2',0)
  test2 <- calculatingforIndividualPatients(pbmc,'TLR2',0)
  test2 <- calculatingforIndividualPatients(pt,'NFKB1',0)
  test2 <- calculatingforIndividualPatients(pt,'MAPK3',0)
  test2 <- calculatingforIndividualPatients(pt,'MAPK3',0)
  COL4A1<- calculatingforIndividualPatients(pt,'COL4A1',0)
  COL4A1_mt<- calculatingforIndividualPatients(mt,'COL4A1',0)
  COL4A2_mt<- t(calculatingforIndividualPatients(mt,'COL4A2',0))
  VlnPlot(mt,features=c('COL4A4','COL4A3','COL4A3BP','COL4A6','COL4A5','COL4A1','COL4A2','COL4A2-AS2','COL4A2-AS1'))
  COL4A2_mt<- calculatingforIndividualPatients(mt,c('COL4A1','COL4A2','COL4A3','COL4A4','COL4A5','COL4A6'),0)
  COL4A2_mt<- calculatingforIndividualPatientsSum(mt,c('COL4A1','COL4A2','COL4A3','COL4A4','COL4A5','COL4A6','COL3A1'),0)
  COL3A1_mt<- calculatingforIndividualPatients(mt,'COL3A1',0)
  
  CD163_mt<-calculatingforIndividualPatientsPerCluster(mt,'CD163',0)
  CD163_pt<-calculatingforIndividualPatientsPerCluster(pt,'CD163',0)
  Spp1_mt<-calculatingforIndividualPatientsPerCluster(mt,'SPP1',0)
  Spp1_pt<-calculatingforIndividualPatientsPerCluster(pt,'SPP1',0)
  Spp1_pt<-calculatingforIndividualPatientsPerCluster(pt,'SPP1',0)
  CD206_mt<-calculatingforIndividualPatientsPerCluster(mt,'MRC1',0)
  CD206_pt<-calculatingforIndividualPatientsPerCluster(pt,'MRC1',0)
  
  CD163_mt_t<-t(CD163_mt$Myeloid)
  CD163_pt<-t(CD163_pt)
  CD206_mt<-t(CD206_mt)
  CD206_pt<-t(CD206_pt)
  
}

#Monocle
{
  
  merged_Myeloid<-readRDS('merged_Myeloid.rds')
  merged_Myeloid_v5<-UpdateSeuratObject(merged_Myeloid)
  Idents(merged_Myeloid_v5)<- merged_Myeloid_v5$celltype_sub
  merged_Myeloid_Filtered<- subset(merged_Myeloid_v5,
                                   idents = c("hC33_cDC-LAMP3","hC34_cDC-MKI67","hC35_cDC1-CLEC9A","hC36_cDC1-CD3E","hC37_cDC2-CD1C","hC38_cDC2-C1QC","hC39_cDC2-FCN1","hC40_cDC2-TIMP1","hC41_cDC2-CD3E","hC41_cDC2-CD3E","hC42_pDC-LILRA4","hC57_Mast-CPA3","hC58_ILC3-KIT"),
                                   invert=TRUE)
  
  merged_Myeloid_Filtered <- FindNeighbors(merged_Myeloid_Filtered, dims = 1:5)
  merged_Myeloid_Filtered <- RunUMAP(merged_Myeloid_Filtered,dims=1:5,n.neighbors=200L,min.dist=0,metric='manhattan')
  f1<-DimPlot(merged_Myeloid_Filtered,group.by='celltype_sub',label=TRUE)+DimPlot(merged_Myeloid_Filtered,group.by='tissue',label=TRUE)+FeaturePlot(merged_Myeloid_Filtered,features=c('VSIG4'))
  saveImagePNG(f1,'Macropahge_Monocyte UMAP Sample and VSIG4 All Tissue',2500,1080)
  #Just Macrophage Myeloid in Liver and Colon
  merged_Myeloid_Filtered_LiverColon<- subset(merged_Myeloid_Filtered,
                                              subset=(tissue==c('metastasis normal','metastasis tumor','primary normal','primary tumor')))
  merged_Myeloid_Filtered_LiverColon <- FindNeighbors(merged_Myeloid_Filtered_LiverColon, dims = 1:10)
  merged_Myeloid_Filtered_LiverColon <- RunUMAP(merged_Myeloid_Filtered_LiverColon,dims=1:10,n.neighbors=75L,min.dist=0.1,metric='manhattan')
  legendsize<-7
  f1<-DimPlot(merged_Myeloid_Filtered_LiverColon,group.by='celltype_sub',label=TRUE,pt.size=3)+theme(legend.title = element_text( size=legendsize), legend.text=element_text(size=legendsize))
  f2<-DimPlot(merged_Myeloid_Filtered_LiverColon,group.by='tissue',label=TRUE,pt.size=3)+theme(legend.title = element_text( size=legendsize), legend.text=element_text(size=legendsize))
  f3<-FeaturePlot(merged_Myeloid_Filtered_LiverColon,features=c('VSIG4'),pt.size=3)+theme(legend.title = element_text( size=legendsize), legend.text=element_text(size=legendsize))
  saveImagePNG(f1+f2+f3,'Macropahge_Monocyte UMAP Sample and VSIG4 All Liver and Colon',1440,720)
  #Now creating seperate umaps for liver and for colon
  merged_Myeloid_Filtered_Liver<- subset(merged_Myeloid_Filtered,
                                         subset=(tissue==c('metastasis normal','metastasis tumor')))
  merged_Myeloid_Filtered_Liver <- FindNeighbors(merged_Myeloid_Filtered_Liver, dims = 1:10)
  merged_Myeloid_Filtered_Liver <- RunUMAP(merged_Myeloid_Filtered_Liver,dims=1:10,n.neighbors=75L,min.dist=0.1,metric='manhattan')
  legendsize<-7
  f1<-DimPlot(merged_Myeloid_Filtered_Liver,group.by='celltype_sub',label=TRUE,pt.size=3)+theme(legend.title = element_text( size=legendsize), legend.text=element_text(size=legendsize))
  f2<-DimPlot(merged_Myeloid_Filtered_Liver,group.by='tissue',label=TRUE,pt.size=3)+theme(legend.title = element_text( size=legendsize), legend.text=element_text(size=legendsize))
  f3<-FeaturePlot(merged_Myeloid_Filtered_Liver,features=c('VSIG4','CLEC4F','TIMD4','CSF1R','CX3CR1'),pt.size=3)+theme(legend.title = element_text( size=legendsize), legend.text=element_text(size=legendsize))
  saveImagePNG(f1+f2+f3,'Macropahge_Monocyte UMAP Sample and VSIG4 Liver',1440,720)
  
  merged_Myeloid_Filtered_Colon<- subset(merged_Myeloid_Filtered,
                                         subset=(tissue==c('primary normal','primary tumor')))
  merged_Myeloid_Filtered_Colon <- FindNeighbors(merged_Myeloid_Filtered_Colon, dims = 1:10)
  merged_Myeloid_Filtered_Colon <- RunUMAP(merged_Myeloid_Filtered_Colon,dims=1:10,n.neighbors=75L,min.dist=0.1,metric='manhattan')
  legendsize<-7
  f1<-DimPlot(merged_Myeloid_Filtered_Colon,group.by='celltype_sub',label=TRUE,pt.size=3)+theme(legend.title = element_text( size=legendsize), legend.text=element_text(size=legendsize))
  f2<-DimPlot(merged_Myeloid_Filtered_Colon,group.by='tissue',label=TRUE,pt.size=3)+theme(legend.title = element_text( size=legendsize), legend.text=element_text(size=legendsize))
  f3<-FeaturePlot(merged_Myeloid_Filtered_Colon,features=c('VSIG4','CLEC4F','TIMD4','CSF1R','CX3CR1'),pt.size=3)+theme(legend.title = element_text( size=legendsize), legend.text=element_text(size=legendsize))
  saveImagePNG(f1+f2+f3,'Macropahge_Monocyte UMAP Sample and VSIG4 Colon',1440,720)
  
  #Wikipedia:
  #TGFβR1 (ALK5) and TGFβR2 have similar ligand-binding affinities and can be distinguished from each other only by peptide mapping. 
  #Both TGFβR1 and TGFβR2 have a high affinity for TGFβ1 and low affinity for TGFβ2. 
  #TGFβR3 (β-glycan) has a high affinity for both homodimeric TGFβ1 and TGFβ2 and in addition the heterodimer TGF-β1.2
  f1<-VlnPlot(merged_Myeloid_Filtered_LiverColon,group.by = 'tissue',features=c('TGFBR1','TGFBR2','IL6R'))
  saveImagePNG(f1,'Macropahge_Monocyte Violin Plot TFGBR1 IL6R Liver Colon',1080,720)
  
  f1<-RidgePlot(merged_Myeloid_Filtered_LiverColon,group.by = 'tissue',features=c('TGFBR1','TGFBR2','IL6R'))
  saveImagePNG(f1,'Macropahge_Monocyte Ridge Plot TFGBR1 IL6R Liver Colon',1080,720)
  
  f1<-DotPlot(merged_Myeloid_Filtered_LiverColon,group.by = 'tissue',features=c('TGFBR1','TGFBR2','IL6R'))
  saveImagePNG(f1,'Macropahge_Monocyte Dot Plot TFGBR1 IL6R Liver Colon',1080,720)
  
  #Monocle
  {
    
    #Doing Monocle analysis
    merged_Myeloid_Filtered_Liver.cs <- as.cell_data_set(merged_Myeloid_Filtered_Liver)
    merged_Myeloid_Filtered_Liver.cs <- estimate_size_factors(merged_Myeloid_Filtered_Liver.cs)
    set.seed(19477)
    #Using their clustering with resolution required to provide cluster enriched with rim tissue #changing 3e-2 then 4.5 e-2
    merged_Myeloid_Filtered_Liver.cs <- cluster_cells(merged_Myeloid_Filtered_Liver.cs,cluster_method = 'leiden' , resolution=7e-4)
    
    p1 <- plot_cells(merged_Myeloid_Filtered_Liver.cs, color_cells_by = "cluster",show_trajectory_graph = FALSE,group_label_size=6,label_roots=TRUE,cell_size = 1)
    p2 <- plot_cells(merged_Myeloid_Filtered_Liver.cs, color_cells_by = "partition", show_trajectory_graph = FALSE,group_label_size=6,label_roots=TRUE,cell_size = 1)
    f1 <- wrap_plots(p1, p2)
    title<-paste0('Monocle 3 Clustering Macrophage Monocyte Liver')
    filename<-paste0('./Analysis Part 7/',title,'.png')
    png(file = filename, width = 1080,height = 720)
    feature_plot<-f1
    print(feature_plot)
    dev.off()
    
    
    #Doing Trajectory Analysis
    merged_Myeloid_Filtered_Liver.cs <- learn_graph(merged_Myeloid_Filtered_Liver.cs)
    #Setting root as cluster in middle of plot though it does not make difference to trajectory
    merged_Myeloid_Filtered_Liver.cs <- order_cells(merged_Myeloid_Filtered_Liver.cs, root_cells = colnames(merged_Myeloid_Filtered_Liver.cs[,monocle3::clusters(merged_Myeloid_Filtered_Liver.cs) == 6]))
    
    
    #Trajectory Plot with labeling and coloring from pseudo
    f1<-plot_cells(merged_Myeloid_Filtered_Liver.cs,
                   color_cells_by = "pseudotime",
                   label_groups_by_cluster=TRUE,
                   group_label_size = 8,
                   trajectory_graph_segment_size = 1.5,
                   label_leaves=FALSE,
                   label_branch_points=FALSE,
                   label_roots = FALSE,
                   #trajectory_graph_color = "grey60",
                   cell_size = 1)
    title<-paste0('Trajectory plot Mph Mono Liver with labeling and coloring from pseudotime')
    filename<-paste0('./Analysis Part 7/',title,'.png')
    png(file = filename, width = 1080,height = 720)
    feature_plot<-f1
    print(feature_plot)
    dev.off()
    
    
    #Trajectory Plot with labeling and coloring from Original 
    f1<-plot_cells(merged_Myeloid_Filtered_Liver.cs,
                   #color_cells_by = "pseudotime",
                   color_cells_by = "ident",
                   group_cells_by = "ident",
                   label_groups_by_cluster=TRUE,
                   group_label_size = 8,
                   trajectory_graph_segment_size = 1.5,
                   label_leaves=FALSE,
                   label_branch_points=FALSE,
                   label_roots = FALSE,
                   #trajectory_graph_color = "grey60",
                   cell_size = 1)
    title<-paste0('Trajectory plot Mph Mono Liver original labelling')
    filename<-paste0('./Analysis Part 7/',title,'.png')
    png(file = filename, width = 1080,height = 720)
    feature_plot<-f1
    print(feature_plot)
    dev.off()
    
    
    merged_Myeloid_Filtered_Liver$pseudotime <- pseudotime(merged_Myeloid_Filtered_Liver.cs)
    merged_Myeloid_Filtered_Liver$pseudotime[!is.finite(merged_Myeloid_Filtered_Liver$pseudotime)]<-0
    f1<-FeaturePlot(merged_Myeloid_Filtered_Liver, "pseudotime")
    title<-paste0('Pseudotime Mph Mono Liver Feature plot')
    filename<-paste0('./Analysis Part 7/',title,'.png')
    png(file = filename, width = 1080,height = 720)
    feature_plot<-f1
    print(feature_plot)
    dev.off()
    
    
    merged_Myeloid_Filtered_Liver.cs_graph_test_results <- graph_test(merged_Myeloid_Filtered_Liver.cs,
                                                                      neighbor_graph = "principal_graph",
                                                                      cores = 16)
    rowData(merged_Myeloid_Filtered_Liver.cs)$gene_short_name <- row.names(rowData(merged_Myeloid_Filtered_Liver.cs))
    merged_Myeloid_Filtered_Liver.cs.deg_ids <- rownames(subset(merged_Myeloid_Filtered_Liver.cs_graph_test_results[order(merged_Myeloid_Filtered_Liver.cs_graph_test_results$morans_I, decreasing = TRUE),], q_value < 0.05))
    
    #Top Genes Along Pseudotime
    f1<-plot_cells(merged_Myeloid_Filtered_Liver.cs,
                   genes=head(merged_Myeloid_Filtered_Liver.cs.deg_ids,n=16),
                   show_trajectory_graph = FALSE,
                   label_cell_groups = FALSE,
                   label_leaves = FALSE,
                   cell_size=1)
    title<-paste0('Top Genes Along Pseudotime Trajectory Mph Mono Liver')
    filename<-paste0('./Analysis Part 7/',title,'.png')
    png(file = filename, width = 1080,height = 720)
    feature_plot<-f1
    print(feature_plot)
    dev.off()
    
    #Choosing top segment
    merged_Myeloid_Filtered_Liver.cs_nearSPP1 <- choose_graph_segments(merged_Myeloid_Filtered_Liver.cs,clear_cds = FALSE)
    merged_Myeloid_Filtered_Liver.cs_nearSPP1 <- cluster_cells(merged_Myeloid_Filtered_Liver.cs_nearSPP1,cluster_method = 'leiden' , resolution=20e-4)
    
    p1 <- plot_cells(merged_Myeloid_Filtered_Liver.cs_nearSPP1, color_cells_by = "cluster",show_trajectory_graph = FALSE,group_label_size=6,label_roots=TRUE,cell_size = 1)
    p2 <- plot_cells(merged_Myeloid_Filtered_Liver.cs_nearSPP1, color_cells_by = "partition", show_trajectory_graph = FALSE,group_label_size=6,label_roots=TRUE,cell_size = 1)
    f1 <- wrap_plots(p1, p2)
    title<-paste0('Monocle 3 Clustering Macrophage Monocyte Liver Near SPP1 Subanalysis')
    filename<-paste0('./Analysis Part 7/',title,'.png')
    png(file = filename, width = 1080,height = 720)
    feature_plot<-f1
    print(feature_plot)
    dev.off()
    
    merged_Myeloid_Filtered_Liver.cs_nearSPP1 <- learn_graph(merged_Myeloid_Filtered_Liver.cs_nearSPP1)
    merged_Myeloid_Filtered_Liver.cs_nearSPP1_graph_test_results <- graph_test(merged_Myeloid_Filtered_Liver.cs_nearSPP1,
                                                                               neighbor_graph = "principal_graph",
                                                                               cores = 16)
    merged_Myeloid_Filtered_Liver.cs_nearSPP1.deg_ids <- rownames(subset(merged_Myeloid_Filtered_Liver.cs_nearSPP1_graph_test_results[order(merged_Myeloid_Filtered_Liver.cs_nearSPP1_graph_test_results$morans_I, decreasing = TRUE),], q_value < 0.05))
    
    # merged_Myeloid_Filtered_Liver.cs_nearSPP1_subset<-merged_Myeloid_Filtered_Liver.cs[c('TGFBR1','IL6R','SPP1')]
    # #merged_Myeloid_Filtered_Liver.cs_nearSPP1_subset<-order_cells(mmerged_Myeloid_Filtered_Liver.cs_nearSPP1_subset)
    # monocle3::plot_genes_in_pseudotime(merged_Myeloid_Filtered_Liver.cs_subset)
    
    f1<-plot_cells(merged_Myeloid_Filtered_Liver.cs_nearSPP1,
                   genes=c('SPP1','FN1'),#'CLEC4F','TIMD4','CX3CR1',
                   show_trajectory_graph = TRUE,
                   label_cell_groups = FALSE,
                   label_leaves = FALSE,
                   cell_size=1)
    title<-paste0('Monocle 3 Clustering Macrophage Monocyte Liver Near SPP1 Subanalysis Psuedotime Plots of SPP1 FN1')
    png(file = paste0('./Analysis Part 7/',title,'.png'), width = 480,height = 360)
    print(f1)
    dev.off()
    
    f1<-plot_cells(merged_Myeloid_Filtered_Liver.cs_nearSPP1,
                   genes=c('CSF1R','VSIG4'),#'CLEC4F','TIMD4','CX3CR1',
                   show_trajectory_graph = TRUE,
                   label_cell_groups = FALSE,
                   label_leaves = FALSE,
                   cell_size=1)
    title<-paste0('Monocle 3 Clustering Macrophage Monocyte Liver Near SPP1 Subanalysis Psuedotime Plots of CSF1R VSIG4')
    png(file = paste0('./Analysis Part 7/',title,'.png'), width = 480,height = 360)
    print(f1)
    dev.off()
    
    f1<-plot_cells(merged_Myeloid_Filtered_Liver.cs_nearSPP1,
                   genes=c('SPP1'),#'CLEC4F','TIMD4','CX3CR1',
                   show_trajectory_graph = TRUE,
                   label_cell_groups = FALSE,
                   label_leaves = FALSE,
                   cell_size=2)
    title<-paste0('Monocle 3 Clustering Macrophage Monocyte Liver Near SPP1 Subanalysis Psuedotime Plots of SPP1')
    png(file = paste0('./Analysis Part 7/',title,'.png'), width = 480,height = 360)
    print(f1)
    dev.off()
    
    
    f1<-plot_cells(merged_Myeloid_Filtered_Liver.cs_nearSPP1,
                   genes=c('TGFBR1','IL6R'),#'CLEC4F','TIMD4','CX3CR1',
                   show_trajectory_graph = TRUE,
                   label_cell_groups = FALSE,
                   label_leaves = FALSE,
                   cell_size=1)
    title<-paste0('Monocle 3 Clustering Macrophage Monocyte Liver Near SPP1 Subanalysis Psuedotime Plots of TGFBR1 IL6R')
    png(file = paste0('./Analysis Part 7/',title,'.png'), width = 480,height = 360)
    print(f1)
    dev.off()
    
    
    f1<-plot_cells(merged_Myeloid_Filtered_Liver.cs_nearSPP1,
                   genes=c('CD163','MRC1'),#'CLEC4F','TIMD4','CX3CR1',
                   show_trajectory_graph = TRUE,
                   label_cell_groups = FALSE,
                   label_leaves = FALSE,
                   cell_size=1)
    title<-paste0('Monocle 3 Clustering Macrophage Monocyte Liver Near SPP1 Subanalysis Psuedotime Plots of CD163 CD206')
    png(file = paste0('./Analysis Part 7/',title,'.png'), width = 480,height = 360)
    print(f1)
    dev.off()
    
    
    saveRDS(merged_Myeloid_Filtered_Liver.cs_nearSPP1,'merged_Myeloid_Filtered_Liver.cs_nearSPP1.RDS')
    saveRDS(merged_Myeloid_Filtered_Liver.cs,'merged_Myeloid_Filtered_Liver.cs.RDS')
    
    #Attempting KC score
    merged_Myeloid_Filtered_Liver_attemptingGSEA<-merged_Myeloid_Filtered_Liver
    merged_Myeloid_Filtered_Liver_attemptingGSEA<-calculate_scGSEA(merged_Myeloid_Filtered_Liver_attemptingGSEA,c('CD5L','VSIG4','SLC1A3','CD163','FOLR2','TIMD4','MARCO','GFRA2','ADRB1','TMEM26','SLC40A1','HMOX1','SLC16A9','VCAM1','SUCNR1'),'KC_score','RNA')
    merged_Myeloid_Filtered_Liver.cs_attemptingGSEA<-merged_Myeloid_Filtered_Liver.cs
    merged_Myeloid_Filtered_Liver.cs_attemptingGSEA@colData$KC_Score<-merged_Myeloid_Filtered_Liver_attemptingGSEA$KC_score
    merged_Myeloid_Filtered_Liver.cs_nearSPP1_attemptingGSEA<-choose_graph_segments(merged_Myeloid_Filtered_Liver.cs_attemptingGSEA,clear_cds = FALSE)
    merged_Myeloid_Filtered_Liver.cs_nearSPP1_attemptingGSEA <- cluster_cells(merged_Myeloid_Filtered_Liver.cs_nearSPP1_attemptingGSEA,cluster_method = 'leiden' , resolution=20e-4)
    merged_Myeloid_Filtered_Liver.cs_nearSPP1_attemptingGSEA <- learn_graph(merged_Myeloid_Filtered_Liver.cs_nearSPP1_attemptingGSEA)
    saveRDS(merged_Myeloid_Filtered_Liver.cs_nearSPP1_attemptingGSEA,'merged_Myeloid_Filtered_Liver.cs_nearSPP1_attemptingGSEA_KCScore.RDS')
    f1<-plot_cells(merged_Myeloid_Filtered_Liver.cs_nearSPP1_attemptingGSEA,
                   #genes=c('SPP1'),#'CLEC4F','TIMD4','CX3CR1',
                   color_cells_by = 'KC_Score',
                   show_trajectory_graph = TRUE,
                   label_cell_groups = FALSE,
                   label_leaves = FALSE,
                   cell_size=1)+scale_colour_gradient2('KC_score',low='Gainsboro',mid='coral',high = 'red',midpoint = 1,transform = 'log10')#,mid='misty rose'
    title<-paste0('Monocle 3 Clustering Macrophage Monocyte Liver Near SPP1 Subanalysis Psuedotime Plots of KC Genes from cell paper rev2')
    png(file = paste0('./Analysis Part 7/',title,'.png'), width = 480,height = 360)
    print(f1)
    dev.off()
    
    
    # #EvenStricter Zoom for genes with pseudotime
    # merged_Myeloid_Filtered_Liver.cs_aplotgenesinpseudotime<-merged_Myeloid_Filtered_Liver.cs
    # merged_Myeloid_Filtered_Liver.cs_aplotgenesinpseudotime<-choose_graph_segments(merged_Myeloid_Filtered_Liver.cs_aplotgenesinpseudotime,clear_cds = FALSE)
    # merged_Myeloid_Filtered_Liver.cs_aplotgenesinpseudotime <- cluster_cells(merged_Myeloid_Filtered_Liver.cs_aplotgenesinpseudotime,cluster_method = 'leiden' , resolution=20e-4)
    # merged_Myeloid_Filtered_Liver.cs_aplotgenesinpseudotime <- learn_graph(merged_Myeloid_Filtered_Liver.cs_aplotgenesinpseudotime)
    # merged_Myeloid_Filtered_Liver.cs_aplotgenesinpseudotime<-order_cells(merged_Myeloid_Filtered_Liver.cs_aplotgenesinpseudotime)
    # merged_Myeloid_Filtered_Liver.cs_aplotgenesinpseudotime_subset<-merged_Myeloid_Filtered_Liver.cs_aplotgenesinpseudotime[c('TGFBR1','IL6R','SPP1','VSIG4','CSF1R')]
    # monocle3::plot_genes_in_pseudotime(merged_Myeloid_Filtered_Liver.cs_aplotgenesinpseudotime_subset,min_expr = -1,trend_formula = "~ splines::ns(pseudotime, df=5)")
    
    #To check if cells express it
    #test<-sum(merged_Myeloid_Filtered_Liver.cs_nearSPP1@assays@data@listData[["counts"]]["VIM",]>0)
    GeneChecker<-function(x,graphtestresults,genelist){
      #genelist<-c('ACVR1B','ARRB2','BMP4','BMP4','CBX8','CCL2','CCN2','CREB3L1','CREB3L1','DDR2','DICER1','ENG','EP300','F2','F2','F2R','F2R','GLI2','IHH','IL18','INHBA','ITGA2','LARP6','LARP6','LTBP1','MKX','MYB','PDGFB','PDGFRB','PRDX5','RGCC','RGCC','SCX','SCX','SERPINB7','SERPINE1','SERPINF2','SERPINF2','SUCO','TGFB1','TGFB1','TGFB1','TGFB3','TGFB3','UCN','UTS2','VIM','VIM','WNT4','WNT4')
      #graphtestresults<-merged_Myeloid_Filtered_Liver.cs_nearSPP1_graph_test_results
      genelist<-unique(genelist)
      genelist<-genelist[genelist %in% rownames(x@assays@data@listData[["counts"]])]
      GeneToKeep<-c()
      for(i in 1:length(genelist)){
        if(sum(x@assays@data@listData[["counts"]][genelist[i],]>0)>3){ #At least 3 cells need to express it
          GeneToKeep<-rbind(GeneToKeep,genelist[i])
        }
      }
      #graphtestresults_subset<-graphtestresults[GeneToKeep,]
      #GeneToKeep <- rownames(subset(graphtestresults_subset[order(graphtestresults_subset$morans_I, decreasing = TRUE),], q_value < 0.05))
      #GeneToKeep <- rownames(subset(graphtestresults_subset[graphtestresults_subset$morans_I,], q_value < 0.05))
      return(GeneToKeep)
    }
    #Fibrosis by KC
    #GO:2000491
    f1<-plot_cells(merged_Myeloid_Filtered_Liver.cs_nearSPP1,
                   genes=c('ACTA2','DDR2','DGAT1','FGFR1','LEP','MYB','PDGFB','PDGFRN','RSP^KA1','SMO'),#'CLEC4F','TIMD4','CX3CR1',
                   show_trajectory_graph = TRUE,
                   label_cell_groups = FALSE,
                   label_leaves = FALSE,
                   cell_size=2)
    title<-paste0('Monocle 3 Clustering Macrophage Monocyte Liver Near SPP1 Subanalysis Fibrotic Genes')
    png(file = paste0('./Analysis Part 7/',title,'.png'), width = 1080,height = 720)
    print(f1)
    dev.off()
    
    merged_Myeloid_Filtered_Liver_attemptingGSEA<-merged_Myeloid_Filtered_Liver
    merged_Myeloid_Filtered_Liver_attemptingGSEA<-calculate_scGSEA(merged_Myeloid_Filtered_Liver_attemptingGSEA,c('ACVR1B','ARRB2','BMP4','BMP4','CBX8','CCL2','CCN2','CREB3L1','CREB3L1','DDR2','DICER1','ENG','EP300','F2','F2','F2R','F2R','GLI2','IHH','IL18','INHBA','ITGA2','LARP6','LARP6','LTBP1','MKX','MYB','PDGFB','PDGFRB','PRDX5','RGCC','RGCC','SCX','SCX','SERPINB7','SERPINE1','SERPINF2','SERPINF2','SUCO','TGFB1','TGFB1','TGFB1','TGFB3','TGFB3','UCN','UTS2','VIM','VIM','WNT4','WNT4'),'Positive_Collagen_Regulation','RNA')
    merged_Myeloid_Filtered_Liver.cs_attemptingGSEA<-merged_Myeloid_Filtered_Liver.cs
    merged_Myeloid_Filtered_Liver.cs_attemptingGSEA@colData$Positive_Collagen_Regulation<-merged_Myeloid_Filtered_Liver_attemptingGSEA$Positive_Collagen_Regulation
    merged_Myeloid_Filtered_Liver.cs_nearSPP1_attemptingGSEA<-choose_graph_segments(merged_Myeloid_Filtered_Liver.cs_attemptingGSEA,clear_cds = FALSE)
    merged_Myeloid_Filtered_Liver.cs_nearSPP1_attemptingGSEA <- cluster_cells(merged_Myeloid_Filtered_Liver.cs_nearSPP1_attemptingGSEA,cluster_method = 'leiden' , resolution=20e-4)
    merged_Myeloid_Filtered_Liver.cs_nearSPP1_attemptingGSEA <- learn_graph(merged_Myeloid_Filtered_Liver.cs_nearSPP1_attemptingGSEA)
    saveRDS(merged_Myeloid_Filtered_Liver.cs_nearSPP1_attemptingGSEA,'merged_Myeloid_Filtered_Liver.cs_nearSPP1_attemptingGSEACollagen.rds')
    #positive regulation of collagen biosynthetic process
    #GO:0032967
    GO0032967<-GeneChecker(merged_Myeloid_Filtered_Liver.cs_nearSPP1,merged_Myeloid_Filtered_Liver.cs_nearSPP1_graph_test_results,c('ACVR1B','ARRB2','BMP4','BMP4','CBX8','CCL2','CCN2','CREB3L1','CREB3L1','DDR2','DICER1','ENG','EP300','F2','F2','F2R','F2R','GLI2','IHH','IL18','INHBA','ITGA2','LARP6','LARP6','LTBP1','MKX','MYB','PDGFB','PDGFRB','PRDX5','RGCC','RGCC','SCX','SCX','SERPINB7','SERPINE1','SERPINF2','SERPINF2','SUCO','TGFB1','TGFB1','TGFB1','TGFB3','TGFB3','UCN','UTS2','VIM','VIM','WNT4','WNT4'))
    f1<-plot_cells(merged_Myeloid_Filtered_Liver.cs_nearSPP1,
                   genes=GO0032967,#'CLEC4F','TIMD4','CX3CR1',
                   show_trajectory_graph = TRUE,
                   label_cell_groups = FALSE,
                   label_leaves = FALSE,
                   cell_size=2)
    title<-paste0('Monocle 3 Clustering Macrophage Monocyte Liver Near SPP1 Subanalysis Collagen Production Genes')
    png(file = paste0('./Analysis Part 7/',title,'.png'), width = 1080,height = 720)
    print(f1)
    dev.off()
    #positive regulation of collagen biosynthetic process
    #GO:0032967
    f1<-plot_cells(merged_Myeloid_Filtered_Liver.cs_nearSPP1_attemptingGSEA,
                   #genes=GO0032967,#'CLEC4F','TIMD4','CX3CR1',
                   color_cells_by = 'Positive_Collagen_Regulation',
                   show_trajectory_graph = TRUE,
                   label_cell_groups = FALSE,
                   label_leaves = FALSE,
                   cell_size=1)+scale_colour_gradient2('Positive_Collagen_Regulation',low='White Smoke',mid='misty rose',high = 'coral',midpoint = 0.25)
    title<-paste0('Monocle 3 Clustering Macrophage Monocyte Liver Near SPP1 Subanalysis Collagen Production scGSEA')
    png(file = paste0('./Analysis Part 7/',title,'.png'), width = 480,height = 360)
    print(f1)
    dev.off()
    
  }
  #Now Trajectory on Colon
  {
    #Doing Monocle analysis
    merged_Myeloid_Filtered_Colon.cs <- as.cell_data_set(merged_Myeloid_Filtered_Colon)
    merged_Myeloid_Filtered_Colon.cs <- estimate_size_factors(merged_Myeloid_Filtered_Colon.cs)
    set.seed(19477)
    #Using their clustering with resolution required to provide cluster enriched with rim tissue #changing 7e-4
    merged_Myeloid_Filtered_Colon.cs <- cluster_cells(merged_Myeloid_Filtered_Colon.cs,cluster_method = 'leiden' , resolution=2.8e-3)
    
    p1 <- plot_cells(merged_Myeloid_Filtered_Colon.cs, color_cells_by = "cluster",show_trajectory_graph = FALSE,group_label_size=6,label_roots=TRUE,cell_size = 1)
    p2 <- plot_cells(merged_Myeloid_Filtered_Colon.cs, color_cells_by = "partition", show_trajectory_graph = FALSE,group_label_size=6,label_roots=TRUE,cell_size = 1)
    f1 <- wrap_plots(p1, p2)
    title<-paste0('Monocle 3 Clustering Macrophage Monocyte Colon')
    filename<-paste0('./Analysis Part 7/',title,'.png')
    png(file = filename, width = 1080,height = 720)
    feature_plot<-f1
    print(feature_plot)
    dev.off()
    
    
    
    #Doing Trajectory Analysis
    merged_Myeloid_Filtered_Colon.cs <- learn_graph(merged_Myeloid_Filtered_Colon.cs,use_partition = FALSE)
    #Setting root as cluster in middle of plot though it does not make difference to trajectory
    merged_Myeloid_Filtered_Colon.cs <- order_cells(merged_Myeloid_Filtered_Colon.cs, root_cells = colnames(merged_Myeloid_Filtered_Colon.cs[,monocle3::clusters(merged_Myeloid_Filtered_Colon.cs) == 5]))
    
    
    #Trajectory Plot with labeling and coloring from pseudo
    f1<-plot_cells(merged_Myeloid_Filtered_Colon.cs,
                   color_cells_by = "pseudotime",
                   label_groups_by_cluster=TRUE,
                   group_label_size = 8,
                   trajectory_graph_segment_size = 1.5,
                   label_leaves=FALSE,
                   label_branch_points=FALSE,
                   label_roots = FALSE,
                   #trajectory_graph_color = "grey60",
                   cell_size = 1)
    title<-paste0('Trajectory plot Mph Mono Colon with labeling and coloring from pseudotime')
    png(file = paste0('./Analysis Part 7/',title,'.png'), width = 1080,height = 720)
    print(f1)
    dev.off()
    
    
    #Trajectory Plot with labeling and coloring from Original 
    f1<-plot_cells(merged_Myeloid_Filtered_Colon.cs,
                   #color_cells_by = "pseudotime",
                   color_cells_by = "ident",
                   group_cells_by = "ident",
                   label_groups_by_cluster=TRUE,
                   group_label_size = 8,
                   trajectory_graph_segment_size = 1.5,
                   label_leaves=FALSE,
                   label_branch_points=FALSE,
                   label_roots = FALSE,
                   #trajectory_graph_color = "grey60",
                   cell_size = 1)
    title<-paste0('Trajectory plot Mph Mono Colon original labelling')
    png(file = paste0('./Analysis Part 7/',title,'.png'), width = 1080,height = 720)
    print(f1)
    dev.off()
    
    
    merged_Myeloid_Filtered_Colon$pseudotime <- pseudotime(merged_Myeloid_Filtered_Colon.cs)
    merged_Myeloid_Filtered_Colon$pseudotime[!is.finite(merged_Myeloid_Filtered_Colon$pseudotime)]<-0
    f1<-FeaturePlot(merged_Myeloid_Filtered_Colon, "pseudotime")
    title<-paste0('Pseudotime Mph Mono Colon Feature plot')
    png(file = paste0('./Analysis Part 7/',title,'.png'), width = 1080,height = 720)
    print(f1)
    dev.off()
    
    
    merged_Myeloid_Filtered_Colon.cs_graph_test_results <- graph_test(merged_Myeloid_Filtered_Colon.cs,
                                                                      neighbor_graph = "principal_graph",
                                                                      cores = 16)
    rowData(merged_Myeloid_Filtered_Colon.cs)$gene_short_name <- row.names(rowData(merged_Myeloid_Filtered_Colon.cs))
    merged_Myeloid_Filtered_Colon.cs.deg_ids <- rownames(subset(merged_Myeloid_Filtered_Colon.cs_graph_test_results[order(merged_Myeloid_Filtered_Colon.cs_graph_test_results$morans_I, decreasing = TRUE),], q_value < 0.05))
    
    
    #Top Genes Along Pseudotime
    f1<-plot_cells(merged_Myeloid_Filtered_Colon.cs,
                   genes=head(merged_Myeloid_Filtered_Colon.cs.deg_ids,n=16),
                   show_trajectory_graph = FALSE,
                   label_cell_groups = FALSE,
                   label_leaves = FALSE,
                   cell_size=1)
    title<-paste0('Top Genes Along Pseudotime Trajectory Mph Mono Colon')
    png(file = paste0('./Analysis Part 7/',title,'.png'), width = 1080,height = 720)
    print(f1)
    dev.off()
    
    
    f1<-plot_cells(merged_Myeloid_Filtered_Colon.cs,
                   genes=c('SPP1'),#'CLEC4F','TIMD4','CX3CR1',
                   show_trajectory_graph = TRUE,
                   label_cell_groups = FALSE,
                   label_leaves = FALSE,
                   cell_size=2)
    title<-paste0('Monocle 3 Clustering Macrophage Monocyte Colon Psuedotime Plots of SPP1')
    png(file = paste0('./Analysis Part 7/',title,'.png'), width = 480,height = 360)
    print(f1)
    dev.off()
    
    f1<-plot_cells(merged_Myeloid_Filtered_Colon.cs,
                   genes=c('CSF1R','VSIG4'),#'CLEC4F','TIMD4','CX3CR1',
                   show_trajectory_graph = TRUE,
                   label_cell_groups = FALSE,
                   label_leaves = FALSE,
                   cell_size=1)
    title<-paste0('Monocle 3 Clustering Macrophage Monocyte Colon Psuedotime Plots of CSF1R VSIG4')
    png(file = paste0('./Analysis Part 7/',title,'.png'), width = 480,height = 360)
    print(f1)
    dev.off()
    
    f1<-plot_cells(merged_Myeloid_Filtered_Colon.cs,
                   genes=c('SPP1','FN1'),#'CLEC4F','TIMD4','CX3CR1',
                   show_trajectory_graph = TRUE,
                   label_cell_groups = FALSE,
                   label_leaves = FALSE,
                   cell_size=1)
    title<-paste0('Monocle 3 Clustering Macrophage Monocyte Colon Psuedotime Plots of SPP1 FN1')
    png(file = paste0('./Analysis Part 7/',title,'.png'), width = 480,height = 360)
    print(f1)
    dev.off()
    
    
    f1<-plot_cells(merged_Myeloid_Filtered_Colon.cs,
                   genes=c('TGFBR1','IL6R'),#'CLEC4F','TIMD4','CX3CR1',
                   show_trajectory_graph = TRUE,
                   label_cell_groups = FALSE,
                   label_leaves = FALSE,
                   cell_size=1)
    title<-paste0('Monocle 3 Clustering Macrophage Monocyte Colon Psuedotime Plots of TGFBR1 IL6R')
    png(file = paste0('./Analysis Part 7/',title,'.png'), width = 480,height = 360)
    print(f1)
    dev.off()
    
    
    f1<-plot_cells(merged_Myeloid_Filtered_Colon.cs,
                   genes=c('CD163','MRC1'),#'CLEC4F','TIMD4','CX3CR1',
                   show_trajectory_graph = TRUE,
                   label_cell_groups = FALSE,
                   label_leaves = FALSE,
                   cell_size=1)
    title<-paste0('Monocle 3 Clustering Macrophage Monocyte LColon Psuedotime Plots of CD163 CD206')
    png(file = paste0('./Analysis Part 7/',title,'.png'), width = 480,height = 360)
    print(f1)
    dev.off()
    
    saveRDS(merged_Myeloid_Filtered_Colon.cs,'merged_Myeloid_Filtered_Colon.cs.RDS')
  }
  #CellChat
  {
    
    ###Read in Files
    {
      setwd("/data/trehanrs/ZhangCRC2022/Zhang 2022 Immune phenotypic linkage/")
      ln <- readRDS( file = "./Analysis Part 1 Output/ln_analysisp1.rds")
      mn <- readRDS(file = "./Analysis Part 1 Output/mn_analysisp1.rds")
      mt <- readRDS(file = "./Analysis Part 1 Output/mt_analysisp1.rds")
      pbmc <- readRDS(file = "./Analysis Part 1 Output/pbmc_analysisp1.rds")
      pn <- readRDS(file = "./Analysis Part 1 Output/pn_analysisp1.rds")
      pt <- readRDS(file = "./Analysis Part 1 Output/pt_analysisp1.rds")
      ln_cd8 <- readRDS(file = "./Analysis Part 2 rev3 Output/ln_cd8_analysisp2r3.rds")
      mn_cd8 <- readRDS(file = "./Analysis Part 2 rev3 Output/mn_cd8_analysisp2r3.rds")
      mt_cd8 <- readRDS(file = "./Analysis Part 2 rev3 Output/mt_cd8_analysisp2r3.rds")
      pbmc_cd8 <- readRDS(file = "./Analysis Part 2 rev3 Output/pbmc_cd8_analysisp2r3.rds")
      pn_cd8 <- readRDS( file = "./Analysis Part 2 rev3 Output/pn_cd8_analysisp2r3.rds")
      pt_cd8 <- readRDS(file = "./Analysis Part 2 rev3 Output/pt_cd8_analysisp2r3.rds")
      ln_cd8_TAS <- readRDS(file = "./Analysis Part 2 rev3 Output/ln_cd8_TAS_analysisp2r3.rds")
      mn_cd8_TAS <- readRDS(file = "./Analysis Part 2 rev3 Output/mn_cd8_TAS_analysisp2r3.rds")
      mt_cd8_TAS <- readRDS(file = "./Analysis Part 2 rev3 Output/mt_cd8_TAS_analysisp2r3.rds")
      mt_cd8_notTAS <- readRDS(file = "./Analysis Part 2 rev3 Output/mt_cd8_notTAS_analysisp2r3.rds")
      pbmc_cd8_TAS <- readRDS(file = "./Analysis Part 2 rev3 Output/pbmc_cd8_TAS_analysisp2r3.rds")
      pn_cd8_TAS <- readRDS(file = "./Analysis Part 2 rev3 Output/pn_cd8_TAS_analysisp2r3.rds")
      pt_cd8_TAS <- readRDS(file = "./Analysis Part 2 rev3 Output/pt_cd8_TAS_analysisp2r3.rds")
    }
    library(CellChat)
    library(patchwork)
    options(stringsAsFactors = FALSE)
    
    
    ###Function to do all the cell Chat Analysis
    doingCellChat <- function(x,y){
      #Lets use seurat object input instead
      #data.input = GetAssayData(object = x,slot= "scale.data")
      #meta = x@meta.data
      cellchat <- createCellChat(object = x, 
                                 
                                 datatype = "RNA",
                                 group.by = deparse(substitute(y)))
      CellChatDB <- CellChatDB.human
      CellChatDB.use <- CellChatDB
      cellchat@DB <- CellChatDB.use
      cellchat <- subsetData(cellchat) # This step is necessary even if using the whole database
      future::plan("multisession", workers = 32) # do parallel
      cellchat <- identifyOverExpressedGenes(cellchat)
      cellchat <- identifyOverExpressedInteractions(cellchat)
      cellchat_Ligand <- computeCommunProb(cellchat,population.size = TRUE) #When analyzing unsorted single-cell transcriptomes, under the assumption that abundant cell populations tend to send collectively stronger signals than the rare cell populations, CellChat can also consider the effect of cell proportion in each cell group in the probability calculation. USER can set population.size = TRUE.
      df_Ligand.net <- subsetCommunication(cellchat) # returns data frame fro ligand/recepter communicaiton, use slot.name = "netP" for signalling pathways
      cellchat_SignalPath <- computeCommunProbPathway(cellchat)
      df_SignalPath.net <- subsetCommunication(cellchat) # returns data frame fro ligand/recepter communicaiton, use slot.name = "netP" for signalling pathways
      cellchat_aggregate <- aggregateNet(cellchat)
      return(list(cellchat,cellchat_Ligand,df_Ligand.net,cellchat_SignalPath,df_SignalPath.net,cellchat_aggregate))
      
    }
    
    ###Changing to create an actual cluster called pTRT (DO NOT SAVE THESE CHANGES - ONLY FOR CELL CHAT)
    mt_cd8_TAS@meta.data$celltype_global <- "pTRT"
    mt_cd8_TAS@meta.data$celltype_major <- "pTRT"
    mt_cd8_TAS@meta.data$celltype_sub <- "pTRT"
    mt@meta.data$celltype_global[colnames(mt) %in% colnames(mt_cd8_TAS)] <- mt_cd8_TAS@meta.data$celltype_global
    mt@meta.data$celltype_major[colnames(mt) %in% colnames(mt_cd8_TAS)] <- mt_cd8_TAS@meta.data$celltype_major
    mt@meta.data$celltype_sub[colnames(mt) %in% colnames(mt_cd8_TAS)] <- mt_cd8_TAS@meta.data$celltype_sub
    
    pt_cd8_TAS@meta.data$celltype_global <- "pTRT"
    pt_cd8_TAS@meta.data$celltype_major <- "pTRT"
    pt_cd8_TAS@meta.data$celltype_sub <- "pTRT"
    pt@meta.data$celltype_global[colnames(pt) %in% colnames(pt_cd8_TAS)] <- pt_cd8_TAS@meta.data$celltype_global
    pt@meta.data$celltype_major[colnames(pt) %in% colnames(pt_cd8_TAS)] <- pt_cd8_TAS@meta.data$celltype_major
    pt@meta.data$celltype_sub[colnames(pt) %in% colnames(pt_cd8_TAS)] <- pt_cd8_TAS@meta.data$celltype_sub
    
    mn_cd8_TAS@meta.data$celltype_global <- "pTRT"
    mn_cd8_TAS@meta.data$celltype_major <- "pTRT"
    mn_cd8_TAS@meta.data$celltype_sub <- "pTRT"
    mn@meta.data$celltype_global[colnames(mn) %in% colnames(mn_cd8_TAS)] <- mn_cd8_TAS@meta.data$celltype_global
    mn@meta.data$celltype_major[colnames(mn) %in% colnames(mn_cd8_TAS)] <- mn_cd8_TAS@meta.data$celltype_major
    mn@meta.data$celltype_sub[colnames(mn) %in% colnames(mn_cd8_TAS)] <- mn_cd8_TAS@meta.data$celltype_sub
    
    pn_cd8_TAS@meta.data$celltype_global <- "pTRT"
    pn_cd8_TAS@meta.data$celltype_major <- "pTRT"
    pn_cd8_TAS@meta.data$celltype_sub <- "pTRT"
    pn@meta.data$celltype_global[colnames(pn) %in% colnames(pn_cd8_TAS)] <- pn_cd8_TAS@meta.data$celltype_global
    pn@meta.data$celltype_major[colnames(pn) %in% colnames(pn_cd8_TAS)] <- pn_cd8_TAS@meta.data$celltype_major
    pn@meta.data$celltype_sub[colnames(pn) %in% colnames(pn_cd8_TAS)] <- pn_cd8_TAS@meta.data$celltype_sub
    
    ln_cd8_TAS@meta.data$celltype_global <- "pTRT"
    ln_cd8_TAS@meta.data$celltype_major <- "pTRT"
    ln_cd8_TAS@meta.data$celltype_sub <- "pTRT"
    ln@meta.data$celltype_global[colnames(ln) %in% colnames(ln_cd8_TAS)] <- ln_cd8_TAS@meta.data$celltype_global
    ln@meta.data$celltype_major[colnames(ln) %in% colnames(ln_cd8_TAS)] <- ln_cd8_TAS@meta.data$celltype_major
    ln@meta.data$celltype_sub[colnames(ln) %in% colnames(ln_cd8_TAS)] <- ln_cd8_TAS@meta.data$celltype_sub
    
    pbmc_cd8_TAS@meta.data$celltype_global <- "pTRT"
    pbmc_cd8_TAS@meta.data$celltype_major <- "pTRT"
    pbmc_cd8_TAS@meta.data$celltype_sub <- "pTRT"
    pbmc@meta.data$celltype_global[colnames(pbmc) %in% colnames(pbmc_cd8_TAS)] <- pbmc_cd8_TAS@meta.data$celltype_global
    pbmc@meta.data$celltype_major[colnames(pbmc) %in% colnames(pbmc_cd8_TAS)] <- pbmc_cd8_TAS@meta.data$celltype_major
    pbmc@meta.data$celltype_sub[colnames(pbmc) %in% colnames(pbmc_cd8_TAS)] <- pbmc_cd8_TAS@meta.data$celltype_sub
    
    
    ###Doing the Cell Chat Analysis 
    celltype_global<-0
    celltype_major <- 0
    celltype_sub <- 0
    
    #celltype_global
    saveReturn1 <- doingCellChat(mt,celltype_global)
    cellchat_mt_global<- saveReturn1[[1]]
    cellchat_Ligand_mt_global<- saveReturn1[[2]]
    df_Ligand.net_mt_global<- saveReturn1[[3]]
    cellchat_SignalPath_mt_global<- saveReturn1[[4]]
    df_SignalPath.net_mt_global<- saveReturn1[[5]]
    cellchat_aggregate_mt_global<- saveReturn1[[6]]
    rm(saveReturn1)
    
    saveReturn <- doingCellChat(pt, celltype_global)
    cellchat_pt_global <- saveReturn[[1]]
    cellchat_Ligand_pt_global <- saveReturn[[2]]
    df_Ligand.net_pt_global <- saveReturn[[3]]
    cellchat_SignalPath_pt_global <- saveReturn[[4]]
    df_SignalPath.net_pt_global <- saveReturn[[5]]
    cellchat_aggregate_pt_global <- saveReturn[[6]]
    
    saveReturn <- doingCellChat(ln, celltype_global)
    cellchat_ln_global <- saveReturn[[1]]
    cellchat_Ligand_ln_global <- saveReturn[[2]]
    df_Ligand.net_ln_global <- saveReturn[[3]]
    cellchat_SignalPath_ln_global <- saveReturn[[4]]
    df_SignalPath.net_ln_global <- saveReturn[[5]]
    cellchat_aggregate_ln_global <- saveReturn[[6]]
    
    saveReturn <- doingCellChat(mn, celltype_global)
    cellchat_mn_global <- saveReturn[[1]]
    cellchat_Ligand_mn_global <- saveReturn[[2]]
    df_Ligand.net_mn_global <- saveReturn[[3]]
    cellchat_SignalPath_mn_global <- saveReturn[[4]]
    df_SignalPath.net_mn_global <- saveReturn[[5]]
    cellchat_aggregate_mn_global <- saveReturn[[6]]
    
    saveReturn <- doingCellChat(pn, celltype_global)
    cellchat_pn_global <- saveReturn[[1]]
    cellchat_Ligand_pn_global <- saveReturn[[2]]
    df_Ligand.net_pn_global <- saveReturn[[3]]
    cellchat_SignalPath_pn_global <- saveReturn[[4]]
    df_SignalPath.net_pn_global <- saveReturn[[5]]
    cellchat_aggregate_pn_global <- saveReturn[[6]]
    
    saveReturn <- doingCellChat(pbmc, celltype_global)
    cellchat_pbmc_global <- saveReturn[[1]]
    cellchat_Ligand_pbmc_global <- saveReturn[[2]]
    df_Ligand.net_pbmc_global <- saveReturn[[3]]
    cellchat_SignalPath_pbmc_global <- saveReturn[[4]]
    df_SignalPath.net_pbmc_global <- saveReturn[[5]]
    cellchat_aggregate_pbmc_global <- saveReturn[[6]]
    
    # 
    # saveReturn <- doingCellChat(pt,celltype_global)
    # list(cellchat_pt_global, cellchat_Ligand_pt_global, df_Ligand.net_pt_global, cellchat_SignalPath_pt_global, df_SignalPath.net_pt_global, cellchat_aggregate_pt_global) <- saveReturn
    # 
    # saveReturn <- doingCellChat(ln,celltype_global)
    # list(cellchat_ln_global, cellchat_Ligand_ln_global, df_Ligand.net_ln_global, cellchat_SignalPath_ln_global, df_SignalPath.net_ln_global, cellchat_aggregate_ln_global) <- saveReturn
    # 
    # saveReturn <- doingCellChat(mn,celltype_global)
    # list(cellchat_mn_global, cellchat_Ligand_mn_global, df_Ligand.net_mn_global, cellchat_SignalPath_mn_global, df_SignalPath.net_mn_global, cellchat_aggregate_mn_global) <- saveReturn
    # 
    # saveReturn <- doingCellChat(pn,celltype_global)
    # list(cellchat_pn_global, cellchat_Ligand_pn_global, df_Ligand.net_pn_global, cellchat_SignalPath_pn_global, df_SignalPath.net_pn_global, cellchat_aggregate_pn_global) <- saveReturn
    # 
    # saveReturn <- doingCellChat(pbmc,celltype_global)
    # list(cellchat_pbmc_global, cellchat_Ligand_pbmc_global, df_Ligand.net_pbmc_global, cellchat_SignalPath_pbmc_global, df_SignalPath.net_pbmc_global, cellchat_aggregate_pbmc_global) <- saveReturn
    
    #celltype_major
    saveReturn <- doingCellChat(mt, celltype_major)
    cellchat_mt_major <- saveReturn[[1]]
    cellchat_Ligand_mt_major <- saveReturn[[2]]
    df_Ligand.net_mt_major <- saveReturn[[3]]
    cellchat_SignalPath_mt_major <- saveReturn[[4]]
    df_SignalPath.net_mt_major <- saveReturn[[5]]
    cellchat_aggregate_mt_major <- saveReturn[[6]]
    
    saveReturn <- doingCellChat(pt, celltype_major)
    cellchat_pt_major <- saveReturn[[1]]
    cellchat_Ligand_pt_major <- saveReturn[[2]]
    df_Ligand.net_pt_major <- saveReturn[[3]]
    cellchat_SignalPath_pt_major <- saveReturn[[4]]
    df_SignalPath.net_pt_major <- saveReturn[[5]]
    cellchat_aggregate_pt_major <- saveReturn[[6]]
    
    saveReturn <- doingCellChat(ln, celltype_major)
    cellchat_ln_major <- saveReturn[[1]]
    cellchat_Ligand_ln_major <- saveReturn[[2]]
    df_Ligand.net_ln_major <- saveReturn[[3]]
    cellchat_SignalPath_ln_major <- saveReturn[[4]]
    df_SignalPath.net_ln_major <- saveReturn[[5]]
    cellchat_aggregate_ln_major <- saveReturn[[6]]
    
    saveReturn <- doingCellChat(mn, celltype_major)
    cellchat_mn_major <- saveReturn[[1]]
    cellchat_Ligand_mn_major <- saveReturn[[2]]
    df_Ligand.net_mn_major <- saveReturn[[3]]
    cellchat_SignalPath_mn_major <- saveReturn[[4]]
    df_SignalPath.net_mn_major <- saveReturn[[5]]
    cellchat_aggregate_mn_major <- saveReturn[[6]]
    
    saveReturn <- doingCellChat(pn, celltype_major)
    cellchat_pn_major <- saveReturn[[1]]
    cellchat_Ligand_pn_major <- saveReturn[[2]]
    df_Ligand.net_pn_major <- saveReturn[[3]]
    cellchat_SignalPath_pn_major <- saveReturn[[4]]
    df_SignalPath.net_pn_major <- saveReturn[[5]]
    cellchat_aggregate_pn_major <- saveReturn[[6]]
    
    saveReturn <- doingCellChat(pbmc, celltype_major)
    cellchat_pbmc_major <- saveReturn[[1]]
    cellchat_Ligand_pbmc_major <- saveReturn[[2]]
    df_Ligand.net_pbmc_major <- saveReturn[[3]]
    cellchat_SignalPath_pbmc_major <- saveReturn[[4]]
    df_SignalPath.net_pbmc_major <- saveReturn[[5]]
    cellchat_aggregate_pbmc_major <- saveReturn[[6]]
    
    # 
    # saveReturn <- doingCellChat(mt,celltype_major)
    # list(cellchat_mt_major, cellchat_Ligand_mt_major, df_Ligand.net_mt_major, cellchat_SignalPath_mt_major, df_SignalPath.net_mt_major, cellchat_aggregate_mt_major) <- saveReturn
    # 
    # saveReturn <- doingCellChat(pt,celltype_major)
    # list(cellchat_pt_major, cellchat_Ligand_pt_major, df_Ligand.net_pt_major, cellchat_SignalPath_pt_major, df_SignalPath.net_pt_major, cellchat_aggregate_pt_major) <- saveReturn
    # 
    # saveReturn <- doingCellChat(ln,celltype_major)
    # list(cellchat_ln_major, cellchat_Ligand_ln_major, df_Ligand.net_ln_major, cellchat_SignalPath_ln_major, df_SignalPath.net_ln_major, cellchat_aggregate_ln_major) <- saveReturn
    # 
    # saveReturn <- doingCellChat(mn,celltype_major)
    # list(cellchat_mn_major, cellchat_Ligand_mn_major, df_Ligand.net_mn_major, cellchat_SignalPath_mn_major, df_SignalPath.net_mn_major, cellchat_aggregate_mn_major) <- saveReturn
    # 
    # saveReturn <- doingCellChat(pn,celltype_major)
    # list(cellchat_pn_major, cellchat_Ligand_pn_major, df_Ligand.net_pn_major, cellchat_SignalPath_pn_major, df_SignalPath.net_pn_major, cellchat_aggregate_pn_major) <- saveReturn
    # 
    # saveReturn <- doingCellChat(pbmc,celltype_major)
    # list(cellchat_pbmc_major, cellchat_Ligand_pbmc_major, df_Ligand.net_pbmc_major, cellchat_SignalPath_pbmc_major, df_SignalPath.net_pbmc_major, cellchat_aggregate_pbmc_major) <- saveReturn
    
    #celltype_sub
    saveReturn <- doingCellChat(mt, celltype_sub)
    cellchat_mt_sub <- saveReturn[[1]]
    cellchat_Ligand_mt_sub <- saveReturn[[2]]
    df_Ligand.net_mt_sub <- saveReturn[[3]]
    cellchat_SignalPath_mt_sub <- saveReturn[[4]]
    df_SignalPath.net_mt_sub <- saveReturn[[5]]
    cellchat_aggregate_mt_sub <- saveReturn[[6]]
    
    saveReturn <- doingCellChat(pt, celltype_sub)
    cellchat_pt_sub <- saveReturn[[1]]
    cellchat_Ligand_pt_sub <- saveReturn[[2]]
    df_Ligand.net_pt_sub <- saveReturn[[3]]
    cellchat_SignalPath_pt_sub <- saveReturn[[4]]
    df_SignalPath.net_pt_sub <- saveReturn[[5]]
    cellchat_aggregate_pt_sub <- saveReturn[[6]]
    
    saveReturn <- doingCellChat(ln, celltype_sub)
    cellchat_ln_sub <- saveReturn[[1]]
    cellchat_Ligand_ln_sub <- saveReturn[[2]]
    df_Ligand.net_ln_sub <- saveReturn[[3]]
    cellchat_SignalPath_ln_sub <- saveReturn[[4]]
    df_SignalPath.net_ln_sub <- saveReturn[[5]]
    cellchat_aggregate_ln_sub <- saveReturn[[6]]
    
    saveReturn <- doingCellChat(mn, celltype_sub)
    cellchat_mn_sub <- saveReturn[[1]]
    cellchat_Ligand_mn_sub <- saveReturn[[2]]
    df_Ligand.net_mn_sub <- saveReturn[[3]]
    cellchat_SignalPath_mn_sub <- saveReturn[[4]]
    df_SignalPath.net_mn_sub <- saveReturn[[5]]
    cellchat_aggregate_mn_sub <- saveReturn[[6]]
    
    saveReturn <- doingCellChat(pn, celltype_sub)
    cellchat_pn_sub <- saveReturn[[1]]
    cellchat_Ligand_pn_sub <- saveReturn[[2]]
    df_Ligand.net_pn_sub <- saveReturn[[3]]
    cellchat_SignalPath_pn_sub <- saveReturn[[4]]
    df_SignalPath.net_pn_sub <- saveReturn[[5]]
    cellchat_aggregate_pn_sub <- saveReturn[[6]]
    
    saveReturn <- doingCellChat(pbmc, celltype_sub)
    cellchat_pbmc_sub <- saveReturn[[1]]
    cellchat_Ligand_pbmc_sub <- saveReturn[[2]]
    df_Ligand.net_pbmc_sub <- saveReturn[[3]]
    cellchat_SignalPath_pbmc_sub <- saveReturn[[4]]
    df_SignalPath.net_pbmc_sub <- saveReturn[[5]]
    cellchat_aggregate_pbmc_sub <- saveReturn[[6]]
    
    
    # saveReturn <- doingCellChat(mt,celltype_sub)
    # list(cellchat_mt_sub, cellchat_Ligand_mt_sub, df_Ligand.net_mt_sub, cellchat_SignalPath_mt_sub, df_SignalPath.net_mt_sub, cellchat_aggregate_mt_sub) <- saveReturn
    # 
    # saveReturn <- doingCellChat(pt,celltype_sub)
    # list(cellchat_pt_sub, cellchat_Ligand_pt_sub, df_Ligand.net_pt_sub, cellchat_SignalPath_pt_sub, df_SignalPath.net_pt_sub, cellchat_aggregate_pt_sub) <- saveReturn
    # 
    # saveReturn <- doingCellChat(ln,celltype_sub)
    # list(cellchat_ln_sub, cellchat_Ligand_ln_sub, df_Ligand.net_ln_sub, cellchat_SignalPath_ln_sub, df_SignalPath.net_ln_sub, cellchat_aggregate_ln_sub) <- saveReturn
    # 
    # saveReturn <- doingCellChat(mn,celltype_sub)
    # list(cellchat_mn_sub, cellchat_Ligand_mn_sub, df_Ligand.net_mn_sub, cellchat_SignalPath_mn_sub, df_SignalPath.net_mn_sub, cellchat_aggregate_mn_sub) <- saveReturn
    # 
    # saveReturn <- doingCellChat(pn,celltype_sub)
    # list(cellchat_pn_sub, cellchat_Ligand_pn_sub, df_Ligand.net_pn_sub, cellchat_SignalPath_pn_sub, df_SignalPath.net_pn_sub, cellchat_aggregate_pn_sub) <- saveReturn
    # 
    # saveReturn <- doingCellChat(pbmc,celltype_sub)
    # list(cellchat_pbmc_sub, cellchat_Ligand_pbmc_sub, df_Ligand.net_pbmc_sub, cellchat_SignalPath_pbmc_sub, df_SignalPath.net_pbmc_sub, cellchat_aggregate_pbmc_sub) <- saveReturn
    
    ###Saving Analysis part 3
    #saveRDS(ln_cd8, file = "./Analysis Part 2 rev3 Output/ln_cd8_analysisp2r3.rds")
    # Define the folder path where you want to save the RDS files
    folder_path <- "/data/trehanrs/ZhangCRC2022/Zhang 2022 Immune phenotypic linkage/Analysis part 4 Output/"
    
    #object_names <- c("cellchat_mt_global", "cellchat_Ligand_mt_global", "df_Ligand.net_mt_global", 
    #                  "cellchat_SignalPath_mt_global", "df_SignalPath.net_mt_global", "cellchat_aggregate_mt_global")
    
    # Create a list of object names and corresponding filenames
    #list(cellchat_mt_global, cellchat_Ligand_mt_global, df_Ligand.net_mt_global, cellchat_SignalPath_mt_global, df_SignalPath.net_mt_global, cellchat_aggregate_mt_global) <- saveReturn
    object_names <- c(
      "cellchat_mt_global", "cellchat_Ligand_mt_global", "df_Ligand.net_mt_global",
      "cellchat_SignalPath_mt_global", "df_SignalPath.net_mt_global", "cellchat_aggregate_mt_global",
      "cellchat_pt_global", "cellchat_Ligand_pt_global", "df_Ligand.net_pt_global",
      "cellchat_SignalPath_pt_global", "df_SignalPath.net_pt_global", "cellchat_aggregate_pt_global",
      "cellchat_ln_global", "cellchat_Ligand_ln_global", "df_Ligand.net_ln_global",
      "cellchat_SignalPath_ln_global", "df_SignalPath.net_ln_global", "cellchat_aggregate_ln_global",
      "cellchat_mn_global", "cellchat_Ligand_mn_global", "df_Ligand.net_mn_global",
      "cellchat_SignalPath_mn_global", "df_SignalPath.net_mn_global", "cellchat_aggregate_mn_global",
      "cellchat_pn_global", "cellchat_Ligand_pn_global", "df_Ligand.net_pn_global",
      "cellchat_SignalPath_pn_global", "df_SignalPath.net_pn_global", "cellchat_aggregate_pn_global",
      "cellchat_pbmc_global", "cellchat_Ligand_pbmc_global", "df_Ligand.net_pbmc_global",
      "cellchat_SignalPath_pbmc_global", "df_SignalPath.net_pbmc_global", "cellchat_aggregate_pbmc_global",
      "cellchat_mt_major", "cellchat_Ligand_mt_major", "df_Ligand.net_mt_major",
      "cellchat_SignalPath_mt_major", "df_SignalPath.net_mt_major", "cellchat_aggregate_mt_major",
      "cellchat_pt_major", "cellchat_Ligand_pt_major", "df_Ligand.net_pt_major",
      "cellchat_SignalPath_pt_major", "df_SignalPath.net_pt_major", "cellchat_aggregate_pt_major",
      "cellchat_ln_major", "cellchat_Ligand_ln_major", "df_Ligand.net_ln_major",
      "cellchat_SignalPath_ln_major", "df_SignalPath.net_ln_major", "cellchat_aggregate_ln_major",
      "cellchat_mn_major", "cellchat_Ligand_mn_major", "df_Ligand.net_mn_major",
      "cellchat_SignalPath_mn_major", "df_SignalPath.net_mn_major", "cellchat_aggregate_mn_major",
      "cellchat_pn_major", "cellchat_Ligand_pn_major", "df_Ligand.net_pn_major",
      "cellchat_SignalPath_pn_major", "df_SignalPath.net_pn_major", "cellchat_aggregate_pn_major",
      "cellchat_pbmc_major", "cellchat_Ligand_pbmc_major", "df_Ligand.net_pbmc_major",
      "cellchat_SignalPath_pbmc_major", "df_SignalPath.net_pbmc_major", "cellchat_aggregate_pbmc_major",
      "cellchat_mt_sub", "cellchat_Ligand_mt_sub", "df_Ligand.net_mt_sub",
      "cellchat_SignalPath_mt_sub", "df_SignalPath.net_mt_sub", "cellchat_aggregate_mt_sub",
      "cellchat_pt_sub", "cellchat_Ligand_pt_sub", "df_Ligand.net_pt_sub",
      "cellchat_SignalPath_pt_sub", "df_SignalPath.net_pt_sub", "cellchat_aggregate_pt_sub",
      "cellchat_ln_sub", "cellchat_Ligand_ln_sub", "df_Ligand.net_ln_sub",
      "cellchat_SignalPath_ln_sub", "df_SignalPath.net_ln_sub", "cellchat_aggregate_ln_sub",
      "cellchat_mn_sub", "cellchat_Ligand_mn_sub", "df_Ligand.net_mn_sub",
      "cellchat_SignalPath_mn_sub", "df_SignalPath.net_mn_sub", "cellchat_aggregate_mn_sub",
      "cellchat_pn_sub", "cellchat_Ligand_pn_sub", "df_Ligand.net_pn_sub",
      "cellchat_SignalPath_pn_sub", "df_SignalPath.net_pn_sub", "cellchat_aggregate_pn_sub",
      "cellchat_pbmc_sub", "cellchat_Ligand_pbmc_sub", "df_Ligand.net_pbmc_sub",
      "cellchat_SignalPath_pbmc_sub", "df_SignalPath.net_pbmc_sub", "cellchat_aggregate_pbmc_sub"
    )
    
    filenames <- paste0(object_names, ".rds")
    
    # Loop through the objects and save them as individual RDS files
    for (i in seq_along(object_names)) {
      object <- get(object_names[i])
      filepath <- file.path(folder_path, filenames[i])
      saveRDS(object, filepath)
    }
  }
}
#CellChat Graphs
{
  ###Setting up Future
  {
    future::plan("multisession", workers = 4) # do parallel processing
    options(future.globals.maxSize = 1000 * 1000 * 1024^2)
    `%<-%` <- future::`%<-%`
  }
  
  
  usingFuture = FALSE
  ###Read in Files
  {
    if(usingFuture){
      cellchat_ln_global %<-% readRDS(file = './Analysis part 4 Output Swarm/cellchat_ln_global.rds')
      cellchat_ln_major %<-% readRDS(file = './Analysis part 4 Output Swarm/cellchat_ln_major.rds')
      cellchat_ln_sub %<-% readRDS(file = './Analysis part 4 Output Swarm/cellchat_ln_sub.rds')
      cellchat_mn_global %<-% readRDS(file = './Analysis part 4 Output Swarm/cellchat_mn_global.rds')
      cellchat_mn_major %<-% readRDS(file = './Analysis part 4 Output Swarm/cellchat_mn_major.rds')
      cellchat_mn_sub %<-% readRDS(file = './Analysis part 4 Output Swarm/cellchat_mn_sub.rds')
      cellchat_mt_global %<-% readRDS(file = './Analysis part 4 Output Swarm/cellchat_mt_global.rds')
      cellchat_mt_major %<-% readRDS(file = './Analysis part 4 Output Swarm/cellchat_mt_major.rds')
      cellchat_mt_sub %<-% readRDS(file = './Analysis part 4 Output Swarm/cellchat_mt_sub.rds')
      cellchat_pbmc_global %<-% readRDS(file = './Analysis part 4 Output Swarm/cellchat_pbmc_global.rds')
      cellchat_pbmc_major %<-% readRDS(file = './Analysis part 4 Output Swarm/cellchat_pbmc_major.rds')
      cellchat_pbmc_sub %<-% readRDS(file = './Analysis part 4 Output Swarm/cellchat_pbmc_sub.rds')
      cellchat_pn_global %<-% readRDS(file = './Analysis part 4 Output Swarm/cellchat_pn_global.rds')
      cellchat_pn_major %<-% readRDS(file = './Analysis part 4 Output Swarm/cellchat_pn_major.rds')
      cellchat_pn_sub %<-% readRDS(file = './Analysis part 4 Output Swarm/cellchat_pn_sub.rds')
      cellchat_pt_global %<-% readRDS(file = './Analysis part 4 Output Swarm/cellchat_pt_global.rds')
      cellchat_pt_major %<-% readRDS(file = './Analysis part 4 Output Swarm/cellchat_pt_major.rds')
      cellchat_pt_sub %<-% readRDS(file = './Analysis part 4 Output Swarm/cellchat_pt_sub.rds')
    }else{
      cellchat_ln_global <- readRDS(file = './Analysis part 4 Output Swarm/cellchat_ln_global.rds')
      cellchat_ln_major <- readRDS(file = './Analysis part 4 Output Swarm/cellchat_ln_major.rds')
      cellchat_ln_sub <- readRDS(file = './Analysis part 4 Output Swarm/cellchat_ln_sub.rds')
      cellchat_mn_global <- readRDS(file = './Analysis part 4 Output Swarm/cellchat_mn_global.rds')
      cellchat_mn_major <- readRDS(file = './Analysis part 4 Output Swarm/cellchat_mn_major.rds')
      cellchat_mn_sub <- readRDS(file = './Analysis part 4 Output Swarm/cellchat_mn_sub.rds')
      cellchat_mt_global <- readRDS(file = './Analysis part 4 Output Swarm/cellchat_mt_global.rds')
      cellchat_mt_major <- readRDS(file = './Analysis part 4 Output Swarm/cellchat_mt_major.rds')
      cellchat_mt_sub <- readRDS(file = './Analysis part 4 Output Swarm/cellchat_mt_sub.rds')
      cellchat_pbmc_global <- readRDS(file = './Analysis part 4 Output Swarm/cellchat_pbmc_global.rds')
      cellchat_pbmc_major <- readRDS(file = './Analysis part 4 Output Swarm/cellchat_pbmc_major.rds')
      cellchat_pbmc_sub <- readRDS(file = './Analysis part 4 Output Swarm/cellchat_pbmc_sub.rds')
      cellchat_pn_global <- readRDS(file = './Analysis part 4 Output Swarm/cellchat_pn_global.rds')
      cellchat_pn_major <- readRDS(file = './Analysis part 4 Output Swarm/cellchat_pn_major.rds')
      cellchat_pn_sub <- readRDS(file = './Analysis part 4 Output Swarm/cellchat_pn_sub.rds')
      cellchat_pt_global <- readRDS(file = './Analysis part 4 Output Swarm/cellchat_pt_global.rds')
      cellchat_pt_major <- readRDS(file = './Analysis part 4 Output Swarm/cellchat_pt_major.rds')
      cellchat_pt_sub <- readRDS(file = './Analysis part 4 Output Swarm/cellchat_pt_sub.rds')
      
    }
    
    cellchat_list <- list(
      cellchat_ln_global = cellchat_ln_global,
      cellchat_ln_major = cellchat_ln_major,
      cellchat_ln_sub = cellchat_ln_sub,
      cellchat_mn_global = cellchat_mn_global,
      cellchat_mn_major = cellchat_mn_major,
      cellchat_mn_sub = cellchat_mn_sub,
      cellchat_mt_global = cellchat_mt_global,
      cellchat_mt_major = cellchat_mt_major,
      cellchat_mt_sub = cellchat_mt_sub,
      cellchat_pbmc_global = cellchat_pbmc_global,
      cellchat_pbmc_major = cellchat_pbmc_major,
      cellchat_pbmc_sub = cellchat_pbmc_sub,
      cellchat_pn_global = cellchat_pn_global,
      cellchat_pn_major = cellchat_pn_major,
      cellchat_pn_sub = cellchat_pn_sub,
      cellchat_pt_global = cellchat_pt_global,
      cellchat_pt_major = cellchat_pt_major,
      cellchat_pt_sub = cellchat_pt_sub
    )
  }
  
  ###Sub-Function for Individual Graphing of Interaction Strength of pTRT only for Sub (if you put it all in one graph it is too small to read)
  individualInteractionStrength_pTRT_Sub<- function(x,z,sampleNam,groupSizeinput){
    groupNames <- c("CD8","CD4_","DC","Mono","Mph","NK","B-","plasma")
    for(indnum in 1:length(groupNames)){
      #x is sender and z is reciever
      #CD8
      #temparray = 1:length(grep("CD8",colnames(x)))
      #Sender
      x_initial <- x[,c(grep(groupNames[indnum],colnames(x)),grep("pTRT",colnames(x)))]
      subset_x <- x_initial[c(grep(groupNames[indnum],rownames(x_initial)),grep("pTRT",rownames(x_initial))),]
      groupSize_x <- groupSizeinput[c(grep(groupNames[indnum],rownames(x_initial)),grep("pTRT",rownames(x_initial)))]
      #Reciever
      z_initial <- z[,c(grep(groupNames[indnum],colnames(z)),grep("pTRT",colnames(z)))]
      subset_z <- z_initial[c(grep(groupNames[indnum],rownames(z_initial)),grep("pTRT",rownames(z_initial))),]
      groupSize_z<-groupSize_x
      #groupSize_z <- as.numeric(table(subset_z))
      #Graphing
      png(filename=paste0("./CellChat Graphs Analysis Part 5/Visualize pTRT Individually ",sampleNam,groupNames[indnum],".png"), width = 1920,height = 1080)
      
      par(mfrow = c(1,2))
      netVisual_circle(subset_x, vertex.weight = groupSize_x, weight.scale = T, title.name = "pTRT Sender Interaction Strength")
      netVisual_circle(subset_z, vertex.weight = groupSize_z, weight.scale = T, title.name = "pTRT Reciever Interaction Strength")
      dev.off()
    }
  }
  
  ###Making Cell Chat Graphs for each Sample Total Number and Aggregate Interactions
  {
    
    #Explnations of edge color/weight, node color/size/shape: In all visualization plots, edge colors are  
    #consistent with the sources as sender, and edge weights are proportional to the interaction strength. 
    #Thicker edge line indicates a stronger signal. In the Hierarchy plot and Circle plot, 
    #circle sizes are proportional to the number of cells in each cell group. In the hierarchy
    #plot, solid and open circles represent source and target, respectively. In the Chord 
    #diagram, the inner thinner bar colors represent the targets that receive signal from
    #the corresponding outer bar. The inner bar size is proportional to the signal strength
    #received by the targets. Such inner bar is helpful for interpreting the complex chord diagram. Note 
    #that there exist some inner bars without any chord for some cell groups, please just igore it because this 
    #is an issue that has not been addressed by circlize package.
    
    #Summary, edge color = sender; edge weight/inner bar size = interaction strength; circle size = # of cells; 
    #solid cricle = source; open circle = target; 
    
    #output_dir <- "./CellChat Graphs Analysis Part 5/"
    
    ##Aggregated Whole All Interactions
    for (i in 1:length(cellchat_list)){
      cellchat <- cellchat_list[[i]]
      groupSize <- as.numeric(table(cellchat@idents))
      
      
      png(filename=paste0("./CellChat Graphs Analysis Part 5/Visualize Aggregated ",names(cellchat_list)[i]," group.png"), width = 1920,height = 1080)
      
      par(mfrow = c(1,2))
      netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
      netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")
      dev.off()
    }
    rm(cellchat)
    rm(i)
    rm(groupSize)
    
    #So the row name is the sender and column name is the reciever
    for (i in 1:length(cellchat_list)){
      #Graphing pTRT as Sender Strength and Reciever Strength
      cellchat <- cellchat_list[[i]]
      groupSize <- as.numeric(table(cellchat@idents))
      mat <- cellchat@net$weight
      #Creating empty matrix
      mat_pTRT_sender<- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
      mat_pTRT_reciever<- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
      #Filling matrixes
      mat_pTRT_sender[which(rownames(mat)=='pTRT'),] = mat[which(rownames(mat)=='pTRT'), ]
      mat_pTRT_reciever[,which(colnames(mat)=='pTRT')] = mat[,which(colnames(mat)=='pTRT')]
      if(i%%3>0)
      {
        #Creating Png Save file
        png(filename=paste0("./CellChat Graphs Analysis Part 5/Visualize Aggregated ",names(cellchat_list)[i]," pTRT.png"), width = 1920,height = 1080)
        
        par(mfrow = c(1,2))
        
        netVisual_circle(mat_pTRT_sender, vertex.weight = groupSize, weight.scale = T, title.name = "pTRT Sender Interaction Strength")
        netVisual_circle(mat_pTRT_reciever, vertex.weight = groupSize, weight.scale = T, title.name = "pTRT Reciever Interaction Strength")
        
        dev.off()
      } else{
        individualInteractionStrength_pTRT_Sub(mat_pTRT_sender,mat_pTRT_reciever,names(cellchat_list)[i],groupSize)
        
      }
      
    }
    
  }
  
  ###Comparing 
  {
    #Creating Difference Objects
    {
      #Major
      object.list_major_diffPtMt <- list(PT = cellchat_pt_major,MT = cellchat_mt_major)
      cellchat_major_diffPtMt <- mergeCellChat(object.list_major_diffPtMt, add.names = names(object.list_major_diffPtMt))
      
      object.list_major_diffPnMn <- list(PN = cellchat_pn_major,MN = cellchat_mn_major)
      cellchat_major_diffPnMn <- mergeCellChat(object.list_major_diffPnMn, add.names = names(object.list_major_diffPnMn))
      
      object.list_major_diffLnMn <- list(LN = cellchat_ln_major,MN = cellchat_mn_major)
      cellchat_major_diffLnMn <- mergeCellChat(object.list_major_diffLnMn, add.names = names(object.list_major_diffLnMn))
      
      object.list_major_diffLnMt <- list(LN = cellchat_ln_major,MT = cellchat_mt_major)
      cellchat_major_diffLnMt <- mergeCellChat(object.list_major_diffLnMt, add.names = names(object.list_major_diffLnMt))
      
      object.list_major_diffMnMt <- list(MN = cellchat_mn_major,MT = cellchat_mt_major)
      cellchat_major_diffMnMt <- mergeCellChat(object.list_major_diffMnMt, add.names = names(object.list_major_diffMnMt))
      
      #Sub
      object.list_sub_diffPtMt <- list(PT = cellchat_pt_sub,MT = cellchat_mt_sub)
      cellchat_sub_diffPtMt <- mergeCellChat(object.list_sub_diffPtMt, add.names = names(object.list_sub_diffPtMt))
      
      object.list_sub_diffPnMn <- list(PN = cellchat_pn_sub,MN = cellchat_mn_sub)
      cellchat_sub_diffPnMn <- mergeCellChat(object.list_sub_diffPnMn, add.names = names(object.list_sub_diffPnMn))
      
      object.list_sub_diffLnMn <- list(LN = cellchat_ln_sub,MN = cellchat_mn_sub)
      cellchat_sub_diffLnMn <- mergeCellChat(object.list_sub_diffLnMn, add.names = names(object.list_sub_diffLnMn))
      
      object.list_sub_diffLnMt <- list(LN = cellchat_ln_sub,MT = cellchat_mt_sub)
      cellchat_sub_diffLnMt <- mergeCellChat(object.list_sub_diffLnMt, add.names = names(object.list_sub_diffLnMt))
      
      object.list_sub_diffMnMt <- list(MN = cellchat_mn_sub,MT = cellchat_mt_sub)
      cellchat_sub_diffMnMt <- mergeCellChat(object.list_sub_diffMnMt, add.names = names(object.list_sub_diffMnMt))
      
      
    }
    ###Overall Number of interactions
    {
      #Starting with Major as it does not matter which sample you use because granularity does not affect overall output
      
      #Compare total number of interactions
      png(filename=paste0("./CellChat Graphs Analysis Part 5 Comparison/MT vs PT Overall (not just pTRT) Interaction Number and Strength.png"), res=120,width = 1920/1,height = 1080/1)
      gg1 <- compareInteractions(cellchat_major_diffPtMt, show.legend = F, group = c(1,2))
      gg2 <- compareInteractions(cellchat_major_diffPtMt, show.legend = F, group = c(1,2), measure = "weight")
      gg1 + gg2
      dev.off()
      
      #Compare total number of interactions
      png(filename=paste0("./CellChat Graphs Analysis Part 5 Comparison/MN vs PN Overall (not just pTRT) Interaction Number and Strength.png"), res=120,width = 1920/1,height = 1080/1)
      gg1 <- compareInteractions(cellchat_major_diffPnMn, show.legend = F, group = c(1,2))
      gg2 <- compareInteractions(cellchat_major_diffPnMn, show.legend = F, group = c(1,2), measure = "weight")
      gg1 + gg2
      dev.off()
      
      #Compare total number of interactions
      png(filename=paste0("./CellChat Graphs Analysis Part 5 Comparison/LN vs MN Overall (not just pTRT) Interaction Number and Strength.png"), res=120,width = 1920/1,height = 1080/1)
      gg1 <- compareInteractions(cellchat_major_diffLnMn, show.legend = F, group = c(1,2))
      gg2 <- compareInteractions(cellchat_major_diffLnMn, show.legend = F, group = c(1,2), measure = "weight")
      gg1 + gg2
      dev.off()
      
      #Compare total number of interactions
      png(filename=paste0("./CellChat Graphs Analysis Part 5 Comparison/LN vs MT Overall (not just pTRT) Interaction Number and Strength.png"), res=120,width = 1920/1,height = 1080/1)
      gg1 <- compareInteractions(cellchat_major_diffLnMt, show.legend = F, group = c(1,2))
      gg2 <- compareInteractions(cellchat_major_diffLnMt, show.legend = F, group = c(1,2), measure = "weight")
      gg1 + gg2
      dev.off()
      
      #Compare total number of interactions
      png(filename=paste0("./CellChat Graphs Analysis Part 5 Comparison/MN vs MT Overall (not just pTRT) Interaction Number and Strength.png"), res=120,width = 1920/1,height = 1080/1)
      gg1 <- compareInteractions(cellchat_major_diffMnMt, show.legend = F, group = c(1,2))
      gg2 <- compareInteractions(cellchat_major_diffMnMt, show.legend = F, group = c(1,2), measure = "weight")
      gg1 + gg2
      dev.off()
    }
    
    
    #Heatmap of difference in interactions for major (as sub would be terrible)
    {
      #MT vs PT
      png(filename=paste0("./CellChat Graphs Analysis Part 5 Comparison/PT vs MT Heatmap.png"), res=120,width = 1920/1,height = 1080/1)
      #Number of interactions
      gg1 <- netVisual_heatmap(cellchat_major_diffPtMt)
      #Strength of interactions
      gg2 <- netVisual_heatmap(cellchat_major_diffPtMt ,measure = "weight")
      mergedgg <- gg1 + gg2
      draw(mergedgg, column_title = "Differential Cell Interaction Hepatic Tumor (red) vs Primary tumor (blue)", column_title_gp = gpar(fontsize = 16))
      dev.off()
      
      #PN VS MN
      png(filename=paste0("./CellChat Graphs Analysis Part 5 Comparison/PN vs MN Heatmap.png"), res=120,width = 1920/1,height = 1080/1)
      #Number of interactions
      gg1 <- netVisual_heatmap(cellchat_major_diffPnMn)
      #Strength of interactions
      gg2 <- netVisual_heatmap(cellchat_major_diffPnMn ,measure = "weight")
      mergedgg <- gg1 + gg2
      draw(mergedgg, column_title = "Differential Cell Interaction Hepatic Normal (red) vs Primary Normal (blue)", column_title_gp = gpar(fontsize = 16))
      dev.off()
      
      #LN vs MN
      png(filename=paste0("./CellChat Graphs Analysis Part 5 Comparison/LN vs MN Heatmap.png"), res=120,width = 1920/1,height = 1080/1)
      #Number of interactions
      gg1 <- netVisual_heatmap(cellchat_major_diffLnMn)
      #Strength of interactions
      gg2 <- netVisual_heatmap(cellchat_major_diffLnMn ,measure = "weight")
      mergedgg <- gg1 + gg2
      draw(mergedgg, column_title = "Differential Cell Interaction Hepatic Normal (red) vs Lymph Node (blue)", column_title_gp = gpar(fontsize = 16))
      dev.off()
      
      #LN vs MT
      png(filename=paste0("./CellChat Graphs Analysis Part 5 Comparison/LN vs MT Heatmap.png"), res=120,width = 1920/1,height = 1080/1)
      #Number of interactions
      gg1 <- netVisual_heatmap(cellchat_major_diffLnMt)
      #Strength of interactions
      gg2 <- netVisual_heatmap(cellchat_major_diffLnMt ,measure = "weight")
      mergedgg <- gg1 + gg2
      draw(mergedgg, column_title = "Differential Cell Interaction Hepatic Tumor (red) vs Lymph Node (blue)", column_title_gp = gpar(fontsize = 16))
      dev.off()
      
      #MN vs MT
      png(filename=paste0("./CellChat Graphs Analysis Part 5 Comparison/MN vs MT Heatmap.png"), res=120,width = 1920/1,height = 1080/1)
      #Number of interactions
      gg1 <- netVisual_heatmap(cellchat_major_diffMnMt)
      #Strength of interactions
      gg2 <- netVisual_heatmap(cellchat_major_diffMnMt ,measure = "weight")
      mergedgg <- gg1 + gg2
      draw(mergedgg, column_title = "Differential Cell Interaction Hepatic Tumor (red) vs Hepatic Normal (blue)", column_title_gp = gpar(fontsize = 16))
      dev.off()
    }
    
    
    #Differential Number of Interactions and Strength for pTRT
    {
      functiontoComparepTRTInteraction <- function(x,z,iterator)
      {
        if(missing(z)){
          png(filename=paste0("./CellChat Graphs Analysis Part 5 Comparison/Redo v2",deparse(substitute(x))," pTRT Sender Reciever.png"),res=120,width = 1920/1,height = 1080/1)
          #xpd lets the title exist outside of the margins
          par(mfrow=c(2,2),mar = c(0,0,2,0),xpd=TRUE)
          {
            netVisual_diffInteraction(x,
                                      sources.use="pTRT", 
                                      weight.scale = T, 
                                      measure = "count", 
                                      label.edge = F,
                                      title.name=paste0("Differential Number of Interactions as pTRT Sender"))
            netVisual_diffInteraction(x, sources.use="pTRT", weight.scale = T, measure = "weight", label.edge = F,
                                      title.name="Differential Interaction Strength as pTRT Sender")                                
            netVisual_diffInteraction(x, targets.use = "pTRT",weight.scale = T, measure = "count", label.edge = F,title.name="Differential Number of Interactions as pTRT Reciever")
            netVisual_diffInteraction(x, targets.use="pTRT", weight.scale = T, measure = "weight", label.edge = F,title.name="Differential Interaction Strength as pTRT Reciever")
          }
          dev.off()
        } else {
          #This cant work because it still shows the rest of the circles 
          png(filename=paste0("./CellChat Graphs Analysis Part 5 Comparison/Redo v2 ",deparse(substitute(x)),"-",z[[iterator]]," pTRT Sender Reciever.png"),res=120,width = 1920/1,height = 1080/1)
          # #xpd lets the title exist outside of the margins
          par(mfrow=c(2,2),mar = c(1,1,2,1),xpd=TRUE)
          {
            netVisual_diffInteraction(x,
                                      sources.use="pTRT", vertex.label.cex = .9,
                                      weight.scale = T, 
                                      measure = "count", 
                                      label.edge = F,
                                      title.name=paste0("Differential Number of Interactions as pTRT Sender"))
            netVisual_diffInteraction(x, sources.use="pTRT", vertex.label.cex = .9,weight.scale = T, measure = "weight", label.edge = F,title.name="Differential Interaction Strength as pTRT Sender")                                
            netVisual_diffInteraction(x, targets.use = "pTRT",vertex.label.cex = .9,weight.scale = T, measure = "count", label.edge = F,title.name="Differential Number of Interactions as pTRT Reciever")
            netVisual_diffInteraction(x, targets.use="pTRT", vertex.label.cex = .9,weight.scale = T, measure = "weight", label.edge = F,title.name="Differential Interaction Strength as pTRT Reciever")
          }
          dev.off()
        }
      }
      #MT vs PT
      #Red is higher signal in MT while Blue is higher Signal in PT
      functiontoComparepTRTInteraction(cellchat_major_diffPtMt)
      functiontoComparepTRTInteraction(cellchat_major_diffPnMn)
      functiontoComparepTRTInteraction(cellchat_major_diffLnMn)
      functiontoComparepTRTInteraction(cellchat_major_diffLnMt)
      functiontoComparepTRTInteraction(cellchat_major_diffMnMt)
      # #For publication
      # {
      #   cellchat_major_diffPtMt@meta$celltype_major[factor(cellchat_major_diffPtMt@meta$celltype_major) =='CD8 T' ] <- 'Nonreactive CD8 T' 
      #   cellchat_major_diffPtMt <- setIdent(cellchat_major_diffPtMt, ident.use = "celltype_major") # set "cell_type" as default cell identity
      # }
      
      clusterListCategories<- list("CD8","CD4 NKT gdT","DC","Mono Mph","ILC NK","B Plasma")
      clusterList <- list(
        c("pTRT",
          "hC01_CD8_Tn-LEF1",
          "hC02_CD8_Tcm-GPR183",
          "hC03_CD8_Temra-CX3CR1",
          "hC04_CD8_Temra-CX3CR1|RGS1",
          "hC05_CD8_Tem-GZMK",
          "hC06_CD8_Tem-CXCR5",
          "hC07_CD8_IEL-CD160",
          "hC08_CD8_Trm-CD6",
          "hC09_CD8_Trm-XCL1",
          "hC10_CD8_Trm-KLRB1",
          "hC11_CD8_Tex-LAYN",
          "hC12_CD8-HSPA1A",
          "hC13_CD8-MKI67",
          "hC14_CD8_MAIT-SLc4A10"),
        c("pTRT","hC5_CD4_Tn-CCR7",
          "hC16_CD4_Tcm-ANXA1",
          "hC17_CD4_Temra-GNLY",
          "hC18_CD4_Temra-GNLY|RGS1",
          "hC19_CD4_Tn-TCF7",
          "hC20_CD4_Tfh-CXCR5",
          "hC21_CD4_Tem-GZMK",
          "hC22_CD4_Trm-CXCR6",
          "hC23_CD4_Th17-IL23R",
          "hC24_CD4_Th1-CXCL13",
          "hC25_CD4_Treg-FOXP3",
          "hC26_CD4_Treg-IL10",
          "hC27_CD4_Treg-CTLA4",
          "hC28_CD4_HSPA1A",
          "hC29_CD4_MKI67",
          "hC30_NKT-FCGR3A",
          "hC31_NKT-CD27",
          "hC32_γδT-TRDV2"),
        c("pTRT","hC33_cDC-LAMP3",
          "hC34_cDC-MKI67",
          "hC35_cDC1-CLEC9A",
          "hC36_cDC1-CD3E",
          "hC37_cDC2-CD1C",
          "hC38_cDC2-C1QC",
          "hC39_cDC2-FCN1",
          "hC40_cDC2-TIMP1",
          "hC41_cDC2-CD3E",
          "hC42_pDC-LILRA4"),
        c("pTRT","hC43_MonoDC-CLEC10A",
          "hC44_Mono-CD14",
          "hC45_Mono-FCGR3A",
          "hC46_Mono-CD14|FCGR3A",
          "hC47_Mono-CD3E",
          "hC48_Mph-FCN1",
          "hC49_Mph-NLRP3",
          "hC50_Mph-IL1B",
          "hC51_Mph-PLTP",
          "hC52_Mph-CXCL12",
          "hC53_Mph-C1QC",
          "hC54_Mph-SPP1",
          "hC55_Mph-MKI67",
          "hC56_Mph-CD3E",
          "hC57_Mast-CPA3"),
        c("pTRT","hC58_ILc3-KIT",
          "hC59_NK-FCGR3A",
          "hC60_NK-SELL",
          "hC61_NK-FCGR3A|IFNG",
          "hC62_NK-DUSP4",
          "hC63_NK-CD160",
          "hC64_NK-CD160|IFNG",
          "hC65_NK-CD160|ICAM1",
          "hC66_NK-MKI67",
          "hC67_NK-IL7R"),
        c("pTRT","hC68_B-CD79B",
          "hC69_B-CD55",
          "hC70_B-GPR183",
          "hC71_B-HSPA1A",
          "hC72_B-IGLc3",
          "hC73_B-RGS13",
          "hC74_B-MKI67",
          "hC75_B-CD3E",
          "hC76_plasmaB-IGHA1",
          "hC77_plasmaB-IGHA1|IGLC2",
          "hC78_plasmaB-IGHG1")
      )
      #Went in and changed subsetCellChat line 140 approximately to add extra comma so now reads values.new <- values[group.existing.index, group.existing.index, ,drop = FALSE]
      for(indent in 1:length(clusterList))
      {
        sub_cellchat_sub_diffPtMt<-subsetCellChat(object=cellchat_sub_diffPtMt,idents.use=clusterList[[indent]])#cells.use=NULL,thresh=0,
        sub_cellchat_sub_diffPnMn<-subsetCellChat(object=cellchat_sub_diffPnMn,idents.use=clusterList[[indent]])#cells.use=NULL,thresh=0,
        #sub_cellchat_sub_diffLnMn<-subsetCellChat(object=cellchat_sub_diffLnMn,idents.use=clusterList[[indent]])#cells.use=NULL,thresh=0,
        #sub_cellchat_sub_diffLnMt<-subsetCellChat(object=cellchat_sub_diffLnMt,idents.use=clusterList[[indent]])#cells.use=NULL,thresh=0,
        sub_cellchat_sub_diffMnMt<-subsetCellChat(object=cellchat_sub_diffMnMt,idents.use=clusterList[[indent]])#cells.use=NULL,thresh=0,
        
        functiontoComparepTRTInteraction(sub_cellchat_sub_diffPtMt, clusterListCategories,indent)
        print(1)
        functiontoComparepTRTInteraction(sub_cellchat_sub_diffPnMn, clusterListCategories,indent)
        print(2)
        #functiontoComparepTRTInteraction(sub_cellchat_sub_diffLnMn, clusterListCategories,indent)
        #print(3)
        #functiontoComparepTRTInteraction(sub_cellchat_sub_diffLnMt, clusterListCategories,indent)
        #print(4)
        functiontoComparepTRTInteraction(sub_cellchat_sub_diffMnMt, clusterListCategories,indent)
        print(5)
        print(indent)
        #}else{print('Could Not Run this One')}
        
        # functiontoComparepTRTInteraction(cellchat_sub_diffPnMn, clusterList[indent], clusterListCategories[indent])
        # functiontoComparepTRTInteraction(cellchat_sub_diffLnMn, clusterList[indent], clusterListCategories[indent])
        # functiontoComparepTRTInteraction(cellchat_sub_diffLnMt, clusterList[indent], clusterListCategories[indent])
        # functiontoComparepTRTInteraction(cellchat_sub_diffMnMt, clusterList[indent], clusterListCategories[indent])
        
      }
      #Just CD8
      
      png(filename=paste0("./CellChat Graphs Analysis Part 5 Comparison/",deparse(substitute(cellchat_major_diffPtMt))," CD8 Sender Reciever.png"),res=120,width = 1920/1,height = 1080/1)
      x<-cellchat_major_diffPtMt
      #xpd lets the title exist outside of the margins
      par(mfrow=c(2,2),mar = c(0,0,2,0),xpd=TRUE)
      {
        netVisual_diffInteraction(x,
                                  sources.use=c("CD8 T"),
                                  targets.use=c("CD45-","CD4 T","Bcell","pTRT","NKT","NK","Monocyte","MAST","DC","gdT",
                                                "ILC3","Macrophage"),
                                  
                                  weight.scale = T,
                                  measure = "count",
                                  label.edge = T,
                                  title.name=paste0("Differential Number of Interactions as CD8T Sender"))
        netVisual_diffInteraction(x, 
                                  sources.use=c("CD8 T"),
                                  targets.use=c("CD45-","CD4 T","Bcell","pTRT","NKT","NK","Monocyte","MAST","DC","gdT",
                                                "ILC3","Macrophage"),
                                  weight.scale = T, measure = "weight", label.edge = F,
                                  title.name="Differential Interaction Strength as  CD8T Sender")
        netVisual_diffInteraction(x, 
                                  targets.use=c("CD8 T"),
                                  sources.use=c("CD45-","CD4 T","Bcell","pTRT","NKT","NK","Monocyte","MAST","DC","gdT",
                                                "ILC3","Macrophage"),
                                  weight.scale = T, measure = "count", label.edge = T,title.name="Differential Number of Interactions CD8T Reciever")
        netVisual_diffInteraction(cellchat_major_diffPtMt, 
                                  targets.use=c("CD8 T"),
                                  sources.use=c("CD45-","CD4 T","Bcell","pTRT","NKT","NK","Monocyte","MAST","DC","gdT",
                                                "ILC3","Macrophage"),
                                  weight.scale = T, measure = "weight", label.edge = F,title.name="Differential Interaction Strength as CD8T Reciever")
      }
      dev.off()
    }
    
    #Differential Pathway Analysis for pTRT
    {
      #Major
      differentialPathwayMajor <-function(x,z){
        filename=paste0("./CellChat Graphs Analysis Part 5 Comparison Pathway/",deparse(substitute(x)),".png")
        
        png(filename=filename, res=120,width = 1920/1,height = 1080/1)
        gg1 <- netVisual_bubble(x, sources.use = "pTRT", comparison = c(1, 2), max.dataset = 2, title.name = paste0("Increased signaling in ", z), angle.x = 45, remove.isolate = T)
        # Comparing communications on a merged object
        gg2 <- netVisual_bubble(x, sources.use = "pTRT", comparison = c(1, 2), max.dataset = 1, title.name = paste0("Increased signaling in ", z), angle.x = 45, remove.isolate = T)
        # Comparing communications on a merged object
        gg1 + gg2
        
      }
      #PT Vs MT
      differentialPathwayMajor(cellchat_major_diffPtMt,"Hepatic Tumor")
      dev.off()
      
      #Pn Vs Mn
      differentialPathwayMajor(cellchat_major_diffPnMn,"Hepatic Normal")
      dev.off()
      
      #LN vs MN
      differentialPathwayMajor(cellchat_major_diffLnMn,"Hepatic Normal")
      dev.off()
      
      #Ln vs MT
      differentialPathwayMajor(cellchat_major_diffLnMt,"Hepatic Tumor")
      dev.off()
      
      
      #Mn vs MT
      differentialPathwayMajor(cellchat_major_diffMnMt,"Hepatic Tumor")
      dev.off()
      
      #Same Thing but with Chord Diagram but to make this reasonable between pTRT and subtype
      
      differentialPathwaySubChord <-function(x,y,z,source,target,targetName,forFilename){
        # define a positive dataset, i.e., the dataset with positive fold change against the other dataset
        pos.dataset = z
        # define a char name used for storing the results of differential expression analysis
        features.name = pos.dataset
        x <- identifyOverExpressedGenes(x, group.dataset = "datasets", pos.dataset = pos.dataset, features.name = features.name, only.pos = FALSE, thresh.pc = 0.1, thresh.fc = 0.1, thresh.p = 1)
        # map the results of differential expression analysis onto the inferred cell-cell communications to easily manage/subset the ligand-receptor pairs of interest
        net <- netMappingDEG(x, features.name = features.name)
        # extract the ligand-receptor pairs with upregulated ligands in z
        net.up <- subsetCommunication(x, net = net, datasets = names(y)[2],ligand.logFC = 0.2, receptor.logFC = NULL)
        # extract the ligand-receptor pairs with upregulated ligands and upregulated recetptors in NL, i.e.,downregulated in LS
        net.down <- subsetCommunication(x, net = net, datasets = names(y)[1],ligand.logFC = -0.1, receptor.logFC = -0.1)
        gene.up <- extractGeneSubsetFromPair(net.up, x)
        gene.down <- extractGeneSubsetFromPair(net.down, x)
        #browser()
        filename=paste0("./CellChat Graphs Analysis Part 5 Comparison Pathway/Pathway Analysis_",forFilename," ",source,"_",targetName, ".png")
        #png(filename=filename, res=300,width = 1920/1,height = 1080/1)
        png(filename=filename, res=142,width = 2560/1,height = 1440/1)#147
        par(mfrow = c(2,2),mar = c(0,1,0,0), xpd=TRUE)
        #y is object.list
        # netVisual_chord_gene(y[[2]], sources.use = source, targets.use = target, slot.name = 'net', net = net.up, lab.cex = 0.8, small.gap = 3.5, legend.pos.x = 10, title.name = paste0("Up-regulated signaling in ", names(y)[2],"-",source," to ",targetName))
        # netVisual_chord_gene(y[[1]], sources.use = source, targets.use = target, slot.name = 'net', net = net.down, lab.cex = 0.8, small.gap = 3.5, legend.pos.x = 10, title.name = paste0("Down-regulated signaling in ", names(y)[2],"-",source," to ",targetName))
        # netVisual_chord_gene(y[[2]], sources.use = target, targets.use = source, slot.name = 'net', net = net.up, lab.cex = 0.8, small.gap = 3.5, legend.pos.x = 10, title.name = paste0("Up-regulated signaling in ", names(y)[2],"-",targetName," to ",source))
        # netVisual_chord_gene(y[[1]], sources.use = target, targets.use = source, slot.name = 'net', net = net.down, lab.cex = 0.8, small.gap = 3.5, legend.pos.x = 10, title.name = paste0("Down-regulated signaling in ", names(y)[2],"-",targetName," to ",source))
        netVisual_chord_gene(y[[2]], sources.use = source, targets.use = target, slot.name = 'net', net = net.up, lab.cex = 1, small.gap = 3.5, legend.pos.x = 230, legend.pos.y = 165, title.name = paste0("Up-regulated signaling in ", names(y)[2],": ",source," to ",targetName))
        netVisual_chord_gene(y[[1]], sources.use = source, targets.use = target, slot.name = 'net', net = net.down, lab.cex = 1, small.gap = 3.5, legend.pos.x = 5, legend.pos.y = 165,title.name = paste0("Down-regulated signaling in ", names(y)[2],": ",source," to ",targetName))
        netVisual_chord_gene(y[[2]], sources.use = target, targets.use = source, slot.name = 'net', net = net.up, lab.cex = 1, small.gap = 3.5, legend.pos.x = 230 ,legend.pos.y = 20, title.name = paste0("Up-regulated signaling in ", names(y)[2],": ",targetName," to ",source))
        netVisual_chord_gene(y[[1]], sources.use = target, targets.use = source, slot.name = 'net', net = net.down, lab.cex = 1, small.gap = 3.5, legend.pos.x = 5,legend.pos.y =10, title.name = paste0("Down-regulated signaling in ", names(y)[2],": ",targetName," to ",source))
        
        dev.off()
      }
      #removing pTRT from cluster list for the subfunction because want targets to not include autocrine but when subsetting want to still include pTRT
      clusterListforChordPathSubfunction <- list(
        c(
          "hC01_CD8_Tn-LEF1",
          "hC02_CD8_Tcm-GPR183",
          "hC03_CD8_Temra-CX3CR1",
          "hC04_CD8_Temra-CX3CR1|RGS1",
          "hC05_CD8_Tem-GZMK",
          "hC06_CD8_Tem-CXCR5",
          "hC07_CD8_IEL-CD160",
          "hC08_CD8_Trm-CD6",
          "hC09_CD8_Trm-XCL1",
          "hC10_CD8_Trm-KLRB1",
          "hC11_CD8_Tex-LAYN",
          "hC12_CD8-HSPA1A",
          "hC13_CD8-MKI67",
          "hC14_CD8_MAIT-SLc4A10"),
        c("hC5_CD4_Tn-CCR7",
          "hC16_CD4_Tcm-ANXA1",
          "hC17_CD4_Temra-GNLY",
          "hC18_CD4_Temra-GNLY|RGS1",
          "hC19_CD4_Tn-TCF7",
          "hC20_CD4_Tfh-CXCR5",
          "hC21_CD4_Tem-GZMK",
          "hC22_CD4_Trm-CXCR6",
          "hC23_CD4_Th17-IL23R",
          "hC24_CD4_Th1-CXCL13",
          "hC25_CD4_Treg-FOXP3",
          "hC26_CD4_Treg-IL10",
          "hC27_CD4_Treg-CTLA4",
          "hC28_CD4_HSPA1A",
          "hC29_CD4_MKI67",
          "hC30_NKT-FCGR3A",
          "hC31_NKT-CD27",
          "hC32_γδT-TRDV2"),
        c("hC33_cDC-LAMP3",
          "hC34_cDC-MKI67",
          "hC35_cDC1-CLEC9A",
          "hC36_cDC1-CD3E",
          "hC37_cDC2-CD1C",
          "hC38_cDC2-C1QC",
          "hC39_cDC2-FCN1",
          "hC40_cDC2-TIMP1",
          "hC41_cDC2-CD3E",
          "hC42_pDC-LILRA4"),
        c("hC43_MonoDC-CLEC10A",
          "hC44_Mono-CD14",
          "hC45_Mono-FCGR3A",
          "hC46_Mono-CD14|FCGR3A",
          "hC47_Mono-CD3E",
          "hC48_Mph-FCN1",
          "hC49_Mph-NLRP3",
          "hC50_Mph-IL1B",
          "hC51_Mph-PLTP",
          "hC52_Mph-CXCL12",
          "hC53_Mph-C1QC",
          "hC54_Mph-SPP1",
          "hC55_Mph-MKI67",
          "hC56_Mph-CD3E",
          "hC57_Mast-CPA3"),
        c("hC58_ILc3-KIT",
          "hC59_NK-FCGR3A",
          "hC60_NK-SELL",
          "hC61_NK-FCGR3A|IFNG",
          "hC62_NK-DUSP4",
          "hC63_NK-CD160",
          "hC64_NK-CD160|IFNG",
          "hC65_NK-CD160|ICAM1",
          "hC66_NK-MKI67",
          "hC67_NK-IL7R"),
        c("hC68_B-CD79B",
          "hC69_B-CD55",
          "hC70_B-GPR183",
          "hC71_B-HSPA1A",
          "hC72_B-IGLc3",
          "hC73_B-RGS13",
          "hC74_B-MKI67",
          "hC75_B-CD3E",
          "hC76_plasmaB-IGHA1",
          "hC77_plasmaB-IGHA1|IGLC2",
          "hC78_plasmaB-IGHG1")
      )
      #Clean up memory
      rm(cellchat_list,cellchat_ln_global.rds,
         cellchat_ln_major.rds,
         cellchat_ln_sub.rds,
         cellchat_mn_global.rds,
         cellchat_mn_major.rds,
         cellchat_mn_sub.rds,
         cellchat_mt_global.rds,
         cellchat_mt_major.rds,
         cellchat_mt_sub.rds,
         cellchat_pbmc_global.rds,
         cellchat_pbmc_major.rds,
         cellchat_pbmc_sub.rds,
         cellchat_pn_global.rds,
         cellchat_pn_major.rds,
         cellchat_pn_sub.rds,
         cellchat_pt_global.rds,
         cellchat_pt_major.rds,
         cellchat_pt_sub.rds)
      #Parallelize loop
      doFuture::registerDoFuture()
      #foreach(indent = 1:length(clusterList)) %dopar%
      for (indent in 1:length(clusterList))
      {
        #try for catching errors but keep going
        sub_cellchat_sub_diffPtMt <- subsetCellChat(object=cellchat_sub_diffPtMt,idents.use=clusterList[[indent]])#cells.use=NULL,thresh=0,
        sub_cellchat_sub_diffPnMn<-subsetCellChat(object=cellchat_sub_diffPnMn,idents.use=clusterList[[indent]])#cells.use=NULL,thresh=0,
        #sub_cellchat_sub_diffLnMn<-subsetCellChat(object=cellchat_sub_diffLnMn,idents.use=clusterList[[indent]])#cells.use=NULL,thresh=0,
        #sub_cellchat_sub_diffLnMt<-subsetCellChat(object=cellchat_sub_diffLnMt,idents.use=clusterList[[indent]])#cells.use=NULL,thresh=0,
        sub_cellchat_sub_diffMnMt<-subsetCellChat(object=cellchat_sub_diffMnMt,idents.use=clusterList[[indent]])#cells.use=NULL,thresh=0,
        
        
        try(differentialPathwaySubChord(sub_cellchat_sub_diffPtMt, object.list_sub_diffPtMt,"MT","pTRT",clusterListforChordPathSubfunction[[indent]],clusterListCategories[[indent]],deparse(substitute(sub_cellchat_sub_diffPtMt))))
        #try(dev.off())
        print(1)
        try(differentialPathwaySubChord(sub_cellchat_sub_diffPnMn, object.list_sub_diffPnMn, "MN", "pTRT",clusterListforChordPathSubfunction[[indent]],clusterListCategories[[indent]],deparse(substitute(sub_cellchat_sub_diffPnMn))))
        #try(dev.off())
        print(2)
        #try(differentialPathwaySubChord(sub_cellchat_sub_diffLnMn, object.list_sub_diffLnMn,"MN", "pTRT",clusterListforChordPathSubfunction[[indent]],clusterListCategories[[indent]],deparse(substitute(sub_cellchat_sub_diffLnMn))))
        print(3)
        #try(differentialPathwaySubChord(sub_cellchat_sub_diffLnMt, object.list_sub_diffLnMt,"MT", "pTRT",clusterListforChordPathSubfunction[[indent]],clusterListCategories[[indent]],deparse(substitute(sub_cellchat_sub_diffLnMt))))
        print(4)
        try(differentialPathwaySubChord(sub_cellchat_sub_diffMnMt, object.list_sub_diffMnMt,"MT", "pTRT",clusterListforChordPathSubfunction[[indent]],clusterListCategories[[indent]],deparse(substitute(sub_cellchat_sub_diffMnMt))))
        #try(dev.off())
        print(5)
        print(indent)
        #}else{print('Could Not Run this One')}
        
        # functiontoComparepTRTInteraction(cellchat_sub_diffPnMn, clusterList[indent], clusterListCategories[indent])
        # functiontoComparepTRTInteraction(cellchat_sub_diffLnMn, clusterList[indent], clusterListCategories[indent])
        # functiontoComparepTRTInteraction(cellchat_sub_diffLnMt, clusterList[indent], clusterListCategories[indent])
        # functiontoComparepTRTInteraction(cellchat_sub_diffMnMt, clusterList[indent], clusterListCategories[indent])
        
      }
      
      #For Major PT vs MT but now just pTRT and macrophages
      subset_cellchat_MAJOR_diffPtMt <- subsetCellChat(object=cellchat_major_diffPtMt,idents.use=c("pTRT","Macrophage","Monocyte"))#cells.use=NULL,thresh=0,
      try(
        differentialPathwaySubChord(subset_cellchat_MAJOR_diffPtMt,object.list_major_diffPtMt,"MT","pTRT",c("pTRT","Macrophage","Monocyte"),"Macrophage and Monocyte",deparse(substitute(cellchat_major_diffPtMt)))
      )
    }
    
    
  }
}
#Publication Figures
{
  calculate_CellProprtion <- function(x,geneList,assay){
    totalNumberofCells<-ncol(x_expr)
    # Extract the expression matrix for the CD8+ T cells
    x_expr <- LayerData(x, assay = assay,search=NA)
    #x_expr <- GetAssayData(x, assay = "RNA")
    
    # Convert the expression matrix to a numeric matrix
    x_expr <- as.matrix(x_expr)
    
    for(geneIterator in geneList){
      x_expr<-x_expr[,x_expr[geneIterator,]>0]
    }
    
    return(ncol(x_expr)/totalNumberofCells)
  }
  calculate_CellProprtionPerPatient <- function(x,geneList,assay){##ChangGroupingvariableName
    listToReturn<-list()
    for(perPatient in levels(factor(x@meta.data[["patient"]]))){
      
      x_subset<-subset(x,subset=(patient==perPatient))
      # Extract the expression matrix for the CD8+ T cells
      x_expr <- LayerData(x_subset, assay = assay,search=NA)
      totalNumberofCells<-ncol(x_expr)
      #x_expr <- GetAssayData(x, assay = "RNA")
      
      # Convert the expression matrix to a numeric matrix
      x_expr <- as.matrix(x_expr)
      x_expr<-x_expr[rowSums(x_expr)>0,] #remove genes that dont have any values
      
      for(geneIterator in geneList){
        if(geneIterator %in% rownames(x_expr)){
          x_expr<-x_expr[,x_expr[geneIterator,]>0]
          if(!is.null(ncol(x_expr))){x_expr<-x_expr[rowSums(x_expr)>0,]}
        }
        else{
          x_expr<-NULL#matrix(, nrow = 1, ncol = 1)
          break
        }
      }
      if(is.null(x_expr)){
        listToReturn[perPatient]<-0
      }else(
        #tryCatch({listToReturn[perPatient]<-(ncol(x_expr)/totalNumberofCells)*100}, error = function(e){listToReturn[perPatient]<-0})
        listToReturn[perPatient]<-(ncol(as.matrix(x_expr))/totalNumberofCells)*100
      )
      
    }
    listToReturn<-as.data.frame(listToReturn)
    return(listToReturn)
  }
  calculate_CellProprtionExclusivePerPatient <- function(x,geneList,geneListExclusive,assay){##ChangGroupingvariableName
    listToReturn<-list()
    for(perPatient in levels(factor(x@meta.data[["patient"]]))){
      
      x_subset<-subset(x,subset=(patient==perPatient))
      # Extract the expression matrix for the CD8+ T cells
      x_expr <- LayerData(x_subset, assay = assay,search=NA)
      totalNumberofCells<-ncol(x_expr)
      #x_expr <- GetAssayData(x, assay = "RNA")
      
      # Convert the expression matrix to a numeric matrix
      x_expr <- as.matrix(x_expr)
      x_expr<-x_expr[rowSums(x_expr)>0,] #remove genes that dont have any values
      
      #To Cacluate for Positive
      for(geneIterator in geneList){
        if(geneIterator %in% rownames(x_expr)){
          x_expr<-x_expr[,x_expr[geneIterator,]>0]
          if(!is.null(ncol(x_expr))){x_expr<-x_expr[rowSums(x_expr)>0,]}
        }
        else{
          x_expr<-NULL#matrix(, nrow = 1, ncol = 1)
          break
        }
      }
      for(geneIterator in geneListExclusive){
        if(geneIterator %in% rownames(x_expr)){
          x_expr<-x_expr[,x_expr[geneIterator,]==0]
          if(!is.null(ncol(x_expr))){x_expr<-x_expr[rowSums(x_expr)>0,]}
        }
      }#no need for else because if it is not there then it does not express it
      
      if(is.null(x_expr)){
        listToReturn[perPatient]<-0
      }else(
        #tryCatch({listToReturn[perPatient]<-(ncol(x_expr)/totalNumberofCells)*100}, error = function(e){listToReturn[perPatient]<-0})
        listToReturn[perPatient]<-(ncol(as.matrix(x_expr))/totalNumberofCells)*100
      )
      
    }
    listToReturn<-as.data.frame(listToReturn)
    return(listToReturn)
  }
  
  #Myeloid Read in Files
  PredictingPancreaticMetastasis_withoutNA<-readRDS('PredictingPancreaticMetastasis_withoutNA.RDS')
  PredictingPancreaticMetastasis_withoutNA_CD68<-readRDS('PredictingPancreaticMetastasis_withoutNA_CD68.rds')
  cluster2<-readRDS('cluster2.RDS')
  cluster2and3<-readRDS('cluster2and3.RDS')
  merged_Myeloid_Filtered_Liver_umap_subset<-readRDS('merged_Myeloid_Filtered_Liver_umap_subset.rds')
  merged_Myeloid_Filtered_Liver <- readRDS('merged_Myeloid_Filtered_Liver_speck.RDS')
  merged_Myeloid_Filtered_Liver@active.assay<-'RNA'
  merged_Myeloid_Filtered_Colon <- readRDS('merged_Myeloid_Filtered_Colon_speck.RDS')
  merged_Myeloid_Filtered_Colon@active.assay<-'RNA'
  merged_Myeloid_Filtered_Liver_umap_subset_higher_resolution<-readRDS('merged_Myeloid_Filtered_Liver_umap_subset_higher resolution.rds') #this one contains the relevant clustering
  merged_Myeloid_Filtered_Liver.cs_nearSPP1<-readRDS('merged_Myeloid_Filtered_Liver.cs_nearSPP1.RDS')
  merged_Myeloid_Filtered_Liver.cs <- readRDS('merged_Myeloid_Filtered_Liver.cs.RDS')
  merged_Myeloid_Filtered_Colon.cs<-readRDS('merged_Myeloid_Filtered_Colon.cs.RDS')
  merged_Myeloid_Filtered_Liver_nearSPP1_seruat<-merged_Myeloid_Filtered_Liver[,colnames(merged_Myeloid_Filtered_Liver) %in% colnames(merged_Myeloid_Filtered_Liver.cs_nearSPP1)]
  
  #CD8T Read in Files
  ln_cd8 <- readRDS(file = "./Analysis part 2 rev2/ln_cd8.rds")
  mn_cd8 <- readRDS(file = "./Analysis part 2 rev2/mn_cd8.rds")
  mt_cd8 <- readRDS(file = "./Analysis part 2 rev2/mt_cd8.rds")
  pbmc_cd8 <- readRDS(file = "./Analysis part 2 rev2/pbmc_cd8.rds")
  pn_cd8 <- readRDS(file = "./Analysis part 2 rev2/pn_cd8.rds")
  pt_cd8 <- readRDS(file = "./Analysis part 2 rev2/pt_cd8.rds")
  
  # pbmc %<-% readRDS(file = "./Analysis Part 1 Output/pbmc_analysisp1.rds")
  # mt %<-% readRDS(file = "./Analysis Part 1 Output/mt_analysisp1.rds")
  # pt %<-% readRDS(file = "./Analysis Part 1 Output/pt_analysisp1.rds")
  # pn %<-% readRDS(file = "./Analysis Part 1 Output/pn_analysisp1.rds")
  # mn %<-% readRDS(file = "./Analysis Part 1 Output/mn_analysisp1.rds")
  # ln %<-% readRDS(file = "./Analysis Part 1 Output/ln_analysisp1.rds")
  
  pbmc_cd8_TAS <- readRDS(file = "./Analysis Part 2 rev3 Output/pbmc_cd8_TAS_analysisp2r3.rds")
  mt_cd8_TAS <- readRDS(file = "./Analysis Part 2 rev3 Output/mt_cd8_TAS_analysisp2r3.rds")
  pt_cd8_TAS <- readRDS(file = "./Analysis Part 2 rev3 Output/pt_cd8_TAS_analysisp2r3.rds")
  pn_cd8_TAS <- readRDS(file = "./Analysis Part 2 rev3 Output/pn_cd8_TAS_analysisp2r3.rds")
  mn_cd8_TAS <- readRDS(file = "./Analysis Part 2 rev3 Output/mn_cd8_TAS_analysisp2r3.rds")
  ln_cd8_TAS <- readRDS(file = "./Analysis Part 2 rev3 Output/ln_cd8_TAS_analysisp2r3.rds")
  
  makingPerPatientExcel <- function(object_TAS,object,objectname){
    #Frequency of CD8 T that is pTRT per Patient
    write.xlsx(cbind(table(object_TAS$patient),table(object$patient),(100*table(object_TAS$patient)/table(object$patient))) %>% as.data.frame.matrix() %>% setNames(c('Number of pTRT Cells', 'Number of CD8 T cells','% pTRT of CD8')), 
               file="./Analysis Part 7/TAS Frequency Analysis Per Patient.xlsx", 
               sheetName=paste0(objectname," Frequency of pTRT"),
               append=TRUE)#This is why you cannot make it a function
    #Frequency of CD8 T that is particular cluster pTRT per Patient
    write.xlsx(as.data.frame.matrix(sweep(table(object_TAS$celltype_sub,object_TAS$patient), 2, table('CD8T Total' = object$patient), `/`) *100),
               file="./Analysis Part 7/TAS Frequency Analysis Per Patient.xlsx",
               sheetName=paste0(objectname," Freq pTRT of CD8 Per Cluster"),
               append=TRUE)
    write.xlsx(as.data.frame.matrix(sweep(table(object_TAS$celltype_sub,object_TAS$patient), 2, table('CD8T Total' = object_TAS$patient), `/`) *100),
               file="./Analysis Part 7/TAS Frequency Analysis Per Patient.xlsx",
               sheetName=paste0(objectname," Freq per Cluster of pTRT"),
               append=TRUE)
  }
  makingPerPatientExcel(pt_cd8_TAS,pt_cd8,'PT')
  makingPerPatientExcel(pn_cd8_TAS,pn_cd8,'PN')
  makingPerPatientExcel(mt_cd8_TAS,mt_cd8,'MT')
  makingPerPatientExcel(mn_cd8_TAS,mn_cd8,'MN')
  makingPerPatientExcel(ln_cd8_TAS,ln_cd8,'LN')
  makingPerPatientExcel(merge(mt_cd8_TAS,mn_cd8_TAS),merge(mt_cd8,mn_cd8),'Liver MtMn')
  # makingPerPatientExcel(pbmc_cd8_TAS,pbmc_cd8,'pbmc')
  
  #Umaps of Macrophages
  saveImagePNG(DimPlot(merged_Myeloid_Filtered_Liver,label = TRUE,label.box = T,repel = T)+ NoLegend(),'Liver Macrophage Monocyte Clusters UMAP',720,480)
  saveImagePNG(DimPlot(merged_Myeloid_Filtered_Colon,label = TRUE,label.box = T,repel = T)+ NoLegend(),'Colon Macrophage Monocyte Clusters UMAP',720,480)
  
  #Nature medicine paper dot plots including split
  cluster2and3<-calculate_scGSEA(cluster2and3,c('SPP1','IL6R','TGFBR1','VSIG4','CSF1R','CD163'),'IC_Population','RNA')
  saveImagePNG(DotPlot(cluster2and3,features=c('IC_Population','VSIG4','CSF1R','SPP1'),group.by='ClinicalCategory'),'DotPlot of Cluster 2 and 3 Positive Cells IC and Relevant Markers For Publication',520,720)
  saveImagePNG(DotPlot(cluster2and3,features=c('IC_Population'),group.by='ClinicalCategory',scale.min = 0,dot.scale=9)+theme(axis.text.x=element_text(size=18),axis.text.y=element_text(size=18)),'DotPlot of Cluster 2 and 3 Positive Cells Just IC Score For Publication',700,520)
  saveImagePNG(DotPlot(cluster2and3,features=c('VSIG4','CSF1R','SPP1'),group.by='ClinicalCategory',scale.min = 0,dot.scale=9)+theme(axis.text.x=element_text(size=18),axis.text.y=element_text(size=18)),'DotPlot of Cluster 2 and 3 Positive Cells Relevant Markers For Publication',520,720)
  
  
  calculate_CellProprtionPerClinicalCategory <- function(x,geneList,assay){##ChangGroupingvariableName
    listToReturn<-list()
    for(perPatient in levels(factor(x@meta.data[["ClinicalCategory"]]))){
      
      x_subset<-subset(x,subset=(ClinicalCategory==perPatient))
      # Extract the expression matrix for the CD8+ T cells
      x_expr <- LayerData(x_subset, assay = assay,search=NA)
      totalNumberofCells<-ncol(x_expr)
      #x_expr <- GetAssayData(x, assay = "RNA")
      
      # Convert the expression matrix to a numeric matrix
      x_expr <- as.matrix(x_expr)
      
      for(geneIterator in geneList){
        if(geneIterator %in% rownames(x_expr)){
          x_expr<-x_expr[,x_expr[geneIterator,]>0]
        }
        else{
          x_expr<-matrix(, nrow = 0, ncol = 0)
          break
        }
      }
      
      listToReturn[perPatient]<-(ncol(x_expr)/totalNumberofCells)*100
    }
    listToReturn<-as.data.frame(listToReturn)
    return(listToReturn)
  }
  
  frequencyCoExprPaC<-calculate_CellProprtionPerGrouping(cluster2and3,c('SPP1','IL6R','TGFBR1','VSIG4','CSF1R','CD163'),"ClinicalCategory",'RNA')
  
  #Making Ridge Plots of Zhemin Zhang Data set
  merged_Myeloid_Filtered_Liver_umap_subset_higher_resolution_renameIdents<-RenameIdents(merged_Myeloid_Filtered_Liver_umap_subset_higher_resolution,'0'='SPP1 Mph','1'='IC','2'='KC')
  Idents(merged_Myeloid_Filtered_Liver_nearSPP1_seruat)<-Idents(merged_Myeloid_Filtered_Liver_umap_subset_higher_resolution_renameIdents)
  GenesForDifferentiationManifold2<-unique(c('CD163','MRC1','VSIG4','CSF1R','SPP1','IL6R','TGFBR1','IL18','pseudotime'))
  merged_Myeloid_Filtered_Liver_umap_subset_higher_resolution_renameIdents<-calculate_scGSEA(merged_Myeloid_Filtered_Liver_umap_subset_higher_resolution_renameIdents,c('SPP1','IL6R','TGFBR1','VSIG4','CSF1R','CD163'),'DP_Score','RNA')
  merged_Myeloid_Filtered_Liver_umap_subset_higher_resolution_renameIdents<-calculate_scGSEA(merged_Myeloid_Filtered_Liver_umap_subset_higher_resolution_renameIdents,c('SPP1','IL6R','TGFBR1','VSIG4','CSF1R','CD163'),'IC_Score','RNA')
  merged_Myeloid_Filtered_Liver_umap_subset_higher_resolution_renameIdents<-calculate_scGSEA(merged_Myeloid_Filtered_Liver_umap_subset_higher_resolution_renameIdents,c('ACVR1B','ARRB2','BMP4','BMP4','CBX8','CCL2','CCN2','CREB3L1','CREB3L1','DDR2','DICER1','ENG','EP300','F2','F2','F2R','F2R','GLI2','IHH','IL18','INHBA','ITGA2','LARP6','LARP6','LTBP1','MKX','MYB','PDGFB','PDGFRB','PRDX5','RGCC','RGCC','SCX','SCX','SERPINB7','SERPINE1','SERPINF2','SERPINF2','SUCO','TGFB1','TGFB1','TGFB1','TGFB3','TGFB3','UCN','UTS2','VIM','VIM','WNT4','WNT4'),'Positive_Collagen_Regulation','RNA')
  merged_Myeloid_Filtered_Liver_umap_subset_higher_resolution_renameIdents<-calculate_scGSEA(merged_Myeloid_Filtered_Liver_umap_subset_higher_resolution_renameIdents,c('CD5L','VSIG4','SLC1A3','CD163','FOLR2','TIMD4','MARCO','GFRA2','ADRB1','TMEM26','SLC40A1','HMOX1','SLC16A9','VCAM1','SUCNR1'),'KC_score','RNA')
  saveImagePNG(RidgePlot(merged_Myeloid_Filtered_Liver_umap_subset_higher_resolution_renameIdents,features=c('CD68',GenesForDifferentiationManifold2,'Positive_Collagen_Regulation','KC_score')),'Reformatted Liver Trajectory Cells Ridge plot with genes from manifold 2 high resolution with scores For Publication')
  saveImagePNG(RidgePlot(merged_Myeloid_Filtered_Liver_umap_subset_higher_resolution_renameIdents,features=c('CD68',GenesForDifferentiationManifold2,'Positive_Collagen_Regulation','KC_score')),'Reformatted Liver Trajectory Cells Ridge plot with genes from manifold 2 high resolution with scores For Publication Shorter',1440,720)
  saveImagePNG(DimPlot(merged_Myeloid_Filtered_Liver_nearSPP1_seruat),'Reformatted Liver Trajectory Cells UMAP Near SPP1 For Publication',480,360)
  saveImagePNG(DimPlot(merged_Myeloid_Filtered_Liver_nearSPP1_seruat,label.box=TRUE,label=TRUE,pt.size=4)+NoLegend(),'Reformatted Liver Trajectory Cells UMAP Near SPP1 For Publication Label On Top',480,360)
  RidgePlot(merged_Myeloid_Filtered_Liver_umap_subset_higher_resolution_renameIdents,features=c())
  #Remake dot plot comparison with new height
  # merged_Myeloid_Filtered_Colon<-calculate_scGSEA(merged_Myeloid_Filtered_Colon,c('SPP1','IL6R','TGFBR1','VSIG4','CSF1R','CD163'),'IC_Population','RNA')
  # frequencyCoExprColon<-calculate_CellProprtion(merged_Myeloid_Filtered_Colon,c('SPP1','IL6R','TGFBR1','VSIG4','CSF1R','CD163'),'RNA')
  # merged_Myeloid_Filtered_Colon$TissueRT<-'Colon'
  # merged_Myeloid_Filtered_Liver<-calculate_scGSEA(merged_Myeloid_Filtered_Liver,c('SPP1','IL6R','TGFBR1','VSIG4','CSF1R','CD163'),'IC_Population','RNA')
  # merged_Myeloid_Filtered_Liver$TissueRT<-'Liver'
  # frequencyCoExprLiver<-calculate_CellProprtion(merged_Myeloid_Filtered_Liver,c('SPP1','IL6R','TGFBR1','VSIG4','CSF1R','CD163'),'RNA')
  # merged_Myeloid_Filtered_LiverColon<-merge(merged_Myeloid_Filtered_Colon,merged_Myeloid_Filtered_Liver)
  # merged_Myeloid_Filtered_LiverColon$TissueRT<-factor(x = merged_Myeloid_Filtered_LiverColon$TissueRT, levels = c("Liver", "Colon"))
  # saveImagePNG(DotPlot(merged_Myeloid_Filtered_LiverColon,feature = 'IC_Population',group.by = 'TissueRT',scale.min=0,scale.max=100),'Double Positive Score Liver and Colon Dot Plot for Publication',480,480)
  frequencyCoExprLiver<-calculate_CellProprtionPerPatient(merged_Myeloid_Filtered_Liver,c('SPP1','IL6R','TGFBR1','VSIG4','CSF1R','CD163'),'RNA')
  frequencyCoExprColon<-calculate_CellProprtionPerPatient(merged_Myeloid_Filtered_Colon,c('SPP1','IL6R','TGFBR1','VSIG4','CSF1R','CD163'),'RNA')
  frequencyCoExprMT<-calculate_CellProprtionPerPatient(subset(merged_Myeloid_Filtered_Liver,subset=(tissue=="metastasis tumor")),c('SPP1','IL6R','TGFBR1','VSIG4','CSF1R','CD163'),'RNA')
  frequencyCoExprMN<-calculate_CellProprtionPerPatient(subset(merged_Myeloid_Filtered_Liver,subset=(tissue=="metastasis normal")),c('SPP1','IL6R','TGFBR1','VSIG4','CSF1R','CD163'),'RNA')
  frequencyCoExprPT<-calculate_CellProprtionPerPatient(subset(merged_Myeloid_Filtered_Colon,subset=(tissue=="primary tumor")),c('SPP1','IL6R','TGFBR1','VSIG4','CSF1R','CD163'),'RNA')
  frequencyCoExprPN<-calculate_CellProprtionPerPatient(subset(merged_Myeloid_Filtered_Colon,subset=(tissue=="primary normal")),c('SPP1','IL6R','TGFBR1','VSIG4','CSF1R','CD163'),'RNA')
  
  frequencyCoExprLiver_withoutCD163<-calculate_CellProprtionPerPatient(merged_Myeloid_Filtered_Liver,c('SPP1','IL6R','TGFBR1','VSIG4','CSF1R'),'RNA')
  frequencyCoExprColon_withoutCD163<-calculate_CellProprtionPerPatient(merged_Myeloid_Filtered_Colon,c('SPP1','IL6R','TGFBR1','VSIG4','CSF1R'),'RNA')
  frequencyCoExprMT_withoutCD163<-calculate_CellProprtionPerPatient(subset(merged_Myeloid_Filtered_Liver,subset=(tissue=="metastasis tumor")),c('SPP1','IL6R','TGFBR1','VSIG4','CSF1R'),'RNA')
  frequencyCoExprMN_withoutCD163<-calculate_CellProprtionPerPatient(subset(merged_Myeloid_Filtered_Liver,subset=(tissue=="metastasis normal")),c('SPP1','IL6R','TGFBR1','VSIG4','CSF1R'),'RNA')
  frequencyCoExprPT_withoutCD163<-calculate_CellProprtionPerPatient(subset(merged_Myeloid_Filtered_Colon,subset=(tissue=="primary tumor")),c('SPP1','IL6R','TGFBR1','VSIG4','CSF1R'),'RNA')
  frequencyCoExprPN_withoutCD163<-calculate_CellProprtionPerPatient(subset(merged_Myeloid_Filtered_Colon,subset=(tissue=="primary normal")),c('SPP1','IL6R','TGFBR1','VSIG4','CSF1R'),'RNA')
  
  frequencyCoExprLiver_just3 <-calculate_CellProprtionPerPatient(merged_Myeloid_Filtered_Liver,c('SPP1','VSIG4','CSF1R'),'RNA')
  frequencyCoExprColon_just3 <-calculate_CellProprtionPerPatient(merged_Myeloid_Filtered_Colon,c('SPP1','VSIG4','CSF1R'),'RNA')
  
  #writing tables for SPP1 cluster percentages
  write.xlsx(as.data.frame.matrix(table(merged_Myeloid_Filtered_Liver$celltype_sub,merged_Myeloid_Filtered_Liver$patient)),
             file="./Analysis Part 7/SPP1 Counts Analysis Per Patient.xlsx",
             sheetName=paste0("Liver Counts of Macrophage Clusters"),
             append=TRUE)
  write.xlsx(as.data.frame.matrix(table(merged_Myeloid_Filtered_Colon$celltype_sub,merged_Myeloid_Filtered_Colon$patient)), 
             file="./Analysis Part 7/SPP1 Counts Analysis Per Patient.xlsx", 
             sheetName=paste0("Colon Counts of Macrophage Clusters"),
             append=TRUE)#This is why you cannot make it a function
  #SPP1 Exclusive
  frequencyOnlySPP1Liver<-calculate_CellProprtionExclusivePerPatient(merged_Myeloid_Filtered_Liver,c('SPP1'),c('IL6R','TGFBR1','VSIG4','CSF1R','CD163'),'RNA')#'IL6R','TGFBR1',,'CSF1R','CD163'),'RNA'
  frequencyOnlySPP1Colon<-calculate_CellProprtionExclusivePerPatient(merged_Myeloid_Filtered_Colon,c('SPP1'),c('IL6R','TGFBR1','VSIG4','CSF1R','CD163'),'RNA')
  frequencyOnlyVSIG4CSF1RCD163Liver<-calculate_CellProprtionExclusivePerPatient(merged_Myeloid_Filtered_Liver,c('VSIG4','CSF1R','CD163'),c('SPP1','IL6R','TGFBR1'),'RNA')#'IL6R','TGFBR1',,'CSF1R','CD163'),'RNA'
  frequencyOnlyVSIG4CSF1RCD163Colon<-calculate_CellProprtionExclusivePerPatient(merged_Myeloid_Filtered_Colon,c('VSIG4','CSF1R','CD163'),c('SPP1','IL6R','TGFBR1'),'RNA')
  
  #Saving excel files of graph teset results
  merged_Myeloid_Filtered_Liver.cs_graph_test_results <- graph_test(merged_Myeloid_Filtered_Liver.cs,
                                                                    neighbor_graph = "principal_graph",
                                                                    cores = 16)
  write.xlsx(merged_Myeloid_Filtered_Liver.cs_graph_test_results,'./Analysis Part 7/merged_Myeloid_Filtered_Liver.cs_graph_test_results.xlsx')
  
  merged_Myeloid_Filtered_Liver.cs_nearSPP1_graph_test_results <- graph_test(merged_Myeloid_Filtered_Liver.cs_nearSPP1,
                                                                             neighbor_graph = "principal_graph",
                                                                             cores = 16)
  write.xlsx(merged_Myeloid_Filtered_Liver.cs_nearSPP1_graph_test_results,'./Analysis Part 7/merged_Myeloid_Filtered_Liver.cs_nearSPP1_graph_test_results.xlsx')
  merged_Myeloid_Filtered_Colon.cs_graph_test_results <- graph_test(merged_Myeloid_Filtered_Colon.cs,
                                                                    neighbor_graph = "principal_graph",
                                                                    cores = 16)
  write.xlsx(merged_Myeloid_Filtered_Colon.cs_graph_test_results,'./Analysis Part 7/merged_Myeloid_Filtered_Colon.cs_graph_test_results.xlsx')
  
  #Changing Ratios on Colon Trajectories now that it is in supplement
  f1<-plot_cells(merged_Myeloid_Filtered_Colon.cs,
                 color_cells_by = "pseudotime",
                 label_groups_by_cluster=TRUE,
                 group_label_size = 8,
                 trajectory_graph_segment_size = 1.5,
                 label_leaves=FALSE,
                 label_branch_points=FALSE,
                 label_roots = FALSE,
                 #trajectory_graph_color = "grey60",
                 cell_size = 1)
  title<-paste0('Trajectory plot Mph Mono Colon with labeling and coloring from pseudotime for publciation')
  png(file = paste0('./Analysis Part 7/',title,'.png'), width = 1080,height = 1080)
  print(f1)
  dev.off()
  
  #changing Ratios for Liver Trajectory again
  f1<-plot_cells(merged_Myeloid_Filtered_Liver.cs_nearSPP1,
                 genes=c('TGFBR1','IL6R'),#'CLEC4F','TIMD4','CX3CR1',
                 show_trajectory_graph = TRUE,
                 label_cell_groups = FALSE,
                 label_leaves = FALSE,
                 cell_size=1.5)
  title<-paste0('Monocle 3 Clustering Macrophage Monocyte Liver Near SPP1 Subanalysis Psuedotime Plots of TGFBR1 IL6R For Publication')
  png(file = paste0('./Analysis Part 7/',title,'.png'), width = 600,height = 360)
  print(f1)
  dev.off()
  
  f1<-plot_cells(merged_Myeloid_Filtered_Liver.cs_nearSPP1,
                 genes=c('CSF1R','VSIG4'),#'CLEC4F','TIMD4','CX3CR1',
                 show_trajectory_graph = TRUE,
                 label_cell_groups = FALSE,
                 label_leaves = FALSE,
                 cell_size=1.5)
  title<-paste0('Monocle 3 Clustering Macrophage Monocyte Liver Near SPP1 Subanalysis Psuedotime Plots of CSF1R VSIG4 For Publication')
  png(file = paste0('./Analysis Part 7/',title,'.png'), width = 600,height = 360)
  print(f1)
  dev.off()
  
  f1<-plot_cells(merged_Myeloid_Filtered_Liver.cs_nearSPP1,
                 genes=c('CD163','MRC1'),#'CLEC4F','TIMD4','CX3CR1',
                 show_trajectory_graph = TRUE,
                 label_cell_groups = FALSE,
                 label_leaves = FALSE,
                 cell_size=1.5)
  title<-paste0('Monocle 3 Clustering Macrophage Monocyte Liver Near SPP1 Subanalysis Psuedotime Plots of CD163 CD206 For Publication')
  png(file = paste0('./Analysis Part 7/',title,'.png'), width = 600,height = 360)
  print(f1)
  dev.off()
  
  f1<-plot_cells(merged_Myeloid_Filtered_Liver.cs_nearSPP1,
                 genes=c('SPP1'),#'CLEC4F','TIMD4','CX3CR1',
                 show_trajectory_graph = TRUE,
                 label_cell_groups = FALSE,
                 label_leaves = FALSE,
                 cell_size=2)
  title<-paste0('Monocle 3 Clustering Macrophage Monocyte Liver Near SPP1 Subanalysis Psuedotime Plots of SPP1 For Publication')
  png(file = paste0('./Analysis Part 7/',title,'.png'), width = 480,height = 360)
  print(f1)
  dev.off()
  
  saveImagePNG(DimPlot(merged_Myeloid_Filtered_Colon,label = TRUE,label.box = T,repel = T)+ NoLegend(),'Colon Macrophage Monocyte Clusters UMAP for publication',720,720)
  
  #Feature Plot of KC markers for liver and colon
  saveImagePNG(FeaturePlot(merged_Myeloid_Filtered_Liver,features=c('VSIG4','CSF1R'))|FeaturePlot(merged_Myeloid_Filtered_Colon,features=c('VSIG4','CSF1R')),'Feature Plot KC Liver and Colon for Publication',1080,260)
}