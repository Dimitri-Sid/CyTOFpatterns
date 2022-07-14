library(readxl)
library(ggplot2)
library(reshape2)
library(dplyr)
library(pheatmap)
library(viridis)
library(tidyr)
library(FlowSOM)
library(flowCore)
library(premessa)

#####################
#### Data Processing
#####################

metaDataFile='ViralHCC_metadata_forDS.xlsx'
panelDataFile='ViralHCC_panel.xlsx'
dataDirectory='Data'
md <- read_excel(metaDataFile)
rownames(md) = md$sample_id
#Load FCS files
fcs_raw <- read.flowSet(md$file_name, path = dataDirectory, transformation = FALSE, truncate_max_range = FALSE)
sample_ids <- rep(md$sample_id, fsApply(fcs_raw, nrow))
panel <- read_excel(panelDataFile)
panel_fcs <- pData(parameters(fcs_raw[[1]]))
# use metal+isotope as mapping between panel from xlsx and panel from the fcs files
rownames(panel_fcs) = panel_fcs$name
pData(parameters(fcs_raw[[1]])) <- panel_fcs

## arcsinh transformation and column subsetting
fcs <- fsApply(fcs_raw, function(x, cofactor = 5){
  colnames(x) <- panel_fcs$name
  expr <- exprs(x)
  expr <- asinh(expr/cofactor)
  exprs(x) <- expr
  x
})
fcs <- list('fcs'=fcs,'sample_ids'=sample_ids,'meta_data'=md)
aggr_flowframe <- fsApply(fcs$fcs, exprs)

markers <- paste0(panel$Metal,panel$Isotope,"Di") # add Di to match metal isotope nomenclature
aggr_flowframe <- aggr_flowframe[,markers] #subset flowframe aggregate to just the variables/antibodies of interest for analysis 
colnames(aggr_flowframe) <- panel$Antigen
barcodes <- which(colnames(aggr_flowframe) == "CD45")
aggr_flowframe <- aggr_flowframe[,-barcodes] #remove CD45 since it is merely used for barcodes
texpr<-t(aggr_flowframe) #transpose
rm(aggr_flowframe)

## Develop complete metadata file 
ids <- data.frame("samples_ids" = fcs$sample_ids)
md <- data.frame(fcs$meta_data)
md <- data.frame(lapply(md, as.character), stringsAsFactors=FALSE)
cells <- length(rownames(ids))
total_md <-data.frame(matrix(nrow=cells, ncol=length(colnames(md))))
colnames(total_md) <- colnames(md)
# fill data in total md df
total_md$sample_id <- ids$samples_ids
for (sample in unique(ids$samples_ids)){
  for (variable in colnames(total_md)){
    total_md[total_md$sample_id==sample,variable] <- md[md$sample_id==sample,variable]
  }
}

## Create cell data set object
# requirements for CDS
fdat<-matrix(rownames(texpr))
rownames(fdat)<-rownames(texpr)
colnames(fdat) <- c('gene_short_name')
pdat<-as.matrix(total_md)
rownames(pdat)<-c(1:(numCellsPerPatient*length(sample_ids)))
# need to detach flowSOM before using monocle3
detach(package:FlowSOM)
library(monocle3)
cds <- new_cell_data_set(texpr,
                         cell_metadata = pdat,
                         gene_metadata = fdat)
# randomly sample cells to reduce run time 
cds <- cds[,sample(colnames(cds),200000)]
cds <- cds[,Matrix::colSums(exprs(cds)) != 0] #use if rows or cols with all zeros
cds <- estimate_size_factors(cds)
cds <- preprocess_cds(cds, num_dim = 10, norm_method = "none", scaling = FALSE) #log is default
plot_pc_variance_explained(cds)
cds <- align_cds(cds,alignment_group = "batch")
cds <- reduce_dimension(cds, max_components=2, reduction_method = c("UMAP")) #make sure it is UMAP
cds = cluster_cells(cds, resolution=1e-5)
pData(cds)$clusters <- ""
pData(cds)$clusters <- clusters(cds)
pData(cds)$partitions <- "" 
pData(cds)$partitions <- partitions(cds)

################
#### Functions
################

plotCells <- function(object,marker,title,save){ # plots UMAPs or saves UMAP file
  if(save=="yes"){
    tiff(paste0(title," ",marker,".tiff"))
    print(ggplot(data.frame(reducedDims(object)$UMAP),  aes(x = X1, y = X2, color = t(exprs(object))[,marker] ))+
            geom_point(size = 1) +
            theme_classic()+
            theme(panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank())+
            scale_color_viridis(discrete = FALSE, option = "D")+
            scale_fill_viridis(discrete = FALSE)+ggtitle(marker)+ylab("UMAP2") + xlab("UMAP1")+labs(color = "Expression"))
    dev.off()
  } else {
    ggplot(data.frame(reducedDims(object)$UMAP),  aes(x = X1, y = X2, color = t(exprs(object))[,marker] ))+
      geom_point(size = 1) +
      theme_classic()+
      theme(panel.grid.major = element_blank(),
            panel.grid.minor = element_blank())+
      scale_color_viridis(discrete = FALSE, option = "D")+
      scale_fill_viridis(discrete = FALSE)+ggtitle(title)+ylab("UMAP2") + xlab("UMAP1")+labs(color = "Expression")
  }
}

trajInfer <- function(object,title){ #learns pseudotime, user picks start node, saves figure
  object <-cluster_cells(object,resolution=1e-5) 
  object<-learn_graph(object)
  object<-order_cells(object)
  pData(object)$pseudotime <- pseudotime(object)
  tiff(paste0(title,".tiff"))
  print(plot_cells(object, color_cells_by = "pseudotime",label_cell_groups=FALSE,label_leaves=FALSE,label_branch_points=FALSE, graph_label_size=1.5))
  dev.off()
  pData(object)$subset <- title
  return(object)
}

#################
#### Plot & save
#################

plotCells(cds,"CD3","PDAC Th","no")
plotCells(cds,"CD4","PDAC Th","no")
plotCells(cds,"CD4","PDAC Th","no")

# plot umap and violins
plot_cells(cds, group_cells_by="partition", color_cells_by="virustype", cell_size = 0.3,labels_per_group = 0,label_cell_groups = FALSE, group_label_size = 6)
plot_cells(cds, group_cells_by="partition", color_cells_by="response", cell_size = 0.3,labels_per_group = 0,label_cell_groups = FALSE, group_label_size = 6)
plot_cells(cds, group_cells_by="partition", color_cells_by="batch", cell_size = 0.3,labels_per_group = 0,label_cell_groups = FALSE, group_label_size = 6)
for (i in 1:length(rownames(cds))){
  print(i)
  png(paste0("./figures/",rownames(cds)[i],".png"))
  print(ggplot(data.frame(reducedDims(cds)$UMAP),  aes(x = X1, y = X2, color = t(exprs(cds))[,rownames(cds[i])] ))+
          geom_point(size = 1) +
          theme_bw()+
          theme(panel.grid.major = element_blank(),
                panel.grid.minor = element_blank())+
          scale_color_viridis(discrete = FALSE, option = "D")+
          scale_fill_viridis(discrete = FALSE)+
          xlab("UMAP1")+ylab("UMAP2")+labs(color = rownames(cds[i])))
  dev.off()
}
dev.off()

save(cds,file="hcc_cds.rda")
