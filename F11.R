

### single_cell_THCA

rm(list=ls())
options(stringsAsFactors = F) 
library(Seurat)
library(ggplot2)
library(clustree)
library(cowplot)
library(dplyr)
library(hdf5r)
library(data.table)
library(Seurat)
library(multtest)
library(dplyr)
library(ggplot2)
library(patchwork)
#library(ifnb.SeuratData)
if(!require(harmony))devtools::install_github("immunogenomics/harmony")
library(harmony)

setwd("~/outspace/4T/THCA/")

library(Seurat)

h5_files <- list.files("data/", pattern = "\\.h5$")
seurat_list <- list()
for (h5_file in h5_files) {
  data.path <- paste0("data/", h5_file)
  seurat_data <- Read10X_h5(filename = data.path)
  sample_name <- tools::file_path_sans_ext(basename(h5_file))
  seurat_obj <- CreateSeuratObject(counts = seurat_data,
                                   project = sample_name,
                                   min.features = 200,
                                   min.cells = 3)
  seurat_list <- append(seurat_list, seurat_obj)
}

sample_names <- sub("_.*", "", h5_files)
seurat_combined <- merge(seurat_list[[1]],
                         y = seurat_list[-1],
                         add.cell.ids = sample_names)
print(seurat_combined)
pbmc <- seurat_combined
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-") 
head(pbmc@meta.data,5)
p1 <- VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3) 
print(p1)
pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 25 )   
p2 <- VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3) 
print(p2)

p1|p2


pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
pbmc <- ScaleData(pbmc, features = rownames(pbmc))
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))

DimPlot(pbmc, reduction = "pca")
DimHeatmap(pbmc, dims = 1:15, cells = 500, balanced = TRUE)
pbmc <- JackStraw(pbmc, num.replicate = 100)
pbmc <- ScoreJackStraw(pbmc, dims = 1:20)
JackStrawPlot(pbmc, dims = 1:10)
ElbowPlot(pbmc) 

pbmc <- FindNeighbors(pbmc, dims = 1:15)
pbmc <- FindClusters(pbmc, resolution = 0.8)
pbmc <- RunUMAP(pbmc, dims = 1:15)
pbmc <- RunTSNE(pbmc, dims = 1:15)

p4 <- DimPlot(pbmc, reduction = "umap", group.by = "ident",   pt.size=0.5, label = TRUE,repel = TRUE)+theme(
  axis.line = element_blank(),
  axis.ticks = element_blank(),axis.text = element_blank()
)

p4
DimPlot(pbmc, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()

library(ggplot2)
p <- DimPlot(pbmc, reduction = "umap", label = TRUE, pt.size = 0.5) + 
  NoLegend() +
  theme(text = element_text(size = 32))  # 设置所有文本的字体大小为 12
print(p)

all_marker_gene  <- c("PAX8","KRT18", ##Follicular cells
                      "CD3D","CD8B", # T cells
                      "HIGD1B","CSRP2", #Pericyte
                      "FCER1G","LYZ", # Myeloid cells
                      "VWF","CLDN5", # Endothelial cells
                      "PDGFRA","COL3A1",# Fibroblast
                      "CD79A","CD79B",# B cells
                      "TPSAB1","KIT"   #Mast cells
)

all_marker_gene  <- intersect(all_marker_gene , rownames(pbmc@assays[["RNA"]]))
p_all_markers <- DotPlot(pbmc, features = all_marker_gene, assay = 'RNA') + 
  coord_flip() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 20))  # 设置字体大小为 10

print(p_all_markers)
new.cluster.ids <- c("Follicular cells","T cells","Follicular cells","T cells","T cells",
                     "Pericyte","T cells","Pericyte","Myeloid cells",
                     "Follicular cells","Endothelial cells","Pericyte","Endothelial cells",
                     "Endothelial cells","T cells","Fibroblast","T cells",
                     "B cells","B cells","Myeloid cells","Follicular cells",
                     "B cells","Fibroblast","B cells","B cells",
                     "Mast cells")

names(new.cluster.ids) <- levels(pbmc)
pbmc <- RenameIdents(pbmc, new.cluster.ids)
Idents(pbmc) <- "celltype"
pbmc@meta.data$celltype

DimPlot(pbmc, reduction = "tsne", label = TRUE, pt.size = 0.5) + NoLegend()
DimPlot(pbmc, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
pbmc@meta.data$celltype <- pbmc@active.ident
saveRDS(pbmc,'deal_data/deal_data_7.9.rds')




VlnPlot(pbmc, features = "ISCA1",pt.size = 0)+NoLegend()
FeaturePlot(object = pbmc, features = "ISCA1",slot = "data",label.size = 6,pt.size = 1.2) 
pbmc <- ScaleData(pbmc, features =  rownames(pbmc)) 
heatmap <- DoHeatmap(pbmc, features = c( "ISCA1"), size = 3)
FeaturePlot(object = pbmc, features = c( "ISCA1"),slot = "data",label.size = 6,pt.size = 1.2) 
ggsave("heatmap.pdf", plot = heatmap,height = 6,width = 18)





setwd("~/outspace/4T/THCA/")
pbmc <- readRDS("deal_data/deal_data_7.9.rds")
library(Seurat)
library(dplyr)
gene_name <- "ISCA1"
gene_expression <- FetchData(pbmc, vars = gene_name)
pbmc <- AddMetaData(pbmc, gene_expression, col.name = gene_name)
cell_types <- unique(pbmc@meta.data$celltype)
pbmc@meta.data$ISCA1_group <- NA
for (cell_type in cell_types) {
  cells <- WhichCells(pbmc, idents = cell_type)
  expression_values <- pbmc@meta.data[cells, gene_name]
  high_group <- expression_values > 0
  pbmc@meta.data[cells, "ISCA1_group"] <- ifelse(high_group, "High", "Low")
}
head(pbmc@meta.data)
Idents(pbmc) <- "ISCA1_group"





#  cellratio
stat<-table(pbmc@meta.data$ISCA1_group,pbmc$celltype)
stat<-as.data.frame(stat)
colnames(stat)<-c('celltype','sample','Freq')

p1<-ggplot(data = stat,aes(x = sample,y = Freq, fill = celltype))+
  geom_bar(stat = 'identity',position = 'stack')
p1

color<-c("#00468B","#925E9F","#759EDD","#0099B4","#76D1B1","#42B540","#B8D24D",
         "#EDE447","#FAB158","#FF7777","#FD0000","#AD002A","#AE8691","#CE9573",
         "#756455")
p2<-ggplot(data = stat,aes(x = sample,y = Freq, fill = celltype))+
  geom_bar(stat = 'identity',position = 'fill',width = 0.5)+
  labs(x = "Sample",y = "precent")+
  scale_fill_manual(values = color)+
  guides(fill = guide_legend(ncol = 1,bycol = T,override.aes = list(size = 5)))
p2

mytheme<-theme(panel.background = element_rect(fill = 'white',color = 'black'),
               panel.grid = element_line(color = 'white'),
               axis.title.y = element_text(face = 'bold',color = 'black',size = 14),
               axis.title.x = element_text(face = 'bold',color = 'black',size = 14),
               axis.text.y = element_text(colour = 'black',size = 12),
               axis.text.x = element_text(color = 'black',size = 12,angle = 45,vjust = 0.6),
               legend.title = element_blank(),
               legend.text = element_text(color = 'black',size = 14),
)
p3<-p2+mytheme
p3
ggsave(filename = "stat.pdf",plot = p3,width = 20, height = 15, units = "cm")








pbmc <- readRDS("deal_data/deal_data_7.9.rds")
names(pbmc@meta.data)
pbmc$celltype.group <- paste(pbmc$celltype, pbmc$ISCA1_group, sep = "_")
Idents(pbmc) <- "celltype.group"
unique(pbmc$celltype)
unique(pbmc$celltype.group)
cellfordeg<-levels(pbmc$celltype)
CELLDEG_ALL <- c()

for(i in 1:length(cellfordeg)){
  CELLDEG <- FindMarkers(pbmc, ident.1 = paste0(cellfordeg[i],"_High"), ident.2 = paste0(cellfordeg[i],"_Low"), verbose = FALSE,min.pct = 0.1,logfc.threshold = 0.25)
  write.csv(CELLDEG,paste0("DEGS/",cellfordeg[i],"-0.25",".CSV"))
  CELLDEG$cluster <- cellfordeg[i]
  CELLDEG$gene <- rownames(CELLDEG)
  CELLDEG_ALL <- rbind(CELLDEG_ALL,CELLDEG)
}

CELLDEG_ALL <- subset(CELLDEG_ALL, p_val_adj < 0.05)
write.csv(CELLDEG_ALL,paste0("DEGS/all-0.25",".CSV"))


library(scRNAtoolVis)
Fe_genes <- read.csv("allgene2.csv")
Fe_genes <- Fe_genes$allgene2
result <- subset(CELLDEG_ALL, rownames(CELLDEG_ALL) %in% Fe_genes)
Fe_genes <- rownames(result)
write.csv(result,paste0("DEGS/FE-0-28-0.25",".CSV"))
volcanal_p <- jjVolcano(diffData = result,
                        log2FC.cutoff = 0.25, 
                        size  = 3, #设置点的大小  #00FFCC
                        fontface = 'italic', #设置字体形式
                        aesCol = c('#66FFFF','#FF6666'), #设置点的颜色
                        tile.col = c("red", "blue", "green", "orange", "purple", 
                                     "yellow", "cyan", "magenta", "brown", "gray", 
                                     "darkgreen", "darkblue"),
                        topGeneN = 3)

ggsave(paste("cell-volca.pdf", sep = ""), 
       plot = volcanal_p, width = 10, height = 8)
Idents(pbmc) <- pbmc@meta.data$ISCA1_group
cells_to_keep <- c("High")
unique(Idents(pbmc))
pbmc_High <- subset(pbmc, idents = cells_to_keep)
Idents(pbmc_High) <- "celltype"
DimPlot(pbmc_High, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()



###cellchat

suppressMessages(if(!require(CellChat))devtools::install_github("sqjin/CellChat"))
suppressMessages(if(!require(ggplot2))install.packages('ggplot2'))
suppressMessages(if(!require(patchwork))install.packages('patchwork') )
suppressMessages(if(!require(ggalluvial))install.packages('ggalluvial'))
suppressMessages(if(!require(igraph))install.packages('igraph'))
suppressMessages(if(!require(dplyr))install.packages('dplyr'))
suppressMessages(options(stringsAsFactors = FALSE))
suppressMessages(options(futrue.globlas.Maxsize=2*1024**3))
suppressWarnings(suppressMessages(future::plan("multiprocess", workers = 8)))


time_select_all <- c("High","Low")
unique(pbmc$ISCA1_group)
Idents(pbmc) <- "celltype"
DimPlot(pbmc, reduction = "umap", group.by = "ident", pt.size = 0.5) + theme_void()
for (time_select in time_select_all) {
  class(pbmc)
  data.input <- pbmc[["RNA"]]@data
  meta <- pbmc@meta.data 
  unique(meta$group)
  cell.use <- rownames(meta)[meta$ISCA1_group == time_select] 
  data.input <- data.input[, cell.use]
  meta = meta[cell.use, ]
  unique(meta$celltype)
  cellchat <- createCellChat(object = data.input, meta = meta, group.by = "celltype")
  cellchat <- addMeta(cellchat, meta = meta)
  cellchat <- setIdent(cellchat, ident.use = "celltype") 
  levels(cellchat@idents) 
  groupSize <- as.numeric(table(cellchat@idents)) 
  groupSize
  CellChatDB <- CellChatDB.human 
  showDatabaseCategory(CellChatDB)
  dplyr::glimpse(CellChatDB$interaction)
  CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling")
  cellchat@DB <- CellChatDB.use
  cellchat <- subsetData(cellchat,features = NULL)
  
  cellchat <- identifyOverExpressedGenes(cellchat)
  cellchat <- identifyOverExpressedInteractions(cellchat)
  cellchat <- projectData(cellchat, PPI.human)
  cellchat <- computeCommunProb(cellchat, raw.use = F)
  cellchat <- filterCommunication(cellchat, min.cells = 10)
  df.net <- subsetCommunication(cellchat)
  DT::datatable(df.net)
  cellchat <- computeCommunProbPathway(cellchat)
  cellchat <- aggregateNet(cellchat)
  saveRDS(cellchat,paste0("communicate_data/",time_select,"-cellchat.rds"))
  
}



library(CellChat)
library(patchwork)
library(cowplot)
cellchat.LS <- readRDS("communicate_data/High-cellchat.rds")
cellchat.NL <- readRDS("communicate_data/Low-cellchat.rds")

object.list <- list(ISCA1_Nonactiavted = cellchat.NL, ISCA1_Actiavted = cellchat.LS)
cellchat <- mergeCellChat(object.list, add.names = names(object.list))
## Merge the following slots: 'data.signaling','net', 'netP','meta', 'idents', 'var.features' , 'DB', and 'LR'.
cellchat
levels(cellchat@idents$joint)

gg1 <- compareInteractions(cellchat, show.legend = F, group = c(1,2))
gg2 <- compareInteractions(cellchat, show.legend = F, group = c(1,2), measure = "weight")
gg1 + gg2

gg1 <- rankNet(cellchat, mode = "comparison", stacked = T, do.stat = TRUE)
gg2 <- rankNet(cellchat, mode = "comparison", stacked = F, do.stat = TRUE)
gg1 + gg2

signal_paths_gg1 <- gg1$data$name 
Fe_genes <- read.csv("allgene2.csv")
Fe_genes <- Fe_genes$allgene2
common_genes <- intersect(signal_paths_gg1,Fe_genes)

par(mfrow = c(1,2), xpd=TRUE)
netVisual_diffInteraction(cellchat, weight.scale = T)
netVisual_diffInteraction(cellchat, weight.scale = T, measure = "weight")



weight.max <- getMaxWeight(object.list, attribute = c("idents","count"))
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_circle(object.list[[i]]@net$count, weight.scale = T, label.edge= F, edge.weight.max = weight.max[2], edge.width.max = 12, title.name = paste0("Number of interactions - ", names(object.list)[i]))
}

weight.max <- getMaxWeight(object.list, attribute = c("idents","weight"))
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_circle(object.list[[i]]@net$weight, weight.scale = T, label.edge= F, edge.weight.max = weight.max[2], edge.width.max = 12, title.name = paste0("Strength of interactions - ", names(object.list)[i]))
}

ggsave(paste0("num_diff.pdf"))

gg1 <- netVisual_heatmap(cellchat)
## Do heatmap based on a merged object
gg2 <- netVisual_heatmap(cellchat, measure = "weight")
## Do heatmap based on a merged object
gg1 + gg2

num.link <- sapply(object.list, function(x) {rowSums(x@net$count) +
    colSums(x@net$count)-diag(x@net$count)})
weight.MinMax <- c(min(num.link), max(num.link)) # control the dot size in the different datasets
gg <- list()
for (i in 1:length(object.list)) {
  gg[[i]] <- netAnalysis_signalingRole_scatter(object.list[[i]],
                                               title = names(object.list)[i], 
                                               weight.MinMax = weight.MinMax)
}
patchwork::wrap_plots(plots = gg)





suppressMessages(library(ComplexHeatmap))
i = 1
object.list[[i]] <- netAnalysis_computeCentrality(object.list[[i]])
object.list[[i+1]] <- netAnalysis_computeCentrality(object.list[[i+1]])

pathway.union <- union(object.list[[i]]@netP$pathways, object.list[[i+1]]@netP$pathways)

ht1 = netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "outgoing", signaling = pathway.union, title = names(object.list)[i], width = 18, height = 18)
ht2 = netAnalysis_signalingRole_heatmap(object.list[[i+1]], pattern = "outgoing", signaling = pathway.union, title = names(object.list)[i+1], width = 18, height = 18)

ht_list = ht1 + ht2

ht_list

#incoming
ht1 = netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "incoming", signaling = pathway.union, title = names(object.list)[i], width = 18, height = 18, color.heatmap = "GnBu")
ht2 = netAnalysis_signalingRole_heatmap(object.list[[i+1]], pattern = "incoming", signaling = pathway.union, title = names(object.list)[i+1], width = 18, height = 18, color.heatmap = "GnBu")

ht_list = ht1 + ht2

ht_list

#overall
ht1 = netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "all", signaling = pathway.union, title = names(object.list)[i], width = 18, height = 18, color.heatmap = "OrRd")
ht2 = netAnalysis_signalingRole_heatmap(object.list[[i+1]], pattern = "all", signaling = pathway.union, title = names(object.list)[i+1], width = 18, height = 18, color.heatmap = "OrRd")
# 合并并绘制热图
ht_list = ht1 + ht2

ht_list
saveRDS(object = cellchat, file = "cellchat.rds")


p1 <- netVisual_bubble(cellchat, sources.use = c(1,2,3,4,5,6,7,8), targets.use = c(1,2,3,4,5,6,7,8),  comparison = c(1, 2), angle.x = 45, signaling = c("MK","VISFATIN","VEGF","MIF","CCL","PTN","TGFb","CXCL"))
p1

p1 <- netVisual_bubble(cellchat, sources.use = c("NAMPT","MIF","MDK"), targets.use = c(1,2,3,4,5,6,7,8),  comparison = c(1, 2), angle.x = 45)
p1

unique(cellchat@meta$celltype)

p2 <- netVisual_bubble(cellchat, sources.use = c(1,2,3,4,6,7,8,9,10,11,12,13,14,15), targets.use = c(5),  comparison = c(1, 2), angle.x = 45)
p2

p_single <- netVisual_bubble(cellchat, sources.use = c(4,5), 
                             targets.use = c(1:15),  comparison = c(1, 2), angle.x = 45)
p_single
p_all <- netVisual_bubble(cellchat, comparison = c(1, 2), angle.x = 45)
p_all
levels(cellchat@idents$joint)

ggsave("Compare_LR_p1.pdf", p1, width = 12, height = 50,limitsize = FALSE)
ggsave("Compare_LR_p2.pdf", p2, width = 12, height = 50,limitsize = FALSE)
ggsave("Compare_single.pdf", p_single, width = 20, height = 50,limitsize = FALSE)
ggsave("Compare_ALL.pdf", p_all, width = 100, height = 100,limitsize = FALSE)

gg1 <- netVisual_bubble(cellchat, sources.use = 4, targets.use = c(5:11),  comparison = c(1, 2), max.dataset = 2, title.name = "Increased signaling in LS", angle.x = 45, remove.isolate = T)
## Comparing communications on a merged object
gg2 <- netVisual_bubble(cellchat, sources.use = 4, targets.use = c(5:11),  comparison = c(1, 2), max.dataset = 1, title.name = "Decreased signaling in LS", angle.x = 45, remove.isolate = T)
## Comparing communications on a merged object
p <- gg1 + gg2

ggsave("up-down_LR_bubble.pdf", p, width = 12, height = 30)
p1 <- netVisual_bubble(cellchat, sources.use = c(4,5), targets.use = c(1,2,3,6), comparison = c(1, 2), 
                       max.dataset = 2, title.name = "Increased signaling in TIL", angle.x = 45, remove.isolate = T)
p2 <- netVisual_bubble(cellchat, sources.use = c(4,5), targets.use = c(1,2,3,6), comparison = c(1, 2), 
                       max.dataset = 1, title.name = "Decreased signaling in TIL", angle.x = 45, remove.isolate = T)
pc <- p1 + p2
ggsave("Compare_LR_regulated.pdf", pc, width = 12, height = 5.5)


