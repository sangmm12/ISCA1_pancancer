
### RNA-immune regulatory genes

setwd("D:/R/Pancancer ISCA1/gene-immune regulation/")
cancers <- c("ACC","BLCA","BRCA","CESC","CHOL","COAD","DLBC","ESCA","GBM","HNSC","KICH","KIRC","KIRP","LAML","LGG","LIHC","LUAD","LUSC","MESO","OV","PAAD","PCPG",
             "PRAD","READ","SARC","SKCM","STAD","TGCT","THCA","THYM","UCEC","UCS","UVM")
dat.r <- c()
allgenename <- c("ISCA1")

for (cancer in cancers) {
  
  library(data.table)
  
  
  Inhibitory <- c("ADORA2A","ARG1","BTLA","CD274","CD276","CTLA4","EDNRB","HAVCR2","IDO1","IL10",
                  "IL13","IL4","KIR2DL1","KIR2DL3","LAG3","PDCD1","SLAMF7","TGFB1","TIGIT","VEGFA",
                  "VEGFB","C10orf54","VTCN1","IL12A") 
  
  Stimulaotry <- c("GZMA","BTN3A1","BTN3A2","CCL5","CD27","CD28","CD40","CD40LG","CD70","CD80",
                   "CX3CL1","CXCL10","CXCL9","ENTPD1","HMGB1","ICAM1","ICOS","ICOSLG","IFNA1",
                   "IFNA2","IFNG","IL1A","IL1B","IL2","IL2RA","ITGB2","PRF1","SELP","TLR4","TNF",
                   "TNFRSF14","TNFRSF18","TNFRSF4","TNFRSF9","TNFSF4","TNFSF9" ) 
  
  chemokine <- c("CCL1","CCL2","CCL3","CCL4","CCL5","CCL7","CCL8","CCL11","CCL13","CCL14","CCL15",
                 "CCL16","CCL17","CCL18","CCL19","CCL20","CCL21","CCL22","CCL23","CCL24","CCL25",
                 "CCL26","CCL27","CCL28","CX3CL1","CXCL1","CXCL2","CXCL3","CXCL5","CXCL6","CXCL8",
                 "CXCL9","CXCL10","CXCL11","CXCL12","CXCL13","CXCL14","CXCL16","CXCL17","XCL1","XCL2")
  
  Immunoinhibitor <- c("ADORA2A","BTLA","CD160","CD244","CD274","CD96","CSF1R","CTLA4","HAVCR2",
                       "IDO1","IL10","IL10RB","KDR","KIR2DL1","KIR2DL3","LAG3","LGALS9","PDCD1",
                       "PDCD1LG2","PVRL2","TGFB1","TGFBR1","TIGIT","VTCN1")
  
  
  Immunostimulator <- c( "BTNL2","C10orf54","CD27","CD276","CD28","CD40","CD40LG","CD48","CD70",     
                         "CD80","CD86","CXCL12","CXCR4","ENTPD1","HHLA2","ICOS","ICOSLG","IL2RA",    
                         "IL6","IL6R","KLRC1","KLRK1","LTA","MICB","NT5E", "PVR","RAET1E",   
                         "TMEM173","TMIGD2","TNFRSF13B","TNFRSF13C","TNFRSF14","TNFRSF17","TNFRSF18","TNFRSF25","TNFRSF4",  
                         "TNFRSF8","TNFRSF9","TNFSF13","TNFSF13B","TNFSF14","TNFSF15","TNFSF18","TNFSF4","TNFSF9",   
                         "ULBP1")
  
  MHC <- c("B2M","HLA-A","HLA-B","HLA-C","HLA-DMA","HLA-DMB","HLA-DOA","HLA-DOB","HLA-DPA1","HLA-DPB1",
           "HLA-DQA1","HLA-DQA2","HLA-DQB1","HLA-DRA","HLA-DRB1","HLA-E","HLA-F","HLA-G","TAP1","TAP2",    
           "TAPBP")
  
  receptor <- c("CCR1","CCR2","CCR3","CCR4","CCR5","CCR6","CCR7","CCR8","CCR9","CCR10","CXCR1","CXCR2","CXCR3",
                "CXCR4","CXCR5","CXCR6","XCR1","CX3CR1")
  
  allgene0 <- c(Inhibitory,Stimulaotry,chemokine,Immunoinhibitor,Immunostimulator,MHC,receptor)
  allgene1 <- unique(allgene0)
  
  library(data.table)
  
  immuneregulaygene <-  fread(paste("D:/R/Pancancer ISCA1/Clinical/TCGA_new/",cancer,".txt",sep=''))
  immuneregulaygene <- as.data.frame(immuneregulaygene)
  
  rownames(immuneregulaygene) <- immuneregulaygene[,1]
  immuneregulaygene <- immuneregulaygene[,-1]
  data_im <- immuneregulaygene[match(allgene1,rownames(immuneregulaygene)),]
  
  data_im <-t(data_im)
  library(data.table)
  
  dat <-  fread(paste("D:/R/Pancancer ISCA1/Clinical/TCGA_new/",cancer,".txt",sep=''))
  dat <- as.data.frame(dat)
  rownames(dat) <- dat[,1]
  dat <- dat[,-1]
  data <- t(dat)
  
  data_gene <- data[,match( allgenename,colnames(data))]
  data_gene <- as.data.frame(data_gene)
  all_name <- names(which(table(c(rownames(data_gene),rownames(data_im) ))==2))
  dat_gene <- data_gene[match(all_name,rownames(data_gene)),]
  
  dat_im <- data_im[match(all_name,rownames(data_im)),]
  dat_gene <- as.data.frame(dat_gene)
  colnames(dat_gene) <- "ISCA1"
  dat_im <-as.data.frame(dat_im)
  
  library(psych)
  dat_im <- dat_im[,!is.na(dat_im[1, ])]
  data.corr <- corr.test(dat_im, dat_gene, method="pearson", adjust="fdr")
  data.r <- data.corr$r  # Correlation coefficient
  data.p <- data.corr$p  # p-value
  
  if(length(dat.r)==0){
    dat.r <- data.r
    dat.p <- data.p
  }else {
    dat.r <- cbind(dat.r,data.r)
    dat.p <- cbind(dat.p,data.p)
  }
  
}

colnames(dat.r) <- c("ACC","BLCA","BRCA","CESC","CHOL","COAD","DLBC","ESCA","GBM","HNSC","KICH","KIRC","KIRP","LAML","LGG","LIHC","LUAD","LUSC","MESO","OV","PAAD","PCPG",
                     "PRAD","READ","SARC","SKCM","STAD","TGCT","THCA","THYM","UCEC","UCS","UVM")

colnames(dat.p) <- c("ACC","BLCA","BRCA","CESC","CHOL","COAD","DLBC","ESCA","GBM","HNSC","KICH","KIRC","KIRP","LAML","LGG","LIHC","LUAD","LUSC","MESO","OV","PAAD","PCPG",
                     "PRAD","READ","SARC","SKCM","STAD","TGCT","THCA","THYM","UCEC","UCS","UVM")


write.csv(dat.r,file= paste("Generation/data.r.csv",sep=''),quote=F)
write.csv(dat.p,file= paste("Generation/data.p.csv",sep=''),quote=F)

data.r=dat.r
data.p=dat.p

# Filter data with P < 0.05 for 20
countp <- apply(dat.p, 1, function(row) sum(row < 0.05))
selected_rows <- names(which(countp > 15))



data.p1 <- data.p[match( selected_rows,rownames(data.p)),]
data.r1 <- data.r[match( selected_rows,rownames(data.r)),]

data.r=data.r1
data.p=data.p1
rows_with_na <- is.na(data.r[, 1])
data.r <- data.r[!rows_with_na, ]
data.p <- data.p[!rows_with_na, ]


write.csv(data.p,file= paste("Generation/gene-all.csv",sep=''),quote=F)

data.p[is.na(data.p)] = 2

library(pheatmap)
getSig <- function(dc) {
  print(dc)
  sc <- ' '
  if (dc < 0.0001) {sc <- '****'}
  else if (dc < 0.001){sc <- '***'}
  else if (dc < 0.01){sc <- '**'}
  else if (dc < 0.05) {sc <- '*'}
  else{sc <- ''
  }
  return(sc)
}
sig.mat <- matrix(sapply(data.p, getSig), nrow=nrow(data.p))
str(sig.mat)

paletteLength <- 1000
myColor <- colorRampPalette(c("#36648b", "white", "#e94644"))(paletteLength)

test <- data.r


test[is.na(test)] = 0


myBreaks <- c(seq(min(test), 0, length.out=ceiling(paletteLength/2) + 1), 
              seq(max(test)/paletteLength, max(test), length.out=floor(paletteLength/2)))
class1 <- rep("Inhibitory", times=6)
class2 <- rep("Stimulatory", times=8)
class3 <- rep("Chemokine", times=2)
class4 <- rep("Immunoinhibitor", times=6)
class5 <- rep("Immunostimulator", times=4)
class6 <- rep("MHC", times=5)
class7 <- rep("Receptor", times=2)

class<- c(class1,class2,class3,class4,class5,class6,class7)
annotation_row <- data.frame(class)
rownames(annotation_row) <- rownames(data.p)
pdf(paste("ISCA1-Immune regulation-RNA.pdf",sep=''),width =10,height =7)
pheatmap(data.r,
         color=myColor,
         breaks=myBreaks,
         clustering_method="average", cluster_rows=F,annotation_row = annotation_row, cluster_cols=F,display_numbers=sig.mat
         
)

dev.off()







### CNV-immune regulatory genes

setwd("D:/R/Pancancer ISCA1/gene-immune regulation/Generation")

cancers <- c("ACC","BLCA","BRCA","CESC","CHOL","COAD","DLBC","ESCA","GBM","HNSC","KICH","KIRC","KIRP","LAML","LGG","LIHC","LUAD","LUSC","MESO","OV","PAAD","PCPG",
             "PRAD","READ","SARC","SKCM","STAD","TGCT","THCA","THYM","UCEC","UCS","UVM")
dat.r <- c()
allgene1 <- c("ISCA1")
library(data.table)
gene <-  fread(paste("D:/R/Pancancer ISCA1/gene-immune regulation/Generation/gene-all.csv",sep=''))
gene <- as.data.frame(gene)
generbind=gene$V1
allgene2 <- unique(generbind)
for (cancer in cancers) {
  library(data.table)
  
  dat <-  fread(paste("D:/R/Pancancer ISCA1/Clinical/TCGA_new/",cancer,".txt",sep=''))
  dat <- as.data.frame(dat)
  rownames(dat) <- dat[,1]
  dat <- dat[,-1]
  data <- dat
  # Normal and tumor numbers, characters 14 and 15, 01-09 are cancers, 10-19 are normal, 20-29 are adjacent to cancer
  group=sapply(strsplit(colnames(data),"\\-"), "[", 4)
  group=sapply(strsplit(group,""), "[", 1)
  # Only retain tumor samples
  data = data[,group == 0]
  data <- t(data)
  
  folder_path <- "D:/R/Pancancer ISCA1/CNV-SNV-甲基化/CNV/data/"
  new_file_path <- file.path(folder_path, paste0(cancer, ".cnv.tsv.gz"))
  gz_file <- gzfile(new_file_path, "rt")  # Open the.gz file
  # Use read.table() function to read the data
  dataMethy <- read.table(gz_file, header = TRUE, sep = "\t")  # Take tab delimiter as an example, you can choose the delimiter according to the actual situation
  
  # Delete rows with duplicate names in the 3rd column
  duplicated_rows <- duplicated(dataMethy[,3])
  dataMethy <- dataMethy[!duplicated_rows, ]
  rownames(dataMethy)=dataMethy$symbol
  dataMethy <-dataMethy[,-1:-3]
  colnames(dataMethy) <- gsub("\\.", "-", colnames(dataMethy))
  dataMethy=t(dataMethy)
  
  original_id <-rownames(dataMethy)
  modified_id <- substr(original_id, 1, 15)
  rownames(dataMethy)=modified_id
  
  all_name <- names(which(table(c(rownames(data),rownames(dataMethy)))==2))
  data <- data[match(all_name,rownames(data)),]
  dataMethy <-  dataMethy[match(all_name,rownames( dataMethy)),]
  dat_Methy <- dataMethy[,match( allgene1,colnames(dataMethy))]
  dat_im <- data[,match(allgene2,colnames(data))]
  dat_im <- dat_im[, colSums(is.na(dat_im)) == 0]
  library(psych)
  data.corr <- corr.test(dat_im, dat_Methy, method="pearson", adjust="fdr")
  data.r <- data.corr$r  # Correlation coefficient
  data.p <- data.corr$p  # p-value
  data.r <- data.r[match(allgene2,rownames(data.r)),]
  data.p <- data.p[match(allgene2,rownames(data.p)),]
  if(length(dat.r)==0){
    dat.r <- data.r
    dat.p <- data.p
  }else {
    dat.r <- cbind(dat.r,data.r)
    dat.p <- cbind(dat.p,data.p)
  }
  
}
colnames(dat.r) <- c("ACC","BLCA","BRCA","CESC","CHOL","COAD","DLBC","ESCA","GBM","HNSC","KICH","KIRC","KIRP","LAML","LGG","LIHC","LUAD","LUSC","MESO","OV","PAAD","PCPG",
                     "PRAD","READ","SARC","SKCM","STAD","TGCT","THCA","THYM","UCEC","UCS","UVM")
colnames(dat.p) <- c("ACC","BLCA","BRCA","CESC","CHOL","COAD","DLBC","ESCA","GBM","HNSC","KICH","KIRC","KIRP","LAML","LGG","LIHC","LUAD","LUSC","MESO","OV","PAAD","PCPG",
                     "PRAD","READ","SARC","SKCM","STAD","TGCT","THCA","THYM","UCEC","UCS","UVM")
dat.r <- dat.r[!is.na(rownames(dat.r)), ]
dat.p <- dat.p[!is.na(rownames(dat.p)), ]
data.r=dat.r
data.p=dat.p
library(pheatmap)
getSig <- function(dc) {
  print(dc)
  sc <- ' '
  if (dc < 0.0001) {sc <- '****'}
  else if (dc < 0.001){sc <- '***'}
  else if (dc < 0.01){sc <- '**'}
  else if (dc < 0.05) {sc <- '*'}
  else{sc <- ''
  }
  return(sc)
}
sig.mat <- matrix(sapply(data.p, getSig), nrow=nrow(data.p))
str(sig.mat)
paletteLength <- 1000
myColor <- colorRampPalette(c("#36648b", "white", "#e94644"))(paletteLength)
test <- data.r
test[is.na(test)] = 0
myBreaks <- c(seq(min(test), 0, length.out=ceiling(paletteLength/2) + 1), 
              seq(max(test)/paletteLength, max(test), length.out=floor(paletteLength/2)))

class1 <- rep("Inhibitory", times=6)
class2 <- rep("Stimulatory", times=8)
class3 <- rep("Chemokine", times=2)
class4 <- rep("Immunoinhibitor", times=6)
class5 <- rep("Immunostimulator", times=4)
class6 <- rep("MHC", times=5)
class7 <- rep("Receptor", times=2)
class<- c(class1,class2,class3,class4,class5,class6,class7)
annotation_row <- data.frame(class)
rownames(annotation_row) <- rownames(data.p)

pdf(paste("cor-immune_CNV_NEW.pdf",sep=''),width =10,height =7)
pheatmap(data.r,
         color=myColor,
         breaks=myBreaks,
         clustering_method="average", cluster_rows=F,annotation_row = annotation_row, cluster_cols=F,display_numbers=sig.mat
)
dev.off()




### Methylation-immune regulatory genes

setwd("D:/R/Pancancer ISCA1/gene-immune regulation/generation")
# 
cancers <- c("ACC","BLCA","BRCA","CESC","CHOL","COAD","DLBC","ESCA","GBM","HNSC","KICH","KIRC","KIRP","LAML","LGG","LIHC","LUAD","LUSC","MESO","OV","PAAD","PCPG",
             "PRAD","READ","SARC","SKCM","STAD","TGCT","THCA","THYM","UCEC","UCS","UVM")
dat.r <- c()

allgene1 <- c("ISCA1")

library(data.table)

Inhibitory <- c("ADORA2A","ARG1","BTLA","CD274","CD276","CTLA4","EDNRB","HAVCR2","IDO1","IL10",
                "IL13","IL4","KIR2DL1","KIR2DL3","LAG3","PDCD1","SLAMF7","TGFB1","TIGIT","VEGFA",
                "VEGFB","C10orf54","VTCN1","IL12A") 

Stimulatory <- c("GZMA","BTN3A1","BTN3A2","CCL5","CD27","CD28","CD40","CD40LG","CD70","CD80",
                 "CX3CL1","CXCL10","CXCL9","ENTPD1","HMGB1","ICAM1","ICOS","ICOSLG","IFNA1",
                 "IFNA2","IFNG","IL1A","IL1B","IL2","IL2RA","ITGB2","PRF1","SELP","TLR4","TNF",
                 "TNFRSF14","TNFRSF18","TNFRSF4","TNFRSF9","TNFSF4","TNFSF9" ) 

Chemokine <- c("CCL1","CCL2","CCL3","CCL4","CCL5","CCL7","CCL8","CCL11","CCL13","CCL14","CCL15",
               "CCL16","CCL17","CCL18","CCL19","CCL20","CCL21","CCL22","CCL23","CCL24","CCL25",
               "CCL26","CCL27","CCL28","CX3CL1","CXCL1","CXCL2","CXCL3","CXCL5","CXCL6","CXCL8",
               "CXCL9","CXCL10","CXCL11","CXCL12","CXCL13","CXCL14","CXCL16","CXCL17","XCL1","XCL2")

Immunoinhibitor <- c("ADORA2A","BTLA","CD160","CD244","CD274","CD96","CSF1R","CTLA4","HAVCR2",
                     "IDO1","IL10","IL10RB","KDR","KIR2DL1","KIR2DL3","LAG3","LGALS9","PDCD1",
                     "PDCD1LG2","PVRL2","TGFB1","TGFBR1","TIGIT","VTCN1")

Immunostimulator <- c( "BTNL2","C10orf54","CD27","CD276","CD28","CD40","CD40LG","CD48","CD70",     
                       "CD80","CD86","CXCL12","CXCR4","ENTPD1","HHLA2","ICOS","ICOSLG","IL2RA",    
                       "IL6","IL6R","KLRC1","KLRK1","LTA","MICB","NT5E", "PVR","RAET1E",   
                       "TMEM173","TMIGD2","TNFRSF13B","TNFRSF13C","TNFRSF14","TNFRSF17","TNFRSF18","TNFRSF25","TNFRSF4",  
                       "TNFRSF8","TNFRSF9","TNFSF13","TNFSF13B","TNFSF14","TNFSF15","TNFSF18","TNFSF4","TNFSF9",   
                       "ULBP1")

MHC <- c("B2M","HLA-A","HLA-B","HLA-C","HLA-DMA","HLA-DMB","HLA-DOA","HLA-DOB","HLA-DPA1","HLA-DPB1",
         "HLA-DQA1","HLA-DQA2","HLA-DQB1","HLA-DRA","HLA-DRB1","HLA-E","HLA-F","HLA-G","TAP1","TAP2",    
         "TAPBP")

Receptor <- c("CCR1","CCR2","CCR3","CCR4","CCR5","CCR6","CCR7","CCR8","CCR9","CCR10","CXCR1","CXCR2","CXCR3",
              "CXCR4","CXCR5","CXCR6","XCR1","CX3CR1")

allgene0 <- c(Inhibitory,Stimulatory,Chemokine,Immunoinhibitor,Immunostimulator,MHC,Receptor)
allgene2 <- unique(allgene0)

for (cancer in cancers) {
  
  library(data.table)
  
  dat <-  fread(paste("D:/R/Pancancer ISCA1/clinical/TCGA_new/",cancer,".txt",sep=''))
  dat <- as.data.frame(dat)
  rownames(dat) <- dat[,1]
  dat <- dat[,-1]
  data <- dat
  #Normal and tumor number, 14th and 15th characters, 01-09 are cancer, 10-19 are normal, 20-29 are adjacent to cancer
  group=sapply(strsplit(colnames(data),"\\-"), "[", 4)
  group=sapply(strsplit(group,""), "[", 1)
  #Only retain tumor samples
  data = data[,group == 0]
  
  methy <- fread(paste("D:/R/Pancancer ISCA1/linear/new-methy/ISCA1-methy/",cancer,"-Methy450.csv",sep=''))
  methy <- data.frame(methy)
  colnames(methy) <- gsub("\\.", "-", colnames(methy))
  
  all_name <- names(which(table(c(colnames(data),colnames(methy)))==2))
  
  data <- data[,match(all_name,colnames(data))]
  dat_Methy <-  methy[,match(all_name,colnames(methy))]
  dat_im <- data[match( allgene2,rownames(data)),]
  cols_with_na <- is.na(dat_im[, 1])
  dat_im <- dat_im[!cols_with_na,]
  
  library(psych)
  dat_Methy <- t(dat_Methy)
  dat_im <- t(dat_im)
  
  dat_Methy=as.numeric(dat_Methy)
  data.corr <- corr.test(dat_im, dat_Methy, method="pearson", adjust="fdr")
  data.r <- data.corr$r  # Correlation coefficient
  data.p <- data.corr$p  # p-value
  
  if(length(dat.r)==0){
    dat.r <- data.r
    dat.p <- data.p
  }else {
    dat.r <- cbind(dat.r,data.r)
    dat.p <- cbind(dat.p,data.p)
  }
  
}
colnames(dat.r) <- c("ACC","BLCA","BRCA","CESC","CHOL","COAD","DLBC","ESCA","GBM","HNSC","KICH","KIRC","KIRP","LAML","LGG","LIHC","LUAD","LUSC","MESO","OV","PAAD","PCPG",
                     "PRAD","READ","SARC","SKCM","STAD","TGCT","THCA","THYM","UCEC","UCS","UVM")
colnames(dat.p) <- c("ACC","BLCA","BRCA","CESC","CHOL","COAD","DLBC","ESCA","GBM","HNSC","KICH","KIRC","KIRP","LAML","LGG","LIHC","LUAD","LUSC","MESO","OV","PAAD","PCPG",
                     "PRAD","READ","SARC","SKCM","STAD","TGCT","THCA","THYM","UCEC","UCS","UVM")

dat.r <- dat.r[!is.na(rownames(dat.r)), ]
dat.p <- dat.p[!is.na(rownames(dat.p)), ]
write.csv(dat.r,file= paste(cancer,"rbin_data.r.csv",sep=''),quote=F)
write.csv(dat.p,file= paste(cancer,"rbin_data.p.csv",sep=''),quote=F)

data.r=dat.r
data.p=dat.p

data.p[is.na(data.p)] = 2
library(pheatmap)
getSig <- function(dc) {
  print(dc)
  sc <- 
    if (dc < 0.0001) {sc <- '****'}
  else if (dc < 0.001){sc <- '***'}
  else if (dc < 0.01){sc <- '**'}
  else if (dc < 0.05) {sc <- '*'}
  else{sc <- ''
  }
  return(sc)
}
sig.mat <- matrix(sapply(data.p, getSig), nrow=nrow(data.p))
str(sig.mat)
paletteLength <- 1000
myColor <- colorRampPalette(c("#36648b", "white", "#e94644"))(paletteLength)
test <- data.r
test[is.na(test)] = 0
myBreaks <- c(seq(min(test), 0, length.out=ceiling(paletteLength/2) + 1), 
              seq(max(test)/paletteLength, max(test), length.out=floor(paletteLength/2)))
class1 <- rep("Inhibitory", times=23)
class2 <- rep("Stimulatory", times=36)
class3 <- rep("Chemokine", times=37)
class4 <- rep("Immunoinhibitor", times=9)
class5 <- rep("Immunostimulator", times=26)
class6 <- rep("MHC", times=21)
class7 <- rep("Receptor", times=17)
class<- c(class1,class2,class3,class4,class5,class6,class7)


annotation_row <- data.frame(class)
rownames(annotation_row) <- rownames(data.p)

pdf(paste("all-immune_methy.pdf",sep=''),width =10,height =40)
pheatmap(data.r,
         color=myColor,
         breaks=myBreaks,
         clustering_method="average", cluster_rows=F,annotation_row = annotation_row, cluster_cols=F,display_numbers=sig.mat
)
dev.off()



### Isoform-immune regulatory genes


setwd("D:/R/Pancancer ISCA1/gene-immune regulation/Generation")

cancers <- c("ACC","BLCA","BRCA","CESC","CHOL","COAD","DLBC","ESCA","GBM","HNSC","KICH","KIRC","KIRP","LAML","LGG","LIHC","LUAD","LUSC","OV","PAAD","PCPG",
             "PRAD","READ","SARC","SKCM","STAD","TGCT","THCA","THYM","UCEC","UCS")

Isoformall <- c("ENST00000311534","ENST00000326094","ENST00000375991")

for (Isoform_single in Isoformall) {
  dat.r <- c()
  allgene1 <- c("ISCA1")
  library(data.table)
  gene <-  fread(paste("D:/R/Pancancer ISCA1/gene-immune regulation/Generation/gene-all.csv",sep=''))
  gene <- as.data.frame(gene)
  generbind=gene$V1
  allgene2 <- unique(generbind)
  for (cancer in cancers) {
    library(data.table)
    dat <-  fread(paste("D:/R/Pancancer ISCA1/Clinical/TCGA_new/",cancer,".txt",sep=''))
    dat <- as.data.frame(dat)
    rownames(dat) <- dat[,1]
    dat <- dat[,-1]
    data <- dat
    # Normal and tumor numbers, characters 14 and 15, 01-09 are cancers, 10-19 are normal, 20-29 are adjacent to cancer
    group=sapply(strsplit(colnames(data),"\\-"), "[", 4)
    group=sapply(strsplit(group,""), "[", 1)
    # Only retain tumor samples
    data = data[,group == 0]
    data <- t(data)
    Isoform <- read.csv(paste0("D:/R/Pancancer ISCA1/Isoform/",Isoform_single,"_transcript_pancan_dist.CSV"))
    Isoform <- subset(Isoform, tissue == cancer)
    rownames(Isoform) <- Isoform$sample
    all_name <- names(which(table(c(rownames(data),Isoform$sample))==2))
    data <- data[match(all_name,rownames(data)),]
    Isoform <-  Isoform[match(all_name,Isoform$sample),]
    dat_im <- data[,match(allgene2,colnames(data))]
    dat_im <- dat_im[, colSums(is.na(dat_im)) == 0]
    Isoform <- Isoform$tpm
    library(psych)
    data.corr <- corr.test(dat_im, Isoform, method="pearson", adjust="fdr")
    data.r <- data.corr$r  # Correlation coefficient
    data.p <- data.corr$p  # p-value
    data.r <- data.r[match(allgene2,rownames(data.r)),]
    data.p <- data.p[match(allgene2,rownames(data.p)),]
    if(length(dat.r)==0){
      dat.r <- data.r
      dat.p <- data.p
    }else {
      dat.r <- cbind(dat.r,data.r)
      dat.p <- cbind(dat.p,data.p)
    }
    
  }
  colnames(dat.r) <- cancers
  colnames(dat.p) <- cancers
  dat.r <- dat.r[!is.na(rownames(dat.r)), ]
  dat.p <- dat.p[!is.na(rownames(dat.p)), ]
  data.r=dat.r
  data.p=dat.p
  library(pheatmap)
  getSig <- function(dc) {
    print(dc)
    sc <- ' '
    if (dc < 0.0001) {sc <- 1}
    else if (dc < 0.001){sc <- 1}
    else if (dc < 0.01){sc <- 1}
    else if (dc < 0.05) {sc <- 1}
    else{sc <- ''
    }
    return(sc)
  }
  sig.mat <- matrix(sapply(data.p, getSig), nrow=nrow(data.p))
  str(sig.mat)
  paletteLength <- 1000
  myColor <- colorRampPalette(c("#36648b", "white", "#e94644"))(paletteLength)
  test <- data.r
  test[is.na(test)] = 0
  myBreaks <- c(seq(min(test), 0, length.out=ceiling(paletteLength/2) + 1), 
                seq(max(test)/paletteLength, max(test), length.out=floor(paletteLength/2)))
  colnames(sig.mat) <- cancers
  rownames(sig.mat) <-  rownames(data.p)
  less_than_0.05 <- colSums(data.p < 0.05)
  less_than_0.05 <- as.data.frame(t(less_than_0.05))
  data.p <- rbind(data.p, less_than_0.05)
  row.names(data.p)[nrow(data.p)] <- "genes_num_0.05"
  write.csv(data.p,file= paste(Isoform_single,"_","data.p.csv",sep=''),quote=F)
  class1 <- rep("Inhibitory", times=6)
  class2 <- rep("Stimulatory", times=8)
  class3 <- rep("Chemokine", times=2)
  class4 <- rep("Immunoinhibitor", times=6)
  class5 <- rep("Immunostimulator", times=4)
  class6 <- rep("MHC", times=5)
  class7 <- rep("Receptor", times=2)
  class<- c(class1,class2,class3,class4,class5,class6,class7)
  annotation_row <- data.frame(class)
  rownames(annotation_row) <- rownames(data.p)
  pdf(paste("cor-immune_isoform_",Isoform_single,".pdf",sep=''),width =10,height =7)
  pheatmap(data.r,
           color=myColor,
           breaks=myBreaks,
           clustering_method="average", cluster_rows=F,annotation_row = annotation_row, cluster_cols=F,display_numbers=sig.mat
  )
  dev.off()
}


