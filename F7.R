
### RNA-immune cell

setwd("D:/R/Pancancer ISCA1/8cell")

cancers <- c("ACC","BLCA","BRCA","CESC","CHOL","COAD","DLBC","ESCA","GBM","HNSC","KICH","KIRC","KIRP","LAML","LGG","LIHC","LUAD","LUSC","MESO","OV","PAAD","PCPG",
             "PRAD","READ","SARC","SKCM","STAD","TGCT","THCA","THYM","UCEC","UVM","UCS")
Isoform_single <- "RNA"
algorithms <- c("cibersort")
for (algorithm in algorithms) {
  dat.r <- c()
  allgene1 <- c("ISCA1")
  for (cancer in cancers) {
    library(data.table)
    CELL <-  fread(paste("8 kinds of algorithms single/TCGA-",cancer,"/im_",algorithm,".csv",sep=''))
    CELL <- as.data.frame(CELL)
    original_id <-CELL[,1]
    modified_id <- substr(original_id, 1, 15)
    CELL[,1] <- modified_id
    
    duplicated_rows <- duplicated(CELL[,1])
    unique_data_im <- CELL[!duplicated_rows, ]
    CELL <-unique_data_im
    CELL <- CELL[, -c((ncol(CELL)-2):ncol(CELL))]
    col_names <- colnames(CELL)
    
    col_names_new <- gsub("[_-]", " ", col_names)
    
    colnames(CELL) <- col_names_new
    rownames(CELL) <- CELL[,1]
    data_im <-CELL[,-1]
    
    library(data.table)
    
    dat <-  fread(paste("D:/R/Pancancer ISCA1/Clinical/TCGA_new/",cancer,".txt",sep=''))
    dat <- as.data.frame(dat)
    rownames(dat) <- dat[,1]
    dat <- dat[,-1]
    data_gene <- t(dat)
    
    data_gene <- as.data.frame(data_gene)
    all_name <- names(which(table(c(rownames(data_gene),rownames(data_im) ))==2))
    data_gene <- data_gene[match(all_name,rownames(data_gene)),]
    data_gene <- data_gene[match(all_name,rownames(data_gene)),]
    
    dat_gene <- as.data.frame(data_gene)
    dat_im <- data_im[match(all_name,rownames(data_im)),]
    dat_gene <- dat_gene[,match(allgene1,colnames(dat_gene))]
    
    dat_gene <- as.numeric(unlist(dat_gene))
    library(psych)
    
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
  colnames(dat.r) <- cancers
  colnames(dat.p) <- cancers
  write.csv(dat.r,file= paste("Generation/","_",algorithm,"_","data.r.csv",sep=''),quote=F)
  write.csv(dat.p,file= paste("Generation/","_",algorithm,"_","data.p.csv",sep=''),quote=F)
  data.r_all=dat.r
  data.p_all=dat.p
  data.p[is.na(data.p)] = 2
}




setwd("D:/R/Pancancer ISCA1/8cell")

cancers <- c("ACC","BLCA","BRCA","CESC","CHOL","COAD","DLBC","ESCA","GBM","HNSC","KICH","KIRC","KIRP","LAML","LGG","LIHC","LUAD","LUSC","MESO","OV","PAAD","PCPG",
             "PRAD","READ","SARC","SKCM","STAD","TGCT","THCA","THYM","UCEC","UVM","UCS")

algorithms <- c("epic","mcpcounter","quantiseq","timer","xcell")
for (algorithm in algorithms) {
  dat.r <- c()
  allgene1 <- c("ISCA1")
  for (cancer in cancers) {
    library(data.table)
    CELL <-  fread(paste("8 kinds of algorithms single/TCGA-",cancer,"/im_",algorithm,".csv",sep=''))
    CELL <- as.data.frame(CELL)
    original_id <-CELL[,1]
    modified_id <- substr(original_id, 1, 15)
    CELL[,1] <- modified_id
    duplicated_rows <- duplicated(CELL[,1])
    unique_data_im <- CELL[!duplicated_rows, ]
    CELL <-unique_data_im
    CELL <- CELL[, -c((ncol(CELL)-2):ncol(CELL))]
    col_names <- colnames(CELL)
    col_names_new <- gsub("[_-]", " ", col_names)
    colnames(CELL) <- col_names_new
    rownames(CELL) <- CELL[,1]
    data_im <-CELL[,-1]
    library(data.table)
    
    dat <-  fread(paste("D:/R/Pancancer ISCA1/Clinical/TCGA_new/",cancer,".txt",sep=''))
    dat <- as.data.frame(dat)
    rownames(dat) <- dat[,1]
    dat <- dat[,-1]
    data_gene <- t(dat)
    data_gene <- as.data.frame(data_gene)
    all_name <- names(which(table(c(rownames(data_gene),rownames(data_im) ))==2))
    data_gene <- data_gene[match(all_name,rownames(data_gene)),]
    dat_gene <- as.data.frame(data_gene)
    dat_im <- data_im[match(all_name,rownames(data_im)),]
    dat_gene <- dat_gene[,match(allgene1,colnames(dat_gene))]
    dat_gene <- as.numeric(unlist(dat_gene))
    
    library(psych)
    
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
  
  colnames(dat.r) <- cancers
  colnames(dat.p) <- cancers
  
  write.csv(dat.r,file= paste("Generation/","_",algorithm,"_","data.r.csv",sep=''),quote=F)
  write.csv(dat.p,file= paste("Generation/","_",algorithm,"_","data.p.csv",sep=''),quote=F)
  
  data.r=dat.r
  data.p=dat.p
  
  data.p[is.na(data.p)] = 2
  
  data.r_all <- rbind(data.r_all,data.r)
  data.p_all <- rbind(data.p_all,data.p)
  
}

data.r<- data.r_all
data.p<- data.p_all

countp <- apply(data.p, 1, function(row) sum(row < 0.05))
selected_rows <- names(which(countp > 16))

data.p1 <- data.p[match( selected_rows,rownames(data.p)),]
data.r1 <- data.r[match( selected_rows,rownames(data.r)),]

write.csv(data.p1,file= paste("Generation/","_",algorithm,"_","filter.csv",sep=''),quote=F)

data.r=data.r1
data.p=data.p1

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
less_than_0.05 <- colSums(data.p < 0.05)
less_than_0.05 <- as.data.frame(t(less_than_0.05))
data.p_out <- rbind(data.p, less_than_0.05)
row.names(data.p_out)[nrow(data.p_out)] <- "cells_num_0.05"
less_than_0.05 <- rowSums(data.p_out < 0.05)
cancer_num_0.05 <- as.data.frame(less_than_0.05)
data.p_out <- cbind(data.p_out, cancer_num_0.05)
write.csv(data.p_out,file= paste('num/',"RNA_","data.p_out.csv",sep=''),quote=F)
###Add annotations
class1 <- rep("CIBERSORT", times=2)
class2 <- rep("EPIC", times=3)
class3 <- rep("MCPcounter", times=2)
class4 <- rep("quantiseq", times=2)
class5 <- rep("TIMER", times=1)
class6 <- rep("xCell", times=18)
class<- c(class1,class2,class3,class4,class5,class6)
annotation_row <- data.frame(class = class)
rownames(annotation_row) <- rownames(data.r)
pdf(paste("Filter cor-",Isoform_single,"-cell-rbind-isoform-16.pdf",sep=''),width =12,height =8)
pheatmap(data.r,
         color=myColor,
         breaks=myBreaks,
         clustering_method="average", cluster_rows=F,annotation_row = annotation_row, cluster_cols=F,display_numbers=sig.mat
)
dev.off()

### RNA-ImmuneScore

setwd('D:/R/Pancancer ISCA1/leidap/')

mean_hot <- c()
mean_cold <- c()

library(data.table)
# Read data R
dat_gene_type <-  fread(paste("D:/R/Pancancer ISCA1/leidap/","ISCA1_data.r.csv",sep=''),header = T)
dat_gene_type <- as.data.frame(dat_gene_type)

gene_type_list <- dat_gene_type[,1]
dat_gene_type <-dat_gene_type[,-1]
rownames(dat_gene_type)=gene_type_list
colnames(dat_gene_type)
# Set radar value range
max_min1 <- data.frame(
  ACC=c(0.8,-0.8),BLCA=c(0.8,-0.8),BRCA=c(0.8,-0.8),CESC=c(0.8,-0.8),CHOL=c(0.8,-0.8),COAD=c(0.8,-0.8),DLBC=c(0.8,-0.8),ESCA=c(0.8,-0.8),GBM=c(0.8,-0.8),HNSC=c(0.8,-0.8),KICH=c(0.8,-0.8),
  KIRC=c(0.8,-0.8),KIRP=c(0.8,-0.8),LAML=c(0.8,-0.8),LGG=c(0.8,-0.8),LIHC=c(0.8,-0.8),LUAD=c(0.8,-0.8),LUSC=c(0.8,-0.8),MESO=c(0.8,-0.8),OV=c(0.8,-0.8),PAAD=c(0.8,-0.8),PCPG=c(0.8,-0.8),
  PRAD=c(0.8,-0.8),READ=c(0.8,-0.8),SARC=c(0.8,-0.8),SKCM=c(0.8,-0.8),STAD=c(0.8,-0.8),TGCT=c(0.8,-0.8),THCA=c(0.8,-0.8),THYM=c(0.8,-0.8),UCEC=c(0.8,-0.8),UCS=c(0.8,-0.8),UVM=c(0.8,-0.8)
)

rownames(max_min1) <- c("Max", "Min")

dat_gene_type=t(dat_gene_type)

dat_c1 <- as.matrix(dat_gene_type[,1])
dat_c1=t(dat_c1)
dat_c2 <- as.matrix(dat_gene_type[,2])
dat_c2=t(dat_c2)
dat_c3 <- as.matrix(dat_gene_type[,3])
dat_c3=t(dat_c3)

dat_c1 <-as.numeric(dat_c1)
dat_c2 <-as.numeric(dat_c2)
dat_c3 <-as.numeric(dat_c3)

df_c1c2 <- data.frame(rbind(max_min1,dat_c1,dat_c2,dat_c3))

# Define row names and column names

rownames(df_c1c2) <- c("Max", "Min", "StromalScore", "ESTIMATEScore", "ImmuneScore"  )

# Original vector
combined_vector <- gene_type_list

library(fmsb)
pdf(paste("infiltration-ISCA-Modified.pdf", sep = ''), width = 10, height = 8)

radarchart(
  df_c1c2, axistype = 1,
  # Customize the polygon
  pcol = c("green3", scales::alpha("red3"), scales::alpha("blue3")), # Use green, red and blue
  pfcol = scales::alpha(c("green3", "red3", "blue3"), 0), # Transparency
  plwd = 2, # Set the width of the polygon line, indirectly affecting the size of the corner ring
  plty = 19, # 19 means solid dots
  # Customize the grid
  cglcol = "grey", cglty = 1, cglwd = 0.8,
  # Customize the axis
  axislabcol = "black", 
  # Variable labels with increased font size
  vlcex = 1,  # Increase the font size to 1.5 times
  vlabels = colnames(df_c1c2),  # Column name labels
  caxislabels = c(-0.8, -0.4, 0, 0.4, 0.8)
)


# Add horizontal legend
legend(
  x = -1,
  legend = rownames(df_c1c2[-c(1, 2),]),
  horiz = FALSE,  # Display the legend vertically
  bty = "n",
  pch = 19,
  col = c("green3", scales::alpha("red3"), scales::alpha("blue3")),
  pt.cex = 2,
  text.col = "black",
  cex = 1
)

dev.off()


### CNV-immune cell



setwd("D:/R/Pancancer ISCA1/8cell")
cancers <- c("ACC","BLCA","BRCA","CESC","CHOL","COAD","DLBC","ESCA","GBM","HNSC","KICH","KIRC","KIRP","LAML","LGG","LIHC","LUAD","LUSC","MESO","OV","PAAD","PCPG",
             "PRAD","READ","SARC","SKCM","STAD","TGCT","THCA","THYM","UCEC","UCS","UVM")
algorithms <- c("cibersort")
for (algorithm in algorithms) {
  dat.r <- c()
  allgene1 <- c("ISCA1")
  for (cancer in cancers) {
    library(data.table)
    CELL <-  fread(paste("8 algorithms individual/TCGA-",cancer,"/im_",algorithm,".csv",sep=''))
    CELL <- as.data.frame(CELL)
    original_id <-CELL[,1]
    modified_id <- substr(original_id, 1, 15)
    CELL[,1] <- modified_id
    duplicated_rows <- duplicated(CELL[,1])
    unique_data_im <- CELL[!duplicated_rows, ]
    CELL <-unique_data_im
    CELL <- CELL[, -c((ncol(CELL)-2):ncol(CELL))]
    col_names <- colnames(CELL)
    col_names_new <- gsub("[_-]", " ", col_names)
    colnames(CELL) <- col_names_new
    rownames(CELL) <- CELL[,1]
    data_im <-CELL[,-1]
    library(data.table)
    folder_path <- "D:/R/Pancancer ISCA1/CNV-SNV-methylation/CNV/data/"
    new_file_path <- file.path(folder_path, paste0(cancer, ".cnv.tsv.gz"))
    gz_file <- gzfile(new_file_path, "rt") 
    dataMethy <- read.table(gz_file, header = TRUE, sep = "\t")
    duplicated_rows <- duplicated(dataMethy[,3])
    dataMethy <- dataMethy[!duplicated_rows, ]
    rownames(dataMethy)=dataMethy$symbol
    dataMethy <-dataMethy[,-1:-3]
    colnames(dataMethy) <- gsub("\\.", "-", colnames(dataMethy))
    dataMethy=t(dataMethy)
    
    original_id <-rownames(dataMethy)
    modified_id <- substr(original_id, 1, 15)
    rownames(dataMethy)=modified_id
    data_gene <- dataMethy[,match( allgene1,colnames(dataMethy))]
    data_gene <- as.data.frame(data_gene)
    all_name <- names(which(table(c(rownames(data_gene),rownames(data_im) ))==2))
    dat_gene <- data_gene[match(all_name,rownames(data_gene)),]
    dat_im <- data_im[match(all_name,rownames(data_im)),]
    dat_gene <- as.data.frame(dat_gene)
    colnames(dat_gene) <- "ISCA1"
    library(psych)
    data.corr <- corr.test(dat_im, dat_gene, method="pearson", adjust="fdr")
    data.r <- data.corr$r  # Correlation coefficient
    data.p <- data.corr$p  # p value
    if(length(data.r)==0){
      data.r <- data.r
      data.p <- data.p
    }else {
      data.r <- cbind(data.r,data.r)
      data.p <- cbind(data.p,data.p)
    }
  }
  colnames(data.r) <- c("ACC","BLCA","BRCA","CESC","CHOL","COAD","DLBC","ESCA","GBM","HNSC","KICH","KIRC","KIRP","LAML","LGG","LIHC","LUAD","LUSC","MESO","OV","PAAD","PCPG",
                        "PRAD","READ","SARC","SKCM","STAD","TGCT","THCA","THYM","UCEC","UCS","UVM")
  colnames(data.p) <- c("ACC","BLCA","BRCA","CESC","CHOL","COAD","DLBC","ESCA","GBM","HNSC","KICH","KIRC","KIRP","LAML","LGG","LIHC","LUAD","LUSC","MESO","OV","PAAD","PCPG",
                        "PRAD","READ","SARC","SKCM","STAD","TGCT","THCA","THYM","UCEC","UCS","UVM")
  
  write.csv(data.r,file= paste("Generated/","_",algorithm,"_","data.r.csv",sep=''),quote=F)
  write.csv(data.p,file= paste("Generated/","_",algorithm,"_","data.p.csv",sep=''),quote=F)
  data.r_all=data.r
  data.p_all=data.p
  
  data.p[is.na(data.p)] = 2
  
  
}



setwd("D:/R/Pancancer ISCA1/8cell")

cancers <- c("ACC","BLCA","BRCA","CESC","CHOL","COAD","DLBC","ESCA","GBM","HNSC","KICH","KIRC","KIRP","LAML","LGG","LIHC","LUAD","LUSC","MESO","OV","PAAD","PCPG",
             "PRAD","READ","SARC","SKCM","STAD","TGCT","THCA","THYM","UCEC","UCS","UVM")

algorithms <- c("epic","mcpcounter","quantiseq","timer","xcell")
for (algorithm in algorithms) {
  data.r <- c()
  allgene1 <- c("ISCA1")
  for (cancer in cancers) {
    
    library(data.table)
    CELL <-  fread(paste("8 algorithms individual/TCGA-",cancer,"/im_",algorithm,".csv",sep=''))
    CELL <- as.data.frame(CELL)
    original_id <-CELL[,1]
    modified_id <- substr(original_id, 1, 15)
    CELL[,1] <- modified_id
    duplicated_rows <- duplicated(CELL[,1])
    unique_data_im <- CELL[!duplicated_rows, ]
    CELL <-unique_data_im
    CELL <- CELL[, -c((ncol(CELL)-2):ncol(CELL))]
    col_names <- colnames(CELL)
    col_names_new <- gsub("[_-]", " ", col_names)
    colnames(CELL) <- col_names_new
    rownames(CELL) <- CELL[,1]
    data_im <-CELL[,-1]
    library(data.table)
    folder_path <- "D:/R/Pancancer ISCA1/CNV-SNV-methylation/CNV/data/"
    # Use file.path() function to construct the file path
    new_file_path <- file.path(folder_path, paste0(cancer, ".cnv.tsv.gz"))
    gz_file <- gzfile(new_file_path, "rt")  # Open.gz file
    # Use read.table() function to read the data
    dataMethy <- read.table(gz_file, header = TRUE, sep = "\t")  # Take the tab delimiter as an example, you can choose the delimiter according to the actual situation
    duplicated_rows <- duplicated(dataMethy[,3])
    dataMethy <- dataMethy[!duplicated_rows, ]
    rownames(dataMethy)=dataMethy$symbol
    dataMethy <-dataMethy[,-1:-3]
    colnames(dataMethy) <- gsub("\\.", "-", colnames(dataMethy))
    dataMethy=t(dataMethy)
    original_id <-rownames(dataMethy)
    modified_id <- substr(original_id, 1, 15)
    rownames(dataMethy)=modified_id
    data_gene <- dataMethy[,match( allgene1,colnames(dataMethy))]
    data_gene <- as.data.frame(data_gene)
    all_name <- names(which(table(c(rownames(data_gene),rownames(data_im) ))==2))
    dat_gene <- data_gene[match(all_name,rownames(data_gene)),]
    dat_im <- data_im[match(all_name,rownames(data_im)),]
    dat_gene <- as.data.frame(dat_gene)
    colnames(dat_gene) <- "ISCA1"
    library(psych)
    data.corr <- corr.test(dat_im, dat_gene, method="pearson", adjust="fdr")
    data.r <- data.corr$r  # Correlation coefficient
    data.p <- data.corr$p  # p value
    
    if(length(data.r)==0){
      data.r <- data.r
      data.p <- data.p
    }else {
      data.r <- cbind(data.r,data.r)
      data.p <- cbind(data.p,data.p)
    }
    
  }
  colnames(data.r) <- c("ACC","BLCA","BRCA","CESC","CHOL","COAD","DLBC","ESCA","GBM","HNSC","KICH","KIRC","KIRP","LAML","LGG","LIHC","LUAD","LUSC","MESO","OV","PAAD","PCPG",
                        "PRAD","READ","SARC","SKCM","STAD","TGCT","THCA","THYM","UCEC","UCS","UVM")
  colnames(data.p) <- c("ACC","BLCA","BRCA","CESC","CHOL","COAD","DLBC","ESCA","GBM","HNSC","KICH","KIRC","KIRP","LAML","LGG","LIHC","LUAD","LUSC","MESO","OV","PAAD","PCPG",
                        "PRAD","READ","SARC","SKCM","STAD","TGCT","THCA","THYM","UCEC","UCS","UVM")
  write.csv(data.r,file= paste("Generated/","_",algorithm,"_","data.r.csv",sep=''),quote=F)
  write.csv(data.p,file= paste("Generated/","_",algorithm,"_","data.p.csv",sep=''),quote=F)
  data.r=data.r
  data.p=data.p
  data.p[is.na(data.p)] = 2
  data.r_all <- rbind(data.r_all,data.r)
  data.p_all <- rbind(data.p_all,data.p)
  
}
data.r<- data.r_all
data.p<- data.p_all
selected_rows <- read.csv("Generated/1 screening.csv")
selected_rows <- selected_rows$X
data.r1 <- data.r[match(selected_rows,rownames(data.r)),]
data.p1 <- data.p[match(selected_rows,rownames(data.p)),]
data.r=data.r1
data.p=data.p1
data.p[is.na(data.p)] = 2
less_than_0.05 <- colSums(data.p < 0.05)
less_than_0.05 <- as.data.frame(t(less_than_0.05))
data.p_out <- rbind(data.p, less_than_0.05)
row.names(data.p_out)[nrow(data.p_out)] <- "cells_num_0.05"
less_than_0.05 <- rowSums(data.p_out < 0.05)
cancer_num_0.05 <- as.data.frame(less_than_0.05)
data.p_out <- cbind(data.p_out, cancer_num_0.05)
write.csv(data.p_out,file= paste('num/',"CNV_","data.p_out.csv",sep=''),quote=F)

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

class1 <- rep("CIBERSORT", times=2)
class2 <- rep("EPIC", times=3)
class3 <- rep("MCPcounter", times=2)
class4 <- rep("quantiseq", times=2)
class5 <- rep("TIMER", times=1)
class6 <- rep("xCell", times=22)
class<- c(class1,class2,class3,class4,class5,class6)

annotation_row <- data.frame(class = class)
rownames(annotation_row) <- rownames(data.r)

pdf(paste("cor-cell-rbind-CNV-screening 16.pdf",sep=''),width =12,height =8)
pheatmap(data.r,
         color=myColor,
         breaks=myBreaks,
         clustering_method="average", cluster_rows=F,annotation_row = annotation_row, cluster_cols=F,display_numbers=sig.mat
)
dev.off()



### CNV-ImmuneScore



setwd('D:/R/Pancancer ISCA1/leidap/')

mean_hot <- c()
mean_cold <- c()

library(data.table)
# Read data R
dat_gene_type <-  fread(paste("D:/R/Pancancer ISCA1/leidap/","CNV_data.r.csv",sep=''),header = T)
dat_gene_type <- as.data.frame(dat_gene_type)

gene_type_list <- dat_gene_type[,1]
dat_gene_type <-dat_gene_type[,-1]
rownames(dat_gene_type)=gene_type_list

colnames(dat_gene_type)

# Set radar value range
max_min1 <- data.frame(
  ACC=c(0.8,-0.8),BLCA=c(0.8,-0.8),BRCA=c(0.8,-0.8),CESC=c(0.8,-0.8),CHOL=c(0.8,-0.8),COAD=c(0.8,-0.8),DLBC=c(0.8,-0.8),ESCA=c(0.8,-0.8),GBM=c(0.8,-0.8),HNSC=c(0.8,-0.8),KICH=c(0.8,-0.8),
  KIRC=c(0.8,-0.8),KIRP=c(0.8,-0.8),LAML=c(0.8,-0.8),LGG=c(0.8,-0.8),LIHC=c(0.8,-0.8),LUAD=c(0.8,-0.8),LUSC=c(0.8,-0.8),MESO=c(0.8,-0.8),OV=c(0.8,-0.8),PAAD=c(0.8,-0.8),PCPG=c(0.8,-0.8),
  PRAD=c(0.8,-0.8),READ=c(0.8,-0.8),SARC=c(0.8,-0.8),SKCM=c(0.8,-0.8),STAD=c(0.8,-0.8),TGCT=c(0.8,-0.8),THCA=c(0.8,-0.8),THYM=c(0.8,-0.8),UCEC=c(0.8,-0.8),UCS=c(0.8,-0.8),UVM=c(0.8,-0.8)
)

rownames(max_min1) <- c("Max", "Min")

dat_gene_type=t(dat_gene_type)

dat_c1 <- as.matrix(dat_gene_type[,1])
dat_c1=t(dat_c1)
dat_c2 <- as.matrix(dat_gene_type[,2])
dat_c2=t(dat_c2)
dat_c3 <- as.matrix(dat_gene_type[,3])
dat_c3=t(dat_c3)

dat_c1 <-as.numeric(dat_c1)
dat_c2 <-as.numeric(dat_c2)
dat_c3 <-as.numeric(dat_c3)

df_c1c2 <- data.frame(rbind(max_min1,dat_c1,dat_c2,dat_c3))
rownames(df_c1c2) <- c("Max", "Min", "StromalScore", "ESTIMATEScore", "ImmuneScore")
# Original vector
combined_vector <- gene_type_list

library(fmsb)
pdf(paste("CNV-ISCA1-Modified.pdf", sep = ''), width = 10, height = 8)

radarchart(
  df_c1c2, axistype = 1,
  # Customize the polygon
  pcol = c("green3", scales::alpha("red3"), scales::alpha("blue3")), # Use green, red and blue
  pfcol = scales::alpha(c("green3", "red3", "blue3"), 0), # Transparency
  plwd = 2, # Set the width of the polygon line, indirectly affecting the size of the corner ring
  plty = 19, # 19 means solid dots
  # Customize the grid
  cglcol = "grey", cglty = 1, cglwd = 0.8,
  # Customize the axis
  axislabcol = "black", 
  # Variable labels with increased font size
  vlcex = 1,  # Increase the font size to 1.5 times
  vlabels = colnames(df_c1c2),  # Column name labels
  caxislabels = c(-0.8, -0.4, 0, 0.4, 0.8)
)

# Add horizontal legend
legend(
  x = -1,
  legend = rownames(df_c1c2[-c(1, 2),]),
  horiz = FALSE,  # Display the legend vertically
  bty = "n",
  pch = 19,
  col = c("green3", scales::alpha("red3"), scales::alpha("blue3")),
  pt.cex = 2,
  text.col = "black",
  cex = 1
)

dev.off()


### Methylation-immune cell


setwd("D:/R/Pancancer ISCA1/8cell")

cancers <- c("ACC","BLCA","BRCA","CESC","CHOL","COAD","DLBC","ESCA","GBM","HNSC","KICH","KIRC","KIRP","LAML","LGG","LIHC","LUAD","LUSC","MESO","OV","PAAD","PCPG",
             "PRAD","READ","SARC","SKCM","STAD","TGCT","THCA","THYM","UCEC","UCS","UVM")

algorithms <- c("cibersort")

for (algorithm in algorithms) {
  dat.r <- c()
  allgene1 <- c("ISCA1")
  for (cancer in cancers) {
    library(data.table)
    CELL <-  fread(paste("8 kinds of algorithms single/TCGA-",cancer,"/im_",algorithm,".csv",sep=''))
    CELL <- as.data.frame(CELL)
    original_id <-CELL[,1]
    modified_id <- substr(original_id, 1, 15)
    CELL[,1] <- modified_id
    duplicated_rows <- duplicated(CELL[,1])
    unique_data_im <- CELL[!duplicated_rows, ]
    CELL <-unique_data_im
    CELL <- CELL[, -c((ncol(CELL)-2):ncol(CELL))]
    col_names <- colnames(CELL)
    col_names_new <- gsub("[_-]", " ", col_names)
    colnames(CELL) <- col_names_new
    rownames(CELL) <- CELL[,1]
    data_im <-CELL[,-1]
    library(data.table)
    methy <- fread(paste("D:/R/Pancancer ISCA1/Linear/new-methy/ISCA1-methy/",cancer,"-Methy450.csv",sep=''))
    methy <- data.frame(methy)
    colnames(methy) <- gsub("\\.", "-", colnames(methy))
    data_gene <- t(methy)
    all_name <- names(which(table(c(rownames(data_gene),rownames(data_im) ))==2))
    data_gene <- data_gene[match(all_name,rownames(data_gene)),]
    dat_gene <- as.data.frame(data_gene)
    
    dat_im <- data_im[match(all_name,rownames(data_im)),]
    colnames(dat_gene) <- "ISCA1"
    dat_gene <- as.numeric(unlist(dat_gene))
    library(psych)
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
  colnames(dat.r) <- cancers
  colnames(dat.p) <- cancers
  
  write.csv(dat.r,file= paste("Generation/","_",algorithm,"_","data.r.csv",sep=''),quote=F)
  write.csv(dat.p,file= paste("Generation/","_",algorithm,"_","data.p.csv",sep=''),quote=F)
  
  data.r_all=dat.r
  data.p_all=dat.p
  
  data.p[is.na(data.p)] = 2
  
}



setwd("D:/R/Pancancer ISCA1/8cell")

cancers <- c("ACC","BLCA","BRCA","CESC","CHOL","COAD","DLBC","ESCA","GBM","HNSC","KICH","KIRC","KIRP","LAML","LGG","LIHC","LUAD","LUSC","MESO","OV","PAAD","PCPG",
             "PRAD","READ","SARC","SKCM","STAD","TGCT","THCA","THYM","UCEC","UCS","UVM")
algorithms <- c("epic","mcpcounter","quantiseq","timer","xcell")
for (algorithm in algorithms) {
  dat.r <- c()
  allgene1 <- c("ISCA1")
  for (cancer in cancers) {
    library(data.table)
    CELL <-  fread(paste("8 kinds of algorithms single/TCGA-",cancer,"/im_",algorithm,".csv",sep=''))
    CELL <- as.data.frame(CELL)
    original_id <-CELL[,1]
    modified_id <- substr(original_id, 1, 15)
    CELL[,1] <- modified_id
    duplicated_rows <- duplicated(CELL[,1])
    unique_data_im <- CELL[!duplicated_rows, ]
    CELL <-unique_data_im
    CELL <- CELL[, -c((ncol(CELL)-2):ncol(CELL))]
    col_names <- colnames(CELL)
    col_names_new <- gsub("[_-]", " ", col_names)
    colnames(CELL) <- col_names_new
    rownames(CELL) <- CELL[,1]
    data_im <-CELL[,-1]
    library(data.table)
    methy <- fread(paste("D:/R/Pancancer ISCA1/Linear/new-methy/ISCA1-methy/",cancer,"-Methy450.csv",sep=''))
    methy <- data.frame(methy)
    colnames(methy) <- gsub("\\.", "-", colnames(methy))
    data_gene <- t(methy)
    all_name <- names(which(table(c(rownames(data_gene),rownames(data_im) ))==2))
    
    data_gene <- data_gene[match(all_name,rownames(data_gene)),]
    dat_gene <- as.data.frame(data_gene)
    
    dat_im <- data_im[match(all_name,rownames(data_im)),]
    colnames(dat_gene) <- "ISCA1"
    dat_gene <- as.numeric(unlist(dat_gene))
    
    library(psych)
    
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
  
  write.csv(dat.r,file= paste("Generation/","_",algorithm,"_","data.r.csv",sep=''),quote=F)
  write.csv(dat.p,file= paste("Generation/","_",algorithm,"_","data.p.csv",sep=''),quote=F)
  
  data.r=dat.r
  data.p=dat.p
  
  data.p[is.na(data.p)] = 2
  
  data.r_all <- rbind(data.r_all,data.r)
  data.p_all <- rbind(data.p_all,data.p)
}

data.r<- data.r_all
data.p<- data.p_all

selected_rows <- read.csv("Generation/1Filter.csv")
selected_rows <- selected_rows$X

data.r1 <- data.r[match(selected_rows,rownames(data.r)),]
data.p1 <- data.p[match(selected_rows,rownames(data.p)),]

data.r=data.r1
data.p=data.p1
data.p[is.na(data.p)] = 2
less_than_0.05 <- colSums(data.p < 0.05)
less_than_0.05 <- as.data.frame(t(less_than_0.05))
data.p_out <- rbind(data.p, less_than_0.05)
row.names(data.p_out)[nrow(data.p_out)] <- "cells_num_0.05"
less_than_0.05 <- rowSums(data.p_out < 0.05)
cancer_num_0.05 <- as.data.frame(less_than_0.05)
data.p_out <- cbind(data.p_out, cancer_num_0.05)
write.csv(data.p_out,file= paste('num/',"methy_","data.p_out.csv",sep=''),quote=F)
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
class1 <- rep("CIBERSORT", times=2)
class2 <- rep("EPIC", times=3)
class3 <- rep("MCPcounter", times=2)
class4 <- rep("quantiseq", times=2)
class5 <- rep("TIMER", times=1)
class6 <- rep("xCell", times=22)
class<- c(class1,class2,class3,class4,class5,class6)
annotation_row <- data.frame(class = class)
rownames(annotation_row) <- rownames(data.r)
pdf(paste("cor-cell-rbind-methy-newFilter16.pdf",sep=''),width =12,height =8)
pheatmap(data.r,
         color=myColor,
         breaks=myBreaks,
         clustering_method="average", cluster_rows=F,annotation_row = annotation_row, cluster_cols=F,display_numbers=sig.mat
)
dev.off()


### Methylation-ImmuneScore

setwd('D:/R/Pancancer ISCA1/leidap/')

mean_hot <- c()
mean_cold <- c()

library(data.table)
# Read data R
data_gene_type <-  fread(paste("D:/R/Pancancer ISCA1/leidap/","methy-data.r.csv",sep=''),header = T)
data_gene_type <- as.data.frame(data_gene_type)


gene_type_list <- data_gene_type[,1]
data_gene_type <-data_gene_type[,-1]
rownames(data_gene_type)=gene_type_list

# Set radar value range
max_min1 <- data.frame(
  ACC=c(0.8,-0.8),BLCA=c(0.8,-0.8),BRCA=c(0.8,-0.8),CESC=c(0.8,-0.8),CHOL=c(0.8,-0.8),COAD=c(0.8,-0.8),DLBC=c(0.8,-0.8),ESCA=c(0.8,-0.8),GBM=c(0.8,-0.8),HNSC=c(0.8,-0.8),KICH=c(0.8,-0.8),
  KIRC=c(0.8,-0.8),KIRP=c(0.8,-0.8),LAML=c(0.8,-0.8),LGG=c(0.8,-0.8),LIHC=c(0.8,-0.8),LUAD=c(0.8,-0.8),LUSC=c(0.8,-0.8),MESO=c(0.8,-0.8),OV=c(0.8,-0.8),PAAD=c(0.8,-0.8),PCPG=c(0.8,-0.8),
  PRAD=c(0.8,-0.8),READ=c(0.8,-0.8),SARC=c(0.8,-0.8),SKCM=c(0.8,-0.8),STAD=c(0.8,-0.8),TGCT=c(0.8,-0.8),THCA=c(0.8,-0.8),THYM=c(0.8,-0.8),UCEC=c(0.8,-0.8),UCS=c(0.8,-0.8),UVM=c(0.8,-0.8)
)

rownames(max_min1) <- c("Max", "Min")


data_gene_type=t(data_gene_type)

data_c1 <- as.matrix(data_gene_type[,1])
data_c1=t(data_c1)
data_c2 <- as.matrix(data_gene_type[,2])
data_c2=t(data_c2)
data_c3 <- as.matrix(data_gene_type[,3])
data_c3=t(data_c3)


data_c1 <-as.numeric(data_c1)
data_c2 <-as.numeric(data_c2)
data_c3 <-as.numeric(data_c3)

df_c1c2 <- data.frame(rbind(max_min1,data_c1,data_c2,data_c3))

rownames(df_c1c2) <- c("Max", "Min", "StromalScore", "ESTIMATEScore" ,"ImmuneScore"  )

# Original vector
combined_vector <- gene_type_list

library(fmsb)
pdf(paste("Methy-ISCA1-new1.pdf", sep = ''), width = 10, height = 8)

radarchart(
  df_c1c2, axistype = 1,
  # Customize the polygon
  pcol = c("green3", scales::alpha("red3"), scales::alpha("blue3")), # Use green, red and blue
  pfcol = scales::alpha(c("green3", "red3", "blue3"), 0), # Transparency
  plwd = 2, # Set the width of the polygon line, indirectly affecting the size of the corner ring
  plty = 19, # 19 means solid dots
  # Customize the grid
  cglcol = "grey", cglty = 1, cglwd = 0.8,
  # Customize the axis
  axislabcol = "black", 
  # Variable labels with increased font size
  vlcex = 1,  # Increase font size to 1.5 times
  vlabels = colnames(df_c1c2),  # Column name labels
  caxislabels = c(-0.8, -0.4, 0, 0.4, 0.8)
)


# Add horizontal legend
legend(
  x = -1,
  legend = rownames(df_c1c2[-c(1, 2),]),
  horiz = FALSE,  # Display the legend vertically
  bty = "n",
  pch = 19,
  col = c("green3", scales::alpha("red3"), scales::alpha("blue3")),
  pt.cex = 2,
  text.col = "black",
  cex = 1
)

dev.off()


### Isoform-immune cell

rm (list = ls ())

setwd("D:/R/Pancancer ISCA1/8cell")

cancers <- c("ACC","BLCA","BRCA","CESC","CHOL","COAD","DLBC","ESCA","GBM","HNSC","KICH","KIRC","KIRP","LAML","LGG","LIHC","LUAD","LUSC","OV","PAAD","PCPG",
             "PRAD","READ","SARC","SKCM","STAD","TGCT","THCA","THYM","UCEC","UCS")

isoform_singles <- c("ENST00000375991","ENST00000311534","ENST00000326094")

for (isoform_single in isoform_singles) {
  algorithms <- c("cibersort")
  for (algorithm in algorithms) {
    
    dat.r <- c()
    allgene1 <- c("ISCA1")
    for (cancer in cancers) {
      library(data.table)
      CELL <-  fread(paste("8 algorithms individual/TCGA-",cancer,"/im_",algorithm,".csv",sep=''))
      CELL <- as.data.frame(CELL)
      original_id <-CELL[,1]
      modified_id <- substr(original_id, 1, 15)
      CELL[,1] <- modified_id
      duplicated_rows <- duplicated(CELL[,1])
      unique_data_im <- CELL[!duplicated_rows, ]
      CELL <-unique_data_im
      CELL <- CELL[, -c((ncol(CELL)-2):ncol(CELL))]
      col_names <- colnames(CELL)
      col_names_new <- gsub("[_-]", " ", col_names)
      colnames(CELL) <- col_names_new
      
      rownames(CELL) <- CELL[,1]
      data_im <-CELL[,-1]
      library(data.table)
      Isoform <- read.csv(paste0("D:/R/Pancancer ISCA1/Isoform/",isoform_single,"_transcript_pancan_dist.CSV"))
      Isoform <- subset(Isoform, tissue == cancer)
      rownames(Isoform) <- Isoform$sample
      data_gene <- Isoform
      all_name <- names(which(table(c(rownames(data_gene),rownames(data_im) ))==2))
      data_gene <- data_gene[match(all_name,rownames(data_gene)),]
      data_gene <- data_gene$tpm
      dat_gene <- as.data.frame(data_gene)
      dat_im <- data_im[match(all_name,rownames(data_im)),]
      colnames(dat_gene) <- "ISCA1"
      
      dat_gene <- as.numeric(unlist(dat_gene))
      library(psych)
      data.corr <- corr.test(dat_im, dat_gene, method="pearson", adjust="fdr")
      data.r <- data.corr$r  # Correlation coefficient
      data.p <- data.corr$p  # p value
      
      if(length(data.r)==0){
        data.r <- data.r
        data.p <- data.p
      }else {
        data.r <- cbind(data.r,data.r)
        data.p <- cbind(data.p,data.p)
      }
      
    }
    colnames(data.r) <- cancers
    colnames(data.p) <- cancers
    write.csv(data.r,file= paste("Generated/","_",algorithm,"_","data.r.csv",sep=''),quote=F)
    write.csv(data.p,file= paste("Generated/","_",algorithm,"_","data.p.csv",sep=''),quote=F)
    data.r_all=data.r
    data.p_all=data.p
    data.p[is.na(data.p)] = 2
  }
  
  
  setwd("D:/R/Pancancer ISCA1/8cell")
  
  cancers <- c("ACC","BLCA","BRCA","CESC","CHOL","COAD","DLBC","ESCA","GBM","HNSC","KICH","KIRC","KIRP","LAML","LGG","LIHC","LUAD","LUSC","OV","PAAD","PCPG",
               "PRAD","READ","SARC","SKCM","STAD","TGCT","THCA","THYM","UCEC","UCS")
  algorithms <- c("epic","mcpcounter","quantiseq","timer","xcell")
  for (algorithm in algorithms) {
    data.r <- c()
    allgene1 <- c("ISCA1")
    for (cancer in cancers) {
      
      library(data.table)
      CELL <-  fread(paste("8 algorithms individual/TCGA-",cancer,"/im_",algorithm,".csv",sep=''))
      CELL <- as.data.frame(CELL)
      original_id <-CELL[,1]
      modified_id <- substr(original_id, 1, 15)
      CELL[,1] <- modified_id
      duplicated_rows <- duplicated(CELL[,1])
      unique_data_im <- CELL[!duplicated_rows, ]
      CELL <-unique_data_im
      CELL <- CELL[, -c((ncol(CELL)-2):ncol(CELL))]
      col_names <- colnames(CELL)
      col_names_new <- gsub("[_-]", " ", col_names)
      colnames(CELL) <- col_names_new
      rownames(CELL) <- CELL[,1]
      data_im <-CELL[,-1]
      library(data.table)
      Isoform <- read.csv(paste0("D:/R/Pancancer ISCA1/Isoform/",isoform_single,"_transcript_pancan_dist.CSV"))
      # Use subset function for screening
      Isoform <- subset(Isoform, tissue == cancer)
      rownames(Isoform) <- Isoform$sample
      data_gene <- Isoform
      all_name <- names(which(table(c(rownames(data_gene),rownames(data_im) ))==2))
      data_gene <- data_gene[match(all_name,rownames(data_gene)),]
      data_gene <- data_gene$tpm
      dat_gene <- as.data.frame(data_gene)
      
      dat_im <- data_im[match(all_name,rownames(data_im)),]
      colnames(dat_gene) <- "ISCA1"
      dat_gene <- as.numeric(unlist(dat_gene))
      library(psych)
      data.corr <- corr.test(dat_im, dat_gene, method="pearson", adjust="fdr")
      data.r <- data.corr$r  # Correlation coefficient
      data.p <- data.corr$p  # p value
      
      if(length(data.r)==0){
        data.r <- data.r
        data.p <- data.p
      }else {
        data.r <- cbind(data.r,data.r)
        data.p <- cbind(data.p,data.p)
      }
    }
    
    colnames(data.r) <- cancers
    colnames(data.p) <- cancers
    
    write.csv(data.r,file= paste("Generated/","_",algorithm,"_","data.r.csv",sep=''),quote=F)
    write.csv(data.p,file= paste("Generated/","_",algorithm,"_","data.p.csv",sep=''),quote=F)
    
    data.r=data.r
    data.p=data.p
    
    data.p[is.na(data.p)] = 2
    data.r_all <- rbind(data.r_all,data.r)
    data.p_all <- rbind(data.p_all,data.p)
  }
  data.r<- data.r_all
  data.p<- data.p_all
  selected_rows <- read.csv("Generated/1 screening.csv")
  selected_rows <- selected_rows$X
  data.r1 <- data.r[match(selected_rows,rownames(data.r)),]
  data.p1 <- data.p[match(selected_rows,rownames(data.p)),]
  data.r=data.r1
  data.p=data.p1
  
  data.p[is.na(data.p)] = 2
  less_than_0.05 <- colSums(data.p < 0.05)
  less_than_0.05 <- as.data.frame(t(less_than_0.05))
  data.p_out <- rbind(data.p, less_than_0.05)
  row.names(data.p_out)[nrow(data.p_out)] <- "cells_num_0.05"
  less_than_0.05 <- rowSums(data.p_out < 0.05)
  # Ensure the result is a numeric vector consistent with the number of columns of the data frame
  cancer_num_0.05 <- as.data.frame(less_than_0.05)
  data.p_out <- cbind(data.p_out, cancer_num_0.05)
  write.csv(data.p_out,file= paste('num/',isoform_single,"_","data.p_out.csv",sep=''),quote=F)
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
  
  class1 <- rep("CIBERSORT", times=2)
  class2 <- rep("EPIC", times=3)
  class3 <- rep("MCPcounter", times=2)
  class4 <- rep("quantiseq", times=2)
  class5 <- rep("TIMER", times=1)
  class6 <- rep("xCell", times=22)
  class<- c(class1,class2,class3,class4,class5,class6)
  annotation_row <- data.frame(class = class)
  rownames(annotation_row) <- rownames(data.r)
  pdf(paste("Screened cor-",isoform_single,"-cell-rbind-isoform.pdf",sep=''),width =12,height =8)
  pheatmap(data.r,
           color=myColor,
           breaks=myBreaks,
           clustering_method="average", cluster_rows=F,annotation_row = annotation_row, cluster_cols=F,display_numbers=sig.mat
  )
  
  dev.off()
}



### Isoform-ImmuneScore


setwd('D:/R/Pancancer ISCA1/leidap/')

mean_hot <- c()
mean_cold <- c()

library(data.table)
# Read data R
data_gene_type <-  fread(paste("D:/R/Pancancer ISCA1/leidap/","ENST00000311534-data.r.csv",sep=''),header = T)
data_gene_type <- as.data.frame(data_gene_type)

gene_type_list <- data_gene_type[,1]
data_gene_type <-data_gene_type[,-1]
rownames(data_gene_type)=gene_type_list

# Set radar value range
max_min1 <- data.frame(
  ACC=c(0.8,-0.8),BLCA=c(0.8,-0.8),BRCA=c(0.8,-0.8),CESC=c(0.8,-0.8),CHOL=c(0.8,-0.8),COAD=c(0.8,-0.8),DLBC=c(0.8,-0.8),ESCA=c(0.8,-0.8),GBM=c(0.8,-0.8),HNSC=c(0.8,-0.8),KICH=c(0.8,-0.8),
  KIRC=c(0.8,-0.8),KIRP=c(0.8,-0.8),LAML=c(0.8,-0.8),LGG=c(0.8,-0.8),LIHC=c(0.8,-0.8),LUAD=c(0.8,-0.8),LUSC=c(0.8,-0.8),OV=c(0.8,-0.8),PAAD=c(0.8,-0.8),PCPG=c(0.8,-0.8),
  PRAD=c(0.8,-0.8),READ=c(0.8,-0.8),SARC=c(0.8,-0.8),SKCM=c(0.8,-0.8),STAD=c(0.8,-0.8),TGCT=c(0.8,-0.8),THCA=c(0.8,-0.8),THYM=c(0.8,-0.8),UCEC=c(0.8,-0.8),UCS=c(0.8,-0.8)
)

rownames(max_min1) <- c("Max", "Min")


data_gene_type=t(data_gene_type)

data_c1 <- as.matrix(data_gene_type[,1])
data_c1=t(data_c1)
data_c2 <- as.matrix(data_gene_type[,2])
data_c2=t(data_c2)
data_c3 <- as.matrix(data_gene_type[,3])
data_c3=t(data_c3)


data_c1 <-as.numeric(data_c1)
data_c2 <-as.numeric(data_c2)
data_c3 <-as.numeric(data_c3)

df_c1c2 <- data.frame(rbind(max_min1,data_c1,data_c2,data_c3))

rownames(df_c1c2) <- c("Max", "Min", "StromalScore", "ESTIMATEScore" ,"ImmuneScore"  )

# Original vector
combined_vector <- gene_type_list

library(fmsb)
pdf(paste("ENST00000311534-data.r-ISCA1.pdf", sep = ''), width = 10, height = 8)

radarchart(
  df_c1c2, axistype = 1,
  # Customize the polygon
  pcol = c("green3", scales::alpha("red3"), scales::alpha("blue3")), # Use green, red and blue
  pfcol = scales::alpha(c("green3", "red3", "blue3"), 0), # Transparency
  plwd = 2, # Set the width of the polygon line, indirectly affecting the size of the corner ring
  plty = 19, # 19 means solid dots
  # Customize the grid
  cglcol = "grey", cglty = 1, cglwd = 0.8,
  # Customize the axis
  axislabcol = "black", 
  # Variable labels with increased font size
  vlcex = 1,  # Increase font size to 1.5 times
  vlabels = colnames(df_c1c2),  # Column name labels
  caxislabels = c(-0.8, -0.4, 0, 0.4, 0.8)
)


# Add horizontal legend
legend(
  x = -1,
  legend = rownames(df_c1c2[-c(1, 2),]),
  horiz = FALSE,  # Display the legend vertically
  bty = "n",
  pch = 19,
  col = c("green3", scales::alpha("red3"), scales::alpha("blue3")),
  pt.cex = 2,
  text.col = "black",
  cex = 1
)

dev.off()

