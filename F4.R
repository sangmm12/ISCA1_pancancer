
### RNA-related analysis of ferroptosis gene associations

setwd("D:/R/Pancancer ISCA1/Linear/")

cancers <- c("ACC","BLCA","BRCA","CESC","CHOL","COAD","DLBC","ESCA","GBM","HNSC","KICH","KIRC","KIRP","LAML","LGG","LIHC","LUAD","LUSC","MESO","OV","PAAD","PCPG",
             "PRAD","READ","SARC","SKCM","STAD","TGCT","THCA","THYM","UCEC","UCS","UVM")

dat.r <- c()

allgene1 <- c("ISCA1")

library(data.table)

gene <-  fread(paste("D:/R/Pancancer ISCA1/Linear/","ferroptosis_driver.csv",sep=''))
gene <- as.data.frame(gene)

gene1 <-  fread(paste("D:/R/Pancancer ISCA1/Linear/","ferroptosis_marker.csv",sep=''))
gene1 <- as.data.frame(gene1)

gene2 <-  fread(paste("D:/R/Pancancer ISCA1/Linear/","ferroptosis_suppressor.csv",sep=''))
gene2 <- as.data.frame(gene2)

gene3 <-  fread(paste("D:/R/Pancancer ISCA1/Linear/","ferroptosis_unclassified.csv",sep=''))
gene3 <- as.data.frame(gene3)


generbind=rbind(gene,gene1,gene2,gene3)

allgenerbind=generbind$symbol
allgenerbind=unique(allgenerbind)


allgene2 <- allgenerbind
allgene2 <- unique(allgene2)


for (cancer in cancers) {
  # cancer <- c("BLCA")
  
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
  
  
  
  dat_gene <- data[,match( allgene1,colnames(data))]
  dat_im <- data[,match(allgene2,colnames(data))]
  
  dat_im <- dat_im[, colSums(is.na(dat_im)) == 0]
  
  library(psych)
  
  
  data.corr <- corr.test(dat_im, dat_gene, method="pearson", adjust="fdr")
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


# Filter data with P < 0.05
countp <- apply(dat.p, 1, function(row) sum(row < 0.05))
selected_rows <- names(which(countp > 28))

data.p1 <- data.p[match( selected_rows,rownames(data.p)),]
data.r1 <- data.r[match( selected_rows,rownames(data.r)),]

data.r=data.r1
data.p=data.p1


# Find rows with NA in the first column
rows_with_na <- is.na(data.r[, 1])

# Delete rows with NA
data.r <- data.r[!rows_with_na, ]
data.p <- data.p[!rows_with_na, ]

data.rWRITE <- data.r
data.pWRITE <- data.p


write.csv(data.rWRITE,file= paste(cancer,"rbin_data.r.csv",sep=''),quote=F)
write.csv(data.pWRITE,file= paste(cancer,"rbin_data.p.csv",sep=''),quote=F)


### Add comments
class1 <- rep("ferroptosis_driver", times=246)
class2 <- rep("ferroptosis_marker", times=6)
class3 <- rep("ferroptosis_suppressor", times=187)
class4 <- rep("ferroptosis_unclassified", times=73)
class<- c(class1,class2,class3,class4)

annotation_row <- data.frame(class = class)
rownames(annotation_row) <- rownames(data.r)




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


pdf(paste("cor-ferroptosis_rbind.pdf",sep=''),width =12,height =120)
pheatmap(data.r,
         color=myColor,
         annotation_row = annotation_row,
         breaks=myBreaks,
         clustering_method="average", cluster_rows=F, cluster_cols=F,display_numbers=sig.mat
)

dev.off()




### CNV-related analysis of ferroptosis gene associations

setwd("D:/R/Pancancer ISCA1/Linear/Fe-die")

cancers <- c("ACC","BLCA","BRCA","CESC","CHOL","COAD","DLBC","ESCA","GBM","HNSC","KICH","KIRC","KIRP","LAML","LGG","LIHC","LUAD","LUSC","MESO","OV","PAAD","PCPG",
             "PRAD","READ","SARC","SKCM","STAD","TGCT","THCA","THYM","UCEC","UCS","UVM")

dat.r <- c()

allgene1 <- c("ISCA1")

library(data.table)

gene <-  fread(paste("D:/R/Pancancer ISCA1/Linear/","ferroptosis_driver.csv",sep=''))
gene <- as.data.frame(gene)

gene1 <-  fread(paste("D:/R/Pancancer ISCA1/Linear/","ferroptosis_marker.csv",sep=''))
gene1 <- as.data.frame(gene1)

gene2 <-  fread(paste("D:/R/Pancancer ISCA1/Linear/","ferroptosis_suppressor.csv",sep=''))
gene2 <- as.data.frame(gene2)

gene3 <-  fread(paste("D:/R/Pancancer ISCA1/Linear/","ferroptosis_unclassified).csv",sep=''))
gene3 <- as.data.frame(gene3)



generbind=rbind(gene,gene1,gene2,gene3)

allgenerbind=generbind$symbol
allgenerbind=unique(allgenerbind)


allgene2 <- allgenerbind
allgene2 <- unique(allgene2)

for (cancer in cancers) {
  # cancer <- c("BLCA")
  
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
  
  # CNV data
  folder_path <- "D:/R/Pancancer ISCA1/CNV-SNV-methy/CNV/data/"
  # cancer <- "ACC"
  # Use file.path() function to construct the file path
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

# Find rows with NA in the first column
rows_with_na <- is.na(data.r[, 1])

# Delete rows with NA
data.r <- data.r[!rows_with_na, ]
data.p <- data.p[!rows_with_na, ]

write.csv(data.r,file= paste(cancer,"rbin_data.r.csv",sep=''),quote=F)
write.csv(data.p,file= paste(cancer,"rbin_data.p.csv",sep=''),quote=F)

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

### Add comments
class1 <- rep("ferroptosis_driver", times=246)
class2 <- rep("ferroptosis_marker", times=6)
class3 <- rep("ferroptosis_suppressor", times=187)
class4 <- rep("ferroptosis_unclassified", times=73)
class<- c(class1,class2,class3,class4)

annotation_row <- data.frame(class = class)
rownames(annotation_row) <- rownames(data.r)


pdf(paste("cor-ferroptosis_CNV-rbind-ALL.pdf",sep=''),width =12,height =120)
pheatmap(data.r,
         color=myColor,
         breaks=myBreaks,
         clustering_method="average", cluster_rows=F,annotation_row = annotation_row, cluster_cols=F,display_numbers=sig.mat
)

dev.off()




### Methylation-related analysis of ferroptosis gene associations

setwd("D:/R/Pancancer ISCA1/Linear/Fe-die")

cancers <- c("ACC","BLCA","BRCA","CESC","CHOL","COAD","DLBC","ESCA","GBM","HNSC","KICH","KIRC","KIRP","LAML","LGG","LIHC","LUAD","LUSC","MESO","OV","PAAD","PCPG",
             "PRAD","READ","SARC","SKCM","STAD","TGCT","THCA","THYM","UCEC","UCS","UVM")

dat.r <- c()

allgene1 <- c("ISCA1")

library(data.table)

gene <-  fread(paste("D:/R/Pancancer ISCA1/Linear/","ferroptosis_driver.csv",sep=''))
gene <- as.data.frame(gene)

gene1 <-  fread(paste("D:/R/Pancancer ISCA1/Linear/","ferroptosis_marker.csv",sep=''))
gene1 <- as.data.frame(gene1)

gene2 <-  fread(paste("D:/R/Pancancer ISCA1/Linear/","ferroptosis_suppressor.csv",sep=''))
gene2 <- as.data.frame(gene2)

gene3 <-  fread(paste("D:/R/Pancancer ISCA1/Linear/","ferroptosis_unclassified).csv",sep=''))
gene3 <- as.data.frame(gene3)



generbind=rbind(gene,gene1,gene2,gene3)

allgenerbind=generbind$symbol
allgenerbind=unique(allgenerbind)


allgene2 <- allgenerbind
allgene2 <- unique(allgene2)

for (cancer in cancers) {
  
  #  cancer <- "ACC"
  #  cg13456106_ISCA1
  
  library(data.table)
  
  dat <-  fread(paste("D:/R/Pancancer ISCA1/Clinical/TCGA_new/",cancer,".txt",sep=''))
  dat <- as.data.frame(dat)
  rownames(dat) <- dat[,1]
  dat <- dat[,-1]
  data <- dat
  group=sapply(strsplit(colnames(data),"\\-"), "[", 4)
  group=sapply(strsplit(group,""), "[", 1)
  data = data[,group == 0]
  
  
  methy <- fread(paste("D:/R/Pancancer ISCA1/Linear/new-methy/ISCA1-methy/",cancer,"-Methy450.csv",sep=''))
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
  data.r <- data.corr$r  # 相关系数
  data.p <- data.corr$p  # p值
  
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

class1 <- rep("ferroptosis driver", times=254)
class2 <- rep("ferroptosis marker", times=6)
class3 <- rep("ferroptosis suppressor", times=194)
class4 <- rep("ferroptosis unclassified", times=75)
class<- c(class1,class2,class3,class4)
annotation_row <- data.frame(class)
rownames(annotation_row) <- rownames(data.p)

pdf(paste("all-Fe_methy.pdf",sep=''),width =10,height =120)
pheatmap(data.r,
         color=myColor,
         breaks=myBreaks,
         clustering_method="average", cluster_rows=F,annotation_row = annotation_row, cluster_cols=F,display_numbers=sig.mat
)

dev.off()



### Isoform-related analysis of ferroptosis gene associations

setwd("D:/R/Pancancer ISCA1/Linear/Isoform")

cancers <- c("ACC","BLCA","BRCA","CESC","CHOL","COAD","DLBC","ESCA","GBM","HNSC","KICH","KIRC","KIRP","LAML","LGG","LIHC","LUAD","LUSC","MESO","OV","PAAD","PCPG",
             "PRAD","READ","SARC","SKCM","STAD","TGCT","THCA","THYM","UCEC","UCS","UVM")
Isoformall <- c("ENST00000311534","ENST00000326094","ENST00000375991")

for (Isoform_single in Isoformall) {
  #  Isoform_single <- "ENST00000311534"
  dat.r <- c()
  allgene1 <- c("ISCA1")
  library(data.table)
  
  gene <-  fread(paste("D:/R/Pancancer ISCA1/Linear/new_CNV/gene-all.csv",sep=''))
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
    # Use subset function for filtering
    Isoform <- subset(Isoform, tissue == cancer)
    rownames(Isoform) <- Isoform$sample
    all_name <- names(which(table(c(rownames(data),Isoform$sample))==2))
    data <- data[match(all_name,rownames(data)),]
    Isoform <-  Isoform[match(all_name,Isoform$sample),]
    dat_im <- data[,match(allgene2,colnames(data))]
    dat_im <- dat_im[, colSums(is.na(dat_im)) == 0]
    Isoform <- Isoform$tpm
    library(psych)
    special_cancers <- c("MESO","UVM")
    if (cancer %in% special_cancers) {
      dat.r <- cbind(dat.r,data.r)
      dat.p <- cbind(dat.p,data.p)
    } else {
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
  class1 <- rep("ferroptosis_driver", times=29)
  class2 <- rep("ferroptosis_suppressor", times=15)
  class3 <- rep("ferroptosis_unclassified", times=2)
  class4 <- rep("ferroptosis_marker", times=3)
  class<- c(class1,class2,class3,class4)
  annotation_row <- data.frame(class)
  rownames(annotation_row) <- rownames(data.p)
  pdf(paste("cor-ferroptosis_isoform_",Isoform_single,".pdf",sep=''),width =10,height =10)
  pheatmap(data.r,
           color=myColor,
           breaks=myBreaks,
           clustering_method="average", cluster_rows=F,annotation_row = annotation_row, cluster_cols=F,display_numbers=sig.mat
  )
  dev.off()
}

