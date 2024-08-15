
### RNA-Stemness

setwd("D:/R/Pancancer ISCA1/Stemness/")

cancers <- c("ACC","BLCA","BRCA","CESC","CHOL","COAD","DLBC","ESCA","GBM","HNSC","KICH","KIRC","KIRP","LAML","LGG","LIHC","LUAD","LUSC","MESO","OV","PAAD","PCPG",
             "PRAD","READ","SARC","SKCM","STAD","TGCT","THCA","THYM","UCEC","UCS","UVM")

hets <- c("DMPss","DNAss","ENHss","EREG.EXPss","EREG.METHss","RNAss")
data.r <- c()
data.p <- c()
for (cancer in cancers) {
  r <- c()
  P <- c()
  for(het in hets){
    library(data.table)
    dat <- fread(paste("D:/R/Pancancer ISCA1/Clinical/TCGA_new/",cancer,".txt",sep=''))
    dat <- as.data.frame(dat)
    rownames(dat) <- dat[,1]
    dat <- dat[,-1]
    dat_gene <- dat
    data =dat_gene
    # Normal and tumor numbers, characters 14 and 15, 01-09 are cancers, 10-19 are normal, 20-29 are adjacent to cancer
    group=sapply(strsplit(colnames(data),"\\-"), "[", 4)
    group=sapply(strsplit(group,""), "[", 1)
    # Only retain tumor samples
    dat_gene = data[,group == 0]
    dat_gene <- t(dat_gene)
    datLOH <- read.csv(paste("D:/R/Pancancer ISCA1/Stemness/",het,".csv",sep=''))
    if(cancer=="GBM"){
      n1 <- grep(cancer,datLOH$X)
      dataLOH <- datLOH[n1,]
      n2 <- grep("GBMLGG",dataLOH$X)
      dataLOH <- dataLOH[-n2,]
    }else 
      if(cancer=="LGG"){
        n1 <- grep(cancer,datLOH$X)
        dataLOH <- datLOH[n1,]
        n2 <- grep("GBMLGG",dataLOH$X)
        dataLOH <- dataLOH[-n2,]
      }else 
        if(cancer=="COAD"){
          
          n1 <- grep(cancer,datLOH$X)
          dataLOH <- datLOH[n1,]
          n2 <- grep("COADREAD",dataLOH$X)
          dataLOH <- dataLOH[-n2,]
          
        }else 
          if(cancer=="READ"){
            n1 <- grep(cancer,datLOH$X)
            dataLOH <- datLOH[n1,]
            n2 <- grep("COADREAD",dataLOH$X)
            dataLOH <- dataLOH[-n2,]
          }else{
            n1 <- grep(cancer,datLOH$X)
            dataLOH <- datLOH[n1,]
          }
    dataLOH <- dataLOH[,-1]
    dataLOH <- data.frame(dataLOH)
    all_name <- names(which(table(c(rownames(dat_gene),dataLOH[,3] ))==2))
    if(length(all_name)==0){
      dat.r <- NA
      dat.p <- NA
    }else {
      data_gene <- dat_gene[match(all_name,rownames(dat_gene)),]
      dataLOH1<- dataLOH[match(all_name,dataLOH[,3]),]
      data<- data.frame(data_gene,Heterogeneity=dataLOH1[,1])
      dat_1 <- data[,grep("\\bISCA1\\b",colnames(data),value = T)]
      dat_2 <- data[,grep("Heterogeneity",colnames(data))]
      dat.cor <- cor.test(dat_1,dat_2, method="pearson", adjust="fdr")
      dat.r <- dat.cor$estimate  # Correlation coefficient
      dat.p <- dat.cor$p.value  # p-value
    }
    if(length(r)==0){
      r <- dat.r
      p <- dat.p
    }else {
      r <- c(r,dat.r)
      p <- c(p,dat.p)
    }
  }
  
  if(length(data.r)==0){
    data.r <- r
    data.p <- p
    
  }else {
    data.r <- cbind(data.r,r)
    data.p <- cbind(data.p,p)
    
  }
  
}
colnames(data.r) <- cancers 
rownames(data.r) <- hets
colnames(data.p) <- cancers 
rownames(data.p) <- hets
write.csv(data.r,file= paste("Generation/data.r.csv",sep=''),quote=F)
write.csv(data.p,file= paste("Generation/data.p.csv",sep=''),quote=F)
data.p[is.na(data.p)] = 2

library(cowplot) # For plotting
library(ggplot2) # For plotting
library(RColorBrewer) # For plotting
library(dplyr) # For data processing
library(tidyr) # For data processing

genes <-colnames(data.r)
tumor_types <-rownames((data.r))
gene_tumor_combinations <- expand.grid(CancerALL = genes, MirnaALL = tumor_types)

library(tidyr)
# Assuming your data frame is named df
df=as.data.frame(data.r)
df=t(df)
df=as.data.frame(df)
dfP=as.data.frame(data.p)
dfP=t(dfP)
dfP=as.data.frame(dfP)


# Use the gather function to combine all row data into one column
gathered_dataR <- gather(df, key = "mRNA", value = "cor")
gathered_dataP <- gather(dfP, key = "mRNA", value = "P-Value")

library(dplyr)

# If gathered_dataR is a data frame, extract one column and then combine
combined_data <- cbind(gene_tumor_combinations, Correlation = gathered_dataR$cor)
combined_data <- cbind(combined_data, Pvalue = gathered_dataP$`P-Value`)
merged_data <-combined_data
colnames(merged_data)=c("Gene","Tumor_Type","Correlation","P_Value")

# Data processing
merged_data$fdr_group <- cut(merged_data$P_Value, breaks = c(0, 0.0001, 0.001, 0.01, 0.05, 0.1, 1), labels = c("<0.0001", "<0.001", "<0.01", "<0.05", "<0.1",">0.1"))
merged_data$fdr_group <- factor(merged_data$fdr_group, levels = c(labels = c("<0.0001","<0.001", "<0.01", "<0.05", "<0.1",">0.1")))
merged_data$fdr_group_2 <- cut(merged_data$P_Value, breaks = c(0, 0.0001, 0.001, 0.01, 0.05, 1), labels = c("<0.0001","<0.001", "<0.01", "<0.05", ">0.05"))
merged_data$fdr_group_2 <- factor(merged_data$fdr_group_2, levels = c(labels = c("<0.0001","<0.001", "<0.01", "<0.05", ">0.05")))

merged_data$Gene <- factor(merged_data$Gene, levels = unique(merged_data$Gene))

ggplot2::ggplot() +
  geom_hline(yintercept = tumor_types, color = "#E8E8E8") +
  geom_point(data = subset(merged_data, Correlation < 0), shape = 19, stroke = 0,
             aes(x = Gene, y = Tumor_Type, size = abs(Correlation), color = fdr_group)) +
  scale_x_discrete(drop = FALSE) +
  geom_point(data = subset(merged_data, Correlation > 0), shape = 21, stroke = 0.1,
             aes(x = Gene, y = Tumor_Type, size = abs(Correlation), fill = fdr_group_2)) +
  scale_x_discrete(drop = FALSE)+
  scale_size_continuous(limits = c(0, 1),breaks = c(0.1, 0.2, 0.4, 0.6, 0.8, 1))+
  scale_fill_manual(values = c(rev(brewer.pal(9,'Reds'))[5:8],"#E8E8E8"))+
  scale_color_manual(values = c(rev(brewer.pal(9,'YlGnBu'))[2:5],"#E8E8E8"))+
  cowplot::theme_cowplot() + 
  theme(panel.grid.major = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
        plot.title = element_text(hjust = 0.5))+
  labs(color = "Negative\ncorrelation\ p-value", size = "R-value",
       fill = "Positive\ncorrelation\ p-value",x = "",y = "",
       title = "Correlation of ISCA1 with Stemness") +
  guides(fill = guide_legend(override.aes = list(size = 4),order = 3),
         color = guide_legend(override.aes = list(size = 4),order = 2),
         size = guide_legend(override.aes = list(size = c(2:7), fill = "white"),order = 1,))

ggsave(paste("1-Stemness-Pancancer-cor-ISCA1-scatter plot.pdf"), plot = last_plot(), device = "pdf", width = 12, height = 6)


### RNA-Heterogeneity


setwd("D:/R/Pancancer ISCA1/Heterogeneity/")

cancers <- c("ACC","BLCA","BRCA","CESC","CHOL","COAD","DLBC","ESCA","GBM","HNSC","KICH","KIRC","KIRP","LAML","LGG","LIHC","LUAD","LUSC","MESO","OV","PAAD","PCPG",
             "PRAD","READ","SARC","SKCM","STAD","TGCT","THCA","THYM","UCEC","UCS","UVM")
hets <- c("LOH","MATH","MSI","Neoantigen","ploidy","purity","TMB","HRD")
data.r <- c()
data.p <- c()
for (cancer in cancers) {
  r <- c()
  P <- c()
  for(het in hets){
    library(data.table)
    dat <- fread(paste("D:/R/Pancancer ISCA1/Clinical/TCGA_new/",cancer,".txt",sep=''))
    #dat <- fread("D:/R/Pancancer ISCA1/Clinical/TCGA_new/",cancer,".txt",sep='')
    dat <- as.data.frame(dat)
    rownames(dat) <- dat[,1]
    dat <- dat[,-1]
    dat_gene <- t(dat)
    datLOH <- read.csv(paste("D:/R/Pancancer ISCA1/Heterogeneity/",het,".csv",sep=''))
    if(cancer=="GBM"){
      
      n1 <- grep(cancer,datLOH$X)
      dataLOH <- datLOH[n1,]
      n2 <- grep("GBMLGG",dataLOH$X)
      dataLOH <- dataLOH[-n2,]
      
    }else 
      #cancer <- "LGG"
      if(cancer=="LGG"){
        
        
        n1 <- grep(cancer,datLOH$X)
        dataLOH <- datLOH[n1,]
        n2 <- grep("GBMLGG",dataLOH$X)
        dataLOH <- dataLOH[-n2,]
        
        
      }else 
        #cancer <- "COAD"
        if(cancer=="COAD"){
          
          n1 <- grep(cancer,datLOH$X)
          dataLOH <- datLOH[n1,]
          n2 <- grep("COADREAD",dataLOH$X)
          dataLOH <- dataLOH[-n2,]
          
        }else 
          #cancer <- "READ"
          if(cancer=="READ"){
            
            n1 <- grep(cancer,datLOH$X)
            dataLOH <- datLOH[n1,]
            n2 <- grep("COADREAD",dataLOH$X)
            dataLOH <- dataLOH[-n2,]
            
          }else{
            n1 <- grep(cancer,datLOH$X)
            dataLOH <- datLOH[n1,]
          }
    dataLOH <- dataLOH[,-1]
    dataLOH <- data.frame(dataLOH)
    all_name <- names(which(table(c(rownames(dat_gene),dataLOH[,3] ))==2))
    if(length(all_name)==0){
      
      dat.r <- NA
      dat.p <- NA
    }else {
      data_gene <- dat_gene[match(all_name,rownames(dat_gene)),]
      dataLOH1<- dataLOH[match(all_name,dataLOH[,3]),]
      dat_1 <- data[,grep("\\bISCA1\\b",colnames(data),value = T)]
      dat_2 <- data[,grep("Heterogeneity",colnames(data))]
      dat.cor <- cor.test(dat_1,dat_2, method="pearson", adjust="fdr")
      dat.r <- dat.cor$estimate  # Correlation coefficient
      dat.p <- dat.cor$p.value  # p-value
    }
    if(length(r)==0){
      
      r <- dat.r
      p <- dat.p
      
      
    }else {
      r <- c(r,dat.r)
      p <- c(p,dat.p)
    }
    
  }
  
  if(length(data.r)==0){
    data.r <- r
    data.p <- p
    
  }else {
    data.r <- cbind(data.r,r)
    data.p <- cbind(data.p,p)
    
  }
  
}

colnames(data.r) <- cancers 
rownames(data.r) <- hets

colnames(data.p) <- cancers 
rownames(data.p) <- hets


write.csv(data.r,file= paste("Generation/data.r.csv",sep=''),quote=F)
write.csv(data.p,file= paste("Generation/data.p.csv",sep=''),quote=F)
data.p[is.na(data.p)] = 2
library(cowplot) # For plotting
library(ggplot2) # For plotting
library(RColorBrewer) # For plotting
library(dplyr) # For data processing
library(tidyr) # For data processing

genes <-colnames(data.r)
tumor_types <-rownames((data.r))
gene_tumor_combinations <- expand.grid(CancerALL = genes, MirnaALL = tumor_types)



library(tidyr)
# Assume your data frame is named df
df=as.data.frame(data.r)
df=t(df)
df=as.data.frame(df)
dfP=as.data.frame(data.p)
dfP=t(dfP)
dfP=as.data.frame(dfP)


# Use the gather function to combine all row data into one column
gathered_dataR <- gather(df, key = "mRNA", value = "cor")
gathered_dataP <- gather(dfP, key = "mRNA", value = "P-Value")

library(dplyr)

# If gathered_dataR is a data frame, extract one column and then combine
combined_data <- cbind(gene_tumor_combinations, Correlation = gathered_dataR$cor)
combined_data <- cbind(combined_data, Pvalue = gathered_dataP$`P-Value`)
merged_data <-combined_data
colnames(merged_data)=c("Gene","Tumor_Type","Correlation","P_Value")



# Data processing
merged_data$fdr_group <- cut(merged_data$P_Value, breaks = c(0, 0.0001, 0.001, 0.01, 0.05,  1), labels = c("<0.0001", "<0.001", "<0.01", "<0.05", ">0.05"))
merged_data$fdr_group <- factor(merged_data$fdr_group, levels = c(labels = c("<0.0001","<0.001", "<0.01", "<0.05", ">0.05")))
merged_data$fdr_group_2 <- cut(merged_data$P_Value, breaks = c(0, 0.0001, 0.001, 0.01, 0.05, 1), labels = c("<0.0001","<0.001", "<0.01", "<0.05", ">0.05"))
merged_data$fdr_group_2 <- factor(merged_data$fdr_group_2, levels = c(labels = c("<0.0001","<0.001", "<0.01", "<0.05", ">0.05")))
merged_data$Gene <- factor(merged_data$Gene, levels = unique(merged_data$Gene))
ggplot2::ggplot() +
  geom_hline(yintercept = tumor_types, color = "#E8E8E8") +
  geom_point(data = subset(merged_data, Correlation < 0), shape = 19, stroke = 0,
             aes(x = Gene, y = Tumor_Type, size = abs(Correlation), color = fdr_group)) +
  scale_x_discrete(drop = FALSE) +
  geom_point(data = subset(merged_data, Correlation > 0), shape = 21, stroke = 0.1,
             aes(x = Gene, y = Tumor_Type, size = abs(Correlation), fill = fdr_group_2)) +
  scale_x_discrete(drop = FALSE)+
  scale_size_continuous(limits = c(0, 1),breaks = c(0.1, 0.2, 0.4, 0.6, 0.8, 1))+
  scale_fill_manual(values = c(rev(brewer.pal(9,'Reds'))[1:4],"#E8E8E8"))+
  scale_color_manual(values = c(rev(brewer.pal(9,'Greens'))[1:4],"#E8E8E8"))+
  cowplot::theme_cowplot() + 
  theme(panel.grid.major = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
        plot.title = element_text(hjust = 0.5))+
  labs(color = "Negative\ncorrelation\ p-value", size = "R-value",
       fill = "Positive\ncorrelation\ p-value",x = "",y = "",
       title = "Correlation of ISCA1 with Heterogeneity") +
  guides(fill = guide_legend(override.aes = list(size = 4),order = 3),
         color = guide_legend(override.aes = list(size = 4),order = 2),
         size = guide_legend(override.aes = list(size = c(2:7), fill = "white"),order = 1,))
ggsave(paste("Heterogeneity-Pancancer-ISCA1-scatterplot.pdf"), plot = last_plot(), device = "pdf", width = 12, height = 6)


### CNV-Stemness

setwd("D:/R/Pancancer ISCA1/Stemness/")

cancers <- c("ACC","BLCA","BRCA","CESC","CHOL","COAD","DLBC","ESCA","GBM","HNSC","KICH","KIRC","KIRP","LAML","LGG","LIHC","LUAD","LUSC","MESO","OV","PAAD","PCPG",
             "PRAD","READ","SARC","SKCM","STAD","TGCT","THCA","THYM","UCEC","UCS","UVM")
hets <- c("DMPss","DNAss","ENHss","EREG.EXPss","EREG.METHss","RNAss")
data.r <- c()
data.p <- c()

for (cancer in cancers) {
  r <- c()
  P <- c()
  for(het in hets){
    # CNV data
    folder_path <- "D:/R/Pancancer ISCA1/CNV-SNV-methy/CNV/data/"
    #cancer <- "ACC"
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
    
    dat_gene <- dataMethy
    datLOH <- read.csv(paste("D:/R/Pancancer ISCA1/Stemness/",het,".csv",sep=''))
    if(cancer=="GBM"){
      
      n1 <- grep(cancer,datLOH$X)
      dataLOH <- datLOH[n1,]
      n2 <- grep("GBMLGG",dataLOH$X)
      dataLOH <- dataLOH[-n2,]
      
    }else 
      #cancer <- "LGG"
      if(cancer=="LGG"){
        
        
        n1 <- grep(cancer,datLOH$X)
        dataLOH <- datLOH[n1,]
        n2 <- grep("GBMLGG",dataLOH$X)
        dataLOH <- dataLOH[-n2,]
        
        
      }else 
        #cancer <- "COAD"
        if(cancer=="COAD"){
          
          n1 <- grep(cancer,datLOH$X)
          dataLOH <- datLOH[n1,]
          n2 <- grep("COADREAD",dataLOH$X)
          dataLOH <- dataLOH[-n2,]
          
        }else 
          #cancer <- "READ"
          if(cancer=="READ"){
            
            n1 <- grep(cancer,datLOH$X)
            dataLOH <- datLOH[n1,]
            n2 <- grep("COADREAD",dataLOH$X)
            dataLOH <- dataLOH[-n2,]
            
          }else{
            n1 <- grep(cancer,datLOH$X)
            dataLOH <- datLOH[n1,]
          }
    dataLOH <- dataLOH[,-1]
    dataLOH <- data.frame(dataLOH)
    all_name <- names(which(table(c(rownames(dat_gene),dataLOH[,3] ))==2))
    if(length(all_name)==0){
      dat.r <- NA
      dat.p <- NA
    }else {
      data_gene <- dat_gene[match(all_name,rownames(dat_gene)),]
      dataLOH1<- dataLOH[match(all_name,dataLOH[,3]),]
      data<- data.frame(data_gene,Heterogeneity=dataLOH1[,1])
      dat_1 <- data[,grep("\\bISCA1\\b",colnames(data),value = T)]
      dat_2 <- data[,grep("Heterogeneity",colnames(data))]
      dat.cor <- cor.test(dat_1,dat_2, method="pearson", adjust="fdr")
      dat.r <- dat.cor$estimate  # Correlation coefficient
      dat.p <- dat.cor$p.value  # p-value
    }
    if(length(r)==0){
      
      r <- dat.r
      p <- dat.p
    }else {
      r <- c(r,dat.r)
      p <- c(p,dat.p)
    }
    
  }
  if(length(data.r)==0){
    data.r <- r
    data.p <- p
    
  }else {
    data.r <- cbind(data.r,r)
    data.p <- cbind(data.p,p)
    
  }
}
colnames(data.r) <- cancers 
rownames(data.r) <- hets

colnames(data.p) <- cancers 
rownames(data.p) <- hets
write.csv(data.r,file= paste("Generation/data.r.csv",sep=''),quote=F)
write.csv(data.p,file= paste("Generation/data.p.csv",sep=''),quote=F)
data.p[is.na(data.p)] = 2

library(cowplot) # For plotting
library(ggplot2) # For plotting
library(RColorBrewer) # For plotting
library(dplyr) # For data processing
library(tidyr) # For data processing

genes <-colnames(data.r)
tumor_types <-rownames((data.r))
gene_tumor_combinations <- expand.grid(CancerALL = genes, MirnaALL = tumor_types)
library(tidyr)
# Assuming your data frame is named df
df=as.data.frame(data.r)
df=t(df)
df=as.data.frame(df)
dfP=as.data.frame(data.p)
dfP=t(dfP)
dfP=as.data.frame(dfP)
# Use the gather function to combine all row data into one column
gathered_dataR <- gather(df, key = "mRNA", value = "cor")
gathered_dataP <- gather(dfP, key = "mRNA", value = "P-Value")

library(dplyr)

# If gathered_dataR is a data frame, extract one column and then combine
combined_data <- cbind(gene_tumor_combinations, Correlation = gathered_dataR$cor)
combined_data <- cbind(combined_data, Pvalue = gathered_dataP$`P-Value`)
merged_data <-combined_data
colnames(merged_data)=c("Gene","Tumor_Type","Correlation","P_Value")
# Data processing
merged_data$fdr_group <- cut(merged_data$P_Value, breaks = c(0, 0.0001, 0.001, 0.01, 0.05,  1), labels = c("<0.0001", "<0.001", "<0.01", "<0.05", ">0.05"))
merged_data$fdr_group <- factor(merged_data$fdr_group, levels = c(labels = c("<0.0001","<0.001", "<0.01", "<0.05", ">0.05")))
merged_data$fdr_group_2 <- cut(merged_data$P_Value, breaks = c(0, 0.0001, 0.001, 0.01, 0.05, 1), labels = c("<0.0001","<0.001", "<0.01", "<0.05", ">0.05"))
merged_data$fdr_group_2 <- factor(merged_data$fdr_group_2, levels = c(labels = c("<0.0001","<0.001", "<0.01", "<0.05", ">0.05")))
ggplot2::ggplot() +
  geom_hline(yintercept = tumor_types, color = "#E8E8E8") +
  geom_point(data = subset(merged_data, Correlation < 0), shape = 19, stroke = 0,
             aes(x = Gene, y = Tumor_Type, size = abs(Correlation), color = fdr_group)) +
  scale_x_discrete(drop = FALSE) +
  geom_point(data = subset(merged_data, Correlation > 0), shape = 21, stroke = 0.1,
             aes(x = Gene, y = Tumor_Type, size = abs(Correlation), fill = fdr_group_2)) +
  scale_x_discrete(drop = FALSE)+
  scale_size_continuous(limits = c(0, 1),breaks = c(0.1, 0.2, 0.4, 0.6, 0.8, 1))+
  scale_fill_manual(values = c(rev(brewer.pal(9,'Reds'))[5:8],"#E8E8E8"))+
  scale_color_manual(values = c(rev(brewer.pal(9,'YlGnBu'))[2:5],"#E8E8E8"))+
  cowplot::theme_cowplot() + 
  theme(panel.grid.major = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
        plot.title = element_text(hjust = 0.5))+
  labs(color = "Negative\ncorrelation\ p-value", size = "R-value",
       fill = "Positive\ncorrelation\ p-value",x = "",y = "",
       title = "Correlation of CNV with Stemness") +
  guides(fill = guide_legend(override.aes = list(size = 4),order = 3),
         color = guide_legend(override.aes = list(size = 4),order = 2),
         size = guide_legend(override.aes = list(size = c(2:7), fill = "white"),order = 1,))
ggsave(paste("1-Stemness-Pancancer-CNV-scatter plot.pdf"), plot = last_plot(), device = "pdf", width = 12, height = 6)


### CNV-Heterogeneity




setwd("D:/R/Pancancer ISCA1/Heterogeneity/")

cancers <- c("ACC","BLCA","BRCA","CESC","CHOL","COAD","DLBC","ESCA","GBM","HNSC","KICH","KIRC","KIRP","LAML","LGG","LIHC","LUAD","LUSC","MESO","OV","PAAD","PCPG",
             "PRAD","READ","SARC","SKCM","STAD","TGCT","THCA","THYM","UCEC","UCS","UVM")

hets <- c("LOH","MATH","MSI","Neoantigen","ploidy","purity","TMB","HRD")
data.r <- c()
data.p <- c()
for (cancer in cancers) {
  r <- c()
  P <- c()
  for(het in hets){
    
    library(data.table)
    # CNV data
    folder_path <- "D:/R/Pancancer ISCA1/CNV-SNV-methylation/CNV/data/"
    #cancer <- "ACC"
    # Use file.path() function to construct the file path
    new_file_path <- file.path(folder_path, paste0(cancer, ".cnv.tsv.gz"))
    gz_file <- gzfile(new_file_path, "rt")  # Open the.gz file
    # Use read.table() function to read the data
    dataMethy <- read.table(gz_file, header = TRUE, sep = "\t")  # Take the tab delimiter as an example, you can choose the delimiter according to the actual situation
    # Delete the row data with duplicate names in the 3rd column
    duplicated_rows <- duplicated(dataMethy[,3])
    dataMethy <- dataMethy[!duplicated_rows, ]
    rownames(dataMethy)=dataMethy$symbol
    dataMethy <-dataMethy[,-1:-3]
    colnames(dataMethy) <- gsub("\\.", "-", colnames(dataMethy))
    dataMethy=t(dataMethy)
    original_id <-rownames(dataMethy)
    modified_id <- substr(original_id, 1, 15)
    rownames(dataMethy)=modified_id
    dat_gene <- dataMethy
    datLOH <- read.csv(paste("D:/R/Pancancer ISCA1/Heterogeneity/",het,".csv",sep=''))
    if(cancer=="GBM"){
      
      n1 <- grep(cancer,datLOH$X)
      dataLOH <- datLOH[n1,]
      n2 <- grep("GBMLGG",dataLOH$X)
      dataLOH <- dataLOH[-n2,]
      
    }else 
      #cancer <- "LGG"
      if(cancer=="LGG"){
        
        
        n1 <- grep(cancer,datLOH$X)
        dataLOH <- datLOH[n1,]
        n2 <- grep("GBMLGG",dataLOH$X)
        dataLOH <- dataLOH[-n2,]
        
        
      }else 
        #cancer <- "COAD"
        if(cancer=="COAD"){
          
          n1 <- grep(cancer,datLOH$X)
          dataLOH <- datLOH[n1,]
          n2 <- grep("COADREAD",dataLOH$X)
          dataLOH <- dataLOH[-n2,]
          
        }else 
          #cancer <- "READ"
          if(cancer=="READ"){
            
            n1 <- grep(cancer,datLOH$X)
            dataLOH <- datLOH[n1,]
            n2 <- grep("COADREAD",dataLOH$X)
            dataLOH <- dataLOH[-n2,]
            
          }else{
            n1 <- grep(cancer,datLOH$X)
            dataLOH <- datLOH[n1,]
          }
    dataLOH <- dataLOH[,-1]
    dataLOH <- data.frame(dataLOH)
    all_name <- names(which(table(c(rownames(dat_gene),dataLOH[,3] ))==2))
    if(length(all_name)==0){
      
      dat.r <- NA
      dat.p <- NA
    }else {
      data_gene <- dat_gene[match(all_name,rownames(dat_gene)),]
      dataLOH1<- dataLOH[match(all_name,dataLOH[,3]),]
      data<- data.frame(data_gene,Heterogeneity=dataLOH1[,1])
      dat_1 <- data[,grep("\\bISCA1\\b",colnames(data),value = T)]
      dat_2 <- data[,grep("Heterogeneity",colnames(data))]
      dat.cor <- cor.test(dat_1,dat_2, method="pearson", adjust="fdr")
      dat.r <- dat.cor$estimate  # Correlation coefficient
      dat.p <- dat.cor$p.value  # p value
    }
    if(length(r)==0){
      
      r <- dat.r
      p <- dat.p
      
      
    }else {
      r <- c(r,dat.r)
      p <- c(p,dat.p)
    }
    
  }
  
  if(length(data.r)==0){
    data.r <- r
    data.p <- p
    
  }else {
    data.r <- cbind(data.r,r)
    data.p <- cbind(data.p,p)
  }
  
}
colnames(data.r) <- cancers 
rownames(data.r) <- hets

colnames(data.p) <- cancers 
rownames(data.p) <- hets

write.csv(data.r,file= paste("Generated/data.r.csv",sep=''),quote=F)
write.csv(data.p,file= paste("Generated/data.p.csv",sep=''),quote=F)
data.p[is.na(data.p)] = 2
library(cowplot) # For plotting
library(ggplot2) # For plotting
library(RColorBrewer) # For plotting
library(dplyr) # For data processing
library(tidyr) # For data processing
genes <-colnames(data.r)
tumor_types <-rownames((data.r))
gene_tumor_combinations <- expand.grid(CancerALL = genes, MirnaALL = tumor_types)
library(tidyr)
# Assume your data frame is named df
df=as.data.frame(data.r)
df=t(df)
df=as.data.frame(df)
dfP=as.data.frame(data.p)
dfP=t(dfP)
dfP=as.data.frame(dfP)
# Use gather function to combine all row data into one column
gathered_dataR <- gather(df, key = "mRNA", value = "cor")
gathered_dataP <- gather(dfP, key = "mRNA", value = "P-Value")
library(dplyr)
# If gathered_dataR is a data frame, extract one column and then combine
combined_data <- cbind(gene_tumor_combinations, Correlation = gathered_dataR$cor)
combined_data <- cbind(combined_data, Pvalue = gathered_dataP$`P-Value`)
merged_data <-combined_data
colnames(merged_data)=c("Gene","Tumor_Type","Correlation","P_Value")
# Data processing
merged_data$fdr_group <- cut(merged_data$P_Value, breaks = c(0, 0.0001, 0.001, 0.01, 0.05,  1), labels = c("<0.0001", "<0.001", "<0.01", "<0.05", ">0.05"))
merged_data$fdr_group <- factor(merged_data$fdr_group, levels = c(labels = c("<0.0001","<0.001", "<0.01", "<0.05", ">0.05")))
merged_data$fdr_group_2 <- cut(merged_data$P_Value, breaks = c(0, 0.0001, 0.001, 0.01, 0.05, 1), labels = c("<0.0001","<0.001", "<0.01", "<0.05", ">0.05"))
merged_data$fdr_group_2 <- factor(merged_data$fdr_group_2, levels = c(labels = c("<0.0001","<0.001", "<0.01", "<0.05", ">0.05")))
merged_data$Gene <- factor(merged_data$Gene, levels = unique(merged_data$Gene))
ggplot2::ggplot() +
  geom_hline(yintercept = tumor_types, color = "#E8E8E8") +
  geom_point(data = subset(merged_data, Correlation < 0), shape = 19, stroke = 0,
             aes(x = Gene, y = Tumor_Type, size = abs(Correlation), color = fdr_group)) +
  scale_x_discrete(drop = FALSE) +
  geom_point(data = subset(merged_data, Correlation > 0), shape = 21, stroke = 0.1,
             aes(x = Gene, y = Tumor_Type, size = abs(Correlation), fill = fdr_group_2)) +
  scale_x_discrete(drop = FALSE)+
  scale_size_continuous(limits = c(0, 1),breaks = c(0.1, 0.2, 0.4, 0.6, 0.8, 1))+
  scale_fill_manual(values = c(rev(brewer.pal(9,'Reds'))[1:4],"#E8E8E8"))+
  scale_color_manual(values = c(rev(brewer.pal(9,'Greens'))[1:4],"#E8E8E8"))+
  cowplot::theme_cowplot() + 
  theme(panel.grid.major = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
        plot.title = element_text(hjust = 0.5))+
  labs(color = "Negative\ncorrelation\ p-value", size = "R-value",
       fill = "Positive\ncorrelation\ p-value",x = "",y = "",
       title = "Correlation of CNV with heterogeneity") +
  guides(fill = guide_legend(override.aes = list(size = 4),order = 3),
         color = guide_legend(override.aes = list(size = 4),order = 2),
         size = guide_legend(override.aes = list(size = c(2:7), fill = "white"),order = 1,))
ggsave(paste("2-Heterogeneity-Pancancer-CNV-scatter plot.pdf"), plot = last_plot(), device = "pdf", width = 12, height = 6)


### Isoform-Stemness


rm (list = ls ())

setwd("D:/R/Pancancer ISCA1/Stemness/")
cancers <- c("ACC","BLCA","BRCA","CESC","CHOL","COAD","DLBC","ESCA","GBM","HNSC","KICH","KIRC","KIRP","LAML","LGG","LIHC","LUAD","LUSC","MESO","OV","PAAD","PCPG",
             "PRAD","READ","SARC","SKCM","STAD","TGCT","THCA","THYM","UCEC","UCS","UVM")
hets <- c("DMPss","DNAss","ENHss","EREG.EXPss","EREG.METHss","RNAss")
Isoform_singles <- c("ENST00000375991","ENST00000311534","ENST00000326094")
for (Isoform_single in Isoform_singles) {
  data.r <- c()
  data.p <- c()
  for (cancer in cancers) {
    r <- c()
    P <- c()
    for(het in hets){
      library(data.table)
      Isoform <- read.csv(paste0("D:/R/Pancancer ISCA1/Isoform/",Isoform_single,"_transcript_pancan_dist.CSV"))
      # Use the subset function for filtering
      Isoform <- subset(Isoform, tissue == cancer)
      rownames(Isoform) <- Isoform$sample
      dat_gene <- Isoform
      # Create a logical condition to select rows where the value in the fifth column is "tumor"
      condition <- dat_gene$type2 == "tumor"
      # Filter the row data that meets the condition
      dat_gene <- dat_gene[condition, ]
      datLOH <- read.csv(paste("D:/R/Pancancer ISCA1/Stemness/",het,".csv",sep=''))
      if(cancer=="GBM"){
        n1 <- grep(cancer,datLOH$X)
        dataLOH <- datLOH[n1,]
        n2 <- grep("GBMLGG",dataLOH$X)
        dataLOH <- dataLOH[-n2,]
      }else 
        #cancer <- "LGG"
        if(cancer=="LGG"){
          n1 <- grep(cancer,datLOH$X)
          dataLOH <- datLOH[n1,]
          n2 <- grep("GBMLGG",dataLOH$X)
          dataLOH <- dataLOH[-n2,]
        }else 
          #cancer <- "COAD"
          if(cancer=="COAD"){
            n1 <- grep(cancer,datLOH$X)
            dataLOH <- datLOH[n1,]
            n2 <- grep("COADREAD",dataLOH$X)
            dataLOH <- dataLOH[-n2,]
            
          }else 
            #cancer <- "READ"
            if(cancer=="READ"){
              
              n1 <- grep(cancer,datLOH$X)
              dataLOH <- datLOH[n1,]
              n2 <- grep("COADREAD",dataLOH$X)
              dataLOH <- dataLOH[-n2,]
              
            }else{
              n1 <- grep(cancer,datLOH$X)
              dataLOH <- datLOH[n1,]
            }
      dataLOH <- dataLOH[,-1]
      dataLOH <- data.frame(dataLOH)
      all_name <- names(which(table(c(rownames(dat_gene),dataLOH[,3] ))==2))
      if(length(all_name)==0){
        dat.r <- NA
        dat.p <- NA
      }else {
        data_gene <- dat_gene[match(all_name,rownames(dat_gene)),]
        dataLOH1<- dataLOH[match(all_name,dataLOH[,3]),]
        data<- data.frame(data_gene,Heterogeneity=dataLOH1[,1])
        dat_1 <- data[,grep("tpm",colnames(data),value = T)]
        dat_2 <- data[,grep("Heterogeneity",colnames(data))]
        dat.cor <- cor.test(dat_1,dat_2, method="pearson", adjust="fdr")
        dat.r <- dat.cor$estimate  # Correlation coefficient
        dat.p <- dat.cor$p.value  # p-value
      }
      if(length(r)==0){
        
        r <- dat.r
        p <- dat.p
        
        
      }else {
        r <- c(r,dat.r)
        p <- c(p,dat.p)
      }
      
    }
    
    if(length(data.r)==0){
      data.r <- r
      data.p <- p
      
    }else {
      data.r <- cbind(data.r,r)
      data.p <- cbind(data.p,p)
      
    }
    
  }
  colnames(data.r) <- cancers 
  rownames(data.r) <- hets
  
  colnames(data.p) <- cancers 
  rownames(data.p) <- hets
  write.csv(data.r,file= paste("Generation/data.r.csv",sep=''),quote=F)
  write.csv(data.p,file= paste("Generation/data.p.csv",sep=''),quote=F)
  data.p[is.na(data.p)] = 2
  library(cowplot) # For plotting
  library(ggplot2) # For plotting
  library(RColorBrewer) # For plotting
  library(dplyr) # For data processing
  library(tidyr) # For data processing
  
  genes <-colnames(data.r)
  tumor_types <-rownames((data.r))
  gene_tumor_combinations <- expand.grid(CancerALL = genes, MirnaALL = tumor_types)
  library(tidyr)
  # Assume your data frame is named df
  df=as.data.frame(data.r)
  df=t(df)
  df=as.data.frame(df)
  dfP=as.data.frame(data.p)
  dfP=t(dfP)
  dfP=as.data.frame(dfP)
  # Use the gather function to combine all row data into one column
  gathered_dataR <- gather(df, key = "mRNA", value = "cor")
  gathered_dataP <- gather(dfP, key = "mRNA", value = "P-Value")
  
  library(dplyr)
  
  # If gathered_dataR is a data frame, extract one column and then combine
  combined_data <- cbind(gene_tumor_combinations, Correlation = gathered_dataR$cor)
  combined_data <- cbind(combined_data, Pvalue = gathered_dataP$`P-Value`)
  merged_data <-combined_data
  colnames(merged_data)=c("Gene","Tumor_Type","Correlation","P_Value")
  # Data processing
  merged_data$fdr_group <- cut(merged_data$P_Value, breaks = c(0, 0.0001, 0.001, 0.01, 0.05,  1), labels = c("<0.0001", "<0.001", "<0.01", "<0.05", ">0.05"))
  merged_data$fdr_group <- factor(merged_data$fdr_group, levels = c(labels = c("<0.0001","<0.001", "<0.01", "<0.05", ">0.05")))
  merged_data$fdr_group_2 <- cut(merged_data$P_Value, breaks = c(0, 0.0001, 0.001, 0.01, 0.05, 1), labels = c("<0.0001","<0.001", "<0.01", "<0.05", ">0.05"))
  merged_data$fdr_group_2 <- factor(merged_data$fdr_group_2, levels = c(labels = c("<0.0001","<0.001", "<0.01", "<0.05", ">0.05")))
  
  ggplot2::ggplot() +
    geom_hline(yintercept = tumor_types, color = "#E8E8E8") +
    geom_point(data = subset(merged_data, Correlation < 0), shape = 19, stroke = 0,
               aes(x = Gene, y = Tumor_Type, size = abs(Correlation), color = fdr_group)) +
    scale_x_discrete(drop = FALSE) +
    geom_point(data = subset(merged_data, Correlation > 0), shape = 21, stroke = 0.1,
               aes(x = Gene, y = Tumor_Type, size = abs(Correlation), fill = fdr_group_2)) +
    scale_x_discrete(drop = FALSE)+
    scale_size_continuous(limits = c(0, 1),breaks = c(0.1, 0.2, 0.4, 0.6, 0.8, 1))+
    scale_fill_manual(values = c(rev(brewer.pal(9,'Reds'))[5:8],"#E8E8E8"))+
    scale_color_manual(values = c(rev(brewer.pal(9,'YlGnBu'))[2:5],"#E8E8E8"))+
    cowplot::theme_cowplot() +
    theme(panel.grid.major = element_blank(),
          axis.text.x = element_text(angle = 45, hjust = 1),
          panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
          plot.title = element_text(hjust = 0.5))+
    labs(color = "Negative\ncorrelation\ p-value", size = "R-value",
         fill = "Positive\ncorrelation\ p-value",x = "",y = "",
         title = paste0("Correlation of ",Isoform_single," with Stemness")) +
    guides(fill = guide_legend(override.aes = list(size = 4),order = 3),
           color = guide_legend(override.aes = list(size = 4),order = 2),
           size = guide_legend(override.aes = list(size = c(2:7), fill = "white"),order = 1,))
  ggsave(paste("1-",Isoform_single,"-Pancancer-Stemness-bubbleplot.pdf",sep=''), plot = last_plot(), device = "pdf", width = 12, height = 6)
}


### Isoform-Heterogeneity


rm (list = ls ())
setwd("D:/R/Pancancer ISCA1/Heterogeneity/")
cancers <- c("ACC","BLCA","BRCA","CESC","CHOL","COAD","DLBC","ESCA","GBM","HNSC","KICH","KIRC","KIRP","LAML","LGG","LIHC","LUAD","LUSC","MESO","OV","PAAD","PCPG",
             "PRAD","READ","SARC","SKCM","STAD","TGCT","THCA","THYM","UCEC","UCS","UVM")
hets <- c("LOH","MATH","MSI","Neoantigen","ploidy","purity","TMB","HRD")
Isoform_singles <- c("ENST00000375991","ENST00000311534","ENST00000326094")
for (Isoform_single in Isoform_singles) {
  data.r <- c()
  data.p <- c()
  for (cancer in cancers) {
    r <- c()
    P <- c()
    for(het in hets){
      library(data.table)
      Isoform <- read.csv(paste0("D:/R/Pancancer ISCA1/Isoform/",Isoform_single,"_transcript_pancan_dist.CSV"))
      # Use the subset function for filtering
      Isoform <- subset(Isoform, tissue == cancer)
      rownames(Isoform) <- Isoform$sample
      dat_gene <- Isoform
      # Create a logical condition to select rows where the value in the fifth column is "tumor"
      condition <- dat_gene$type2 == "tumor"
      # Filter the row data that meets the condition
      dat_gene <- dat_gene[condition, ]
      datLOH <- read.csv(paste("D:/R/Pancancer ISCA1/Heterogeneity/",het,".csv",sep=''))
      if(cancer=="GBM"){
        n1 <- grep(cancer,datLOH$X)
        dataLOH <- datLOH[n1,]
        n2 <- grep("GBMLGG",dataLOH$X)
        dataLOH <- dataLOH[-n2,]
      }else 
        if(cancer=="LGG"){
          n1 <- grep(cancer,datLOH$X)
          dataLOH <- datLOH[n1,]
          n2 <- grep("GBMLGG",dataLOH$X)
          dataLOH <- dataLOH[-n2,]
        }else 
          if(cancer=="COAD"){
            n1 <- grep(cancer,datLOH$X)
            dataLOH <- datLOH[n1,]
            n2 <- grep("COADREAD",dataLOH$X)
            dataLOH <- dataLOH[-n2,]
          }else 
            if(cancer=="READ"){
              n1 <- grep(cancer,datLOH$X)
              dataLOH <- datLOH[n1,]
              n2 <- grep("COADREAD",dataLOH$X)
              dataLOH <- dataLOH[-n2,]
            }else{
              n1 <- grep(cancer,datLOH$X)
              dataLOH <- datLOH[n1,]
            }
      dataLOH <- dataLOH[,-1]
      dataLOH <- data.frame(dataLOH)
      all_name <- names(which(table(c(rownames(dat_gene),dataLOH[,3] ))==2))
      if(length(all_name)==0){
        dat.r <- NA
        dat.p <- NA
      }else {
        data_gene <- dat_gene[match(all_name,rownames(dat_gene)),]
        dataLOH1<- dataLOH[match(all_name,dataLOH[,3]),]
        data<- data.frame(data_gene,Heterogeneity=dataLOH1[,1])
        dat_1 <- data[,grep("tpm",colnames(data),value = T)]
        dat_2 <- data[,grep("Heterogeneity",colnames(data))]
        dat.cor <- cor.test(dat_1,dat_2, method="pearson", adjust="fdr")
        dat.r <- dat.cor$estimate  # Correlation coefficient
        dat.p <- dat.cor$p.value  # p-value
      }
      if(length(r)==0){
        
        r <- dat.r
        p <- dat.p
        
        
      }else {
        r <- c(r,dat.r)
        p <- c(p,dat.p)
      }
      
    }
    
    if(length(data.r)==0){
      data.r <- r
      data.p <- p
      
    }else {
      data.r <- cbind(data.r,r)
      data.p <- cbind(data.p,p)
      
    }
    
  }
  colnames(data.r) <- cancers 
  rownames(data.r) <- hets
  colnames(data.p) <- cancers 
  rownames(data.p) <- hets
  write.csv(data.r,file= paste("Generation/data.r.csv",sep=''),quote=F)
  write.csv(data.p,file= paste("Generation/data.p.csv",sep=''),quote=F)
  data.p[is.na(data.p)] = 2
  #######Bubble correlation plot conversion
  library(cowplot) # For plotting
  library(ggplot2) # For plotting
  library(RColorBrewer) # For plotting
  library(dplyr) # For data processing
  library(tidyr) # For data processing
  
  genes <-colnames(data.r)
  tumor_types <-rownames((data.r))
  gene_tumor_combinations <- expand.grid(CancerALL = genes, MirnaALL = tumor_types)
  library(tidyr)
  # Assume your data frame is named df
  df=as.data.frame(data.r)
  df=t(df)
  df=as.data.frame(df)
  dfP=as.data.frame(data.p)
  dfP=t(dfP)
  dfP=as.data.frame(dfP)
  # Use the gather function to combine all row data into one column
  gathered_dataR <- gather(df, key = "mRNA", value = "cor")
  gathered_dataP <- gather(dfP, key = "mRNA", value = "P-Value")
  library(dplyr)
  # If gathered_dataR is a data frame, extract one column and then combine
  combined_data <- cbind(gene_tumor_combinations, Correlation = gathered_dataR$cor)
  combined_data <- cbind(combined_data, Pvalue = gathered_dataP$`P-Value`)
  merged_data <-combined_data
  colnames(merged_data)=c("Gene","Tumor_Type","Correlation","P_Value")
  # Data processing
  merged_data$fdr_group <- cut(merged_data$P_Value, breaks = c(0, 0.0001, 0.001, 0.01, 0.05,  1), labels = c("<0.0001", "<0.001", "<0.01", "<0.05", ">0.05"))
  merged_data$fdr_group <- factor(merged_data$fdr_group, levels = c(labels = c("<0.0001","<0.001", "<0.01", "<0.05", ">0.05")))
  merged_data$fdr_group_2 <- cut(merged_data$P_Value, breaks = c(0, 0.0001, 0.001, 0.01, 0.05, 1), labels = c("<0.0001","<0.001", "<0.01", "<0.05", ">0.05"))
  merged_data$fdr_group_2 <- factor(merged_data$fdr_group_2, levels = c(labels = c("<0.0001","<0.001", "<0.01", "<0.05", ">0.05")))
  merged_data$Gene <- factor(merged_data$Gene, levels = unique(merged_data$Gene))
  ggplot2::ggplot() +
    geom_hline(yintercept = tumor_types, color = "#E8E8E8") +
    geom_point(data = subset(merged_data, Correlation < 0), shape = 19, stroke = 0,
               aes(x = Gene, y = Tumor_Type, size = abs(Correlation), color = fdr_group)) +
    scale_x_discrete(drop = FALSE) +
    geom_point(data = subset(merged_data, Correlation > 0), shape = 21, stroke = 0.1,
               aes(x = Gene, y = Tumor_Type, size = abs(Correlation), fill = fdr_group_2)) +
    scale_x_discrete(drop = FALSE)+
    scale_size_continuous(limits = c(0, 1),breaks = c(0.1, 0.2, 0.4, 0.6, 0.8, 1))+
    scale_fill_manual(values = c(rev(brewer.pal(9,'Reds'))[1:4],"#E8E8E8"))+
    scale_color_manual(values = c(rev(brewer.pal(9,'Greens'))[1:4],"#E8E8E8"))+
    cowplot::theme_cowplot() +
    theme(panel.grid.major = element_blank(),
          axis.text.x = element_text(angle = 45, hjust = 1),
          panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
          plot.title = element_text(hjust = 0.5))+
    labs(color = "Negative\ncorrelation\ p-value", size = "R-value",
         fill = "Positive\ncorrelation\ p-value",x = "",y = "",
         title = paste0("Correlation of ",Isoform_single," with Heterogeneity")) +
    guides(fill = guide_legend(override.aes = list(size = 4),order = 3),
           color = guide_legend(override.aes = list(size = 4),order = 2),
           size = guide_legend(override.aes = list(size = c(2:7), fill = "white"),order = 1,))
  ggsave(paste("2-Heterogeneity-",Isoform_single,"-Pancancer-bubbleplot.pdf",sep=''), plot = last_plot(), device = "pdf", width = 12, height = 6)
  
}