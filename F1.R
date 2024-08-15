

                ### Clinical differential analysis


setwd("D:/R/Pancancer ISCA1/linhcuang")
fenqi <- c("Age")

cancers <- c("ACC","BLCA","BRCA","CESC","CHOL","COAD",
             "DLBC","ESCA","GBM","HNSC","KICH","KIRC",
             "KIRP","LAML","LGG","LIHC","LUAD","LUSC",
             "MESO","OV","PAAD","PCPG","PRAD","READ",
             "SARC","SKCM","STAD","TGCT","THCA","THYM",
             "UCEC","UCS","UVM")

dat <- c()

fenqilist <- c("Age")

for(fenqi in fenqilist){
  
  cancers <- c("ACC","BLCA","BRCA","CESC","CHOL","COAD",
               "DLBC","ESCA","GBM","HNSC","KICH","KIRC",
               "KIRP","LAML","LGG","LIHC","LUAD","LUSC",
               "MESO","OV","PAAD","PCPG","PRAD","READ",
               "SARC","SKCM","STAD","TGCT","THCA","THYM",
               "UCEC","UCS","UVM")
  
  dat <- c()
  
  for (cancer in cancers) {
    
    cluster <- read.csv(paste("D:/R/Pancancer ISCA1/linhcuang/gene/fenqi/",fenqi,".csv",sep=''),header = F)
    
    if(cancer=="GBM"){
      
      n1 <- grep(cancer,cluster$V1)
      cluster <- cluster[n1,]
      n2 <- grep("GBMLGG",cluster$V1)
      cluster <- cluster[-n2,]
      
    }else 
      
      if(cancer=="LGG"){
        
        n1 <- grep(cancer,cluster$V1)
        cluster <- cluster[n1,]
        n2 <- grep("GBMLGG",cluster$V1)
        cluster <- cluster[-n2,]
        
      }else 
        
        if(cancer=="COAD"){
          
          n1 <- grep(cancer,cluster$V1)
          cluster <- cluster[n1,]
          n2 <- grep("COADREAD",cluster$V1)
          cluster <- cluster[-n2,]
          
        }else 
          
          if(cancer=="READ"){
            
            n1 <- grep(cancer,cluster$V1)
            cluster <- cluster[n1,]
            n2 <- grep("COADREAD",cluster$V1)
            cluster <- cluster[-n2,]
            
            
          }else{
            n1 <- grep(cancer,cluster$V1)
            cluster <- cluster[n1,]
          }
    
    cluster <- cluster[,-1]
    rownames(cluster) <- cluster$V4
    
    
    
    
    library(data.table)
    data <-  fread(paste("D:/R/Pancancer ISCA1/linhcuang/TCGA_new/",cancer,".txt",sep=''))
    data <- as.data.frame(data)
    rownames(data) <- data[,1]
    data <- data[,-1]
    
    data <- t(data)
    
    
    all_name <- names(which(table(c(rownames(data),rownames(cluster) ))==2))
    
    cluster1 <- cluster[match(all_name,rownames(cluster)),]
    
    data1 <- data[match(all_name,rownames(data)),]
    
    aa1 <- rep(cancer, times=length(rownames(data1)))
    
    exp <- unlist(data1[,match("ISCA1",colnames(data))])
    
    dat1 <- data.frame(aa1,cluster1[,1],exp)
    
    if(length(dat)==0){
      dat <- dat1
      
    }else {
      dat <- rbind(dat,dat1)
    }
  }
  
  colnames(dat) <- c("Cancer","Group","value")
  xxxx <- matrix(unique(dat[,1])) #colname
  temp_length <- length(xxxx)*2
  library(ggpubr)
  df <- compare_means(value ~ Group, data = dat, group.by = "Cancer",method = "anova")
  
  pdf(paste("fig2_",fenqi,"_ns.pdf",sep=''),width=temp_length/2.5,height = 8)
  p <- ggboxplot(dat, x = "Cancer", y = "value",
                 color = "Group", palette = "FA7F6F", 
                 add = "jitter") 
  p <- p+xlab("Cancer")+ylab("log(TPM+1)")
  p <- p + stat_compare_means(aes(group = Group), label = "p.signif",method = "anova")
  p <- p + theme(axis.text.x = element_text(size = 20, angle = 45, hjust = 1)) 
  p <- p + theme(axis.text.y = element_text(size = 20)) 
  p <- p + theme(axis.text.y = element_text(size = 20,  hjust = 1))  
  p <- p + labs(x = "Cancer Type", y = "Gene Expression (TPM)") +
    theme(axis.title.x = element_text(size = 20),  
          axis.title.y = element_text(size = 20))  
  p
  print(p)
  dev.off()
  
  
  
  xx <- compare_means(value ~ Group, data = dat, group.by = "Cancer",method = "anova")
  p_value <- as.matrix(xx$p)
  
  p_value[is.na(p_value)] <- 0.5
  
  for(j in 1:length(p_value))
  {
    if (p_value[j] > 0.05)
    {
      dat <- dat[dat[,"Cancer"]!=xxxx[j],]
      temp_length <- temp_length - 2
    }
    
  }  
  
  pdf(paste("fig2_",fenqi,".pdf",sep=''),width=temp_length/1.5,height = 8)
  p <- ggboxplot(dat, x = "Cancer", y = "value",
                 color = "Group", palette = "FA7F6F", 
                 add = "jitter") 
  p <- p+xlab("Cancer")+ylab("log(TPM+1)")
  p <- p + stat_compare_means(aes(group = Group), label = "p.signif",method = "anova")
  p <- p + theme(axis.text.x = element_text(size = 20, angle = 45, hjust = 1)) 
  p <- p + theme(axis.text.y = element_text(size = 20))  
  p <- p + theme(axis.text.y = element_text(size = 20,  hjust = 1)) 
  p <- p + labs(x = "Cancer Type", y = "Gene Expression (TPM)") +
    theme(axis.title.x = element_text(size = 20),  
          axis.title.y = element_text(size = 20))  
  p
  print(p)
  dev.off()
  
  
}










fenqilist <- c("Gender")

for(fenqi in fenqilist){
  
  cancers <- c("ACC","BLCA","BRCA","CESC","CHOL","COAD",
               "DLBC","ESCA","GBM","HNSC","KICH","KIRC",
               "KIRP","LAML","LGG","LIHC","LUAD","LUSC",
               "MESO","OV","PAAD","PCPG","PRAD","READ",
               "SARC","SKCM","STAD","TGCT","THCA","THYM",
               "UCEC","UCS","UVM")
  
  dat <- c()
  
  for (cancer in cancers) {
    
    cluster <- read.csv(paste("D:/R/Pancancer ISCA1/linhcuang/gene/fenqi/",fenqi,".csv",sep=''),header = F)
    
    
    #cancer <- "GBM"
    if(cancer=="GBM"){
      
      n1 <- grep(cancer,cluster$V1)
      cluster <- cluster[n1,]
      n2 <- grep("GBMLGG",cluster$V1)
      cluster <- cluster[-n2,]
      
    }else 
      #cancer <- "LGG"
      if(cancer=="LGG"){
        
        n1 <- grep(cancer,cluster$V1)
        cluster <- cluster[n1,]
        n2 <- grep("GBMLGG",cluster$V1)
        cluster <- cluster[-n2,]
        
      }else 
        #cancer <- "COAD"
        if(cancer=="COAD"){
          
          n1 <- grep(cancer,cluster$V1)
          cluster <- cluster[n1,]
          n2 <- grep("COADREAD",cluster$V1)
          cluster <- cluster[-n2,]
          
        }else 
          #cancer <- "READ"
          if(cancer=="READ"){
            
            n1 <- grep(cancer,cluster$V1)
            cluster <- cluster[n1,]
            n2 <- grep("COADREAD",cluster$V1)
            cluster <- cluster[-n2,]
            
            
          }else{
            n1 <- grep(cancer,cluster$V1)
            cluster <- cluster[n1,]
          }
    
    cluster <- cluster[,-1]
    rownames(cluster) <- cluster$V4
    
    
    
    
    library(data.table)
    data <-  fread(paste("D:/R/Pancancer ISCA1/linchuang/TCGA_new/",cancer,".txt",sep=''))
    data <- as.data.frame(data)
    rownames(data) <- data[,1]
    data <- data[,-1]
    
    data <- t(data)
    
    
    all_name <- names(which(table(c(rownames(data),rownames(cluster) ))==2))
    
    
    cluster1 <- cluster[match(all_name,rownames(cluster)),]
    
    data1 <- data[match(all_name,rownames(data)),]
    
    aa1 <- rep(cancer, times=length(rownames(data1)))
    
    exp <- unlist(data1[,match("ISCA1",colnames(data))])
    
    dat1 <- data.frame(aa1,cluster1[,1],exp)
    
    if(length(dat)==0){
      dat <- dat1
      
    }else {
      dat <- rbind(dat,dat1)
    }
  }
  
  colnames(dat) <- c("Cancer","Group","value")
  
  xxxx <- matrix(unique(dat[,1])) #colname
  
  temp_length <- length(xxxx)*2
  
  library(ggpubr)
  df <- compare_means(value ~ Group, data = dat, group.by = "Cancer",method = "anova")
  
  pdf(paste("fig2_",fenqi,"_ns.pdf",sep=''),width=temp_length/2.5,height = 8)
  p <- ggboxplot(dat, x = "Cancer", y = "value",
                 color = "Group", palette = "FA7F6F", 
                 add = "jitter")   
  p <- p+xlab("Cancer")+ylab("log(TPM+1)")
  p <- p + stat_compare_means(aes(group = Group), label = "p.signif",method = "anova")
  p <- p + theme(axis.text.x = element_text(size = 20, angle = 45, hjust = 1))    
  p <- p + theme(axis.text.y = element_text(size = 20))   
  p <- p + theme(axis.text.y = element_text(size = 20,  hjust = 1))    
  p <- p + labs(x = "Cancer Type", y = "Gene Expression (TPM)") +
    theme(axis.title.x = element_text(size = 20),    
          axis.title.y = element_text(size = 20))    
  p
  print(p)
  dev.off()
  
  
  
  xx <- compare_means(value ~ Group, data = dat, group.by = "Cancer",method = "anova")
  p_value <- as.matrix(xx$p)
  
  p_value[is.na(p_value)] <- 0.5
  
  for(j in 1:length(p_value))
  {
    if (p_value[j] > 0.05)
    {
      dat <- dat[dat[,"Cancer"]!=xxxx[j],]
      temp_length <- temp_length - 2
    }
    
  }  
  
  
  pdf(paste("fig2_",fenqi,".pdf",sep=''),width=temp_length/0.6,height = 8)
  p <- ggboxplot(dat, x = "Cancer", y = "value",
                 color = "Group", palette = "FA7F6F", 
                 add = "jitter")   
  p <- p+xlab("Cancer")+ylab("log(TPM+1)")
  p <- p + stat_compare_means(aes(group = Group), label = "p.signif",method = "anova")
  p <- p + theme(axis.text.x = element_text(size = 20, angle = 45, hjust = 1))    
  p <- p + theme(axis.text.y = element_text(size = 20))   
  p <- p + theme(axis.text.y = element_text(size = 20,  hjust = 1))    
  p <- p + labs(x = "Cancer Type", y = "Gene Expression (TPM)") +
    theme(axis.title.x = element_text(size = 20),    
          axis.title.y = element_text(size = 20))    
  p
  print(p)
  dev.off()
  
  
}































fenqilist <- c("M")

for(fenqi in fenqilist){
  
  cancers <- c("ACC","BLCA","BRCA","CESC","CHOL","COAD","DLBC","ESCA","GBM","HNSC","KICH","KIRC","KIRP","LAML","LGG","LIHC","LUAD","LUSC","MESO","OV","PAAD","PCPG",
               "PRAD","READ","SARC","SKCM","STAD","TGCT","THCA","THYM","UCEC","UCS","UVM")
  
  dat <- c()
  
  for (cancer in cancers) {
    
    cluster <- read.csv(paste("D:/R/Pancancer ISCA1/linchuang/gene/fenqi/",fenqi,".csv",sep=''),header = F)
    
    
    #cancer <- "GBM"
    if(cancer=="GBM"){
      
      n1 <- grep(cancer,cluster$V1)
      cluster <- cluster[n1,]
      n2 <- grep("GBMLGG",cluster$V1)
      cluster <- cluster[-n2,]
      
    }else 
      #cancer <- "LGG"
      if(cancer=="LGG"){
        
        n1 <- grep(cancer,cluster$V1)
        cluster <- cluster[n1,]
        n2 <- grep("GBMLGG",cluster$V1)
        cluster <- cluster[-n2,]
        
      }else 
        #cancer <- "COAD"
        if(cancer=="COAD"){
          
          n1 <- grep(cancer,cluster$V1)
          cluster <- cluster[n1,]
          n2 <- grep("COADREAD",cluster$V1)
          cluster <- cluster[-n2,]
          
        }else 
          #cancer <- "READ"
          if(cancer=="READ"){
            
            n1 <- grep(cancer,cluster$V1)
            cluster <- cluster[n1,]
            n2 <- grep("COADREAD",cluster$V1)
            cluster <- cluster[-n2,]
            
            
          }else{
            n1 <- grep(cancer,cluster$V1)
            cluster <- cluster[n1,]
          }
    
    cluster <- cluster[,-1]
    rownames(cluster) <- cluster$V4
    
    
    
    
    library(data.table)
    data <-  fread(paste("D:/R/Pancancer ISCA1/linchuang/TCGA_NEW/",cancer,".txt",sep=''))
    data <- as.data.frame(data)
    rownames(data) <- data[,1]
    data <- data[,-1]
    
    data <- t(data)
    
    
    all_name <- names(which(table(c(rownames(data),rownames(cluster) ))==2))
    
    
    cluster1 <- cluster[match(all_name,rownames(cluster)),]
    
    data1 <- data[match(all_name,rownames(data)),]
    
    aa1 <- rep(cancer, times=length(rownames(data1)))
    
    exp <- unlist(data1[,match("ISCA1",colnames(data))])
    
    dat1 <- data.frame(aa1,cluster1[,1],exp)
    
    if(length(dat)==0){
      dat <- dat1
      
    }else {
      dat <- rbind(dat,dat1)
    }
  }
  
  colnames(dat) <- c("Cancer","Group","value")
  
  xxxx <- matrix(unique(dat[,1])) #colname
  
  temp_length <- length(xxxx)*2
  
  library(ggpubr)
  df <- compare_means(value ~ Group, data = dat, group.by = "Cancer",method = "anova")
  
  pdf(paste("fig2_",fenqi,"_ns.pdf",sep=''),width=temp_length/2.5,height = 8)
  p <- ggboxplot(dat, x = "Cancer", y = "value",
                 color = "Group", palette = "FA7F6F", 
                 add = "jitter")   
  p <- p+xlab("Cancer")+ylab("log(TPM+1)")
  p <- p + stat_compare_means(aes(group = Group), label = "p.signif",method = "anova")
  p <- p + theme(axis.text.x = element_text(size = 20, angle = 45, hjust = 1))    
  p <- p + theme(axis.text.y = element_text(size = 20))   
  p <- p + theme(axis.text.y = element_text(size = 20,  hjust = 1))    
  p <- p + labs(x = "Cancer Type", y = "Gene Expression (TPM)") +
    theme(axis.title.x = element_text(size = 20),    
          axis.title.y = element_text(size = 20))    
  p
  print(p)
  dev.off()
  
  
  
  xx <- compare_means(value ~ Group, data = dat, group.by = "Cancer",method = "anova")
  p_value <- as.matrix(xx$p)
  
  p_value[is.na(p_value)] <- 0.5
  
  for(j in 1:length(p_value))
  {
    if (p_value[j] > 0.05)
    {
      dat <- dat[dat[,"Cancer"]!=xxxx[j],]
      temp_length <- temp_length - 2
    }
    
  }  
  
  
  pdf(paste("fig2_",fenqi,".pdf",sep=''),width=temp_length/1.5,height = 8)
  p <- ggboxplot(dat, x = "Cancer", y = "value",
                 color = "Group", palette = "FA7F6F", 
                 add = "jitter")   
  p <- p+xlab("Cancer")+ylab("log(TPM+1)")
  p <- p + stat_compare_means(aes(group = Group), label = "p.signif",method = "anova")
  p <- p + theme(axis.text.x = element_text(size = 20, angle = 45, hjust = 1))    
  p <- p + theme(axis.text.y = element_text(size = 20))   
  p <- p + theme(axis.text.y = element_text(size = 20,  hjust = 1))    
  p <- p + labs(x = "Cancer Type", y = "Gene Expression (TPM)") +
    theme(axis.title.x = element_text(size = 20),    
          axis.title.y = element_text(size = 20))    
  p
  print(p)
  dev.off()
  
  
}











fenqilist <- c("N")

for(fenqi in fenqilist){
  
  cancers <- c("ACC","BLCA","BRCA","CESC","CHOL","COAD","DLBC","ESCA","GBM","HNSC","KICH","KIRC","KIRP","LAML","LGG","LIHC","LUAD","LUSC","MESO","OV","PAAD","PCPG",
               "PRAD","READ","SARC","SKCM","STAD","TGCT","THCA","THYM","UCEC","UCS","UVM")
  
  dat <- c()
  
  for (cancer in cancers) {
    
    cluster <- read.csv(paste("D:/R/Pancancer ISCA1/linchuang/gene/fenqi/",fenqi,".csv",sep=''),header = F)
    
    
    #cancer <- "GBM"
    if(cancer=="GBM"){
      
      n1 <- grep(cancer,cluster$V1)
      cluster <- cluster[n1,]
      n2 <- grep("GBMLGG",cluster$V1)
      cluster <- cluster[-n2,]
      
    }else 
      #cancer <- "LGG"
      if(cancer=="LGG"){
        
        n1 <- grep(cancer,cluster$V1)
        cluster <- cluster[n1,]
        n2 <- grep("GBMLGG",cluster$V1)
        cluster <- cluster[-n2,]
        
      }else 
        #cancer <- "COAD"
        if(cancer=="COAD"){
          
          n1 <- grep(cancer,cluster$V1)
          cluster <- cluster[n1,]
          n2 <- grep("COADREAD",cluster$V1)
          cluster <- cluster[-n2,]
          
        }else 
          #cancer <- "READ"
          if(cancer=="READ"){
            
            n1 <- grep(cancer,cluster$V1)
            cluster <- cluster[n1,]
            n2 <- grep("COADREAD",cluster$V1)
            cluster <- cluster[-n2,]
            
            
          }else{
            n1 <- grep(cancer,cluster$V1)
            cluster <- cluster[n1,]
          }
    
    cluster <- cluster[,-1]
    rownames(cluster) <- cluster$V4
    
    
    
    
    library(data.table)
    data <-  fread(paste("D:/R/Pancancer ISCA1/linchuang/TCGA_NEW/",cancer,".txt",sep=''))
    data <- as.data.frame(data)
    rownames(data) <- data[,1]
    data <- data[,-1]
    
    data <- t(data)
    
    
    all_name <- names(which(table(c(rownames(data),rownames(cluster) ))==2))
    
    
    cluster1 <- cluster[match(all_name,rownames(cluster)),]
    
    data1 <- data[match(all_name,rownames(data)),]
    
    aa1 <- rep(cancer, times=length(rownames(data1)))
    
    exp <- unlist(data1[,match("ISCA1",colnames(data))])
    
    dat1 <- data.frame(aa1,cluster1[,1],exp)
    
    if(length(dat)==0){
      dat <- dat1
      
    }else {
      dat <- rbind(dat,dat1)
    }
  }
  
  colnames(dat) <- c("Cancer","Group","value")
  
  xxxx <- matrix(unique(dat[,1])) #colname
  
  temp_length <- length(xxxx)*2
  
  library(ggpubr)
  df <- compare_means(value ~ Group, data = dat, group.by = "Cancer",method = "anova")
  
  pdf(paste("fig2_",fenqi,"_ns.pdf",sep=''),width=temp_length/2.5,height = 8)
  p <- ggboxplot(dat, x = "Cancer", y = "value",
                 color = "Group", palette = "FA7F6F", 
                 add = "jitter")   
  p <- p+xlab("Cancer")+ylab("log(TPM+1)")
  p <- p + stat_compare_means(aes(group = Group), label = "p.signif",method = "anova")
  p <- p + theme(axis.text.x = element_text(size = 20, angle = 45, hjust = 1))    
  p <- p + theme(axis.text.y = element_text(size = 20))   
  p <- p + theme(axis.text.y = element_text(size = 20,  hjust = 1))    
  p <- p + labs(x = "Cancer Type", y = "Gene Expression (TPM)") +
    theme(axis.title.x = element_text(size = 20),    
          axis.title.y = element_text(size = 20))    
  p
  print(p)
  dev.off()
  
  
  
  xx <- compare_means(value ~ Group, data = dat, group.by = "Cancer",method = "anova")
  p_value <- as.matrix(xx$p)
  
  p_value[is.na(p_value)] <- 0.5
  
  for(j in 1:length(p_value))
  {
    if (p_value[j] > 0.05)
    {
      dat <- dat[dat[,"Cancer"]!=xxxx[j],]
      temp_length <- temp_length - 2
    }
    
  }  
  
  
  pdf(paste("fig2_",fenqi,".pdf",sep=''),width=temp_length/1.5,height = 8)
  p <- ggboxplot(dat, x = "Cancer", y = "value",
                 color = "Group", palette = "FA7F6F", 
                 add = "jitter")   
  p <- p+xlab("Cancer")+ylab("log(TPM+1)")
  p <- p + stat_compare_means(aes(group = Group), label = "p.signif",method = "anova")
  p <- p + theme(axis.text.x = element_text(size = 20, angle = 45, hjust = 1))    
  p <- p + theme(axis.text.y = element_text(size = 20))   
  p <- p + theme(axis.text.y = element_text(size = 20,  hjust = 1))    
  p <- p + labs(x = "Cancer Type", y = "Gene Expression (TPM)") +
    theme(axis.title.x = element_text(size = 20),    
          axis.title.y = element_text(size = 20))    
  p
  print(p)
  dev.off()
  
  
}















fenqilist <- c("Grade")

for(fenqi in fenqilist){
  
  cancers <- c("ACC","BLCA","BRCA","CESC","CHOL","COAD","DLBC","ESCA","GBM","HNSC","KICH","KIRC","KIRP","LAML","LGG","LIHC","LUAD","LUSC","MESO","OV","PAAD","PCPG",
               "PRAD","READ","SARC","SKCM","STAD","TGCT","THCA","THYM","UCEC","UCS","UVM")
  
  dat <- c()
  
  for (cancer in cancers) {
    
    cluster <- read.csv(paste("D:/R/Pancancer ISCA1/linchuang/gene/fenqi/",fenqi,".csv",sep=''),header = F)
    
    
    #cancer <- "GBM"
    if(cancer=="GBM"){
      
      n1 <- grep(cancer,cluster$V1)
      cluster <- cluster[n1,]
      n2 <- grep("GBMLGG",cluster$V1)
      cluster <- cluster[-n2,]
      
    }else 
      #cancer <- "LGG"
      if(cancer=="LGG"){
        
        n1 <- grep(cancer,cluster$V1)
        cluster <- cluster[n1,]
        n2 <- grep("GBMLGG",cluster$V1)
        cluster <- cluster[-n2,]
        
      }else 
        #cancer <- "COAD"
        if(cancer=="COAD"){
          
          n1 <- grep(cancer,cluster$V1)
          cluster <- cluster[n1,]
          n2 <- grep("COADREAD",cluster$V1)
          cluster <- cluster[-n2,]
          
        }else 
          #cancer <- "READ"
          if(cancer=="READ"){
            
            n1 <- grep(cancer,cluster$V1)
            cluster <- cluster[n1,]
            n2 <- grep("COADREAD",cluster$V1)
            cluster <- cluster[-n2,]
            
            
          }else{
            n1 <- grep(cancer,cluster$V1)
            cluster <- cluster[n1,]
          }
    
    cluster <- cluster[,-1]
    rownames(cluster) <- cluster$V4
    
    
    
    
    library(data.table)
    data <-  fread(paste("D:/R/Pancancer ISCA1/linchuang/TCGA_NEW/",cancer,".txt",sep=''))
    data <- as.data.frame(data)
    rownames(data) <- data[,1]
    data <- data[,-1]
    
    data <- t(data)
    
    
    all_name <- names(which(table(c(rownames(data),rownames(cluster) ))==2))
    
    
    cluster1 <- cluster[match(all_name,rownames(cluster)),]
    
    data1 <- data[match(all_name,rownames(data)),]
    
    aa1 <- rep(cancer, times=length(rownames(data1)))
    
    exp <- unlist(data1[,match("ISCA1",colnames(data))])
    
    dat1 <- data.frame(aa1,cluster1[,1],exp)
    
    if(length(dat)==0){
      dat <- dat1
      
    }else {
      dat <- rbind(dat,dat1)
    }
  }
  
  colnames(dat) <- c("Cancer","Group","value")
  
  xxxx <- matrix(unique(dat[,1])) #colname
  
  temp_length <- length(xxxx)*2
  
  library(ggpubr)
  df <- compare_means(value ~ Group, data = dat, group.by = "Cancer",method = "anova")
  
  pdf(paste("fig2_",fenqi,"_ns.pdf",sep=''),width=temp_length/2.5,height = 8)
  p <- ggboxplot(dat, x = "Cancer", y = "value",
                 color = "Group", palette = "FA7F6F", 
                 add = "jitter")   
  p <- p+xlab("Cancer")+ylab("log(TPM+1)")
  p <- p + stat_compare_means(aes(group = Group), label = "p.signif",method = "anova")
  p <- p + theme(axis.text.x = element_text(size = 20, angle = 45, hjust = 1))    
  p <- p + theme(axis.text.y = element_text(size = 20))   
  p <- p + theme(axis.text.y = element_text(size = 20,  hjust = 1))    
  p <- p + labs(x = "Cancer Type", y = "Gene Expression (TPM)") +
    theme(axis.title.x = element_text(size = 20),    
          axis.title.y = element_text(size = 20))    
  p
  print(p)
  dev.off()
  
  
  
  xx <- compare_means(value ~ Group, data = dat, group.by = "Cancer",method = "anova")
  p_value <- as.matrix(xx$p)
  
  p_value[is.na(p_value)] <- 0.5
  
  for(j in 1:length(p_value))
  {
    if (p_value[j] > 0.05)
    {
      dat <- dat[dat[,"Cancer"]!=xxxx[j],]
      temp_length <- temp_length - 2
    }
    
  }  
  
  
  pdf(paste("fig2_",fenqi,".pdf",sep=''),width=temp_length/1.5,height = 8)
  p <- ggboxplot(dat, x = "Cancer", y = "value",
                 color = "Group", palette = "FA7F6F", 
                 add = "jitter")   
  p <- p+xlab("Cancer")+ylab("log(TPM+1)")
  p <- p + stat_compare_means(aes(group = Group), label = "p.signif",method = "anova")
  p <- p + theme(axis.text.x = element_text(size = 20, angle = 45, hjust = 1))    
  p <- p + theme(axis.text.y = element_text(size = 20))   
  p <- p + theme(axis.text.y = element_text(size = 20,  hjust = 1))    
  p <- p + labs(x = "Cancer Type", y = "Gene Expression (TPM)") +
    theme(axis.title.x = element_text(size = 20),    
          axis.title.y = element_text(size = 20))    
  p
  print(p)
  dev.off()
  
  
}









fenqilist <- c("Stage")
for(fenqi in fenqilist){
  
  cancers <- c("ACC","BLCA","BRCA","CESC","CHOL","COAD","DLBC","ESCA","GBM","HNSC","KICH","KIRC","KIRP","LAML","LGG","LIHC","LUAD","LUSC","MESO","OV","PAAD","PCPG",
               "PRAD","READ","SARC","SKCM","STAD","TGCT","THCA","THYM","UCEC","UCS","UVM")
  
  dat <- c()
  
  for (cancer in cancers) {
    
    cluster <- read.csv(paste("D:/R/Pancancer ISCA1/linchuang/gene/fenqi/",fenqi,".csv",sep=''),header = F)
    
    
    #cancer <- "GBM"
    if(cancer=="GBM"){
      
      n1 <- grep(cancer,cluster$V1)
      cluster <- cluster[n1,]
      n2 <- grep("GBMLGG",cluster$V1)
      cluster <- cluster[-n2,]
      
    }else 
      #cancer <- "LGG"
      if(cancer=="LGG"){
        
        n1 <- grep(cancer,cluster$V1)
        cluster <- cluster[n1,]
        n2 <- grep("GBMLGG",cluster$V1)
        cluster <- cluster[-n2,]
        
      }else 
        #cancer <- "COAD"
        if(cancer=="COAD"){
          
          n1 <- grep(cancer,cluster$V1)
          cluster <- cluster[n1,]
          n2 <- grep("COADREAD",cluster$V1)
          cluster <- cluster[-n2,]
          
        }else 
          #cancer <- "READ"
          if(cancer=="READ"){
            
            n1 <- grep(cancer,cluster$V1)
            cluster <- cluster[n1,]
            n2 <- grep("COADREAD",cluster$V1)
            cluster <- cluster[-n2,]
            
            
          }else{
            n1 <- grep(cancer,cluster$V1)
            cluster <- cluster[n1,]
          }
    
    cluster <- cluster[,-1]
    rownames(cluster) <- cluster$V4
    
    
    
    
    library(data.table)
    data <-  fread(paste("D:/R/Pancancer ISCA1/linchuang/TCGA_NEW/",cancer,".txt",sep=''))
    data <- as.data.frame(data)
    rownames(data) <- data[,1]
    data <- data[,-1]
    
    data <- t(data)
    
    
    all_name <- names(which(table(c(rownames(data),rownames(cluster) ))==2))
    
    
    cluster1 <- cluster[match(all_name,rownames(cluster)),]
    
    data1 <- data[match(all_name,rownames(data)),]
    
    aa1 <- rep(cancer, times=length(rownames(data1)))
    
    exp <- unlist(data1[,match("ISCA1",colnames(data))])
    
    dat1 <- data.frame(aa1,cluster1[,1],exp)
    
    if(length(dat)==0){
      dat <- dat1
      
    }else {
      dat <- rbind(dat,dat1)
    }
  }
  
  colnames(dat) <- c("Cancer","Group","value")
  
  xxxx <- matrix(unique(dat[,1])) #colname
  
  temp_length <- length(xxxx)*2
  
  library(ggpubr)
  df <- compare_means(value ~ Group, data = dat, group.by = "Cancer",method = "anova")
  
  pdf(paste("fig2_",fenqi,"_ns.pdf",sep=''),width=temp_length/2.5,height = 8)
  p <- ggboxplot(dat, x = "Cancer", y = "value",
                 color = "Group", palette = "FA7F6F", 
                 add = "jitter")   
  p <- p+xlab("Cancer")+ylab("log(TPM+1)")
  p <- p + stat_compare_means(aes(group = Group), label = "p.signif",method = "anova")
  p <- p + theme(axis.text.x = element_text(size = 20, angle = 45, hjust = 1))    
  p <- p + theme(axis.text.y = element_text(size = 20))   
  p <- p + theme(axis.text.y = element_text(size = 20,  hjust = 1))    
  p <- p + labs(x = "Cancer Type", y = "Gene Expression (TPM)") +
    theme(axis.title.x = element_text(size = 20),    
          axis.title.y = element_text(size = 20))    
  p
  print(p)
  dev.off()
  
  
  
  xx <- compare_means(value ~ Group, data = dat, group.by = "Cancer",method = "anova")
  p_value <- as.matrix(xx$p)
  
  p_value[is.na(p_value)] <- 0.5
  
  for(j in 1:length(p_value))
  {
    if (p_value[j] > 0.05)
    {
      dat <- dat[dat[,"Cancer"]!=xxxx[j],]
      temp_length <- temp_length - 2
    }
    
  }  
  
  
  pdf(paste("fig2_",fenqi,".pdf",sep=''),width=temp_length/1,height = 8)
  p <- ggboxplot(dat, x = "Cancer", y = "value",
                 color = "Group", palette = "FA7F6F", 
                 add = "jitter")   
  p <- p+xlab("Cancer")+ylab("log(TPM+1)")
  p <- p + stat_compare_means(aes(group = Group), label = "p.signif",method = "anova")
  p <- p + theme(axis.text.x = element_text(size = 20, angle = 45, hjust = 1))    
  p <- p + theme(axis.text.y = element_text(size = 20))   
  p <- p + theme(axis.text.y = element_text(size = 20,  hjust = 1))    
  p <- p + labs(x = "Cancer Type", y = "Gene Expression (TPM)") +
    theme(axis.title.x = element_text(size = 20),    
          axis.title.y = element_text(size = 20))    
  p
  print(p)
  dev.off()
  
  
}












fenqilist <- c("T")
for(fenqi in fenqilist){
  
  cancers <- c("ACC","BLCA","BRCA","CESC","CHOL","COAD","DLBC","ESCA","GBM","HNSC","KICH","KIRC","KIRP","LAML","LGG","LIHC","LUAD","LUSC","MESO","OV","PAAD","PCPG",
               "PRAD","READ","SARC","SKCM","STAD","TGCT","THCA","THYM","UCEC","UCS","UVM")
  
  dat <- c()
  
  for (cancer in cancers) {
    
    cluster <- read.csv(paste("D:/R/Pancancer ISCA1/linchuang/gene/fenqi/",fenqi,".csv",sep=''),header = F)
    
    
    #cancer <- "GBM"
    if(cancer=="GBM"){
      
      n1 <- grep(cancer,cluster$V1)
      cluster <- cluster[n1,]
      n2 <- grep("GBMLGG",cluster$V1)
      cluster <- cluster[-n2,]
      
    }else 
      #cancer <- "LGG"
      if(cancer=="LGG"){
        
        n1 <- grep(cancer,cluster$V1)
        cluster <- cluster[n1,]
        n2 <- grep("GBMLGG",cluster$V1)
        cluster <- cluster[-n2,]
        
      }else 
        #cancer <- "COAD"
        if(cancer=="COAD"){
          
          n1 <- grep(cancer,cluster$V1)
          cluster <- cluster[n1,]
          n2 <- grep("COADREAD",cluster$V1)
          cluster <- cluster[-n2,]
          
        }else 
          #cancer <- "READ"
          if(cancer=="READ"){
            
            n1 <- grep(cancer,cluster$V1)
            cluster <- cluster[n1,]
            n2 <- grep("COADREAD",cluster$V1)
            cluster <- cluster[-n2,]
            
            
          }else{
            n1 <- grep(cancer,cluster$V1)
            cluster <- cluster[n1,]
          }
    
    cluster <- cluster[,-1]
    rownames(cluster) <- cluster$V4
    
    
    
    
    library(data.table)
    data <-  fread(paste("D:/R/Pancancer ISCA1/linchuang/TCGA_NEW/",cancer,".txt",sep=''))
    data <- as.data.frame(data)
    rownames(data) <- data[,1]
    data <- data[,-1]
    
    data <- t(data)
    
    
    all_name <- names(which(table(c(rownames(data),rownames(cluster) ))==2))
    
    
    cluster1 <- cluster[match(all_name,rownames(cluster)),]
    
    data1 <- data[match(all_name,rownames(data)),]
    
    aa1 <- rep(cancer, times=length(rownames(data1)))
    
    exp <- unlist(data1[,match("ISCA1",colnames(data))])
    
    dat1 <- data.frame(aa1,cluster1[,1],exp)
    
    if(length(dat)==0){
      dat <- dat1
      
    }else {
      dat <- rbind(dat,dat1)
    }
  }
  
  colnames(dat) <- c("Cancer","Group","value")
  
  xxxx <- matrix(unique(dat[,1])) #colname
  
  temp_length <- length(xxxx)*2
  
  library(ggpubr)
  df <- compare_means(value ~ Group, data = dat, group.by = "Cancer",method = "anova")
  
  pdf(paste("fig2_",fenqi,"_ns.pdf",sep=''),width=temp_length/2.5,height = 8)
  p <- ggboxplot(dat, x = "Cancer", y = "value",
                 color = "Group", palette = "FA7F6F", 
                 add = "jitter")   
  p <- p+xlab("Cancer")+ylab("log(TPM+1)")
  p <- p + stat_compare_means(aes(group = Group), label = "p.signif",method = "anova")
  p <- p + theme(axis.text.x = element_text(size = 20, angle = 45, hjust = 1))    
  p <- p + theme(axis.text.y = element_text(size = 20))   
  p <- p + theme(axis.text.y = element_text(size = 20,  hjust = 1))    
  p <- p + labs(x = "Cancer Type", y = "Gene Expression (TPM)") +
    theme(axis.title.x = element_text(size = 20),    
          axis.title.y = element_text(size = 20))    
  p
  print(p)
  dev.off()
  
  
  
  xx <- compare_means(value ~ Group, data = dat, group.by = "Cancer",method = "anova")
  p_value <- as.matrix(xx$p)
  
  p_value[is.na(p_value)] <- 0.5
  
  for(j in 1:length(p_value))
  {
    if (p_value[j] > 0.05)
    {
      dat <- dat[dat[,"Cancer"]!=xxxx[j],]
      temp_length <- temp_length - 2
    }
    
  }  
  
  
  pdf(paste("fig2_",fenqi,".pdf",sep=''),width=temp_length/1.5,height = 8)
  p <- ggboxplot(dat, x = "Cancer", y = "value",
                 color = "Group", palette = "FA7F6F", 
                 add = "jitter")   
  p <- p+xlab("Cancer")+ylab("log(TPM+1)")
  p <- p + stat_compare_means(aes(group = Group), label = "p.signif",method = "anova")
  p <- p + theme(axis.text.x = element_text(size = 20, angle = 45, hjust = 1))    
  p <- p + theme(axis.text.y = element_text(size = 20))   
  p <- p + theme(axis.text.y = element_text(size = 20,  hjust = 1))    
  p <- p + labs(x = "Cancer Type", y = "Gene Expression (TPM)") +
    theme(axis.title.x = element_text(size = 20),    
          axis.title.y = element_text(size = 20))    
  p
  print(p)
  dev.off()
  
}




                            ### Pan-cancer differential analysis

setwd("D:/R/Pancancer ISCA1/fig1")

dat <- c()
cancers <- c("ACC","BLCA","BRCA","CESC","CHOL","COAD",
             "DLBC","ESCA","GBM","HNSC","KICH","KIRC",
             "KIRP","LAML","LGG","LIHC","LUAD","LUSC",
             "MESO","OV","PAAD","PCPG","PRAD","READ",
             "SARC","SKCM","STAD","TGCT","THCA","THYM",
             "UCEC","UCS","UVM")

for (cancer in cancers) {
  #cancer="ACC"
  library(data.table)
  data1 <- fread("D:/R/Pancancer ISCA1/fig1/EPCAM.csv")
  data1 <- as.data.frame(data1)
  cluster <- data1
  cluster[,1] <- gsub("\\([^\\)]+\\)", "", cluster[,1])
  #i=1
  nn <- which(cluster$V1==cancer)
  dd <- cluster[nn,]
  rownames(dd) <- dd$SampleName
  exp <- dd$Expression
  
  dat1 <- data.frame(dd[,-4])
  if(length(dat)==0){
    dat <- dat1
    
  }else {
    dat <- rbind(dat,dat1)
  }
  
}

colnames(dat) <- c("Cancer","Group","value")
unique(cluster$V1)
dat <- dat[dat$Cancer != "GBMLGG", ] 
dat <- dat[dat$Cancer != "COADREAD", ] 
dat <- dat[dat$Cancer != "ALL", ] 
dat <- dat[dat$Cancer != "KIPAN", ] 
xxxx <- matrix(unique(dat[,1])) #colname
temp_length <- length(xxxx)*2

library(ggpubr)
df <- compare_means(value ~ Group, data = dat, group.by = "Cancer")

pdf("fig1_ns.pdf",width=temp_length/2.5,height = 8)
p <- ggboxplot(dat, x = "Cancer", y = "value",
               color = "Group", palette = "FA7F6F", 
               add = "jitter") 
p <- p+xlab("Cancer")+ylab("log(TPM+0.001)")
p <- p + stat_compare_means(aes(group = Group), label = "p.signif")
p <- p + theme(axis.text.x = element_text(size = 20, angle = 45, hjust = 1))  
p <- p + theme(axis.text.y = element_text(size = 20,  hjust = 1)) 
p <- p + labs(x = "Cancer Type", y = "Gene Expression (log(TPM+0.001))") +
  theme(axis.title.x = element_text(size = 20),   
        axis.title.y = element_text(size = 20))  
p
#print(p)
dev.off()


xx <- compare_means(value ~ Group, data = dat, group.by = "Cancer")
p_value <- as.matrix(xx$p)

p_value[is.na(p_value)] <- 0.5

for(j in 1:length(p_value))
{
  if (p_value[j] > 0.05)
  {
    dat <- dat[dat[,"Cancer"]!=xxxx[j],]
    temp_length <- temp_length - 2
  }
}
pdf("fig1_ISCA1.pdf",width=temp_length/2.5,height = 8)
p <- ggboxplot(dat, x = "Cancer", y = "value",
               color = "Group", palette = "FA7F6F", 
               add = "jitter") 
p <- p + labs(x = "Cancer Type", y = "Gene Expression (log(TPM+1))") 
p <- p + stat_compare_means(aes(group = Group), label = "p.signif")
p <- p + theme(axis.text.x = element_text(size = 20, angle = 45, hjust = 1)) 
p <- p + theme(axis.text.y = element_text(size = 20)) 
p <- p + theme(axis.text.y = element_text(size = 20,  hjust = 1)) 
p <- p + labs(x = "Cancer Type", y = "Gene Expression (TPM)") +
  theme(axis.title.x = element_text(size = 20),   
        axis.title.y = element_text(size = 20))   
p#print(p)
dev.off()











