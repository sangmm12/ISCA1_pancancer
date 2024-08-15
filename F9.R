
### drug-predict

setwd('D:/R/Pancancer ISCA1/drug-new/')

# 2. Preparation of drug data
# 
# 2.1 Reading drug-related data

library(readxl)
rt1 <- read_excel("DTP_NCI60_ZSCORE.xlsx", skip = 7)
#Firstly, the read_excel() function in the readxl package is used to read the drug-related data. As the first 7 rows are annotation information, the parameter skip is used to skip the first 7 rows.

colnames(rt1) <- rt1[1,]
rt1 <- rt1[-1,-c(67,68)]
#Meanwhile, the first row is set as the column name, and the last two columns and other information are removed.

#2.2 Screening criteria for drugs

table(rt1$`FDA status`)
#The table() function is used for the drug status, and the result shows that among them, 75 have undergone clinical trials and 188 have been approved by the FDA.

rt1 <- rt1[rt1$`FDA status` %in% c("FDA approved", "Clinical trial"),]
rt1 <- rt1[,-c(1, 3:6)]
write.table(rt1, file = "drug.txt",sep = "\t",row.names = F,quote = F)
#To ensure the reliability of the analysis results, the drug results that have undergone clinical trials (Clinical trial) and FDA approval (FDA approved) are selected.
#Of course, you can also choose to retain all the drugs.
#Ultimately, a total of 263 drug results are obtained and saved as a txt file for subsequent analysis.

#3. Preparation of gene expression data

rt2 <- read_excel(path = "RNA__RNA_seq_composite_expression.xls", skip = 9)
colnames(rt2) <- rt2[1,]
rt2 <- rt2[-1,-c(2:6)]
write.table(rt2, file = "geneExp.txt",sep = "\t",row.names = F,quote = F)

#Similarly, the expression data is read, sorted and saved for subsequent analysis.


#4. Drug sensitivity analysis

rm(list = ls())
#4.1 Referencing packages

library(impute)
library(limma)

#Firstly, the R packages used for the analysis are loaded, including the impute package and the limma package.

#4.2 Reading the drug input file

# rt <- read.table("drug.txt",sep="\t",header=T,check.names=F)

library(data.table)
rt <- fread("drug.txt")

rt <- as.matrix(rt)
rownames(rt) <- rt[,1]
drug <- rt[,2:ncol(rt)]
dimnames <- list(rownames(drug),colnames(drug))
data <- matrix(as.numeric(as.matrix(drug)),nrow=nrow(drug),dimnames=dimnames)
#Firstly, the previously saved drug sensitivity results are read, the corresponding row names are set, and it is converted into a matrix form.

mat <- impute.knn(data)
drug <- mat$data
drug <- avereps(drug)
#Considering the presence of some NA missing values in the drug sensitivity data, the impute.knn() function is used to evaluate and complete the drug data.
#Among them, the impute.knn() function is a function that uses the average of the nearest neighbors to estimate the missing expression data.


#4.3 Reading the expression input file

exp <- read.table("geneExp.txt", sep="\t", header=T, row.names = 1, check.names=F)
dim(exp)
exp[1:4, 1:4]
#Meanwhile, the gene expression situation in the NCI-60 cell line after sorting is read.
#The result shows: It contains the expression situation of 60 different tumor cell lines and 23,805 genes.

#4.4 Extracting the expression of specific genes

# gene <- read.table("gene.txt",sep="\t",header=F,check.names=F)
# genelist <- as.vector(gene[,1])

genelist=c('ISCA1')
genelist
#The pre-prepared target gene list is read; the result shows that it includes genes such as FANCD2, BRCA1, ABCC1, TP53, and EGFR.

genelist <- gsub(" ","",genelist)
genelist <- intersect(genelist,row.names(exp))
exp <- exp[genelist,]

#4.5 Drug sensitivity calculation
#Firstly, a new empty data frame is created to save the subsequent analysis results.
outTab <- data.frame()
for(Gene in row.names(exp)){ 
  x <- as.numeric(exp[Gene,]) 
  #Looping through drugs 
  for(Drug in row.names(drug)){   
    y <- as.numeric(drug[Drug,])   
    corT <- cor.test(x,y,method="pearson")   
    cor <- corT$estimate   
    pvalue <- corT$p.value   
    if(pvalue < 0.001){    
      outVector <- cbind(Gene,Drug,cor,pvalue)     
      outTab <- rbind(outTab,outVector)   
    } 
  }
}
#Subsequently, using for loops, the Pearson correlation coefficient between each gene expression and different drugs is calculated respectively. Based on the threshold of P value < 0.01, the analysis results are screened, and the results are saved to the variable outTab. The result shows: Ultimately, 63 correlation analysis results are obtained.
outTab <- outTab[order(as.numeric(as.vector(outTab$pvalue))),]
write.table(outTab, file="drugCor.txt", sep="\t", row.names=F, quote=F)
# Finally, the correlation analysis results are output and saved for subsequent visualization
library(ggplot2)
library(ggpubr)

plotlist_1 <- list()
corPlotNum <- 1000
if(nrow(outTab)<corPlotNum){
  corPlotNum=nrow(outTab)
}

for (i in 1:corPlotNum) {
  Gene <- outTab[i,1]
  Drug <- outTab[i,2]
  x <- as.numeric(exp[Gene,])
  y <- as.numeric(drug[Drug,])
  cor <- sprintf("%.03f",as.numeric(outTab[i,3]))
  pvalue=0
  if(as.numeric(outTab[i,4]<0.001)){
    pvalue="p<0.001"
  }else{
    pvalue=paste0("p=",sprintf("%.03f",as.numeric(outTab[i,4])))
  }
  df1 <- as.data.frame(cbind(x,y))
  p1=ggplot(data = df1,aes(x=x,y=y))+
    geom_point(size=1)+
    stat_smooth(method = "lm",se=FALSE,formula = y~x)+
    labs(x="Expression",y="IC50",title = paste0(Gene,",",Drug),subtitle = paste0("Cor=",cor,",",pvalue))
  theme(axis.ticks = element_blank(),axis.text.y = element_blank(),axis.text.x = element_blank())
  theme_bw()
  plotlist_1[[i]]=p1
  
  
}

nrow <- ceiling(sqrt(corPlotNum))
ncol <- ceiling(corPlotNum/nrow)
drugpdf <- ggarrange(plotlist = plotlist_1,nrow = nrow,ncol = ncol)

for (i in 1:corPlotNum) {
  Gene <- outTab[i,1]
  Drug <- outTab[i,2]
  x <- as.numeric(exp[Gene,])
  y <- as.numeric(drug[Drug,])
  cor <- sprintf("%.03f",as.numeric(outTab[i,3]))
  pvalue=0
  if(as.numeric(outTab[i,4]<0.001)){
    pvalue="p<0.001"
  }else{
    pvalue=paste0("p=",sprintf("%.03f",as.numeric(outTab[i,4])))
  }
  df1 <- as.data.frame(cbind(x,y))
  p1 <- ggplot(data = df1,aes(x=x,y=y)) +
    geom_point(size=1) +
    stat_smooth(method = "lm",se=FALSE,formula = y~x) +
    labs(x="Expression",y="IC50",title = paste0(Gene,",",Drug),subtitle = paste0("Cor=",cor,",",pvalue)) +
    theme(axis.ticks = element_blank(),
          axis.text.y = element_blank(),
          axis.text.x = element_blank()) +
    theme_bw() +
    theme(panel.background = element_rect(fill = "white"),
          axis.text = element_text(color = "black"))
  
  plotlist_1[[i]] <- p1
}

nrow <- ceiling(sqrt(corPlotNum))
ncol <- ceiling(corPlotNum/nrow)
ncol <- 2
nrow <- 5
drugpdf <- ggarrange(plotlist = plotlist_1, nrow = nrow, ncol = ncol)
print(drugpdf)

ggsave("drugpdf.pdf", plot = drugpdf, height = 50, width = 20)