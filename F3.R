
                          ### Clinical prognostic relevance

rm(list=ls())
setwd("D:/R/Pancancer ISCA1/live/result")

## RNA

OS <- read.csv("All_cox_result.CSV")

# Assign the value 2 to all entries in the HR column that are greater than 1

OS$Age[OS$Age > 1] <- 2
OS$Gender[OS$Gender > 1] <- 2
OS$Grade[OS$Grade > 1] <- 2
OS$M[OS$M > 1] <- 2
OS$N[OS$N > 1] <- 2
OS$T[OS$T > 1] <- 2
OS$Stage[OS$Stage > 1] <- 2

# Assign the value 1 to all entries in the HR column that are less than 1

OS$Age[OS$Age < 1] <- 1
OS$Gender[OS$Gender < 1] <- 1
OS$Grade[OS$Grade < 1] <- 1
OS$M[OS$M < 1] <- 1
OS$N[OS$N < 1] <- 1
OS$T[OS$T < 1] <- 1
OS$Stage[OS$Stage < 1] <- 1



OS$Age[OS$AgeP >= 0.05] <- 0
OS$Gender[OS$GenderP >= 0.05] <- 0
OS$Grade[OS$GradeP >= 0.05] <- 0
OS$M[OS$MP >= 0.05] <- 0
OS$N[OS$NP >= 0.05] <- 0
OS$Stage[OS$StageP >= 0.05] <- 0
OS$T[OS$TP >= 0.05] <- 0


OS$CancerCode <- gsub("\\(.*\\)", "", OS$CancerCode)
rownames(OS) <- OS$CancerCode
OS <- OS[,-1]
OS <- OS[,c(-1,-3,-5,-7,-9,-11,-13)]
rownames(OS) <- gsub('TCGA-', '', rownames(OS))
colnames(OS) <- sub("P$", "", colnames(OS))
final_sheet_temp <- OS
final_sheet <- OS
final_sheet <- as.matrix(final_sheet)
final_sheet_temp <- as.matrix(final_sheet_temp)

library(ComplexHeatmap)
colors = structure(c("white","#00A087B2","#BB362F"), names = c(0,1,2))
row_boxplot <- rowAnnotation(
  . = anno_barplot(rowSums(final_sheet_temp)[order(-rowSums(final_sheet_temp))]/2,
                   bar_width = 0.8,
                   axis = F,
                   border=F,
                   gp = gpar(fill = 1:10),
                   add_numbers=T,
                   numbers_gp = gpar(fontsize = 12),
                   numbers_offset = unit(4, "mm"),
                   width = unit(20, "mm"))
)




pdf('complexheatmap-Fenqi.pdf',width=24,height=40)
Heatmap(final_sheet,
        na_col='grey',
        col=colors,
        cluster_columns = F,
        cluster_rows = F,
        height=nrow(final_sheet)*unit(8, "mm"),
        width=ncol(final_sheet)*unit(16, "mm"),
        right_annotation = row_boxplot,
        rect_gp = gpar(col = "black", lwd = 1),
        row_names_side = 'left',
        column_gap = unit(5, "mm"),
        column_names_gp = gpar(col = c("#E64B35B2")),
        column_title_gp = gpar(fill = c("#E64B35B2")))
dev.off()


                            ### Isoform prognostic relevance
rm(list=ls())
setwd("D:/R/Pancancer ISCA1/Survival/Isoform/")


### ENST00000311534

OS <- read.csv("ENST00000311534-OS-Isoform Results.csv")
DFS <- read.csv("ENST00000311534-DFS-Isoform Results.csv")
DSS <- read.csv("ENST00000311534-DSS-Isoform Results.csv")
PFS <- read.csv("ENST00000311534-PFS-Isoform Results.csv")


# Assign 2 to values in the HR column that are greater than 1
OS$HR[OS$HR > 1] <- 2
OS$HR[OS$HR < 1] <- 1
DFS$HR[DFS$HR > 1] <- 2
DFS$HR[DFS$HR < 1] <- 1
DSS$HR[DSS$HR > 1] <- 2
DSS$HR[DSS$HR < 1] <- 1
PFS$HR[PFS$HR > 1] <- 2
PFS$HR[PFS$HR < 1] <- 1

OS$HR[OS$pvalue >= 0.05] <- 0
DFS$HR[DFS$pvalue >= 0.05] <- 0
DSS$HR[DSS$pvalue >= 0.05] <- 0
PFS$HR[PFS$pvalue >= 0.05] <- 0

OS$CancerCode <- gsub("\\(.*\\)", "", OS$CancerCode)
DSS$CancerCode <- gsub("\\(.*\\)", "", DSS$CancerCode)
DFS$CancerCode <- gsub("\\(.*\\)", "", DFS$CancerCode)
PFS$CancerCode <- gsub("\\(.*\\)", "", PFS$CancerCode)


rownames(OS) <- OS$CancerCode
rownames(DSS) <- DSS$CancerCode
rownames(DFS) <- DFS$CancerCode
rownames(PFS) <- PFS$CancerCode



# Merge the "HR" columns of the data frames based on the values of the "CancerCode" column
merged_df <- merge(merge(merge(OS, DSS, by = "CancerCode", all = TRUE), DFS, by = "CancerCode", all = TRUE), PFS, by = "CancerCode", all = TRUE)

# Select the "HR" column of the merged data frame
final_merged_HR <- merged_df[, grepl("HR", names(merged_df))]


rownames(final_merged_HR) <- merged_df$CancerCode

colnames(final_merged_HR) <- c("OS","DSS","DFS","PFS")

final_merged_HR$gene <- rownames(final_merged_HR)

colnames(final_merged_HR)[1:4] = paste0('ENST00000311534 ',colnames(final_merged_HR)[1:4])

final_sheet_temp <- final_merged_HR



### ENST00000375991

OS <- read.csv("ENST00000375991-OS-Isoform Results.csv")
DFS <- read.csv("ENST00000375991-DFS-Isoform Results.csv")
DSS <- read.csv("ENST00000375991-DSS-Isoform Results.csv")
PFS <- read.csv("ENST00000375991-PFS-Isoform Results.csv")


# Assign 2 to values in the HR column that are greater than 1
OS$HR[OS$HR > 1] <- 2
OS$HR[OS$HR < 1] <- 1
DFS$HR[DFS$HR > 1] <- 2
DFS$HR[DFS$HR < 1] <- 1
DSS$HR[DSS$HR > 1] <- 2
DSS$HR[DSS$HR < 1] <- 1
PFS$HR[PFS$HR > 1] <- 2
PFS$HR[PFS$HR < 1] <- 1

OS$HR[OS$pvalue >= 0.05] <- 0
DFS$HR[DFS$pvalue >= 0.05] <- 0
DSS$HR[DSS$pvalue >= 0.05] <- 0
PFS$HR[PFS$pvalue >= 0.05] <- 0

OS$CancerCode <- gsub("\\(.*\\)", "", OS$CancerCode)
DSS$CancerCode <- gsub("\\(.*\\)", "", DSS$CancerCode)
DFS$CancerCode <- gsub("\\(.*\\)", "", DFS$CancerCode)
PFS$CancerCode <- gsub("\\(.*\\)", "", PFS$CancerCode)


rownames(OS) <- OS$CancerCode
rownames(DSS) <- DSS$CancerCode
rownames(DFS) <- DFS$CancerCode
rownames(PFS) <- PFS$CancerCode



# Merge the "HR" columns of the data frames based on the values of the "CancerCode" column
merged_df <- merge(merge(merge(OS, DSS, by = "CancerCode", all = TRUE), DFS, by = "CancerCode", all = TRUE), PFS, by = "CancerCode", all = TRUE)

# Select the "HR" column of the merged data frame
final_merged_HR <- merged_df[, grepl("HR", names(merged_df))]


rownames(final_merged_HR) <- merged_df$CancerCode

colnames(final_merged_HR) <- c("OS","DSS","DFS","PFS")

final_merged_HR$gene <- rownames(final_merged_HR)

colnames(final_merged_HR)[1:4] = paste0('ENST00000375991 ',colnames(final_merged_HR)[1:4])



# Load the dplyr package
library(dplyr)

# Merge the data frames using full_join() based on the "gene" column
final_sheet_temp <- full_join(final_sheet_temp, final_merged_HR, by = "gene")





### ENST00000326094

OS <- read.csv("ENST00000326094-OS-Isoform Results.csv")
DFS <- read.csv("ENST00000326094-DFS-Isoform Results.csv")
DSS <- read.csv("ENST00000326094-DSS-Isoform Results.csv")
PFS <- read.csv("ENST00000326094-PFS-Isoform Results.csv")


# Assign 2 to values in the HR column that are greater than 1
OS$HR[OS$HR > 1] <- 2
OS$HR[OS$HR < 1] <- 1
DFS$HR[DFS$HR > 1] <- 2
DFS$HR[DFS$HR < 1] <- 1
DSS$HR[DSS$HR > 1] <- 2
DSS$HR[DSS$HR < 1] <- 1
PFS$HR[PFS$HR > 1] <- 2
PFS$HR[PFS$HR < 1] <- 1

OS$HR[OS$pvalue >= 0.05] <- 0
DFS$HR[DFS$pvalue >= 0.05] <- 0
DSS$HR[DSS$pvalue >= 0.05] <- 0
PFS$HR[PFS$pvalue >= 0.05] <- 0

OS$CancerCode <- gsub("\\(.*\\)", "", OS$CancerCode)
DSS$CancerCode <- gsub("\\(.*\\)", "", DSS$CancerCode)
DFS$CancerCode <- gsub("\\(.*\\)", "", DFS$CancerCode)
PFS$CancerCode <- gsub("\\(.*\\)", "", PFS$CancerCode)


rownames(OS) <- OS$CancerCode
rownames(DSS) <- DSS$CancerCode
rownames(DFS) <- DFS$CancerCode
rownames(PFS) <- PFS$CancerCode



# Merge the "HR" columns of the data frames based on the values of the "CancerCode" column
merged_df <- merge(merge(merge(OS, DSS, by = "CancerCode", all = TRUE), DFS, by = "CancerCode", all = TRUE), PFS, by = "CancerCode", all = TRUE)

# Select the "HR" column of the merged data frame
final_merged_HR <- merged_df[, grepl("HR", names(merged_df))]


rownames(final_merged_HR) <- merged_df$CancerCode

colnames(final_merged_HR) <- c("OS","DSS","DFS","PFS")

final_merged_HR$gene <- rownames(final_merged_HR)

colnames(final_merged_HR)[1:4] = paste0('ENST00000326094 ',colnames(final_merged_HR)[1:4])



# Load the dplyr package
library(dplyr)

# Merge the data frames using full_join() based on the "gene" column
final_sheet_temp <- full_join(final_sheet_temp, final_merged_HR, by = "gene")



# Delete the row data with ID "name1" and "name2" in the first column
final_sheet_temp<- subset(final_sheet_temp,!(gene %in% c("TCGA-COADREAD", "TCGA-GBMLGG","TCGA-KIPAN","TCGA-SKCM-M",
                                                         "TCGA-SKCM-P","TCGA-STES")))
# Remove the "TCGA-" prefix
final_sheet_temp$gene <- gsub('TCGA-', '', final_sheet_temp$gene)

rownames(final_sheet_temp) <- final_sheet_temp$gene

final_sheet_temp <- final_sheet_temp[,-5]

final_sheet <- final_sheet_temp

final_sheet <- as.matrix(final_sheet)
final_sheet_temp <- as.matrix(final_sheet_temp)





library(ComplexHeatmap)
colors <- structure(c("white", "#00A087B2", "#BB362F"), names = c(0, 1, 2))

row_sums <- rowSums(final_sheet_temp, na.rm = TRUE)
row_boxplot <- rowAnnotation(
  . = anno_barplot(row_sums[order(-row_sums)] / 2,
                   bar_width = 0.8,
                   axis = FALSE,
                   border = FALSE,
                   gp = gpar(fill = 1:10),
                   add_numbers = row_sums,
                   numbers_gp = gpar(fontsize = 12),
                   numbers_offset = unit(4, "mm"),
                   width = unit(20, "mm"))
)



pdf('complexheatmap-Isoform.pdf',width=24,height=40)
Heatmap(final_sheet,
        na_col='grey',
        col=colors,
        cluster_columns = F,
        cluster_rows = F,
        height=nrow(final_sheet)*unit(8, "mm"),
        width=ncol(final_sheet)*unit(12, "mm"),
        right_annotation = row_boxplot,
        rect_gp = gpar(col = "black", lwd = 1),
        column_split = c(rep("ENST00000311534",4),rep("ENST00000375991",4),rep("ENST00000326094",4)),
        row_names_side = 'left',
        column_gap = unit(5, "mm"),
        column_names_gp = gpar(col = c("#E64B35B2", "#4DBBD5B2" ,"#F5AE6B")),
        column_title_gp = gpar(fill = c("#E64B35B2", "#4DBBD5B2" ,"#F5AE6B")))
dev.off()



                                  ### Clinical prognostic KM

rm(list = ls())

setwd("D:/R/Pancancer ISCA1/Survival/MutationSurvival")

library(data.table)
library("maxstat")
library("survival")

classification <- c("Age")

cancers <- c("ACC","BLCA","BRCA","CESC","CHOL","COAD","DLBC","ESCA","GBM","HNSC","KICH",
             "KIRC","KIRP","LAML","LGG","LIHC","LUAD","LUSC","MESO","OV","PAAD","PCPG",
             "PRAD","READ","SARC","SKCM","STAD","TGCT","THCA","THYM","UCEC","UCS","UVM")


HR_value <- c()

for (cancerType in cancers) {
  
  cancerType <-  paste0("TCGA-",cancerType,sep='')
  
  dataset0 <-  read.csv("OS.CSV")
  
  dataset0 <- as.data.frame(dataset0)
  
  df1=c()
  
  # Use subset function to filter rows in df2 with row name "cancerType"
  subset_df <- subset(dataset0, cancer == cancerType)
  
  # Merge subset_df to df1
  df1 <- rbind(df1, subset_df)
  
  dataset0=df1
  
  rownames(dataset0)=dataset0$samplename
  dataset0=dataset0[,-1:-2]
  
  connam <- colnames(dataset0)
  i=1
  colnames(dataset0)=c("ISCA1","Time","Status")
  dataset <- data.frame(dataset0[,2:3],dataset0[,i])
  library("glmnet")
  library("survival")
  library("survminer")
  rt <- dataset
  rt$risk <- c("Low")
  cluster <- read.csv(paste("D:/R/Pancancer ISCA1/Clinical/gene/classification/",classification,".csv",sep=''),header = T)
  
  colnames(cluster) <- c("X","Range","Age","SampleName")
  
  matching_rows <- intersect(cluster$SampleName, rownames(rt))
  
  cluster <- cluster[match(matching_rows,cluster$SampleName),]
  
  rt <- rt[match(matching_rows,rownames(rt)),]
  match_rows <- match(rownames(rt), cluster$SampleName)
  rt$risk[match_rows] <- cluster$Range
  write.csv(rt,file= paste(cancerType,"-","rtp.csv",sep=''),quote=F)
  rt$Time  =rt$Time
  #rt$Time  =rt$Time /365
  if(length(rownames(rt))==0){
    hr_coef <- NA
  }
  else{
    cox_model <- coxph(Surv(Time, Status) ~ risk, data = rt)
    hr_coef <- round(exp(cox_model$coefficients), digits = 2)
    fit <- survfit(Surv(Time, Status) ~ risk, data = rt)
    ggsurvplot(fit, data = rt)
    library(survival)
    library(survminer)
    a=3
    fit <- survfit(Surv(Time, Status) ~ risk, data = rt)
    write.csv(rt,file= paste("rt.csv",sep=''),quote=F)
    pdf(file=paste0("classification/",classification,"-",cancerType,".pdf",sep=''),onefile = FALSE,width = 6,height =7.5)
    surPlot=ggsurvplot(fit,
                       data=rt,
                       #font.title = paste(connam[i],sep=''),
                       #ggtitle = paste(connam[i],sep=''),
                       #conf.int=TRUE,
                       legend = "top",
                       legend.title="Risk",
                       pval=TRUE,
                       #pval.method = TRUE,
                       pval.size = 5,
                       xlab="Time",
                       break.time.by = ceiling((max(rt$Time))/4),
                       risk.table.title="",
                       palette=c("#009E73","red"),
                       risk.table=T,
                       risk.table.height=.25,)
    surPlot <- surPlot +
      labs(title = paste(" (HR = ", hr_coef, ")", sep=''))
    print(surPlot)
    dev.off()
  }
  hr_coef <- data.frame(hr_coef)
  rownames(hr_coef) <- cancerType
  if(length(HR_value)==0){
    HR_value <- hr_coef
  }else {
    HR_value <- rbind(HR_value,hr_coef)
  }
}
write.csv(HR_value,paste("classification/",classification,".csv"))




### age
rm (list = ls ())
setwd("D:/R/Pancancer ISCA1/Survival/MutationSurvival")
library(data.table)
library ("maxstat")
library("survival")
fenqi <- c("Age")
cancers <- c("ACC","BLCA","BRCA","CESC","CHOL","COAD","DLBC","ESCA","GBM","HNSC","KICH",
             "KIRC","KIRP","LAML","LGG","LIHC","LUAD","LUSC","MESO","OV","PAAD","PCPG",
             "PRAD","READ","SARC","SKCM","STAD","TGCT","THCA","THYM","UCEC","UCS","UVM")


HR_value <- c()

for (cancertype in cancers) {
  
  cancertype <-  paste0("TCGA-",cancertype,sep='')
  
  
  dataset0 <-  read.csv("OS.CSV")
  
  dataset0 <- as.data.frame(dataset0)
  
  df1=c()
  
  
  subset_df <- subset(dataset0, cancer == cancertype)
  
  
  df1 <- rbind(df1, subset_df)
  
  dataset0=df1
  
  rownames(dataset0)=dataset0$samplename
  dataset0=dataset0[,-1:-2]
  
  connam <- colnames(dataset0)
  i=1
  colnames(dataset0)=c(" ISCA1","Time","Status")
  dataset <- data.frame(dataset0[,2:3],dataset0[,i])
  
  library("glmnet")
  library("survival")
  library("survminer")
  
  
  rt <- dataset
  rt$risk <- c("Low")
  
  
  cluster <- read.csv(paste("D:/R/Pancancer ISCA1/linchuang/gene/fenqi/",fenqi,".csv",sep=''),header = T)
  colnames(cluster) <- c("X","Range","Age","SampleName")
  
  matching_rows <- intersect(cluster$SampleName, rownames(rt))
  cluster <- cluster[match(matching_rows,cluster$SampleName),]
  rt <- rt[match(matching_rows,rownames(rt)),]
  match_rows <- match(rownames(rt), cluster$SampleName)
  rt$risk[match_rows] <- cluster$Range
  write.csv(rt,file= paste(cancertype,"-","rtp.csv",sep=''),quote=F)
  rt$Time  =rt$Time
  #rt$Time  =rt$Time /365
  if(length(rownames(rt))==0){
    hr_coef <- NA
  }
  else{
    cox_model <- coxph(Surv(Time, Status) ~ risk, data = rt)
    hr_coef <- round(exp(cox_model$coefficients), digits = 2)
    fit <- survfit(Surv(Time, Status) ~ risk, data = rt)
    ggsurvplot(fit, data = rt)
    library(survival)
    library(survminer)
    a=3
    fit <- survfit(Surv(Time, Status) ~ risk, data = rt)
    write.csv(rt,file= paste("rt.csv",sep=''),quote=F)
    pdf(file=paste0("fenqi/",fenqi,"-",cancertype,".pdf",sep=''),onefile = FALSE,width = 6,height =7.5)
    surPlot=ggsurvplot(fit,
                       data=rt,
                       #font.title = paste(connam[i],sep=''),
                       #ggtitle = paste(connam[i],sep=''),
                       #conf.int=TRUE,
                       legend = "top",
                       legend.title="Risk",
                       pval=TRUE,
                       #pval.method = TRUE,
                       pval.size = 5 ,
                       xlab="Time",
                       break.time.by = ceiling((max(rt$Time))/4),
                       risk.table.title="",
                       palette=c("#009E73","red"),
                       risk.table=T,
                       risk.table.height=.25,)
    surPlot <- surPlot +
      labs(title = paste(" (HR = ", hr_coef, ")", sep=''))
    print(surPlot)
    dev.off()
  }
  hr_coef <- data.frame(hr_coef)
  rownames(hr_coef) <- cancertype
  
  if(length(HR_value)==0){
    HR_value <- hr_coef
  }else {
    HR_value <- rbind(HR_value,hr_coef)
  }
  
}
write.csv(HR_value,paste("fenqi/",fenqi,".csv"))



### Gender
rm (list = ls ())
setwd(“D:/R/Pancancer ISCA1/Survival/MutationSurvival”)
library(data.table)
library ("maxstat")
library("survival")
fenqi <- c("Gender")
cancers <- c("ACC","BLCA","BRCA","CESC","CHOL","COAD","DLBC","ESCA","GBM","HNSC","KICH",
             "KIRC","KIRP","LAML","LGG","LIHC","LUAD","LUSC","MESO","OV","PAAD","PCPG",
             "PRAD","READ","SARC","SKCM","STAD","TGCT","THCA","THYM","UCEC","UCS","UVM")
HR_value <- c()

for (cancertype in cancers) {
  cancertype <-  paste0("TCGA-",cancertype,sep='')
  dataset0 <-  read.csv("OS.CSV")
  dataset0 <- as.data.frame(dataset0)
  df1=c()
  subset_df <- subset(dataset0, cancer == cancertype)
  df1 <- rbind(df1, subset_df)
  dataset0=df1
  rownames(dataset0)=dataset0$samplename
  dataset0=dataset0[,-1:-2]
  connam <- colnames(dataset0)
  i=1
  colnames(dataset0)=c("ISCA1","Time","Status")
  dataset <- data.frame(dataset0[,2:3],dataset0[,i])
  library("glmnet")
  library("survival")
  library("survminer")
  rt <- dataset
  rt$risk <- c("Low")
  cluster <- read.csv(paste("D:/R/Pancancer ISCA1/linchuang/gene/fenqi/",fenqi,".csv",sep=''),header = T)
  colnames(cluster) <- c("X","Range","Age","SampleName")
  matching_rows <- intersect(cluster$SampleName, rownames(rt))
  cluster <- cluster[match(matching_rows,cluster$SampleName),]
  rt <- rt[match(matching_rows,rownames(rt)),]
  match_rows <- match(rownames(rt), cluster$SampleName)
  rt$risk[match_rows] <- cluster$Range
  write.csv(rt,file= paste(cancertype,"-","rtp.csv",sep=''),quote=F)
  rt$Time  =rt$Time
  #rt$Time  =rt$Time /365
  if(length(rownames(rt))==0){
    hr_coef <- NA
  }
  else{
    cox_model <- coxph(Surv(Time, Status) ~ risk, data = rt)
    hr_coef <- round(exp(cox_model$coefficients), digits = 2)
    fit <- survfit(Surv(Time, Status) ~ risk, data = rt)
    ggsurvplot(fit, data = rt)
    library(survival)
    library(survminer)
    a=3
    fit <- survfit(Surv(Time, Status) ~ risk, data = rt)
    write.csv(rt,file= paste("rt.csv",sep=''),quote=F)
    pdf(file=paste0("fenqi/",fenqi,"-",cancertype,".pdf",sep=''),onefile = FALSE,width = 6,height =7.5)
    surPlot=ggsurvplot(fit,
                       data=rt,
                       #font.title = paste(connam[i],sep=''),
                       #ggtitle = paste(connam[i],sep=''),
                       #conf.int=TRUE,
                       legend = "top",
                       legend.title="Risk",
                       pval=TRUE,
                       #pval.method = TRUE,
                       pval.size = 5 ,
                       xlab="Time",
                       break.time.by = ceiling((max(rt$Time))/4),
                       risk.table.title="",
                       palette=c("#009E73","red"),
                       risk.table=T,
                       risk.table.height=.25,)
    
    surPlot <- surPlot +
      labs(title = paste(" (HR = ", hr_coef, ")", sep=''))
    
    print(surPlot)
    dev.off()
  }
  hr_coef <- data.frame(hr_coef)
  rownames(hr_coef) <- cancertype
  if(length(HR_value)==0){
    HR_value <- hr_coef
  }else {
    HR_value <- rbind(HR_value,hr_coef)
  }
  
}
write.csv(HR_value,paste("fenqi/",fenqi,".csv"))







### Gender

rm (list = ls ())
setwd(“D:/R/Pancancer ISCA1/Survival/MutationSurvival”)
library(data.table)
library ("maxstat")
library("survival")
fenqi <- c("Grade")
cancers <- c("ACC","BLCA","BRCA","CESC","CHOL","COAD","DLBC","ESCA","GBM","HNSC","KICH",
             "KIRC","KIRP","LAML","LGG","LIHC","LUAD","LUSC","MESO","OV","PAAD","PCPG",
             "PRAD","READ","SARC","SKCM","STAD","TGCT","THCA","THYM","UCEC","UCS","UVM")

HR_value <- c()
for (cancertype in cancers) {
  cancertype <-  paste0("TCGA-",cancertype,sep='')
  dataset0 <-  read.csv("OS.CSV")
  dataset0 <- as.data.frame(dataset0)
  df1=c()
  subset_df <- subset(dataset0, cancer == cancertype)
  df1 <- rbind(df1, subset_df)
  dataset0=df1
  rownames(dataset0)=dataset0$samplename
  dataset0=dataset0[,-1:-2]
  connam <- colnames(dataset0)
  i=1
  colnames(dataset0)=c("ISCA1","Time","Status")
  dataset <- data.frame(dataset0[,2:3],dataset0[,i])
  library("glmnet")
  library("survival")
  library("survminer")
  rt <- dataset
  rt$risk <- c("Low")
  cluster <- read.csv(paste("D:/R/Pancancer ISCA1/linchuang/gene/fenqi/",fenqi,".csv",sep=''),header = T)
  colnames(cluster) <- c("X","Range","Age","SampleName")
  cluster$Range <- gsub("G1", "G1~2", cluster$Range)
  cluster$Range <- gsub("G2", "G1~2", cluster$Range)
  cluster$Range <- gsub("G3", "G3~4", cluster$Range)
  cluster$Range <- gsub("G4", "G3~4", cluster$Range)
  matching_rows <- intersect(cluster$SampleName, rownames(rt))
  cluster <- cluster[match(matching_rows,cluster$SampleName),]
  rt <- rt[match(matching_rows,rownames(rt)),]
  match_rows <- match(rownames(rt), cluster$SampleName)
  rt$risk[match_rows] <- cluster$Range
  write.csv(rt,file= paste(cancertype,"-","rtp.csv",sep=''),quote=F)
  rt$Time  =rt$Time
  #rt$Time  =rt$Time /365
  if(length(rownames(rt))==0){
    hr_coef <- NA
  }
  else{
    cox_model <- coxph(Surv(Time, Status) ~ risk, data = rt)
    hr_coef <- round(exp(cox_model$coefficients), digits = 2)
    fit <- survfit(Surv(Time, Status) ~ risk, data = rt)
    ggsurvplot(fit, data = rt)
    library(survival)
    library(survminer)
    a=3
    fit <- survfit(Surv(Time, Status) ~ risk, data = rt)
    write.csv(rt,file= paste("rt.csv",sep=''),quote=F)
    pdf(file=paste0("fenqi/",fenqi,"-",cancertype,".pdf",sep=''),onefile = FALSE,width = 6,height =7.5)
    surPlot=ggsurvplot(fit,
                       data=rt,
                       #font.title = paste(connam[i],sep=''),
                       #ggtitle = paste(connam[i],sep=''),
                       #conf.int=TRUE,
                       legend = "top",
                       legend.title="Risk",
                       pval=TRUE,
                       #pval.method = TRUE,
                       pval.size = 5 ,
                       xlab="Time",
                       break.time.by = ceiling((max(rt$Time))/4),
                       risk.table.title="",
                       palette=c("#009E73","red"),
                       risk.table=T,
                       risk.table.height=.25,)
    
    surPlot <- surPlot +
      labs(title = paste(" (HR = ", hr_coef, ")", sep=''))
    
    print(surPlot)
    dev.off()
  }
  
  hr_coef <- data.frame(hr_coef)
  rownames(hr_coef) <- cancertype
  
  if(length(HR_value)==0){
    HR_value <- hr_coef
  }else {
    HR_value <- rbind(HR_value,hr_coef)
  }
  
}
write.csv(HR_value,paste("fenqi/",fenqi,".csv"))







### M

rm (list = ls ())
setwd(“D:/R/Pancancer ISCA1/Survival/MutationSurvival”)
library(data.table)
library ("maxstat")
library("survival")
fenqi <- c("M")
cancers <- c("ACC","BLCA","BRCA","CESC","CHOL","COAD","DLBC","ESCA","GBM","HNSC","KICH",
             "KIRC","KIRP","LAML","LGG","LIHC","LUAD","LUSC","MESO","OV","PAAD","PCPG",
             "PRAD","READ","SARC","SKCM","STAD","TGCT","THCA","THYM","UCEC","UCS","UVM")
HR_value <- c()
for (cancertype in cancers) {
  cancertype <-  paste0("TCGA-",cancertype,sep='')
  dataset0 <-  read.csv("OS.CSV")
  dataset0 <- as.data.frame(dataset0)
  df1=c()
  subset_df <- subset(dataset0, cancer == cancertype)
  df1 <- rbind(df1, subset_df)
  dataset0=df1
  rownames(dataset0)=dataset0$samplename
  dataset0=dataset0[,-1:-2]
  connam <- colnames(dataset0)
  i=1
  colnames(dataset0)=c("ISCA1","Time","Status")
  dataset <- data.frame(dataset0[,2:3],dataset0[,i])
  library("glmnet")
  library("survival")
  library("survminer")
  rt <- dataset
  rt$risk <- c("Low")
  cluster <- read.csv(paste("D:/R/Pancancer ISCA1/linchuang/gene/fenqi/",fenqi,".csv",sep=''),header = T)
  colnames(cluster) <- c("X","Range","Age","SampleName")
  matching_rows <- intersect(cluster$SampleName, rownames(rt))
  cluster <- cluster[match(matching_rows,cluster$SampleName),]
  rt <- rt[match(matching_rows,rownames(rt)),]
  match_rows <- match(rownames(rt), cluster$SampleName)
  rt$risk[match_rows] <- cluster$Range
  write.csv(rt,file= paste(cancertype,"-","rtp.csv",sep=''),quote=F)
  rt$Time  =rt$Time
  #rt$Time  =rt$Time /365
  if(length(rownames(rt))==0){
    hr_coef <- NA
    
  }
  else{
    cox_model <- coxph(Surv(Time, Status) ~ risk, data = rt)
    hr_coef <- round(exp(cox_model$coefficients), digits = 2)
    fit <- survfit(Surv(Time, Status) ~ risk, data = rt)
    ggsurvplot(fit, data = rt)
    library(survival)
    library(survminer)
    a=3
    fit <- survfit(Surv(Time, Status) ~ risk, data = rt)
    write.csv(rt,file= paste("rt.csv",sep=''),quote=F)
    pdf(file=paste0("fenqi/",fenqi,"-",cancertype,".pdf",sep=''),onefile = FALSE,width = 6,height =7.5)
    surPlot=ggsurvplot(fit,
                       data=rt,
                       #font.title = paste(connam[i],sep=''),
                       #ggtitle = paste(connam[i],sep=''),
                       #conf.int=TRUE,
                       legend = "top",
                       legend.title="Risk",
                       pval=TRUE,
                       #pval.method = TRUE,
                       pval.size = 5 ,
                       xlab="Time",
                       break.time.by = ceiling((max(rt$Time))/4),
                       risk.table.title="",
                       palette=c("#009E73","red"),
                       risk.table=T,
                       risk.table.height=.25,)
    
    surPlot <- surPlot +
      labs(title = paste(" (HR = ", hr_coef, ")", sep=''))
    
    print(surPlot)
    dev.off()
  }
  
  hr_coef <- data.frame(hr_coef)
  rownames(hr_coef) <- cancertype
  
  if(length(HR_value)==0){
    HR_value <- hr_coef
  }else {
    HR_value <- rbind(HR_value,hr_coef)
  }
  
}

write.csv(HR_value,paste("fenqi/",fenqi,".csv"))







### N
rm (list = ls ())
setwd(“D:/R/Pancancer ISCA1/Survival/MutationSurvival”)
library(data.table)
library ("maxstat")
library("survival")
fenqi <- c("N")
cancers <- c("ACC","BLCA","BRCA","CESC","CHOL","COAD","DLBC","ESCA","GBM","HNSC","KICH",
             "KIRC","KIRP","LAML","LGG","LIHC","LUAD","LUSC","MESO","OV","PAAD","PCPG",
             "PRAD","READ","SARC","SKCM","STAD","TGCT","THCA","THYM","UCEC","UCS","UVM")
HR_value <- c()
for (cancertype in cancers) {
  cancertype <-  paste0("TCGA-",cancertype,sep='')
  dataset0 <-  read.csv("OS.CSV")
  dataset0 <- as.data.frame(dataset0)
  df1=c()
  subset_df <- subset(dataset0, cancer == cancertype)
  df1 <- rbind(df1, subset_df)
  dataset0=df1
  rownames(dataset0)=dataset0$samplename
  dataset0=dataset0[,-1:-2]
  connam <- colnames(dataset0)
  i=1
  colnames(dataset0)=c(" ISCA1","Time","Status")
  dataset <- data.frame(dataset0[,2:3],dataset0[,i])
  library("glmnet")
  library("survival")
  library("survminer")
  rt <- dataset
  rt$risk <- c("Low")
  cluster <- read.csv(paste("D:/R/Pancancer ISCA1/linchuang/gene/fenqi/",fenqi,".csv",sep=''),header = T)
  colnames(cluster) <- c("X","Range","Age","SampleName")
  cluster$Range <- gsub("N0", "N0", cluster$Range)
  cluster$Range <- gsub("N1", "N1~3", cluster$Range)
  cluster$Range <- gsub("N2", "N1~3", cluster$Range)
  cluster$Range <- gsub("N3", "N1~3", cluster$Range)
  matching_rows <- intersect(cluster$SampleName, rownames(rt))
  cluster <- cluster[match(matching_rows,cluster$SampleName),]
  rt <- rt[match(matching_rows,rownames(rt)),]
  match_rows <- match(rownames(rt), cluster$SampleName)
  rt$risk[match_rows] <- cluster$Range
  write.csv(rt,file= paste(cancertype,"-","rtp.csv",sep=''),quote=F)
  rt$Time  =rt$Time
  #rt$Time  =rt$Time /365
  if(length(rownames(rt))==0){
    hr_coef <- NA
  }
  else{
    cox_model <- coxph(Surv(Time, Status) ~ risk, data = rt)
    hr_coef <- round(exp(cox_model$coefficients), digits = 2)
    fit <- survfit(Surv(Time, Status) ~ risk, data = rt)
    ggsurvplot(fit, data = rt)
    library(survival)
    library(survminer)
    a=3
    fit <- survfit(Surv(Time, Status) ~ risk, data = rt)
    write.csv(rt,file= paste("rt.csv",sep=''),quote=F)
    pdf(file=paste0("fenqi/",fenqi,"-",cancertype,".pdf",sep=''),onefile = FALSE,width = 6,height =7.5)
    surPlot=ggsurvplot(fit,
                       data=rt,
                       #font.title = paste(connam[i],sep=''),
                       #ggtitle = paste(connam[i],sep=''),
                       #conf.int=TRUE,
                       legend = "top",
                       legend.title="Risk",
                       pval=TRUE,
                       #pval.method = TRUE,
                       pval.size = 5 ,
                       xlab="Time",
                       break.time.by = ceiling((max(rt$Time))/4),
                       risk.table.title="",
                       palette=c("#009E73","red"),
                       risk.table=T,
                       risk.table.height=.25,)
    
    surPlot <- surPlot +
      labs(title = paste(" (HR = ", hr_coef, ")", sep=''))
    
    print(surPlot)
    dev.off()
  }
  
  hr_coef <- data.frame(hr_coef)
  rownames(hr_coef) <- cancertype
  
  if(length(HR_value)==0){
    HR_value <- hr_coef
  }else {
    HR_value <- rbind(HR_value,hr_coef)
  }
}
write.csv(HR_value,paste("fenqi/",fenqi,".csv"))






### T

rm (list = ls ())
setwd(“D:/R/Pancancer ISCA1/Survival/MutationSurvival”)
library(data.table)
library ("maxstat")
library("survival")
fenqi <- c("T")

cancers <- c("ACC","BLCA","BRCA","CESC","CHOL","COAD","DLBC","ESCA","GBM","HNSC","KICH",
             "KIRC","KIRP","LAML","LGG","LIHC","LUAD","LUSC","MESO","OV","PAAD","PCPG",
             "PRAD","READ","SARC","SKCM","STAD","TGCT","THCA","THYM","UCEC","UCS","UVM")
HR_value <- c()

for (cancertype in cancers) {
  cancertype <-  paste0("TCGA-",cancertype,sep='')
  dataset0 <-  read.csv("OS.CSV")
  dataset0 <- as.data.frame(dataset0)
  df1=c()
  subset_df <- subset(dataset0, cancer == cancertype)
  df1 <- rbind(df1, subset_df)
  dataset0=df1
  rownames(dataset0)=dataset0$samplename
  dataset0=dataset0[,-1:-2]
  connam <- colnames(dataset0)
  i=1
  colnames(dataset0)=c(" ISCA1","Time","Status")
  dataset <- data.frame(dataset0[,2:3],dataset0[,i])
  library("glmnet")
  library("survival")
  library("survminer")
  rt <- dataset
  rt$risk <- c("Low")
  cluster <- read.csv(paste("D:/R/Pancancer ISCA1/linchuang/gene/fenqi/",fenqi,".csv",sep=''),header = T)
  colnames(cluster) <- c("X","Range","Age","SampleName")
  cluster$Range <- gsub("T1", "T1~2", cluster$Range)
  cluster$Range <- gsub("T2", "T1~2", cluster$Range)
  cluster$Range <- gsub("T3", "T3~4", cluster$Range)
  cluster$Range <- gsub("T4", "T3~4", cluster$Range)
  matching_rows <- intersect(cluster$SampleName, rownames(rt))
  cluster <- cluster[match(matching_rows,cluster$SampleName),]
  rt <- rt[match(matching_rows,rownames(rt)),]
  match_rows <- match(rownames(rt), cluster$SampleName)
  rt$risk[match_rows] <- cluster$Range
  write.csv(rt,file= paste(cancertype,"-","rtp.csv",sep=''),quote=F)
  rt$Time  =rt$Time
  #rt$Time  =rt$Time /365
  if(length(rownames(rt))==0){
    hr_coef <- NA
  }
  else{
    cox_model <- coxph(Surv(Time, Status) ~ risk, data = rt)
    hr_coef <- round(exp(cox_model$coefficients), digits = 2)
    fit <- survfit(Surv(Time, Status) ~ risk, data = rt)
    ggsurvplot(fit, data = rt)
    library(survival)
    library(survminer)
    a=3
    fit <- survfit(Surv(Time, Status) ~ risk, data = rt)
    write.csv(rt,file= paste("rt.csv",sep=''),quote=F)
    pdf(file=paste0("fenqi/",fenqi,"-",cancertype,".pdf",sep=''),onefile = FALSE,width = 6,height =7.5)
    surPlot=ggsurvplot(fit,
                       data=rt,
                       #font.title = paste(connam[i],sep=''),
                       #ggtitle = paste(connam[i],sep=''),
                       #conf.int=TRUE,
                       legend = "top",
                       legend.title="Risk",
                       pval=TRUE,
                       #pval.method = TRUE,
                       pval.size = 5 ,
                       xlab="Time",
                       break.time.by = ceiling((max(rt$Time))/4),
                       risk.table.title="",
                       palette=c("#009E73","red"),
                       risk.table=T,
                       risk.table.height=.25,)
    
    surPlot <- surPlot +
      labs(title = paste(" (HR = ", hr_coef, ")", sep=''))
    
    print(surPlot)
    dev.off()
  }
  
  hr_coef <- data.frame(hr_coef)
  rownames(hr_coef) <- cancertype
  
  if(length(HR_value)==0){
    HR_value <- hr_coef
  }else {
    HR_value <- rbind(HR_value,hr_coef)
  }
  
}
write.csv(HR_value,paste("fenqi/",fenqi,".csv"))



