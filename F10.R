
### Volcano-GO-KEGG

setwd("D:/R/THCA ISCA1/difference gene1/limma/")

library(data.table)
gene <-  fread(paste("D:/R/Pancancer ISCA1/linear/","ferroptosis_driver.csv",sep=''))
gene <- as.data.frame(gene)
allgene <- gene$symbol
allgene <- unique(allgene)
gene1 <-  fread(paste("D:/R/Pancancer ISCA1/linear/","ferroptosis_marker.csv",sep=''))
gene1 <- as.data.frame(gene1)
allgene1 <- gene1$symbol
allgene1 <- unique(allgene1)
gene2 <-  fread(paste("D:/R/Pancancer ISCA1/linear/","ferroptosis_suppressor.csv",sep=''))
gene2 <- as.data.frame(gene2)
allgene2 <- gene2$symbol
allgene2 <- unique(allgene2)
gene3 <-  fread(paste("D:/R/Pancancer ISCA1/linear/","ferroptosis_unclassified.csv",sep=''))
gene3 <- as.data.frame(gene3)
allgene3 <- gene3$symbol
geneall <- unique(c(allgene,allgene1,allgene2,allgene3))
geneall <- as.character(geneall)
data <- read.table("D:/R/THCA ISCA1/difference gene1/hallmark/THCA-count.txt",row.names = 1)
colnames(data) <- data[1,]
data <- data[-1,]
data <- data[complete.cases(data[,1]), ]

# Normal and tumor numbers, 14th, 15th characters, 01-09 is cancer, 10-19 is normal, 20-29 is adjacent to cancer
group=sapply(strsplit(colnames(data),"\\-"), "[", 4)
group=sapply(strsplit(group,""), "[", 1)
# Only retain tumor samples
data = data[,group == 0]

cluster <- read.csv("D:/R/THCA ISCA1/difference gene1/KM divide into high and low/HIGH-LOW.csv",header = F)
cluster <-cluster[-1,-1]
colnames(cluster)=c("V1","V2")

c1 <- cluster[which(cluster$V2=="Low"),]
c2 <- cluster[which(cluster$V2=="High"),]

data_c1 <- as.matrix(data)[,na.omit(match(c1[,1],colnames(data)))]
data_c2 <- as.matrix(data)[,na.omit(match(c2[,1],colnames(data)))]
data_merge <- cbind(data_c1,data_c2)
data_merge=apply(data_merge,2,as.numeric)
rownames(data_merge) <- rownames(data)
library(limma)
library(edgeR)
data_merge <- data_merge[complete.cases(data_merge[,2]), ]

d0 <- DGEList(data_merge)
d0 <- calcNormFactors(d0, method="TMM")
d <- d0
group_list <- factor(c(rep("c1",dim(data_c1)[2]),
                       rep("c2",dim(data_c2)[2])), levels = c("c1","c2"))
design <- model.matrix(~group_list)
colnames(design) <- levels(group_list)
rownames(design) <- colnames(data)
v <- voom(d, design, plot=F)
fit <- lmFit(v, design)
fit <- eBayes(fit, trend=F)
output <- topTable(fit, coef=2,n=Inf)
res<-subset(output, adj.P.Val<0.05 & abs(logFC) >0.25 )
write.csv(res,file= "DEG_voom-limma.csv")
library(ggpubr)
library(ggthemes)
library(ggrepel)
deg.data <- output
deg.data <-as.data.frame(deg.data)
deg.data$log10P <-  (-log10(deg.data$adj.P.Val))
deg.data$Symbol <- rownames(deg.data)
deg.data$Group <- "not-significant"
deg.data$Group[which((deg.data$adj.P.Val<0.05) & (deg.data$logFC > 0.5) )]="up-Regulated"
deg.data$Group[which((deg.data$adj.P.Val<0.05) & (deg.data$logFC < -0.5) )]="down-Regulated"
table(deg.data$Group)
deg.data$Label=""
deg.data<-deg.data[order(deg.data$adj.P.Val),]

degallgenes <- names(which(table(c(geneall,rownames(deg.data)))==2))
deg.data <- deg.data[match(degallgenes,rownames(deg.data)),]
up.genes<-head(deg.data$Symbol[which(deg.data$Group=="up-Regulated")],10)
down.genes<-head(deg.data$Symbol[which(deg.data$Group=="down-Regulated")],10)
deg.top10.genes<-c(as.character(up.genes),as.character(down.genes))
deg.data$Label[match(deg.top10.genes,deg.data$Symbol)]<-deg.top10.genes
write.csv(deg.data,file= "THCA_DEG_voom-limma.csv")

library(clusterProfiler)
library(org.Hs.eg.db)

gene_name<-deg.data$Symbol[which(deg.data$Group=="up-Regulated")]
write.csv(gene_name,file="THCA_DEGs_up_Fe.csv")
gene_name<-deg.data$Symbol[which(deg.data$Group=="down-Regulated")]
write.csv(gene_name,file="THCA_DEGs_down_Fe.csv")

pdf("volcano-THCA-Fe-count final 1.pdf", width = 10)
p <- ggscatter(deg.data, x = "logFC", y = "log10P",
               color = "Group",
               palette = c("#2f5688", "#BBBBBB", "#CC0000"),
               alpha = 0.8, size = 2.5,
               label = deg.data$Label,
               font.label = 8, repel = TRUE,
               xlab = "logFC", ylab = "-log10(padj)") +
  theme_base() +
  geom_hline(yintercept = 1.30, linetype = "dashed") +
  geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed") +
  xlim(c(-3, 3)) +
  labs(x = "logFC", y = "-log10(padj)") +
  scale_x_continuous(breaks = c(-3, -2, -1,-0.5,0,0.5, 1, 2, 3))  # Add 0.5 scale line
p
dev.off()
p



library(clusterProfiler)
library(org.Hs.eg.db)
gene_name<-deg.data$Symbol[which(deg.data$Group=="up-Regulated")]
write.csv(gene_name,file="DEGs_up.csv")
gene.df <- bitr(gene_name, fromType = "SYMBOL",
                toType = c("ENSEMBL","ENTREZID"),
                OrgDb = org.Hs.eg.db)
ego_plot<-enrichGO(gene=gene.df$ENTREZID,OrgDb=org.Hs.eg.db,ont = "all",
                   pAdjustMethod ="BH", pvalueCutoff =0.05,qvalueCutoff=1,
                   readable= TRUE)#"none"
write.csv(ego_plot,file="GO_up_plot.csv")
pdf("GO_up-dotplot.pdf",height = 15,width=8)
dotplot(ego_plot, split="ONTOLOGY")+ facet_grid(ONTOLOGY~.,scale="free")
dev.off()
#BP#CC#MF
ego_plot<-enrichGO(gene=gene.df$ENTREZID,OrgDb=org.Hs.eg.db,ont = "BP",
                   pAdjustMethod ="BH", pvalueCutoff =0.05,qvalueCutoff=1,
                   readable= TRUE)#"none"
write.csv(ego_plot,file="GO_up_BP_plot.csv")
pdf("GO_up_BP-dotplot.pdf",height = 7,width=8)
dotplot(ego_plot)
dev.off()
#CC#MF
ego_plot<-enrichGO(gene=gene.df$ENTREZID,OrgDb=org.Hs.eg.db,ont = "CC",
                   pAdjustMethod ="BH", pvalueCutoff =0.05,qvalueCutoff=1,
                   readable= TRUE)#"none"
write.csv(ego_plot,file="GO_up_CC_dotplot.csv")
pdf("GO_up_CC-dotplot.pdf",height = 7,width=8)
dotplot(ego_plot)
dev.off()
#MF
ego_plot<-enrichGO(gene=gene.df$ENTREZID,OrgDb=org.Hs.eg.db,ont = "MF",
                   pAdjustMethod ="BH", pvalueCutoff =0.05,qvalueCutoff=1,
                   readable= TRUE)#"none"
write.csv(ego_plot,file="GO_up_MF_dotplot.csv")
pdf("GO_up_MF-dotplot.pdf",height = 7,width=8)
dotplot(ego_plot)
dev.off()

ekk_plot<-enrichKEGG(gene=gene.df$ENTREZID,organism="hsa",
                     keyType = "kegg",pAdjustMethod="none",pvalueCutoff =0.05,qvalueCutoff=1)#"none"

ekk2 <- setReadable(ekk_plot, 'org.Hs.eg.db', 'ENTREZID')
write.csv(ekk2,file="KEGG_up.csv")
pdf('KEGG_up_dotplot.pdf',height = 6,width=8)
dotplot(ekk2)
dev.off()



gene_name<-deg.data$Symbol[which(deg.data$Group=="down-Regulated")]
write.csv(gene_name,file="DEGs_down.csv")
gene.df <- bitr(gene_name, fromType = "SYMBOL",
                toType = c("ENSEMBL","ENTREZID"),
                OrgDb = org.Hs.eg.db)
ego_plot<-enrichGO(gene=gene.df$ENTREZID,OrgDb=org.Hs.eg.db,ont = "all",
                   pAdjustMethod ="BH", pvalueCutoff =0.05,qvalueCutoff=1,
                   readable= TRUE)#"none"
write.csv(ego_plot,file="GO_down_plot.csv")
pdf("GO_down-dotplot.pdf",height = 15,width=8)
dotplot(ego_plot, split="ONTOLOGY")+ facet_grid(ONTOLOGY~.,scale="free")
dev.off()
#BP#CC#MF
ego_plot<-enrichGO(gene=gene.df$ENTREZID,OrgDb=org.Hs.eg.db,ont = "BP",
                   pAdjustMethod ="BH", pvalueCutoff =0.05,qvalueCutoff=1,
                   readable= TRUE)#"none"
write.csv(ego_plot,file="GO_down_BP_plot.csv")
pdf("GO_down_BP-dotplot.pdf",height = 7,width=8)
dotplot(ego_plot)
dev.off()
#CC#MF
ego_plot<-enrichGO(gene=gene.df$ENTREZID,OrgDb=org.Hs.eg.db,ont = "CC",
                   pAdjustMethod ="BH", pvalueCutoff =0.05,qvalueCutoff=1,
                   readable= TRUE)#"none"
write.csv(ego_plot,file="GO_down_CC_plot.csv")
pdf("GO_down_CC-dotplot.pdf",height = 7,width=8)
dotplot(ego_plot)
dev.off()
#MF
ego_plot<-enrichGO(gene=gene.df$ENTREZID,OrgDb=org.Hs.eg.db,ont = "MF",
                   pAdjustMethod ="BH", pvalueCutoff =0.05,qvalueCutoff=1,
                   readable= TRUE)#"none"
write.csv(ego_plot,file="GO_down_MF_plot.csv")
pdf("GO_down_MF-dotplot.pdf",height = 7,width=8)
dotplot(ego_plot)
dev.off()
ekk_plot<-enrichKEGG(gene=gene.df$ENTREZID,organism="hsa",
                     keyType = "kegg",pAdjustMethod="none",pvalueCutoff =0.05,qvalueCutoff=1)#"none"

ekk2 <- setReadable(ekk_plot, 'org.Hs.eg.db', 'ENTREZID')
write.csv(ekk2,file="KEGG_down.csv")
pdf('KEGG_down_dotplot.pdf',height = 6,width=8)
dotplot(ekk2)
dev.off()


### Hallmark



setwd("D:/R/THCA ISCA1/difference gene1/hallmark/")

library(GSVA)
library(GSEABase)
library(limma)
library(Seurat)
library(msigdbr)

human <- msigdbr(species = "Homo sapiens")
human_GO_bp = msigdbr(species = "human",
                      category = "H") %>% 
  dplyr::select(gs_name,gene_symbol)
human_GO_bp_Set = human_GO_bp %>% split(x =.$gene_symbol, f =.$gs_name)
s.sets <- human_GO_bp_Set
library(data.table)
data <- read.table("D:/R/THCA ISCA1/difference gene1/hallmark/THCA-count.txt",row.names = 1)

colnames(data) <- data[1,]
data <- data[-1,]

library(data.table)

group=sapply(strsplit(colnames(data),"\\-"), "[", 4)
group=sapply(strsplit(group,""), "[", 1)
data = data[,group == 0]

cluster <- read.csv("D:/R/THCA ISCA1/difference gene1/KM divide into high and low/HIGH-LOW.csv",header = F)
cluster <-cluster[-1,-1]
colnames(cluster)=c("V1","V2")

c1 <- cluster[which(cluster$V2=="Low"),]
c2 <- cluster[which(cluster$V2=="High"),]

data_c1 <- as.matrix(data)[,na.omit(match(c1[,1],colnames(data)))]
data_c2 <- as.matrix(data)[,na.omit(match(c2[,1],colnames(data)))]

data_merge <- cbind(data_c1,data_c2)
data_merge=apply(data_merge,2,as.numeric)
rownames(data_merge) <- rownames(data)
expr1 <- as.matrix(data_merge)
expr <- expr1

es.matrix <- gsva(
  expr,
  s.sets,
  min.sz = 10,
  max.sz = Inf,
  tau = 1,
  method = "gsva",
  abs.ranking = FALSE,
  verbose = TRUE,
  parallel.sz = 1
)

saveRDS(es.matrix,file="hallmark.rds")
n1 <- 1:dim(data_c1)[2]
n2 <- dim(data_c1)[2]+1:dim(data_c2)[2]
es.matrix.1 <-
  as.data.frame(es.matrix[, n1],
                row.names = row.names(es.matrix))
es.matrix.2 <-
  as.data.frame(es.matrix[, n2],
                row.names = row.names(es.matrix))
es.matrix.f <- cbind(es.matrix.1, es.matrix.2)
grouP <-
  c(rep("case", dim(es.matrix.1)[2]),
    rep("control", dim(es.matrix.2)[2]))

grouP <- as.factor( grouP)
design <- model.matrix(~ grouP + 0)
row.names(design)<-c(colnames(es.matrix.1), colnames(es.matrix.2))
comparE <-
  makeContrasts(grouPcase - grouPcontrol, levels = design)
fit <- lmFit(es.matrix, design)
fit2 <- contrasts.fit(fit, comparE)
fit3 <- eBayes(fit2)
diff <- topTable(fit3, coef = 1, number = dim(es.matrix)[1])
t_results <-
  as.data.frame(diff$t, row.names = rownames(es.matrix))
head(t_results)
colnames(t_results) <- c("t_value")
saveRDS(t_results,file="t_results_hallmark.rds")
write.csv(t_results,file="t_results_hallmark.csv")

library(ggplot2)

library(pheatmap)
rownames(t_results) <- gsub("HALLMARK_", "", rownames(t_results))
rownames(t_results) <- gsub("_", " ", rownames(t_results))
focus.cluster <- "t_value"
sub_t_results <- as.data.frame(t_results[, focus.cluster],
                               row.names = rownames(t_results))
sub_t_results$hallmark <- rownames(sub_t_results)
colnames(sub_t_results) <- c("t", "hallmark")
sub_t_results$hallmark = with(sub_t_results, reorder(hallmark, t))
sub_t_results$fill <- ""
sub_t_results[sub_t_results$t >= 2.58,]$fill <-
  "up"
sub_t_results[sub_t_results$t <= -2.58,]$fill <-
  "down"
sub_t_results[abs(sub_t_results$t) < 2.58,]$fill <-
  "no"
sub_t_results$color <- ""
sub_t_results[abs(sub_t_results$t) < 2.58,]$color <-
  "n"
sub_t_results[abs(sub_t_results$t) >= 2.58,]$color <-
  "y"

sub_t_results <- sub_t_results[c(1:50),]
summary(sub_t_results$t)

p <-
  ggplot(sub_t_results, aes(x = hallmark, y = t, fill = fill)) +
  geom_bar(stat = "identity", width =.8) +
  scale_fill_manual(
    values = c(
      "down" = "#36648b",
      "up" = "#e94644",
      "no" = "#cccccc"
    ),
    guide = F
  ) + ylim(-10,12)+
  geom_hline(
    yintercept = c(-2.58, 2.58),
    color = "white",
    linetype = "dotted",
    size = 0.5
  ) +
  coord_flip() +
  xlab("") +
  geom_text(
    data = subset(sub_t_results, t < 0),
    aes(
      x = hallmark,
      y = 0.1,
      label = paste0(" ", hallmark),
      color = color
    ),
    size = 1.8,
    hjust ="inward"
  ) +geom_text(
    data = subset(sub_t_results, t > 0),
    aes(
      x = hallmark,
      y = -0.1,
      label = paste0(" ", hallmark),
      color = color
    ),
    size = 1.8,
    hjust =  "outward"
  ) +
  scale_colour_manual(values = c("y" = "black", "n" = "#cccccc"),
                      guide = FALSE) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.ticks.y = element_blank(),
    plot.title = element_blank(),
    axis.text.y = element_blank(),
    panel.border = element_rect(
      colour = "black",
      fill = NA,
      size = 0.5
    ),
    panel.background = element_blank(),
    axis.text.x = element_text(colour = "black"),
    axis.ticks.x = element_line(colour = "black", size = 0.5),
  )
#p
ggsave(
  filename = "hallmark-THCA-count-3.27.pdf",
  plot = p,
  height = 6.5,
  width =5
)



### Differential analysis of TIPs


setwd("D:/R/Pancancer ISCA1/TIP/TIP(1)/TIP/")

cancer <- "THCA"

library(data.table)
TIP <-  fread(paste("D:/R/Pancancer ISCA1/TIP/data/",cancer,".txt",sep=''))
TIP  <- as.data.frame(TIP)
rownames(TIP) <- TIP[,1]
TIP <- TIP[,-1]
dat_immu <- t(TIP)
aa <- substr(rownames(dat_immu), 1, 15)
rownames(dat_immu) <- aa

rf_risk <- read.csv("D:/R/THCA ISCA1/difference gene1/KM divide into high and low/HIGH-LOW.csv",header = F)
rf_risk <- rf_risk[,c(2,3)]

colnames(rf_risk) <- rf_risk[1,]
rf_risk <- rf_risk[-1,]

cluster <- subset(rf_risk, risk %in% c("High", "Low"))
cluster <- as.data.frame(cluster)
rownames(cluster) <- cluster[,1]
colnames(cluster) <- c("V1","V2") 

all_name <- names(which(table(c(rownames(dat_immu),cluster[,1] ))==2))

dat_cluster <- cluster[match(all_name,cluster[,1]),]
dat_imm <- dat_immu[match(all_name,rownames(dat_immu)),]
dat_im <- as.data.frame(dat_imm)
dat_cluster <- na.omit(dat_cluster)
all_name <- names(which(table(c(rownames(dat_im),dat_cluster[,1] ))==2))
dat_cluster <- dat_cluster[match(all_name,dat_cluster[,1]),]
dat_im <- dat_im[match(all_name,rownames(dat_im)),]
dat_cluster <- dat_cluster[c('V2')]
dat <- data.frame()
for(coln in colnames(dat_im))
{
  for(i in all_name)
  {
    dat <- rbind(dat,c(coln,dat_cluster[match(i,rownames(dat_cluster)),],dat_im[c(coln)][match(i,rownames(dat_im)),]))
  }
}

dat[,3] = as.numeric(dat[,3])
library(ggpubr)
colnames(dat) <- c("Gene","Group","value")
aa <- compare_means(value ~ Group, data = dat, group.by = "Gene",method = "anova")

pdf(paste("cluster_THCA_TIP.pdf",sep=''),width=16,height = 6)
p <- ggboxplot(dat, x = "Gene", y = "value",
               color = "Group", palette = c("#f8ac8c","#2878b5"), 
               add = "median_q1q3",x.text.angle=60) # palette can select the corresponding color scheme according to the journal, such as "npg", etc.
p <- p+xlab("Step")+ylab("Expression of TIP")
p <- p + stat_compare_means(aes(group = Group), label = "p.signif",method = "anova")
p
print(p)
dev.off()



### TIPS - Lollipop plot


rm (list = ls ())

setwd("D:/R/Pancancer ISCA1/TIP/")
library(data.table)
cancer <- "THCA"
high_low <- "Low"
data <- read.csv(paste(high_low,"-",cancer,"-data-all.csv",sep=''))
colnames(data)=c("Cell","cor","pvalue")

p.col = c('gold','pink','orange','LimeGreen','darkgreen')
fcolor = function(x,p.col){
  color = ifelse(x>0.8,p.col[1],ifelse(x>0.6,p.col[2],ifelse(x>0.4,p.col[3],
                                                             ifelse(x>0.2,p.col[4], p.col[5])
  )))
  return(color)
}
p.cex = seq(2.5, 5.5, length=5)
fcex = function(x){
  x=abs(x)
  cex = ifelse(x<0.1,p.cex[1],ifelse(x<0.2,p.cex[2],ifelse(x<0.3,p.cex[3],
                                                           ifelse(x<0.4,p.cex[4],p.cex[5]))))
  return(cex)
}

points.color = fcolor(x=data$pvalue,p.col=p.col)
data$points.color = points.color
points.cex = fcex(x=data$cor)
data$points.cex = points.cex
library(dplyr)  # Load the dplyr package
data <- arrange(data, desc(Cell))
xlim = ceiling(max(abs(data$cor))*10)/10
pdf(file=paste(high_low,"-",cancer,"-EPIC-Lollipop Plot.pdf",sep=''), width=9, height=7)
layout(mat=matrix(c(1,1,1,1,1,0,2,0,3,0),nc=2),width=c(8,2.2),heights=c(1,2,1,2,1))
par(bg="white",las=1,mar=c(5,18,2,4),cex.axis=1.5,cex.lab=2)
plot(1,type="n",xlim=c(-xlim,xlim),ylim=c(0.5,nrow(data)+0.5),xlab="Correlation Coefficient",ylab="",yaxt="n",yaxs="i",axes=F)
rect(par('usr')[1],par('usr')[3],par('usr')[2],par('usr')[4],col="#F5F5F5",border="#F5F5F5")
grid(ny=nrow(data),col="white",lty=1,lwd=2)
segments(x0=data$cor,y0=1:nrow(data),x1=0,y1=1:nrow(data),lwd=4)
points(x=data$cor,y = 1:nrow(data),col = data$points.color,pch=16,cex=data$points.cex)
text(par('usr')[1],1:nrow(data),data$Cell,adj=1,xpd=T,cex=1.5)
pvalue.text=ifelse(data$pvalue<0.001,'<0.001',sprintf("%.03f",data$pvalue))
redcutoff_cor=0
redcutoff_pvalue=0.05
text(par('usr')[2],1:nrow(data),pvalue.text,adj=0,xpd=T,col=ifelse(abs(data$cor)>redcutoff_cor & data$pvalue<redcutoff_pvalue,"red","black"),cex=1.5)
axis(1,tick=F)
par(mar=c(0,4,3,4))
plot(1,type="n",axes=F,xlab="",ylab="")
legend("left",legend=c(0.1,0.2,0.3,0.4,0.5),col="black",pt.cex=p.cex,pch=16,bty="n",cex=2,title="abs(cor)")
par(mar=c(0,6,4,6),cex.axis=1.5,cex.main=2)
barplot(rep(1,5),horiz=T,space=0,border=NA,col=p.col,xaxt="n",yaxt="n",xlab="",ylab="",main="pvalue")
axis(4,at=0:5,c(1,0.8,0.6,0.4,0.2,0),tick=F)
dev.off() 



### High_Low CIBERSORT heatmap analysis

rm (list = ls ())

setwd("D:/R/Pancancer ISCA1/8cell")

library(linkET)
#devtools::install_github ("Hy4m/linkET", force = TRUE)
library(ggplot2)
library(dplyr)
library(GSVA)
library(data.table)
# high
cancer <- "THCA"
algorithm <- c("cibersort")
high_low <- "Low"
CELL <-  fread(paste("8 algorithms individual/TCGA-",cancer,"/im_",algorithm,".csv",sep=''))
CELL <- as.data.frame(CELL)
original_id <-CELL[,1]
modified_id <- substr(original_id, 1, 15)
CELL[,1] <- modified_id
# Find duplicate row data
duplicated_rows <- duplicated(CELL[,1])
# Extract non-duplicated row data
unique_data_im <- CELL[!duplicated_rows, ]
CELL <-unique_data_im
CELL <- CELL[, -c((ncol(CELL)-2):ncol(CELL))]
# Assume the data frame name is df
# Get column names
col_names <- colnames(CELL)
# Process column names: remove the last _ and the content after it, replace all _ and - with spaces
col_names_new <- gsub("[_-]", " ", col_names)
# Update the column names of the data frame
colnames(CELL) <- col_names_new
rownames(CELL) <- CELL[,1]
data_im <-CELL[,-1]
CIBER=data_im
cluster <- read.csv(paste("D:/R/",cancer," ISCA1/difference gene1/KM divide into high and low/HIGH-LOW.csv",sep=''),header = F)
cluster <-cluster[-1,-1]
cluster <- data.frame(cluster)
rownames(cluster) <- cluster$V2
colnames(cluster)=c("V1","V2")
cluster <- cluster[which(cluster$V2==high_low),]
all_name <- names(which(table(c(rownames(cluster),rownames(CIBER) ))==2))
selected_rows <- CIBER[as.character(all_name), ] 
CIBER <- selected_rows
riskscore <- fread(paste('D:/R/Pancancer ISCA1/clinical/TCGA new/',cancer,".txt",sep=''), header = T)
riskscore <- as.data.frame(riskscore)
data_gene <- riskscore[match( "ISCA1",riskscore$Tag),]
riskscore <-data_gene
riskscore <-rbind(riskscore,data_gene)
riskscore <-t(riskscore)
riskscore[,2]=rownames(riskscore)
colnames(riskscore)=c("ISCA1","SampleName")
riskscore <-riskscore[-1,]
riskscore <- as.data.frame(riskscore)
CIBER$SampleName <- rownames(CIBER)
out <- merge(CIBER,riskscore,by='SampleName')
riskscore <- data.frame(ISCA1 = out$ISCA1)
rownames(out) <- out$SampleName
riskscore <- as.numeric(riskscore$ISCA1)
out <- subset(out,select=-c(SampleName,ISCA1))
data.corr <- qcorrplot(correlate(out), type = "lower", diag = FALSE)
data_cor <- data.corr$data  # Correlation coefficient
write.csv(data_cor,file= paste0(high_low,"_cell.CSV"),quote=F)
mantel2 <- mantel_test(riskscore, out,
                       spec_select = list(riskscore = 1:1)) %>% 
  mutate(rd = cut(r, breaks = c(-Inf, 0.1, 0.3, Inf),
                  labels = c("< 0.1", "0.1 - 0.3", ">= 0.3")),#Divide the correlation coefficient for mapping size
         pd = cut(p, breaks = c(-Inf, 0.01, 0.05, Inf),
                  labels = c("< 0.01", "0.01 - 0.05", ">= 0.05")))#Divide the P value for mapping color
mantel <- mantel2
write.csv(mantel,paste0(high_low,"_CIBERSORT.CSV"))
pdf(paste("linkET_",cancer,"_",algorithm,"_",high_low,".pdf",sep=''),width=length(colnames(CIBER))/2.6,height=length(colnames(CIBER))/2.6)
qcorrplot(correlate(out), type = "lower", diag = FALSE) +#Heatmap drawing
  geom_square() +#Heatmap drawing
  geom_couple(aes(colour = pd, size = rd),data = mantel,curvature = nice_curvature()) +#aes inside is the line format, data corresponds to the mantel test calculation result, curvature controls the line curvature
  scale_fill_gradientn(colours = rev(RColorBrewer::brewer.pal(11, "RdBu")), limits = c(-0.5, 0.5)) +
  scale_size_manual(values = c(0.5, 1, 2)) +
  scale_colour_manual(values = color_pal(3)) +
  guides(size = guide_legend(title = "r",##guides() function adjusts label style
                             override.aes = list(colour = "grey35"), 
                             order = 2),
         colour = guide_legend(title = "p", 
                               override.aes = list(size = 3), 
                               order = 1),
         fill = guide_colorbar(title = "R", order = 3))
dev.off()







setwd("D:/R/Pancancer ISCA1/8cell")

library(linkET)
#devtools::install_github ("Hy4m/linkET", force = TRUE)
library(ggplot2)
library(dplyr)
library(GSVA)
library(data.table)
# high
algorithms <- c("epic","mcpcounter","quantiseq","timer","xcell")
algorithm <- "xcell"
high_low <- "High"
CELL <-  fread(paste("8 algorithms individual/TCGA-HNSC/im_",algorithm,".csv",sep=''))
CELL <- as.data.frame(CELL)
original_id <-CELL[,1]
modified_id <- substr(original_id, 1, 15)
CELL[,1] <- modified_id
# Find duplicate row data
duplicated_rows <- duplicated(CELL[,1])
# Extract non-duplicated row data
unique_data_im <- CELL[!duplicated_rows, ]
CELL <-unique_data_im
CELL <- CELL[, -c((ncol(CELL)-2):ncol(CELL))]
# Assume the data frame name is df
# Get column names
col_names <- colnames(CELL)
# Process column names: remove the last _ and the content after it, replace all _ and - with spaces
col_names_new <- gsub("[_-]", " ", col_names)
# Update the column names of the data frame
colnames(CELL) <- col_names_new
rownames(CELL) <- CELL[,1]
data_im <-CELL[,-1]
CIBER=data_im
# Read the all-high-low table
cluster <- read.csv("D:/R/HNSC ISCA1/difference gene1/KM divide into high and low/HIGH-LOW.csv",header = F)
cluster <-cluster[-1,-1]
cluster <- data.frame(cluster)
rownames(cluster) <- cluster$V2
colnames(cluster)=c("V1","V2")
cluster <- cluster[which(cluster$V2==high_low),]
all_name <- names(which(table(c(rownames(cluster),rownames(CIBER) ))==2))
selected_rows <- CIBER[as.character(all_name), ] 
CIBER <- selected_rows
riskscore <- fread(paste('D:/R/Pancancer ISCA1/clinical/TCGA new/',"HNSC.txt",sep=''), header = T)
riskscore <- as.data.frame(riskscore)
data_gene <- riskscore[match( "ISCA1",riskscore$Tag),]
riskscore <-data_gene
riskscore <-rbind(riskscore,data_gene)
riskscore <-t(riskscore)
riskscore[,2]=rownames(riskscore)
colnames(riskscore)=c("ISCA1","SampleName")
riskscore <-riskscore[-1,]
riskscore <- as.data.frame(riskscore)
CIBER$SampleName <- rownames(CIBER)
out <- merge(CIBER,riskscore,by='SampleName')
riskscore <- data.frame(ISCA1 = out$ISCA1)
rownames(out) <- out$SampleName
riskscore <- as.numeric(riskscore$ISCA1)
out <- subset(out,select=-c(SampleName,ISCA1))
data.corr <- qcorrplot(correlate(out), type = "lower", diag = FALSE)
data_cor <- data.corr$data  # Correlation coefficient
mantel2 <- mantel_test(riskscore, out,
                       spec_select = list(riskscore = 1:1)) %>% 
  mutate(rd = cut(r, breaks = c(-Inf, 0.1, 0.3, Inf),
                  labels = c("< 0.1", "0.1 - 0.3", ">= 0.3")),#Divide the correlation coefficient for mapping size
         pd = cut(p, breaks = c(-Inf, 0.01, 0.05, Inf),
                  labels = c("< 0.01", "0.01 - 0.05", ">= 0.05")))#Divide the P value for mapping color
mantel <- mantel2
pdf(paste("linkET_",algorithm,"_",high_low,".pdf",sep=''),width=length(colnames(CIBER))/2.6,height=length(colnames(CIBER))/2.6)
qcorrplot(correlate(out), type = "lower", diag = FALSE) +#Heatmap drawing
  geom_square() +#Heatmap drawing
  geom_couple(aes(colour = pd, size = rd),data = mantel,curvature = nice_curvature()) +#aes inside is the line format, data corresponds to the mantel test calculation result, curvature controls the line curvature
  scale_fill_gradientn(colours = rev(RColorBrewer::brewer.pal(11, "RdBu")), limits = c(-0.5, 0.5)) +
  scale_size_manual(values = c(0.5, 1, 2)) +
  scale_colour_manual(values = color_pal(3)) +
  guides(size = guide_legend(title = "r",##guides() function Adjusts label style
                             override.aes = list(colour = "grey35"), 
                             order = 2),
         colour = guide_legend(title = "p", 
                               override.aes = list(size = 3), 
                               order = 1),
         fill = guide_colorbar(title = "R", order = 3))
dev.off()