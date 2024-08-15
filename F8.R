
### single cell pancancer

setwd('D:/R/Pancancer ISCA1/Single-cell sequencing')

library(tidyHeatmap)
library(tidyverse)
library(data.table)

# Read data
cell_data <-  fread(paste("D:/R/Pancancer ISCA1/Single-cell sequencing/","TISCH_ISCA1_heatmap改.csv",sep=''))
cell_data <- as.data.frame(cell_data)

# Assign 0 to NA
cell_data[is.na(cell_data)] <- -0
# Change column names
colnames(cell_data)=c('cell','data','value')

# Load reshape2 package
library(reshape2)

# Use acast function to convert the data frame to a matrix
mat <- acast(cell_data, cell ~ data, value.var = "value")

mat=t(mat)

# Delete columns that are all NA
mat <- mat[, colSums(is.na(mat))!= nrow(mat)]

# Assume your matrix is named mat
library(tidyr)

# Convert the matrix to a data frame
df <- as.data.frame(mat)

# Delete columns whose sum is less than or equal to 0.1
df <- subset(df, select = colSums(df) > 0.1)

# Add row name column
df$rowname <- rownames(df)

# Use colSums and logical operation to calculate the number of non-zero elements in each column
non_zero_counts <- colSums(df!= 0)

# Add the result to the last row of the data frame
df_with_counts <- rbind(df, non_zero_counts)

df_with_counts <- df_with_counts[,-30]

# Use gather function to convert the matrix to a data frame with three columns
df_long <- gather(df, key = "property", value = "value", -rowname)

# To be consistent with your example, you may need to rename the columns
df_long <- df_long %>%
  rename(`cell name` = rowname)

cell_data <- df_long

# Rename columns
colnames(cell_data)=c('data','cellname','expression')

# Add a column, this column takes the value before the first _ in the first column.
cell_data$cancer <- sapply(cell_data[, 1], function(x) unlist(strsplit(x, "_"))[1])

# Delete cells with low expression
cell_data <- cell_data[!(cell_data$cellname %in% c("ILC", "Ductal", "Acinar")), ]

# Assume your data frame is named df
library(tibble)
# Convert the data frame to a tibble
df_tibble <- as_tibble(cell_data)

# First reorder the data frame by the cancer column
df_tibble <- df_tibble[order(df_tibble$cancer), ]

# Draw the heatmap
mtcars_heatmap <- df_tibble %>% 
  heatmap(
    .column = cellname, 
    .row = data, 
    .value = expression,
    cluster_rows = FALSE,
    col = c("white","red"),
    rect_gp = grid::gpar(col = "#161616", lwd = 0.3)
  )%>%
  annotation_tile(cancer)

pdf(paste("Single-cell.pdf",sep=''),width =10,height = 12)

mtcars_heatmap

dev.off()