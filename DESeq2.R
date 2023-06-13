library(tidyverse) #view command active

counts_ <- read.csv("/Users/khadija/Desktop/old/cb/count_bc.csv", sep=";",row.names = 1, stringsAsFactors=FALSE)
view(counts_)
coldata_ <- read.csv("/Users/khadija/Desktop/old/cb/col_bc.csv", sep=";", row.names = 1)
view(coldata_)

de = function(counts, colData){
  
  library(DESeq2)
  library(devtools)
  library(BiocParallel)
  library(dplyr)
  library(data.table)
  library(ggplot2)
  
  register(MulticoreParam())
  
  dds = DESeqDataSetFromMatrix(countData = data.matrix(counts),
                               colData = colData,
                               design = ~ Group)
  
  #dds$Group <- relevel(dds$Group, ref = "N") #control
  keep = rowSums(counts(dds)) >= 10
  dds = dds[keep,]
  
  dds <- DESeq(dds, parallel=TRUE, fitType='local')#perform tasks parallel by allocating cores to R 
  
  ## Plot dispersion estimates
  plotDispEsts(dds)
  
  head(dds)
  
  a = DESeq2::vst(dds, blind = TRUE, fitType = "local")
  print(plotPCA(a, intgroup=c("Group"), returnData=FALSE))
  view(plotPCA(a, intgroup=c("Group"), returnData=TRUE))
  
  #dds = DESeq2::DESeq(dds, parallel=TRUE,fitType = "local")
  return(dds)
  
}


De_ <- de(counts_, coldata_)


##### USE IF FILTERING IS REQUIRED 
de_ <- function(counts, colData){
  
  library(DESeq2)
  library(devtools)
  library(BiocParallel)
  library(dplyr)
  library(data.table)
  library(ggplot2)
  
  register(MulticoreParam())
  
  dds <- DESeqDataSetFromMatrix(countData = data.matrix(counts),
                                colData = colData,
                                design = ~ Group)
  
  keep <- rowSums(counts(dds)) >= 10
  dds <- dds[keep,]
  
  dds <- DESeq(dds, parallel=TRUE, fitType='local')#perform tasks parallel by allocating cores to R 
  
  a <- DESeq2::vst(dds, blind = TRUE, fitType = "local")
  print(plotPCA(a, intgroup=c("Group"), returnData=FALSE))
  
  pcaData <- DESeq2::plotPCA(a, intgroup=c("Group"), returnData=TRUE)
  
  pcaData_L <- as.data.frame(pcaData[pcaData$Group %like% "L", ])
  pcaData_sensitive <- as.data.frame(pcaData[pcaData$Group %like% 'C', ])
  
  pcaData_C <- as.data.frame(subset(pcaData_C, PC1 > 0 & PC2 > 0 ))
  pcaData_L <- as.data.frame(subset(pcaData_L, PC1 < 0 & PC2 < 0))
  
  pcaData <- rbind(pcaData_L, pcaData_C)
  sampls <- counts(dds)
  newData <- sampls[ , rownames(pcaData)]
  coldata <- colData[colnames(newData),]
  
  return(list(newData, coldata))
  
}

new_data <- de_(counts_, coldata_)
view(new_data[[1]])
view(new_data[[2]])
write.csv(new_data[[1]],'/Users/khadija/Desktop/c3/new.csv', row.names = TRUE)
write.csv(new_data[[2]],'/Users/khadija/Desktop/c3/new_col.csv', row.names = TRUE)

de_seq <- function(stage, opt){
  
  design_matrx <- stage[[2]]
  
  register(MulticoreParam())
  
  dds <- DESeq2::DESeqDataSetFromMatrix(countData = as.matrix(stage[[1]]),
                                        colData = stage[[2]],
                                        design = ~ Group)
  
  a <- DESeq2::vst(dds, blind = TRUE, fitType = "local")
  print(DESeq2::plotPCA(a, intgroup=c("Group")))
  
  a <- assay(a)
  keep <- rowSums(counts(dds)) >= 10
  dds <- dds[keep,]
  
  dds <- DESeq2::DESeq(dds, parallel=TRUE)
  
  return(dds)
  
}


DE_de <- de_seq(new_data)

###########

PDAC_cancer <- function(dds){
  
  dds$Group <- relevel(dds$Group, ref = "C") #control
  res_ <- DESeq2::results(dds, contrast=c("Group", "L", "C"), alpha=0.05, parallel=TRUE)
  resOrdered <- res_[which(res_$padj < 0.05), ]
  
  #resultsNames(dds)
  resLFC <- lfcShrink(dds, coef="Group_L_vs_C", type="apeglm")
  plotMA(resLFC, ylim=c(-20,20))
  
  up <- resOrdered[which(resOrdered$log2FoldChange > 1), ]
  dwn <- resOrdered[which(resOrdered$log2FoldChange < -1), ]
  
  resNew <- rbind(up, dwn)
  
  write.csv(up,'/Users/khadija/Desktop/graph/d3_type/up.csv', row.names = TRUE)
  write.csv(dwn,'/Users/khadija/Desktop/graph/d3_type/down.csv', row.names = TRUE)
  
  #return(list(resNew))
  return(list(resNew, up, dwn))
  
}

P = PDAC_cancer(DE_de)
view(P)

up = as.data.frame(P[[2]]) #upregulated genes
down = as.data.frame(P[[3]]) #dwnregulated genes

write.csv(up,'/Users/khadija/Desktop/d1/d3_subtypes/up.csv', row.names = TRUE)
write.csv(dwn,'/Users/khadija/Desktop/d1/d3_subtypes/dwn.csv', row.names = TRUE)

PDAC_C <- function(dds){
  
  #dds$Group <- relevel(dds$Group, ref = "P") #control
  res_ <- DESeq2::results(dds, contrast=c("Group", "C", "L"), alpha=0.05, parallel=TRUE)
  resOrdered <- res_[which(res_$padj < 0.05), ]
  
  
  up <- resOrdered[which(resOrdered$log2FoldChange > 1), ]
  dwn <- resOrdered[which(resOrdered$log2FoldChange < -1), ]
  
  resNew <- rbind(up, dwn)
  
  #return(list(resNew))
  #write.csv(up,'/Users/khadija/Desktop/c3/wg/classical/up.csv', row.names = TRUE)
  #write.csv(dwn,'/Users/khadija/Desktop/c3/wg/classical/dwn.csv', row.names = TRUE)
  view(up)
  return(list(resNew, up, dwn))

}

C = PDAC_C(DE_de)
