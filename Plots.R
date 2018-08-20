#Plots
library(SummarizedExperiment)
library(DESeq2)
library(Biobase)
library(DEFormats)
library(DESeq2)
library(magrittr)
library(RColorBrewer)
library(pheatmap)
library(ggplot2)
library(ggrepel)
library(dplyr)
library(BiocGenerics)

base_csv_Tcell <- read.table( "~/Desktop/Tcell_only/Tables/Th1/Tcell_master_p1.txt", header = TRUE)

#PCA - to do beofre replicates are combined
#for PCA

setwd("~/Desktop/Tcell_only/24hr_files/")
rdl <-read.csv("All_24_timepts_fpkms.csv", row.names=NULL)
View(rdl)

#if I want genes names as row.names
rdl2 <- rdl[,-1]
rownames(rdl2) <- rdl[,1]
rdl <- rdl2

#for PCA MDS
SummarizedExperiment()

DESeqDataSetFromMatrix(rdl)                   
colData(rdl)
countdata <- assay(counts_as_RSE)
head(countdata,4)
coldata <- colData(counts_as_RSE)
head(coldata)

condition=as.factor(c(rep("24hr-in",3),rep("0hr",4),rep("2wk",4),rep("72hr",4)))
mycols = data.frame(row.names = TRUE, condition)

dds <- DESeqDataSetFromMatrix(countData = countdata, colData = mycols, design = ~condition)

rld <- rlog(dds, blind = FALSE)

plotPCA(rdl)
rdl.sub.exon <- rdl[,c(2,5:10)]
rdl.sub.intron <- rdl[,c(3,11:16)]
head(rdl.sub.exon)

plotPCA(rdl.sub.exon) + geom_label_repel(aes(label =colnames(rdl.sub.exon)), size = 2.5, label.size = 0, fill = "white") + 
  ggtitle("Exon PCA with labels")
plotPCA(rdl.sub.intron) + geom_label_repel(aes(label =colnames(rdl.sub.intron)), size = 2.5, label.size = 0, fill = "white") + 
  ggtitle("Intron PCA with labels")

#Avg replicates

#if I want genes names as row.names
bsc2 <- base_csv_Tcell[,-1]
rownames(bsc2) <- base_csv_Tcell[,1]
base_csv_Tcell <- bsc2 

base_csv_Tcell$Mean_0hr_exon_fpkm <- rowMeans(base_csv_Tcell[,1:2])
base_csv_Tcell$Mean_2w_exon_fpkm <- rowMeans(base_csv_Tcell[,3:4])
base_csv_Tcell$Mean_72hr_exon_fpkm <- rowMeans(base_csv_Tcell[,5:6])
base_csv_Tcell$Mean_0hr_intron_fpkm <- rowMeans(base_csv_Tcell[,7:8])
base_csv_Tcell$Mean_2w_intron_fpkm <- rowMeans(base_csv_Tcell[,9:10])
base_csv_Tcell$Mean_72hr_intron_fpkm <- rowMeans(base_csv_Tcell[,11:12])

base_csv_Tcell$Mean_0hr_ratio <- rowMeans(base_csv_Tcell[,13:14])
base_csv_Tcell$Mean_2w_ratio <- rowMeans(base_csv_Tcell[,15:16])
base_csv_Tcell$Mean_72hr_ratio <- rowMeans(base_csv_Tcell[,17:18])


mean_base_csv_Tcell <- (base_csv_Tcell[19:27])

log_mean_basecsv_Tcell <- log2(mean_base_csv_Tcell) 
log_mean_basecsv_Tcell <- round(log_mean_basecsv_Tcell, 3)
Genes <- row.names(base_csv_Tcell)

SE_genes <- SE_loost[,1]

#call in SE list
SE_list2 <-read.csv("~/Desktop/Tcell_only/SE_lists/Th1_SE_list.csv")
SE_loost <- SE_list$X.2
SE_loost <- as.data.frame(SE_loost)

library(plyr)
lmb <- log_mean_basecsv_Tcell[,c(11,1:10)]
merged_SEeeee<- merge(lmb, SE_loost,  by = 1, all = F)


#df for SE data
log_mean_basecsv_Tcell$SE_loost <- rownames(log_mean_basecsv_Tcell)
SE_list_final <- merge(SE_list, log_mean_basecsv_Tcell, by = "SE_loost", all= F)
write.csv(SE_list_final, "~/Desktop/Tcell_only/SE_list_plus_my_values.csv")


#plot SE w/ expression data
SExperessionGenes <- SE_list_final$SE_loost
ggplot(SE_list_final, aes(Mean_72hr_ratio, Mean_72hr_exon_fpkm, label=SExperessionGenes)) + geom_point(alpha=0.5, color = "steelblue4") + 
  geom_text_repel(aes(label=ifelse(Mean_72hr_exon_fpkm>5.7,as.character(SExperessionGenes),'')),hjust=0,vjust=0) + 
  ggtitle("SE 72 hr exon/intron ratio to exon fpkm", subtitle = "applied min exon cutoff 2, log transformed, replicates combined")

  #se interactive
library(plotly)

se_p <- plot_ly(data = SE_list_final, x = SE_list_final$Mean_72hr_ratio, y = SE_list_final$Mean_72hr_exon_fpkm, 
             alpha = 0.5, symbol = SExperessionGenes)

#with log transform
ggplot(log_mean_basecsv_Tcell, aes(Mean_72hr_ratio, Mean_72hr_exon_fpkm, label=Genes)) + geom_point(alpha=0.5) + 
  geom_text_repel(aes(label=ifelse(Mean_72hr_exon_fpkm>9.5,as.character(Genes),'')),hjust=0,vjust=0) + 
  ggtitle("T cell 72 hr exon/intron ratio to exon fpkm", subtitle = "applied min exon cutoff 2, log transformed, replicates combined")

#with highlights
#5 ##THIS ONE!!!!!
ggplot(data = log_mean_basecsv_Tcell, aes(log_mean_basecsv_Tcell$Mean_72hr_ratio, log_mean_basecsv_Tcell$Mean_72hr_exon_fpkm)) +  
  geom_point(alpha=0.5, colour="darkgrey") +
  geom_point(data=merged_SEeeee, aes(merged_SEeeee$Mean_72hr_ratio,merged_SEeeee$Mean_72hr_exon_fpkm), alpha=0.5, colour="blue") + 
  ggtitle("T cell 72 hr exon/intron ratio to exon fpkm - SE genes blue", subtitle = "applied min exon cutoff 2, log transformed, replicates combined")

#Boxplot
sorted_merged_SE_bploot <- order(merged_SE_boxploot, na.last = TRUE)  
sorted_merged_SE_bploot <- merged_SE_boxploot[order(merged_SE_boxploot$`CORE ENRICHMENT`),]  

final_SE_w_vals <- sorted_merged_SE_bploot[complete.cases(sorted_merged_SE_bploot), ]
#a simpler check??
ggplot(data.frame(x=df1$x, y=df2$x), aes(x,y)) + geom_point(alpha=0.5, color="cadetblue4")


#interactive plot..crashes comp in browser, very large
library(plotly)

p <- plot_ly(data = log_mean_basecsv_Tcell, x = log_mean_basecsv_Tcell$Mean_72hr_ratio, y = log_mean_basecsv_Tcell$Mean_72hr_exon_fpkm, 
             alpha = 0.5, symbol = row.names(log_mean_basecsv_Tcell))

#exon vs exon FC
#first add columnn w/ FC eq
#aa$z <- with(aa, x + y - 2) OR
#need to make x-axis = FC = ((B-A)/A) = (72-0)/0) = ((avg exon fpkm 72 - avg exon fpk 0) / (avg exon fpkm 0)...will need to create a variable..just new column for this..shouldn't be too hard... 

log_mean_basecsv_Tcell$FC_72hr_exon <- ((log_mean_basecsv_Tcell$Mean_72hr_exon_fpkm - log_mean_basecsv_Tcell$Mean_0hr_exon_fpkm) / log_mean_basecsv_Tcell$Mean_0hr_exon_fpkm)

ggplot(data = log_mean_basecsv_Tcell, aes(log_mean_basecsv_Tcell$Mean_72hr_exon_fpkm, log_mean_basecsv_Tcell$FC_72hr_exon)) +  
  geom_point(alpha=0.5, colour="darkgrey") +
  geom_point(data=merged_SEeeee, aes(merged_SEeeee$Mean_72hr_exon_fpkm,merged_SEeeee$FC_72hr_exon), alpha=0.5, colour="blue") + 
  ggtitle("T cell 72 hr exon fpkm to exon FC - SE genes blue", subtitle = "applied min exon cutoff 2, log transformed, replicates combined")

  #se FC interactive
library(plotly)

se_p_exonFC <- plot_ly(data = merged_SEeeee, x = merged_SEeeee$Mean_72hr_exon_fpkm, y = merged_SEeeee$FC_72hr_exon, 
                alpha = 0.5, symbol = merged_SEeeee$SE_loost)
#FOR 2WKS
#need to do FC calculation:

log_mean_basecsv_Tcell$FC_2wk_exon <- ((log_mean_basecsv_Tcell$Mean_2w_exon_fpkm - log_mean_basecsv_Tcell$Mean_0hr_exon_fpkm) / log_mean_basecsv_Tcell$Mean_0hr_exon_fpkm)
merged_SEeeee$FC_2wk_exon <- ((merged_SEeeee$Mean_2w_exon_fpkm - merged_SEeeee$Mean_0hr_exon_fpkm) / merged_SEeeee$Mean_0hr_exon_fpkm)

ggplot(data = log_mean_basecsv_Tcell, aes(log_mean_basecsv_Tcell$Mean_2w_exon_fpkm, log_mean_basecsv_Tcell$FC_2wk_exon)) +  
  geom_point(alpha=0.5, colour="darkgrey") +
  geom_point(data=merged_SEeeee, aes(merged_SEeeee$Mean_2w_exon_fpkm,merged_SEeeee$FC_2wk_exon), alpha=0.5, colour="blue") + 
  ggtitle("T cell 2wk exon fpkm to exon FC - SE genes blue", subtitle = "applied min exon cutoff 2, log transformed, replicates combined")

  #se FC interactive
library(plotly)

se_p_exonFC_2wks <- plot_ly(data = merged_SEeeee, x = merged_SEeeee$Mean_2w_exon_fpkm, y = merged_SEeeee$FC_2wk_exon, 
                       alpha = 0.5, symbol = merged_SEeeee$SE_loost)

#boxplot 
library(readxl)
SEs <- read_xlsx("~/Desktop/Tcell_only/SE_enriched_forR.xlsx")
avg_exon_fpkms <- read_xlsx("~/Desktop/Tcell_only/avg_exon_fpkm_forR.xlsx")

merged_SE_boxploot <- merge(SEs, avg_exon_fpkms, by = 1, all = TRUE)

sorted_merged_SE_bploot <- order(merged_SE_boxploot, na.last = TRUE)  
sorted_merged_SE_bploot <- merged_SE_boxploot[order(merged_SE_boxploot$`CORE ENRICHMENT`),]  

final_SE_w_vals <- sorted_merged_SE_bploot[complete.cases(sorted_merged_SE_bploot), ]
 
log272 <- log2(final_SE_w_vals$`72h_r-avg_exon_fpkm`) 
log272 <- as.data.frame(log272)
final_SE_w_vals <- cbind(final_SE_w_vals, log272)

ggplot(final_SE_w_vals, aes(x = final_SE_w_vals$`CORE ENRICHMENT`, y = final_SE_w_vals$log272)) + geom_boxplot() + 
  ggtitle("Boxplot of SE genes", subtitle = "log transformed, rep combined 72hr fpkm")

FC_fpkm <- read_xlsx("~/Desktop/Tcell_only/FC_fpkm_7r_forR.xlsx")
FC_ratio <- read_xlsx("~/Desktop/Tcell_only/FC_ratio_72_forR.xlsx")
all_ratio <- read_xlsx("~/Desktop/Tcell_only/ratio_enriched_forR.xlsx")

mergedagain <- merge(final_SE_w_vals, FC_fpkm, by = 1, all = TRUE)
mergedagain <- merge(mergedagain, all_ratio, by = 1, all = TRUE)
mergedagain <- merge(mergedagain, FC_ratio, by = 1, all = TRUE)

mergedagainn <- mergedagain[order(mergedagain$`CORE ENRICHMENT`),]  
mergedagainn <- as.data.frame(mergedagainn)

log272_ratio <- log2(mergedagainn$)
