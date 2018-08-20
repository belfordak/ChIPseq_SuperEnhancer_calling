library(readxl)
library(ggplot2)
library(reshape2)
library(RColorBrewer)

setwd("/Users/belfordak/Desktop/ChIP_comp/")

main_df <- read_xlsx("Tcell_SE_genes_forR.xlsx")
#reshape the data
mdata_FC <- main_df[,c(1,2,4,5,7,8,10,11,13,14)]
mdata_fpkm <- main_df[,c(1,3,4,6,7,9,10,12,13,15)]

mdata2_FC <- melt(mdata_FC)
mdata2_fpkm <-melt(mdata_fpkm)

#make the boxplots
FC_boxplot_all_data <- ggplot(mdata2_FC, aes(x=variable, y=value)) + geom_boxplot(aes(color = variable)) + ggtitle("T cell 0-72hr fold change SE epxression") +
  scale_color_brewer(palette="Set2") + theme(axis.text.x=element_text(color = "black", size=8, angle=30, vjust=.8, hjust=0.8)) +
  xlab("SE combination") + ylab("0-72hr Fold Change")

FC_boxplot_without_il2rb <- ggplot(mdata2_FC, aes(x=variable, y=value)) + geom_boxplot(aes(color = variable)) + ggtitle("T cell 0-72hr fold change SE epxression", subtitle = "Il2rb outlier cut off") +
  scale_color_brewer(palette="Set2") + ylim(-1,25) + theme(axis.text.x=element_text(color = "black", size=8, angle=30, vjust=.8, hjust=0.8)) +
  xlab("SE combination") + ylab("0-72hr Fold Change")


fpkm_boxplot_all_data <- ggplot(mdata2_fpkm, aes(x=variable, y=value)) + geom_boxplot(aes(color = variable)) + ggtitle("T cell 72hr fpkm SE epxression") +
  scale_color_brewer(palette="Set1") + theme(axis.text.x=element_text(color = "black", size=8, angle=30, vjust=.8, hjust=0.8)) +
  xlab("SE combination") + ylab("72hr fpkm")

#log transform
logdmdat_FC <- mdata2_FC
logdmdat_FC$value <- log(mdata2_FC$value, 2)

logdmdat_fpkm <- mdata2_fpkm
logdmdat_fpkm$value <- log(mdata2_fpkm$value, 2)

#plot logs
logd_FC_boxplot_all_data <- ggplot(logdmdat_FC, aes(x=variable, y=value)) + geom_boxplot(aes(color = variable)) + ggtitle("T cell 0-72hr log fold change SE epxression") +
  scale_color_brewer(palette="Set2") + theme(axis.text.x=element_text(color = "black", size=8, angle=30, vjust=.8, hjust=0.8)) +
  xlab("SE combination") + ylab("0-72hr Fold Change")

logd_fpkm_boxplot_all_data <- ggplot(logdmdat_fpkm, aes(x=variable, y=value)) + geom_boxplot(aes(color = variable)) + ggtitle("T cell 72hr log fpkm SE epxression") +
  scale_color_brewer(palette="Set1") + theme(axis.text.x=element_text(color = "black", size=8, angle=30, vjust=.8, hjust=0.8)) +
  xlab("SE combination") + ylab("72hr fpkm")

#do statistic between columns, wilcox-ranking test
#for fold change
wilcox.test(x=mdata_FC$`Expression FC_1plus`, y=mdata_FC$`Expression FC_3TF`)
wilcox.test(x=mdata_FC$`Expression FC_1plus`, y=mdata_FC$`Expression FC_3plus`)
wilcox.test(x=mdata_FC$`Expression FC_1plus`, y=mdata_FC$`Expression FC_justH3K27`)
wilcox.test(x=mdata_FC$`Expression FC_1plus`, y=mdata_FC$`Expression FC_97shared`)

wilcox.test(x=mdata_FC$`Expression FC_3TF`, y=mdata_FC$`Expression FC_3plus`)
wilcox.test(x=mdata_FC$`Expression FC_3TF`, y=mdata_FC$`Expression FC_justH3K27`)
wilcox.test(x=mdata_FC$`Expression FC_3TF`, y=mdata_FC$`Expression FC_97shared`)

wilcox.test(x=mdata_FC$`Expression FC_3plus`, y=mdata_FC$`Expression FC_justH3K27`)
wilcox.test(x=mdata_FC$`Expression FC_3plus`, y=mdata_FC$`Expression FC_97shared`)

wilcox.test(x=mdata_FC$`Expression FC_justH3K27`, y=mdata_FC$`Expression FC_97shared`)

#for fpkm
wilcox.test(x=mdata_fpkm$`Exon fpkm_1plus`, y=mdata_fpkm$`Exon fpkm_3TF`)
wilcox.test(x=mdata_fpkm$`Exon fpkm_1plus`, y=mdata_fpkm$`Exon fpkm_3plus`)
wilcox.test(x=mdata_fpkm$`Exon fpkm_1plus`, y=mdata_fpkm$`Exon fpkm_justH3K27`)
wilcox.test(x=mdata_fpkm$`Exon fpkm_1plus`, y=mdata_fpkm$`Exon fpkm_97shared`)

wilcox.test(x=mdata_fpkm$`Exon fpkm_3TF`, y=mdata_fpkm$`Exon fpkm_3plus`)
wilcox.test(x=mdata_fpkm$`Exon fpkm_3TF`, y=mdata_fpkm$`Exon fpkm_justH3K27`)
wilcox.test(x=mdata_fpkm$`Exon fpkm_3TF`, y=mdata_fpkm$`Exon fpkm_97shared`)

wilcox.test(x=mdata_fpkm$`Exon fpkm_3plus`, y=mdata_fpkm$`Exon fpkm_justH3K27`)
wilcox.test(x=mdata_fpkm$`Exon fpkm_3plus`, y=mdata_fpkm$`Exon fpkm_97shared`)

wilcox.test(x=mdata_fpkm$`Exon fpkm_justH3K27`, y=mdata_fpkm$`Exon fpkm_97shared`)

              
#make the pie charts
#instead of subsetting, just figure out counts then manually build plot

#FC 1TF+ac
fc1_orless <- length(which(main_df$`Expression FC_1plus` < 0)) #85
fc0_2 <- length(which(main_df$`Expression FC_1plus` < 2 & main_df$`Expression FC_1plus` > 0)) #279
fc2_5 <- length(which(main_df$`Expression FC_1plus` > 2 & main_df$`Expression FC_1plus` < 5)) #29
fc5_10 <- length(which(main_df$`Expression FC_1plus` > 5 & main_df$`Expression FC_1plus` < 10)) #13
fc10_50 <- length(which(main_df$`Expression FC_1plus` > 10 & main_df$`Expression FC_1plus` < 50)) #6
fc50_plus <- length(which(main_df$`Expression FC_1plus` > 50)) #1

slices_1 <- c(85,279,29,13,6,1)
lbels <- c("1 or less FC ","0-2 FC","2-5 FC","5-10 FC","10-50 FC","50 + FC")
pie_1 <- pie(slices_1, labels = lbels, main = "1TF+H3K27ac Fold Change 0-72hr expression")

#FC 3TF
fc1_orless <- length(which(main_df$`Expression FC_3TF` < 0)) #39
fc0_2FC <- length(which(main_df$`Expression FC_3TF` < 2 & main_df$`Expression FC_3TF` > 0)) #116
fc2_5fc <- length(which(main_df$`Expression FC_3TF` > 2 & main_df$`Expression FC_3TF` < 5)) #19
fc5_10 <- length(which(main_df$`Expression FC_3TF` > 5 & main_df$`Expression FC_3TF` < 10)) #9
fc10_50 <- length(which(main_df$`Expression FC_3TF` > 10 & main_df$`Expression FC_3TF` < 50)) #4
fc50_plus <- length(which(main_df$`Expression FC_3TF` > 50)) #1

slices_2 <- c(39,116,19,9,4,1)
pie_2 <- pie(slices_2, labels = lbels, main = "3TF Fold Change 0-72hr expression")

#FC 3TF+ac
fc1_orless <- length(which(main_df$`Expression FC_3plus` < 0)) #46
fc0_2FC <- length(which(main_df$`Expression FC_3plus` < 2 & main_df$`Expression FC_3plus` > 0)) #167
fc2_5fc <- length(which(main_df$`Expression FC_3plus` > 2 & main_df$`Expression FC_3plus` < 5)) #17
fc5_10 <- length(which(main_df$`Expression FC_3plus` > 5 & main_df$`Expression FC_3plus` < 10)) #10
fc10_50 <- length(which(main_df$`Expression FC_3plus` > 10 & main_df$`Expression FC_3plus` < 50)) #5
fc50_plus <- length(which(main_df$`Expression FC_3plus` > 50)) #1

slices_3 <- c(46,167,17,10,5,1)
pie_3 <- pie(slices_3, labels = lbels, main = "3TF+H3K27ac Fold Change 0-72hr expression")

#FC just ac
fc1_orless <- length(which(main_df$`Expression FC_justH3K27` < 0)) #110
fc0_2FC <- length(which(main_df$`Expression FC_justH3K27` < 2 & main_df$`Expression FC_justH3K27` > 0)) #271
fc2_5fc <- length(which(main_df$`Expression FC_justH3K27` > 2 & main_df$`Expression FC_justH3K27` < 5)) #32
fc5_10 <- length(which(main_df$`Expression FC_justH3K27` > 5 & main_df$`Expression FC_justH3K27` < 10)) #10
fc10_50 <- length(which(main_df$`Expression FC_justH3K27` > 10 & main_df$`Expression FC_justH3K27` < 50)) #6
fc50_plus <- length(which(main_df$`Expression FC_justH3K27` > 50)) #1

slices_4 <- c(110,271,32,10,6,1)
pie_4 <- pie(slices_4, labels = lbels, main = "Just H3K27ac Fold Change 0-72hr expression")

#FC shared vann group
fc1_orless <- length(which(main_df$`Expression FC_97shared` < 0)) #22
fc0_2FC <- length(which(main_df$`Expression FC_97shared` < 2 & main_df$`Expression FC_97shared` > 0)) #57
fc2_5fc <- length(which(main_df$`Expression FC_97shared` > 2 & main_df$`Expression FC_97shared` < 5)) #11
fc5_10 <- length(which(main_df$`Expression FC_97shared` > 5 & main_df$`Expression FC_97shared` < 10)) #3
fc10_50 <- length(which(main_df$`Expression FC_97shared` > 10 & main_df$`Expression FC_97shared` < 50)) #3
fc50_plus <- length(which(main_df$`Expression FC_97shared` > 50)) #1

slices_5 <- c(22,57,11,3,3,1)
pie_5 <- pie(slices_5, labels = lbels, main = "97 Venn Shared Group Fold Change 0-72hr expression")

#for fpkm
fpkm0_10
fpkm10_50
fpkm50_100
fpkm_100_500
fpkm_500_plus


ip <- as.data.frame(installed.packages()[,c(1,3:4)])
rownames(ip) <- NULL

write.csv(ip, "/Users/belfordak/Desktop/1-June-last-transfers/Rpackages_work_comp")
