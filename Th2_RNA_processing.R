#Th2 and Th17

#1st: General intron/exon processing
library(DESeq2)
library(GenomicAlignments)
library(GenomicFeatures)

#function to get the intervals of Grange object
interval <- function(gr){
  require(GenomicRanges)
  require(GenomicFeatures)
  len.gr <- length(gr)
  if (len.gr <= 1 ) {
    intron <- GRanges()}
  else {
    gap <- gaps(gr)
    intron <- gap[2:length(gap)]
  }
  return(intron)
}

#function to get the only exonic ranges, 
#query is a disjoined Grange object of a gene, exByTxUnlist is the GRange object constains all the exons of transcripts
exonOnlyRange <- function(query, exByTxUnlist){
  exonCounts <- countOverlaps(query, exByTxUnlist)
  return (query[exonCounts == max(exonCounts)])
}

#GrangeList to bedfile
grlToBed <- function(grl){
  gr <- unlist(grl)
  df <- data.frame(seqnames=seqnames(gr),
                   starts=start(gr)-1,
                   ends=end(gr),
                   names=names(gr),
                   scores=c(rep(".", length(gr))),
                   strands=strand(gr))
  
  write.table(df, file=paste0(substitute(grl),".bed"), quote=F, sep="\t", row.names=F, col.names=F)  
}

#GrangeList to GTF
grlToGTF <- function (grl,feature=c("exon","transcript","gene","intron")){
  gr <- unlist(grl)
  gene_id <- paste0('"',names(gr),'";')
  df <- data.frame(seqnames=seqnames(gr),
                   source=c(rep("NIDDK_LGP", length(gr))),
                   feature=feature,
                   starts=start(gr)-1,
                   ends=end(gr),
                   scores=c(rep(".", length(gr))),
                   strands=strand(gr),
                   frame=c(rep(".", length(gr))),
                   attribute=paste("gene_id",gene_id, "transcript_id", gene_id, "gene_name", gene_id, sep = " "))
  write.table(df, file=paste0(substitute(grl),".gtf"), quote=F, sep="\t", row.names = F, col.names = F)
}
#awk -F "\t" '$3 == "exon" { print $9}' Mus_musculus.GRCm38.90.gtf | tr -d ";\"" | awk -F " " '$16 == "protein_coding" {print $2}' | uniq > protein.coding.geneid
#selecting for only "protein coding" elimintes the calling of incorrect genes 
txdb <- makeTxDbFromGFF("/Users/belfordak/Data/Genome_files/UCSC_genes_mm10.gtf")
PC.geneid <- read.table("/Users/belfordak/Data/Genome_files/protein.coding.geneid", sep = '\t', header = F)$V1
library(mygene)
PC.genename <- queryMany(PC.geneid, scopes = "ensembl.gene", fields = "symbol")
PC.symbol <- unique(PC.genename$symbol)

#find out all the protein coding genes
exondb <- exonsBy(txdb, by= "gene", use.names=F) #GRangeList object of genes with overlapping exons
exondb.gr <- unlist(exondb)
exondb.pc <- subset(exondb.gr, names(exondb.gr) %in% PC.symbol)
exondb.pc <- split(exondb.pc, names(exondb.pc))
grlToBed(exondb.pc)

#filter out genes mapped to multiple places
exondb.filtered <- GRangesList()
for (i in 1:length(exondb.pc)) {
  if (length(unique(as.vector(strand(exondb.pc[[i]])))) > 1) next
  if (length(unique(as.vector(seqnames(exondb.pc[[i]])))) > 1) next
  exondb.filtered <- c(exondb.filtered, exondb.pc[i])
}


exonCollapsedb.filtered <- GenomicRanges::reduce(exondb.filtered) #GRangeList object of genes with non-overlapping exons
grlToGTF(exonCollapsedb.filtered,"exon")
exonCollapseTxdb <-makeTxDbFromGFF("exonCollapsedb.filtered.gtf")
introndb.filtered <- intronsByTranscript(exonCollapseTxdb, use.names=T)


names(introndb.filtered) <- paste0("Intron_",names(introndb.filtered))
exonPlusIntrondb.filtered <- c(exonCollapsedb.filtered, introndb.filtered) ##Here's the file for better library normalisation (accounts for both intron and exon, and within the context of the gene, for each file (intron and exon file))

##Here's the point that differs between different samples/conditions. The above code only need be run once (per enviornment)

setwd("/Users/belfordak/Desktop/Tcell_only/bams/Th2/")

bmfls <- list.files(pattern="bam$",full.names=T) #I added in another sample (0hr tt) b/c deseq needs replicates and 24hr doesn't have a rep
bmfls_list <- BamFileList(bmfls)
geneCnt <- summarizeOverlaps(exonPlusIntrondb.filtered,bmfls_list,ignore.strand=T,singleEnd=TRUE,fragments=FALSE)
samples <- colnames(geneCnt)
condition <- c("0h", "0hr", "2w", "2w", "24h","72h", "72h")
colData(geneCnt) <- DataFrame(samples, condition, row.names = samples)
dds <- DESeqDataSet(geneCnt, design = ~condition) 
dds <- DESeq(dds)
fpkm.all <- fpkm(dds)
intron.names <- grep("Intron_", row.names(fpkm.all))
exon.names <- grep("Intron_", row.names(fpkm.all), invert = T)
fpkm.intron <- fpkm.all[intron.names,]
fpkm.exon <- fpkm.all[exon.names,]
rownames(fpkm.intron) <- gsub("Intron_","",rownames(fpkm.intron))
fpkm.transform <- merge(fpkm.exon, fpkm.intron, by=0)
View(fpkm.transform)
names(fpkm.transform) <- c("Gene", "ttN_0hr_r1_exon", "ttN_0hr_r2_exon", "Th2_24_r1_exon","Th2_2w_r1_exon", "Th2_2w_r2_exon","Th2_72h_r1_exon", "Th2_72h_r2_exon",
                           "ttN_0hr_r1_intron", "ttN_0hr_r2_intron","Th2_24_r1_intron","Th2_2w_r1_intron", "Th2_2w_r2_intron","Th2_72h_r1_intron", "Th2_72h_r2_intron")

write.csv(fpkm.transform, file = "~/Desktop/Tcell_only/Tables/Th2/Th2_exon_intron_fpkm.csv", quote = F, row.names = F)

Th2_base <- fpkm.transform
Th2_base <- read.csv("~/Desktop/Tcell_only/Tables/Th2/Th2_exon_intron_fpkm.csv")
#Data cleaning

library(magrittr)
library(dplyr)
library(DESeq2)
library(DEFormats)
library(RColorBrewer)
library(ggplot2)

#Avg replicates

Th2_base$Mean_0hr_exon_fpkm <- rowMeans(Th2_base[,2:3])
Th2_base$Mean_2w_exon_fpkm <- rowMeans(Th2_base[,5:6])
Th2_base$Mean_72hr_exon_fpkm <- rowMeans(Th2_base[,7:8])
Th2_base$Mean_0hr_intron_fpkm <- rowMeans(Th2_base[,9:10])
Th2_base$Mean_2w_intron_fpkm <- rowMeans(Th2_base[,12:13])
Th2_base$Mean_72hr_intron_fpkm <- rowMeans(Th2_base[,14:15])

View(Th2_base)
Th2_base <- mutate(Th2_base, avg_ratio_0hr = Mean_0hr_exon_fpkm / Mean_0hr_intron_fpkm, avg_ratio_2w = Mean_2w_exon_fpkm / Mean_2w_intron_fpkm, avg_ratio_72hr = Mean_72hr_exon_fpkm / Mean_72hr_intron_fpkm)

# Round
Th2_base2 <- Th2_base[,-1]
rownames(Th2_base2) <- Th2_base[,1]

Th2_base2 <-round(Th2_base2[1:23], 3) 

#1. Apply min fpkm value of 2  

trim_test <- Th2_base2[rowSums(Th2_base2[,1:7] < 2) <=1 , , drop = FALSE] #applies min cutoff for atleast one cell of exons to be 2, otherwise delete row
#trim_test2 <- trim_test[apply(trim_test, 1, function(y) !all(is.na(y))),] #would be for romoving only if who row was NAs, this doesn't work for us b/c of the above function
trim_3 <- trim_test[apply(trim_test[,21:23], 1, function(y) !all(is.na(y))),] #removes rows where all ratio counts are NAs

trim_4 <- trim_3[complete.cases(trim_3),] #...just to check                

write.csv(trim_4, "~/Desktop/Tcell_only/Tables/Th2/Actual_Th2_full_avgs_noNA.csv")

#Log transform

log_mean_Th2 <- log2(trim_4) 
log_mean_Th2 <- round(log_mean_Th2, 3) #need to do again

View(log_mean_Th2)

#Plooots!

#can do ratios

#SE ploots!

SE_list_Th2 <- read.csv("~/Desktop/Tcell_only/SE_lists/Th2_SEs.csv")
SE_list_Th2.1 <- as.data.frame(SE_list_Th2[,4])

library(plyr)
library(ggrepel)
library(dplyr)

lmb <- cbind(rownames(log_mean_Th2), data.frame(log_mean_Th2, row.names=NULL))
merged_SEeeee_2 <- merge(lmb, SE_list_Th2.1,  by = 1, all = F)

#ex/int ratio to fpkm, with SEs highlighted and labels
ggplot(data = log_mean_Th2, aes(log_mean_Th2$avg_ratio_72hr, log_mean_Th2$Mean_72hr_exon_fpkm)) +  
  geom_point(alpha=0.5, colour="darkgrey") +
  geom_point(data=merged_SEeeee, aes(merged_SEeeee$avg_ratio_72hr,merged_SEeeee$Mean_72hr_exon_fpkm), alpha=0.5, colour="blue") + 
  ggtitle("Th2 72 hr exon/intron ratio to exon fpkm - SE genes blue", subtitle = "applied min exon cutoff 2, log transformed, replicates combined") +
  geom_text_repel(aes(label=ifelse(log_mean_Th2$avg_ratio_72hr>9.5, as.character(merged_SEeeee$`rownames(log_mean_Th2)`),'')),hjust=0,vjust=0)


#exon vs exon FC
#first add columnn w/ FC eq
#aa$z <- with(aa, x + y - 2) OR
#need to make x-axis = FC = ((B-A)/A) = (72-0)/0) = ((avg exon fpkm 72 - avg exon fpk 0) / (avg exon fpkm 0)...will need to create a variable..just new column for this..shouldn't be too hard... 

log_mean_Th2$FC_72hr_exon <- ((log_mean_Th2$Mean_72hr_exon_fpkm - log_mean_Th2$Mean_0hr_exon_fpkm) / log_mean_Th2$Mean_0hr_exon_fpkm)

ggplot(data = log_mean_Th2, aes(log_mean_Th2$Mean_72hr_exon_fpkm, log_mean_Th2$FC_72hr_exon)) +  
  geom_point(alpha=0.5, colour="darkgrey") +
  geom_point(data=merged_SEeeee, aes(merged_SEeeee$Mean_72hr_exon_fpkm,merged_SEeeee$FC_72hr_exon), alpha=0.5, colour="blue") + 
  ggtitle("Th2 72 hr exon fpkm to exon FC - SE genes blue", subtitle = "applied min exon cutoff 2, log transformed, replicates combined") +
  geom_text_repel(aes(label=ifelse(FC_72hr_exon>5, as.character(merged_SEeeee$`rownames(log_mean_Th2)`),'')),hjust=0,vjust=0)

#above but with other SE lists highlighted too..not as easy as just plopping that SEeee in (that would have the other Th expression appended)
#(SE_list_Th17.1 loaded in enviornment, from Th17_proc.. script, Th1 was not, so shown below)
merged_SE_17_for_Th2 <- merge(lmb, SE_list_Th17.1, by = 1, all = F)

SE_list_Th1 <- read.csv("~/Desktop/Tcell_only/SE_lists/Th1_SE_list.csv")
SE_list_Th1 <- as.data.frame(SE_list_Th1[,4])
merged_SE_1_for_Th2 <- merge(lmb, SE_list_Th1, by = 1, all = F)

library(ggpubr)

ggplot(data = log_mean_Th2, aes(log_mean_Th2$Mean_72hr_exon_fpkm, log_mean_Th2$FC_72hr_exon)) +  
  geom_point(alpha=0.5, colour="darkgrey") +
  geom_point(data=merged_SEeeee, aes(merged_SEeeee$Mean_72hr_exon_fpkm,merged_SEeeee$FC_72hr_exon), alpha=0.5, colour="dodgerblue") + 
  geom_point(data=merged_SE_17_for_Th2, aes(merged_SE_17_for_Th2$Mean_72hr_exon_fpkm,merged_SE_17_for_Th2$FC_72hr_exon), alpha=0.5, colour="gold1") +
  geom_point(data=merged_SE_1_for_Th2, aes(merged_SE_1_for_Th2$Mean_72hr_exon_fpkm,merged_SE_1_for_Th2$FC_72hr_exon), alpha=0.5, colour="firebrick1") +
  ggtitle("Th2 72 hr exon fpkm to exon FC", subtitle = "applied min exon cutoff 2, log transformed, replicates combined") +
 coord_cartesian(xlim = 0:11, ylim = -0.7:7.5) #this line trims the plot, so the one high FC outlier doesn't smoosh everything
  
all3 + facet_grid(Mean_72hr_exon_fpkm ~ FC_72hr_exon)

ggarrange(all3, common.legend = TRUE, legend = "bottom")  
  scale_colour_manual(name="SE genes", labels = c("Th2" = "blue", "Th17" = "red", "Th1" = "green")) 
  
  
#se FC interactive
library(plotly)

se_p_exonFC <- plot_ly(data = merged_SEeeee, x = merged_SEeeee$Mean_72hr_exon_fpkm, y = merged_SEeeee$FC_72hr_exon, 
                       alpha = 0.5, symbol = merged_SEeeee$SE_loost)

