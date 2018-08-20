#straight expression difference (for mRNA decay)

library(ggplot2)
library(ggrepel)
library(ggridges)
library(GGally)
library(magrittr)
library(dplyr)
library(DESeq2)
library(DEFormats)
library(RColorBrewer)

#pull in expression datasets

#Th1
Th1_exp <- read.csv("~/Desktop/Tcell_only/Tables/Th1/noNA_fpkm_int_exon.csv")
#log transform
Th1_exp <- base_csv_Tcell
Th1_exp_log <- log2(Th1_exp) 
Th1_exp_log <- round(Th1_exp_log, 3)

#Th2
Th2_exp <- read.csv("~/Desktop/Tcell_only/Tables/Th2/Actual_Th2_full_avgs_noNA.csv")
#log transform
Th2_exp_log <- log2(Th2_exp[,2:24]) 
Th2_exp_log <- cbind(Th2_exp$X, Th2_exp_log)
Th2_exp_log <- round(Th2_exp_log, 3)

#Th17
Th17_exp <- read.csv("~/Desktop/Tcell_only/Tables/Th17/Actual_Th17_full_avgs_noNA.csv")
#log transform
Th17_exp_log <- log2(Th17_exp[,2:24])
Th17_exp_logg <- cbind(Th17_exp$X, Th17_exp_log)
Th17_exp_log <- round(Th17_exp_log, 3)

#pull in SE lists with values
SE_list <-read.csv("~/Desktop/Tcell_only/Tables/Th1/SE_list_plus_my_values")
#SE_loost <- SE_list$X.2

SE_list_Th2 <- read.csv("~/Desktop/Tcell_only/Tables/Th2/full_avgs_noNA.csv")
#SE_list_Th2.1 <- as.data.frame(SE_list_Th2[,4])

SE_list_Th17 <- read.csv("~/Desktop/Tcell_only/Tables/Th17/Th17_full_avgs_noNA.csv")
SE_list_Th17.1 <- merge(SE_list_Th17, Th17_exp_log, by = 1)

#SE_list_Th17.1 <- as.data.frame(SE_list_Th17[,4])
#i do also have SEs with expression for Th2 and 17 - they're the "full", but no "Accurate-" tables in the /tables directory


#parallel plot, one for each Th plu combined

filtered_5TH2 <- filter(Th2_exp, Th2_exp$Mean_0hr_exon_fpkm > 0,  Th2_exp$Mean_72hr_exon_fpkm > 0,  Th2_exp$Mean_2w_exon_fpkm > 0)
filtered_5TH23 <- filter(Th2_exp, Th2_exp$Mean_0hr_exon_fpkm > 100,  Th2_exp$Mean_72hr_exon_fpkm > 100,  Th2_exp$Mean_2w_exon_fpkm > 100)

filtered_Th1 <- filter(Th1_exp, Th1_exp$Mean_0hr_exon_fpkm > 200,  Th1_exp$Mean_72hr_exon_fpkm > 200,  Th1_exp$Mean_2w_exon_fpkm > 200)
  

ggparcoord(data = filtered_Th1[,c(19,21,20)], scale = "globalminmax", groupColumn = 2 ) + scale_x_discrete(expand = c(0.02, 0.02)) +
  ggtitle("Gene expression at key differentiation time points, Th1") + labs(x="Time Point", y="Expression (fpkm)") 

#sorted df for poster: 
sort_Fig1 <- Th1_exp[order(-Th1_exp$Mean_72hr_exon_fpkm),]
sort_Fig1 <- sort_Fig1[c(19,21,20)]
write.csv(sort_Fig1, "~/Desktop/Tcell_only/Poster_Fig_tables/Fig1.csv") 
  
#just SEs  
ggparcoord(data = SE_list_final[,c(11,13,12)], scale = "globalminmax", groupColumn = 1) + scale_x_discrete(expand = c(0.02, 0.02)) +
  ggtitle("Gene expression at key differentiation time points", subtitle = "Super-Enhancer Genes") + labs(x="Time Point", y="Expression (fpkm)")+
  scale_color_manual(values = c( "1" = "#ef8a62" )) 

#sorted df for poster: 
sort_Fig2 <- SE_list_final[order(-SE_list_final$Mean_72hr_exon_fpkm.y),]
sort_Fig2 <- sort_Fig2[c(1,11,13,12)]
write.csv(sort_Fig2, "~/Desktop/Tcell_only/Poster_Fig_tables/Fig2.csv")

#SE subsets based on slope (not using)
#up-up
#A<B, B<C
up_up <- subset(SE_list_final, SE_list_final$Mean_0hr_exon_fpkm < SE_list_final$Mean_72hr_exon_fpkm.x)
up_up_final <- subset(up_up, up_up$Mean_72hr_exon_fpkm.x < up_up$Mean_2w_exon_fpkm.x) 
#up-down
#A<B, B>C
up_down <- subset(SE_list_final, SE_list_final$Mean_0hr_exon_fpkm < SE_list_final$Mean_72hr_exon_fpkm.x)
up_down_final <- subset(up_down, up_down$Mean_72hr_exon_fpkm.x > up_down$Mean_2w_exon_fpkm.x) 
#down-up
#A>B, B<C
down_up <- subset(SE_list_final, SE_list_final$Mean_0hr_exon_fpkm > SE_list_final$Mean_72hr_exon_fpkm.x)
down_up_final <- subset(down_up, down_up$Mean_72hr_exon_fpkm.x < down_up$Mean_2w_exon_fpkm.x) 
#down-down
#A>B, B>C
down_down <- subset(SE_list_final, SE_list_final$Mean_0hr_exon_fpkm > SE_list_final$Mean_72hr_exon_fpkm.x)
down_down_final <- subset(down_down, down_down$Mean_72hr_exon_fpkm.x > down_down$Mean_2w_exon_fpkm.x) 

ggparcoord(data = down_down_final[,c(11,13,12)], scale = "globalminmax", groupColumn = 2) + scale_x_discrete(expand = c(0.02, 0.02)) 


#non SE Th1 slope subsets ******* all Th1, 1.5x

#up-up
#A<B, B<C
ALL_up_up <- subset(Th1_exp, Th1_exp$Mean_0hr_exon_fpkm *1.5 < Th1_exp$Mean_72hr_exon_fpkm)
ALL_up_up_final <- subset(ALL_up_up, ALL_up_up$Mean_72hr_exon_fpkm *1.5 < ALL_up_up$Mean_2w_exon_fpkm) 

ggparcoord(data = ALL_up_up_final[,c(19,21,20)], scale = "globalminmax", groupColumn = 2) + scale_x_discrete(expand = c(0.02, 0.02)) +
  ggtitle("Gene expression at key differentiation time points, Th1", subtitle = "positive-positive slope, 1.5x") + labs(x="Time Point", y="Expression (fpkm)")

#here's the ratio:
ggparcoord(data = ALL_up_up_final[,c(25,27,26)], scale = "globalminmax", groupColumn = 2) + scale_x_discrete(expand = c(0.02, 0.02)) +
  ggtitle("Exon/intron ratio of genes at key differentiation time points, Th1", subtitle = "positive-positive slope, 1.5x") + labs(x="Time Point", y="Ratio (exon/intron)")

cor5 <- ALL_up_up_final[,c(19,21,20)]
cor6 <- ALL_up_up_final[,c(25,27,26)]
cor(cor5, cor6)

#sorted df for poster: 
sort_Fig3a <- ALL_up_up_final[order(-ALL_up_up_final$Mean_72hr_exon_fpkm),]
sort_Fig3a <- sort_Fig3a[c(19,21,20)]
write.csv(sort_Fig3a, "~/Desktop/Tcell_only/Poster_Fig_tables/Fig3a.csv")

#sorted df for poster (ratio): 
sort_Fig3d <- ALL_up_up_final[order(-ALL_up_up_final$Mean_72hr_ratio),]
sort_Fig3d <- sort_Fig3d[c(25,27,26)]
write.csv(sort_Fig3d, "~/Desktop/Tcell_only/Poster_Fig_tables/Fig3d.csv")


#up-down
#A<B, B>C
ALL_up_down <- subset(Th1_exp, Th1_exp$Mean_0hr_exon_fpkm *1.5 < Th1_exp$Mean_72hr_exon_fpkm)
ALL_up_down_final <- subset(ALL_up_down, ALL_up_down$Mean_72hr_exon_fpkm > 1.5 * ALL_up_down$Mean_2w_exon_fpkm) 

ggparcoord(data = ALL_up_down_final[,c(19,21,20)], scale = "globalminmax", groupColumn = 2) + scale_x_discrete(expand = c(0.02, 0.02)) +
  ggtitle("Gene expression at key differentiation time points, Th1", subtitle = "positive-negative slope, 1.5x") + labs(x="Time Point", y="Expression (fpkm)")

#here's the ratio:
ggparcoord(data = ALL_up_down_final[,c(25,27,26)], scale = "globalminmax", groupColumn = 2) + scale_x_discrete(expand = c(0.02, 0.02)) +
  ggtitle("Exon/intron ratio of genes at key differentiation time points, Th1", subtitle = "positive-negative slope, 1.5x") + labs(x="Time Point", y="Ratio (exon/intron)")

cor3 <- ALL_up_down_final[,c(19,21,20)]
cor4 <- ALL_up_down_final[,c(25,27,26)]
cor(cor3, cor4)

#sorted df for poster: 
sort_Fig3b <- ALL_up_down_final[order(-ALL_up_down_final$Mean_72hr_exon_fpkm),]
sort_Fig3b <- sort_Fig3b[c(19,21,20)]
write.csv(sort_Fig3b, "~/Desktop/Tcell_only/Poster_Fig_tables/Fig3b.csv")

#sorted df for poster (ratio): 
sort_Fig3e <- ALL_up_down_final[order(-ALL_up_down_final$Mean_72hr_ratio),]
sort_Fig3e <- sort_Fig3e[c(25,27,26)]
write.csv(sort_Fig3e, "~/Desktop/Tcell_only/Poster_Fig_tables/Fig3e.csv")

#down-up
#A>B, B<C
ALL_down_up <- subset(Th1_exp, Th1_exp$Mean_0hr_exon_fpkm > 1.5* Th1_exp$Mean_72hr_exon_fpkm)
ALL_down_up_final <- subset(ALL_down_up, ALL_down_up$Mean_72hr_exon_fpkm *1.5 < ALL_down_up$Mean_2w_exon_fpkm) 

ggparcoord(data = ALL_down_up_final[,c(19,21,20)], scale = "globalminmax", groupColumn = 2) + scale_x_discrete(expand = c(0.02, 0.02)) +
  ggtitle("Gene expression at key differentiation time points, Th1", subtitle = "negative-positive slope, 1.5x") + labs(x="Time Point", y="Expression (fpkm)")

#here's the ratio:
ggparcoord(data = ALL_down_up_final[,c(25,27,26)], scale = "globalminmax", groupColumn = 2) + scale_x_discrete(expand = c(0.02, 0.02)) +
  ggtitle("Exon/intron ratio of genes at key differentiation time points, Th1", subtitle = "negative-positive slope, 1.5x") + labs(x="Time Point", y="Ratio (exon/intron)")

cor1 <- ALL_down_up_final[,c(19,21,20)]
cor2 <- ALL_down_up_final[,c(25,27,26)]
cor(cor1, cor2)

#sorted df for poster: 
sort_Fig3c <- ALL_down_up_final[order(-ALL_down_up_final$Mean_72hr_exon_fpkm),]
sort_Fig3c <- sort_Fig3c[c(19,21,20)]
write.csv(sort_Fig3c, "~/Desktop/Tcell_only/Poster_Fig_tables/Fig3c.csv")

#sorted df for poster (ratio): 
sort_Fig3f <- ALL_down_up_final[order(-ALL_down_up_final$Mean_72hr_ratio),]
sort_Fig3f <- sort_Fig3f[c(25,27,26)]
write.csv(sort_Fig3f, "~/Desktop/Tcell_only/Poster_Fig_tables/Fig3f.csv")

#down-down  there's none
#A>B, B>C
ALL_down_down <- subset(Th1_exp, Th1_exp$Mean_0hr_exon_fpkm > 0.5 * Th1_exp$Mean_72hr_exon_fpkm)
ALL_down_down_final <- subset(ALL_down_down, ALL_down_down$Mean_72hr_exon_fpkm.x > 0.5 * ALL_down_down$Mean_2w_exon_fpkm) 

ggparcoord(data = ALL_down_down_final[,c(19,21,20)], scale = "globalminmax", groupColumn = 2) + scale_x_discrete(expand = c(0.02, 0.02)) +
  ggtitle("Gene expression at key differentiation time points, Th1", subtitle = "negative-negative slope, 1.5x") + labs(x="Time Point", y="Expression (fpkm)")


#donut charts!

rn_A_u_u_f <- as.data.frame(row.names(ALL_up_up_final))
rn_A_u_d_f <- as.data.frame(row.names(ALL_up_down_final))
rn_A_d_u_f <- as.data.frame(row.names(ALL_down_up_final))

SE_list$SE_loost %in% log_mean_basecsv_Tcell$SE_loost #test
SE_list$SE_loost[SE_list$SE_loost %in% rn_A_u_u_f$`row.names(ALL_up_up_final)`] #7/96
SE_list$SE_loost[SE_list$SE_loost %in% rn_A_u_d_f$`row.names(ALL_up_down_final)`] #9/630
SE_list$SE_loost[SE_list$SE_loost %in% rn_A_d_u_f$`row.names(ALL_down_up_final)`] #3/41

#make dfs out of just these SE genes
up_up_SEs <- c("Fosl2","Jarid2", "Mrpl33", "Rbpj", "Rbpj", "Tax1bp1", "Traf1")
up_down_SEs <- c("Bcl2", "Bcl3", "Cish", "Rara", "Rps6ka1", "Rps6ka1", "Sell", "Tnfrsf18", "Tomm20")
down_up_SEs <- c("Btg1", "Fasl", "Il7r")

#actually make it below

dat1 = data.frame(count=c(7, 89), category=c("SE", "non-SE"))
dat2 = data.frame(count=c(9, 621), category=c("SE", "non-SE"))
dat3 = data.frame(count=c(3, 38), category=c("SE", "non-SE"))

# Add addition columns, needed for drawing with geom_rect.
dat3$fraction = dat3$count / sum(dat3$count)
dat3 = dat3[order(dat3$fraction), ]
dat3$ymax = cumsum(dat3$fraction)
dat3$ymin = c(0, head(dat3$ymax, n=-1))

# Make the plot
p3 = ggplot(dat3, aes(fill=category, ymax=ymax, ymin=ymin, xmax=4, xmin=3)) +
  geom_rect() +
  coord_polar(theta="y") +
  xlim(c(0, 4)) +
  theme(panel.grid=element_blank()) +
  theme(axis.text=element_blank()) +
  theme(axis.ticks=element_blank()) +
  annotate("text", x = 0, y = 0, label = "SE= 7.3%") +
  labs(title="Super-enhancer and non-super-enhancer genes in the dataset", subtitle = "negative-positive slope subset")

p1 
p2
p3


#boxplot of FC

SE_list_final
bp_Th1SE_in_Th2exp <- merge(SE_list_final$SE_loost, Th2_exp_log, by = 1)
bp_Th1SE_in_Th17exp <- merge(SE_list_final$SE_loost, Th17_exp_logg, by = 1)

#FC calc
bp_Th1SE_in_Th2exp$FC_72hr_exon <- ((bp_Th1SE_in_Th2exp$Mean_72hr_exon_fpkm - bp_Th1SE_in_Th2exp$Mean_0hr_exon_fpkm) / bp_Th1SE_in_Th2exp$Mean_0hr_exon_fpkm)
bp_Th1SE_in_Th17exp$FC_72hr_exon <- ((bp_Th1SE_in_Th17exp$Mean_72hr_exon_fpkm - bp_Th1SE_in_Th17exp$Mean_0hr_exon_fpkm) / bp_Th1SE_in_Th17exp$Mean_0hr_exon_fpkm)

all_Tcells_FC72_for_Th1_SE 
maybe <- cbind( Th2 = bp_Th1SE_in_Th2exp$FC_72hr_exon, Th17 = bp_Th1SE_in_Th17exp$FC_72hr_exon)
maybe2 <- cbind(Th1 = SE_list_final$FC_72hr_exon, maybe)
maybe3 <- cbind(SE_list_final$SE_loost, maybe2)
#I manually put in the classes (excel) because that was easier...
write.csv(bp_Th1SE_in_Th2exp, "~/Desktop/iiintermed_FC_allcells.csv")

write.csv(maybe3, "~/Desktop/intermed_FC_allcells.csv")

all_Tcells_FC72_for_Th1_SE <- read.csv("~/Desktop/intermed_FC_allcells.csv")

bp <- boxplot(all_Tcells_FC72_for_Th1_SE[,2:5], main="Th1 Super-enhancer FC expression in different cell types",xlab="T cell", ylab="Fold Change 0-72hrs (fpkm)", 
        col = c("coral1","darkseagreen", "lightblue2"))

#with ggplot
library(reshape2)

f <- factor((all_Tcells_FC72_for_Th1_SE[,2:4]), label = c("Th1", "Th2", "Th17"))
dat.m <- melt(all_Tcells_FC72_for_Th1_SE,id.vars='Gene', measure.vars=c('Th1','Th2','Th17'))

dat.mm <- merge(dat.m , all_Tcells_FC72_for_Th1_SE, by = 1)

qplot(x =variable, y=value, data = dat.mm, color= variable, xlab = "Cell Type", ylab = "Expression (log fpkm)") + geom_boxplot(lwd=0.4, alpha = 0.5, shape = 1) + theme_bw() + 
  scale_color_manual(values =c("coral1","darkseagreen", "lightblue2")) + ggtitle("Th1 super-enhancer gene expression in different cell types")

qplot(x =variable, y=value, data = dat.mm, color= variable, facets = ~Class, xlab = "Cell Type", ylab = "Expression (log fpkm)") + geom_boxplot(lwd=0.4, alpha = 0.5, shape = 1) + theme_bw() + 
  scale_color_manual(values =c("coral1","darkseagreen", "lightblue2")) + ggtitle("Th1 super-enhancer gene expression in different cell types", 
  subtitle ="As present in slope subset groups")





all_SE_FC_bplot <- rbind(SE_list, SE_list_Th2)

SE_Th1_72hr <- rnorm(SE_list$Mean_72hr_exon_fpkm)
SE_Th2_72hr <- rnorm(SE_list_Th2$)
C <- rnorm(12)
mydf <- data.frame(y=c(A,B,C),x=c(rep(1,length(A)),rep(2,length(B)),rep(3,length(C))))
with(mydf, boxplot(y~x))

#mRNA halflife calculation? it needs a slope



#Alternative SE list on geom plot

SE_alt <- read.csv("~/Desktop/Tcell_only/SE_lists/SE_alt.csv")
#need to get out basecsvtcell row.names as column 1

SE_alt_merged <- merge(SE_alt[,5], base_csv_Tcell, by = 1, all = F)
  
  ggplot(data = Th1_exp_log , aes(Th1_exp_log$Mean_72hr_exon_fpkm, Th1_exp_log$FC_72hr_exon)) +  
  geom_point(alpha=0.5, colour="darkgrey") +
  geom_point(data=SE_alt_merged, aes(SE_alt_merged$Mean_72hr_exon_fpkm,SE_alt_merged$FC_72hr_exon), alpha=0.5, colour="firebrick1") +
  ggtitle("Th1 72 hr exon fpkm to exon FC - Th1 SE ALTERNATE genes", subtitle = "applied min exon cutoff 2, log transformed, replicates combined"))





  
   
#lollipop SE chart
#need FC 
#need to make x-axis = FC = ((B-A)/A) = (72-0)/0) = ((avg exon fpkm 72 - avg exon fpk 0) / (avg exon fpkm 0)...will need to create a variable..just new column for this..shouldn't be too hard... 
SE_list$FC_72hr_exon <- ((SE_list$Mean_72hr_exon_fpkm - SE_list$Mean_0hr_exon_fpkm) / SE_list$Mean_0hr_exon_fpkm)
#except SE_list_final doesnt need it

Above = ifelse(SE_list_final$Mean_72hr_ratio > 0, TRUE, FALSE)
SE_list_final$Mean_72hr_ratio=as.numeric(SE_list_final$Mean_72hr_ratio)

ggplot(SE_list_final, aes(SE_list_final$FC_72hr_exon, SE_list_final$SE_loost, color = Above)) +
  geom_segment(aes(x = FC_72hr_exon, y = SE_loost, xend = Above, yend = SE_loost), color = "grey50") +
  geom_point()


