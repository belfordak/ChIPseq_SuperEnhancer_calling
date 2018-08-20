#Script for each cell type with each cell types SE list (3 x 3 array)

#everything is already loaded into the env., see the ThX_processing scripts for that
#Th1 SE = firebrick1 / Th2_SE = dodgerblue / Th17_SE = gold1

#making the se lists

#Th1
#merged_SEeeeeee works for merged_SE_1_for_Th1
#need to fix log_mean_base....put row.names as column 1

log_mean_basecsv_Tcell$gene_names <- rownames(log_mean_basecsv_Tcell)
View(log_mean_basecsv_Tcell)
log_mean_basecsv_Tcell2 <- log_mean_basecsv_Tcell[,c(11,1:10)]

merged_SE_1_for_Th1 <- merged_SEeeee #just for clarity here
merged_SE_17_for_Th1 <- merge(log_mean_basecsv_Tcell2, SE_list_Th17.1,  by = 1, all = F)
merged_SE_2_for_Th1 <- merge(log_mean_basecsv_Tcell2, SE_list_Th2.1,  by = 1, all = F)

#Th17
merged_SE_17_for_Th17 <- merge(lmb_17, SE_list_Th17.1,  by = 1, all = F)
merged_SE_1_for_Th17 <- merge(lmb_17, SE_list_Th1,  by = 1, all = F)
merged_SE_2_for_Th17 <- merge(lmb_17, SE_list_Th2.1,  by = 1, all = F)

#don't need to trim Th1, plots don't look so bad
#Th1 - Th1_SE #WORKS!
ggplot(data = log_mean_basecsv_Tcell, aes(log_mean_basecsv_Tcell$Mean_72hr_exon_fpkm, log_mean_basecsv_Tcell$FC_72hr_exon)) +  
  geom_point(alpha=0.5, colour="darkgrey") +
  geom_point(data=merged_SEeeee, aes(merged_SEeeee$Mean_72hr_exon_fpkm,merged_SEeeee$FC_72hr_exon), alpha=0.5, colour="firebrick1") +
  ggtitle("Th1 72 hr exon fpkm to exon FC - Th1 SE genes", subtitle = "applied min exon cutoff 2, log transformed, replicates combined") +
  coord_cartesian(xlim = 0:11, ylim = -0.7:6)
  
#a little different for these ones
#Th1 - Th2_SE
ggplot(data = log_mean_basecsv_Tcell2, aes(log_mean_basecsv_Tcell2$Mean_72hr_exon_fpkm, log_mean_basecsv_Tcell2$FC_72hr_exon)) +  
  geom_point(alpha=0.5, colour="darkgrey") +
  geom_point(data=merged_SE_2_for_Th1, aes(merged_SE_2_for_Th1$Mean_72hr_exon_fpkm,merged_SE_2_for_Th1$FC_72hr_exon), alpha=0.5, colour="dodgerblue") + 
  ggtitle("Th1 72 hr exon fpkm to exon FC - Th2 SE genes", subtitle = "applied min exon cutoff 2, log transformed, replicates combined")+
  

#Th1 - Th17_SE
ggplot(data = log_mean_basecsv_Tcell, aes(log_mean_basecsv_Tcell$Mean_72hr_exon_fpkm, log_mean_basecsv_Tcell$FC_72hr_exon)) +  
  geom_point(alpha=0.5, colour="darkgrey") +
  geom_point(data=merged_SE_17_for_Th1, aes(merged_SE_17_for_Th1$Mean_72hr_exon_fpkm,merged_SE_17_for_Th1$FC_72hr_exon), alpha=0.5, colour="gold1") +
  ggtitle("Th1 72 hr exon fpkm to exon FC - Th17 SE genes", subtitle = "applied min exon cutoff 2, log transformed, replicates combined") 

#---------

#Th2 - Th1_SE
ggplot(data = log_mean_Th2, aes(log_mean_Th2$Mean_72hr_exon_fpkm, log_mean_Th2$FC_72hr_exon)) +  
  geom_point(alpha=0.5, colour="darkgrey") +
  geom_point(data=merged_SE_1_for_Th2, aes(merged_SE_1_for_Th2$Mean_72hr_exon_fpkm,merged_SE_1_for_Th2$FC_72hr_exon), alpha=0.5, colour="firebrick1") +
  ggtitle("Th2 72 hr exon fpkm to exon FC - Th1 SE genes", subtitle = "applied min exon cutoff 2, log transformed, replicates combined") +
  coord_cartesian(xlim = 0:11, ylim = -0.7:7.5) #this line trims the plot, so the one high FC outlier doesn't smoosh the plot

#Th2 - Th2_SE
ggplot(data = log_mean_Th2, aes(log_mean_Th2$Mean_72hr_exon_fpkm, log_mean_Th2$FC_72hr_exon)) +  
  geom_point(alpha=0.5, colour="darkgrey") +
  geom_point(data=merged_SEeeee, aes(merged_SEeeee_2$Mean_72hr_exon_fpkm,merged_SEeeee_2$FC_72hr_exon), alpha=0.5, colour="dodgerblue") + 
  ggtitle("Th2 72 hr exon fpkm to exon FC - Th2 SE genes", subtitle = "applied min exon cutoff 2, log transformed, replicates combined") +
  coord_cartesian(xlim = 0:11, ylim = -0.7:7.5) #this line trims the plot, so the one high FC outlier doesn't smoosh the plot

#Th2 - Th17_SE
ggplot(data = log_mean_Th2, aes(log_mean_Th2$Mean_72hr_exon_fpkm, log_mean_Th2$FC_72hr_exon)) +  
  geom_point(alpha=0.5, colour="darkgrey") +
  geom_point(data=merged_SE_17_for_Th2, aes(merged_SE_17_for_Th2$Mean_72hr_exon_fpkm,merged_SE_17_for_Th2$FC_72hr_exon), alpha=0.5, colour="gold1") +
  ggtitle("Th2 72 hr exon fpkm to exon FC - Th17 SE genes", subtitle = "applied min exon cutoff 2, log transformed, replicates combined") +
  coord_cartesian(xlim = 0:11, ylim = -0.7:7.5) #this line trims the plot, so the one high FC outlier doesn't smoosh everything

#---------
  
#Th17 - Th1_SE
ggplot(data = lmb_17, aes(lmb_17$Mean_72hr_exon_fpkm, lmb_17$FC_72hr_exon)) +
  geom_point(alpha=0.5, colour="darkgrey") +
  geom_point(data=merged_SE_1_for_Th17, aes(merged_SE_1_for_Th17$Mean_72hr_exon_fpkm, merged_SE_1_for_Th17$FC_72hr_exon), alpha=0.5, colour="firebrick1") +
  ggtitle("Th17 72 hr exon fpkm to exon FC - Th1 SE genes", subtitle = "applied min exon cutoff 2, log transformed, replicates combined") +
  coord_cartesian(xlim = 0:11, ylim = -0.7:6) #this line trims the plot, so the one high FC outlier doesn't smoosh the plot

#Th17 - Th2_SE
ggplot(data = lmb_17, aes(lmb_17$Mean_72hr_exon_fpkm, lmb_17$FC_72hr_exon)) +  
  geom_point(alpha=0.5, colour="darkgrey") +
  geom_point(data=merged_SE_2_for_Th17, aes(merged_SE_2_for_Th17$Mean_72hr_exon_fpkm, merged_SE_2_for_Th17$FC_72hr_exon), alpha=0.5, colour="dodgerblue") + 
  ggtitle("Th17 72 hr exon fpkm to exon FC - Th2 SE genes", subtitle = "applied min exon cutoff 2, log transformed, replicates combined") +
  coord_cartesian(xlim = 0:11, ylim = -0.7:6) #this line trims the plot, so the one high FC outlier doesn't smoosh the plot

#Th17 - Th17_SE
ggplot(data = lmb_17, aes(lmb_17$Mean_72hr_exon_fpkm, lmb_17$FC_72hr_exon)) +  
  geom_point(alpha=0.5, colour="darkgrey") +
  geom_point(data=merged_SE_17_for_Th17, aes(merged_SE_17_for_Th17$Mean_72hr_exon_fpkm, merged_SE_17_for_Th17$FC_72hr_exon), alpha=0.5, colour="gold1") +
  ggtitle("Th17 72 hr exon fpkm to exon FC - Th17 SE genes", subtitle = "applied min exon cutoff 2, log transformed, replicates combined") +
  coord_cartesian(xlim = 0:11, ylim = -0.7:6) #this line trims the plot, so the one high FC outlier doesn't smoosh everything

#-------------

#Plots with all SEs for each Th1 dataset (these have no legends)

#Th1
ggplot(data = log_mean_basecsv_Tcell2, aes(log_mean_basecsv_Tcell2$Mean_72hr_exon_fpkm,log_mean_basecsv_Tcell2$FC_72hr_exon)) +  
  geom_point(alpha=0.5, colour="darkgrey") +
  geom_point(data=merged_SE_2_for_Th1, aes(merged_SE_2_for_Th1$Mean_72hr_exon_fpkm,merged_SE_2_for_Th1$FC_72hr_exon), alpha=0.5, colour="dodgerblue") + 
  geom_point(data=merged_SE_17_for_Th1, aes(merged_SE_17_for_Th1$Mean_72hr_exon_fpkm,merged_SE_17_for_Th1$FC_72hr_exon), alpha=0.5, colour="gold1") +
  geom_point(data=merged_SE_1_for_Th1, aes(merged_SE_1_for_Th1$Mean_72hr_exon_fpkm,merged_SE_1_for_Th1$FC_72hr_exon), alpha=0.5, colour="firebrick1") +
  ggtitle("Th1 72 hr exon fpkm to exon FC", subtitle = "applied min exon cutoff 2, log transformed, replicates combined") +
  coord_cartesian(xlim = 0:11, ylim = -0.7:7.5) #this line trims the plot, so the one high FC outlier doesn't smoosh everything

#Th2 (this one is also in the Th2 processing file)
ggplot(data = log_mean_Th2, aes(log_mean_Th2$Mean_72hr_exon_fpkm, log_mean_Th2$FC_72hr_exon)) +  
  geom_point(alpha=0.5, colour="darkgrey") +
  geom_point(data=merged_SEeeee_2, aes(merged_SEeeee_2$Mean_72hr_exon_fpkm,merged_SEeeee_2$FC_72hr_exon), alpha=0.5, colour="dodgerblue") + 
  geom_point(data=merged_SE_17_for_Th2, aes(merged_SE_17_for_Th2$Mean_72hr_exon_fpkm,merged_SE_17_for_Th2$FC_72hr_exon), alpha=0.5, colour="gold1") +
  geom_point(data=merged_SE_1_for_Th2, aes(merged_SE_1_for_Th2$Mean_72hr_exon_fpkm,merged_SE_1_for_Th2$FC_72hr_exon), alpha=0.5, colour="firebrick1") +
  ggtitle("Th2 72 hr exon fpkm to exon FC", subtitle = "applied min exon cutoff 2, log transformed, replicates combined") +
  coord_cartesian(xlim = 0:11, ylim = -0.7:7.5) #this line trims the plot, so the one high FC outlier doesn't smoosh everything

#this one is for full + names of top 
#have to do the row.names as column thing for labeling
log_mean_Th2 <- cbind(Row.Names = rownames(log_mean_Th2), log_mean_Th2)

ggplot(data = log_mean_Th2, aes(log_mean_Th2$Mean_72hr_exon_fpkm, log_mean_Th2$FC_72hr_exon)) +  
  geom_point(alpha=0.5, colour="darkgrey") +
  geom_point(data=merged_SEeeee_2, aes(merged_SEeeee_2$Mean_72hr_exon_fpkm,merged_SEeeee_2$FC_72hr_exon), alpha=0.5, colour="dodgerblue") + 
  geom_point(data=merged_SE_17_for_Th2, aes(merged_SE_17_for_Th2$Mean_72hr_exon_fpkm,merged_SE_17_for_Th2$FC_72hr_exon), alpha=0.5, colour="gold1") +
  geom_point(data=merged_SE_1_for_Th2, aes(merged_SE_1_for_Th2$Mean_72hr_exon_fpkm,merged_SE_1_for_Th2$FC_72hr_exon), alpha=0.5, colour="firebrick1") +
  ggtitle("Th2 72 hr exon fpkm to exon FC", subtitle = "applied min exon cutoff 2, log transformed, replicates combined") +
  geom_text_repel(aes(label=ifelse(log_mean_Th2$FC_72hr_exon>4, as.character(log_mean_Th2$Row.Names),'')),hjust=0,vjust=0)

#Th17
ggplot(data = lmb_17, aes(lmb_17$Mean_72hr_exon_fpkm,lmb_17$FC_72hr_exon)) +  
  geom_point(alpha=0.5, colour="darkgrey") +
  geom_point(data=merged_SE_2_for_Th17, aes(merged_SE_2_for_Th17$Mean_72hr_exon_fpkm,merged_SE_2_for_Th17$FC_72hr_exon), alpha=0.5, colour="dodgerblue") + 
  geom_point(data=merged_SE_17_for_Th17, aes(merged_SE_17_for_Th17$Mean_72hr_exon_fpkm,merged_SE_17_for_Th17$FC_72hr_exon), alpha=0.5, colour="gold1") +
  geom_point(data=merged_SE_1_for_Th17, aes(merged_SE_1_for_Th17$Mean_72hr_exon_fpkm,merged_SE_1_for_Th17$FC_72hr_exon), alpha=0.5, colour="firebrick1") +
  ggtitle("Th17 72 hr exon fpkm to exon FC", subtitle = "applied min exon cutoff 2, log transformed, replicates combined") +
  coord_cartesian(xlim = 0:11, ylim = -0.7:6) #this line trims the plot, so the one high FC outlier doesn't smoosh everything
  
#this one is for full + names of top 
ggplot(data = lmb_17, aes(lmb_17$Mean_72hr_exon_fpkm,lmb_17$FC_72hr_exon)) +  
  geom_point(alpha=0.5, colour="darkgrey") +
  geom_point(data=merged_SE_2_for_Th17, aes(merged_SE_2_for_Th17$Mean_72hr_exon_fpkm,merged_SE_2_for_Th17$FC_72hr_exon), alpha=0.5, colour="dodgerblue") + 
  geom_point(data=merged_SE_17_for_Th17, aes(merged_SE_17_for_Th17$Mean_72hr_exon_fpkm,merged_SE_17_for_Th17$FC_72hr_exon), alpha=0.5, colour="gold1") +
  geom_point(data=merged_SE_1_for_Th17, aes(merged_SE_1_for_Th17$Mean_72hr_exon_fpkm,merged_SE_1_for_Th17$FC_72hr_exon), alpha=0.5, colour="firebrick1") +
  ggtitle("Th17 72 hr exon fpkm to exon FC", subtitle = "applied min exon cutoff 2, log transformed, replicates combined") + 
  geom_text_repel(aes(label=ifelse(lmb_17$FC_72hr_exon>4, as.character(lmb_17$`rownames(log_mean_Th17)`),'')),hjust=0,vjust=0)
