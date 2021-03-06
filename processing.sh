#This doc is made of examples from individual scripts used for processing and peakcalling in different dirs...definitely could have done this more gracefully...but here it is so I remember what I did...
#!/bin/bash

##In dir /base_bams

###merge bams 
###(just a crude merge of everything here) 

module load samtools

samtools merge 1TF_plusH3K27.bam TBET_ChIP.bam ChIP_Tcell_WT_H3K27ac_nc_trim.bam
samtools merge 3TF_plusH3K27.bam TBET_ChIP.bam STAT1_ChIP.bam STAT4_ChIP.bam ChIP_Tcell_WT_H3K27ac_nc_trim.bam
samtools merge 3TF_only.bam TBET_ChIP.bam STAT1_ChIP.bam STAT4_ChIP.bam

###sort bams 

module load samtools

samtools sort 1TF_plusH3K27.bam  > sorted_1TF_plusH3K27.bam
samtools sort 3TF_plusH3K27.bam > sorted_3TF_plusH3K27.bam
samtools sort  3TF_only.bam > sorted_3TF_only.bam

##In dir /Combinations

###this is an important step - here we're selecting only for regions (peaks) that overlap atleast 50%, and pulling down that overlap

module load bedtools

bedtools intersect -a sorted_cp_peaks_TBET.bed \
         -b sorted_regions_H3K27.bed -f 0.5 > NEW_1TF_plus_H3K27_sorted_comb.bed
bedtools intersect -a sorted_cp_peaks_TBET.bed \
         -b sorted_cp_peaks_STAT4.bed  sorted_cp_peaks_STAT1.bed sorted_regions_H3K27.bed -f 0.5 > 3TF_plus_H3K27_sorted_comb.bed

###(sort again w/ i.e. bedtools sort -i regions_H3K27.bed > sorted_regions_H3K27.bed)

##In dir /Homer

###make tag directories first

module load homer
module load samtools

makeTagDirectory TBET_dir/ /data/belfordak/ChIP_seq/T_cells/Bowtie_mapped/TBET_ChIP.bam
makeTagDirectory GATA3_dir/ /data/belfordak/ChIP_seq/T_cells/Bowtie_mapped/GATA3_ChIP.bam
makeTagDirectory STAT1_dir/ /data/belfordak/ChIP_seq/T_cells/Bowtie_mapped/STAT1_ChIP.bam
makeTagDirectory STAT4_dir/ /data/belfordak/ChIP_seq/T_cells/Bowtie_mapped/STAT4_ChIP.bam
makeTagDirectory STAT6_dir/ /data/belfordak/ChIP_seq/T_cells/Bowtie_mapped/STAT6_ChIP.bam

###peakcalling

module load homer

findPeaks TBET_dir -style factor -o auto -i inputTcell_dir
findPeaks STAT1_dir -style factor -o auto -i inputTcell_dir
findPeaks STAT4_dir -style factor -o auto -i inputTcell_dir
findPeaks H3K27acTcell_dir -style factor -o auto -i inputTcell_dir
findPeaks STAT6_dir -style factor -o auto -i inputTcell_dir

##.txt to .bed 

module load homer

pos2bed.pl /data/belfordak/ChIP_seq/T_cells/Homer/H3K27acTcell_dir2/regions.txt > /data/belfordak/ChIP_seq/T_cells/Homer/H3K27acTcell_dir2/regions_H3K27.bed

##make UCSC visualisation file to look at tracks

module load homer
module load bedtools

makeUCSCfile H3K27acTcell_dir2/ -o auto
