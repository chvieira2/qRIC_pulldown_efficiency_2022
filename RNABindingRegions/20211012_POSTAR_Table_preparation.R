setwd("K:/Results/20210903_Phospho_qRIC_revisionPaper/RNABindingRegions")

library(data.table)
library(limma)



#### preparation of PORSTAR_2019 overlaps with genes ####
## "Human RBP binding sites" was downloaded from http://lulab.life.tsinghua.edu.cn/postar/download.php
# POSTAR_2019 <- fread(file = "K:/Datasets/202012_POSTAR_human_RBP_binding_sites.txt", sep = "\t", na = "", quote = F)



## Run only once to generate the minimal BED file containing only essential columns
#fwrite(POSTAR_2019[,1:6], file = "POSTAR_2019_essential.bed", sep = "\t", na = "", quote = F, row.names = F, col.names = F)

# library(rtracklayer)
# twist.summits <- import.bed("POSTAR_2019_essential.bed")
# 
# # Change the Chrmossome annotation to match the table below
# seqlevelsStyle(twist.summits)<-"NCBI"
# 
# # Human hg19 annotation file
# anno<-import("K:/Datasets/Homo_sapiens.GRCh37.87.gtf")
# 
# # Select only gene annotations
# anno<-anno[anno$type == "gene"]
# 
# # overlap with my BED file
# olaps <- findOverlaps(twist.summits, anno)
# olaps <- cbind(data.frame(twist.summits)[queryHits(olaps),],data.frame(anno)[subjectHits(olaps),])
# 
# olaps <- subset(olaps, select = c("seqnames", "start", "end","width","strand", "name","score","gene_id", "gene_name", "gene_biotype"))
# 
#fwrite(olaps, file = "Overlap_peaks_Genes.bed", sep = "\t", na = "", quote = F, row.names = F, col.names = T)






#### Add gene info to POSTAR ####
## "Human RBP binding sites" was downloaded from http://lulab.life.tsinghua.edu.cn/postar/download.php
# POSTAR_2019 <- fread(file = "K:/Datasets/202012_POSTAR_human_RBP_binding_sites.txt", sep = "\t", na = "", quote = F)
# 
# 
# # Load overlaps
# olaps <- fread("Overlap_peaks_Genes.bed", sep = "\t", stringsAsFactors = FALSE, data.table = F, check.names = T)
# olaps <- subset(olaps, select = c("name", "gene_id", "gene_name", "gene_biotype"))
# 
# 
# # Set col names for POSTAR
# colnames(POSTAR_2019) <- c("Chromosome", "Peak_start", "Peak_end", "name", "score", "Strand", "RBP_name", "Method", "Cell_line", "Source_ID", "width")
# POSTAR_2019 <- subset(POSTAR_2019, select = c("name", "RBP_name", "Cell_line", "Method"))
# 
# # Merge overlaps and POSTAR
# POSTAR_2019_merged <- merge(POSTAR_2019, olaps, by = "name", all.x = T)
# 
# ## There are duplicated peaks in the POSTAR database. Remove duplicates from POSTAR_2019 database
# POSTAR_2019_merged <- distinct(POSTAR_2019_merged)
# 
# 
# fwrite(POSTAR_2019_merged, file = "POSTAR_2019_merged.txt", sep = "\t", na = "", quote = F, row.names = F)


#### Counts total ####
# POSTAR_2019_merged <- fread(file = "POSTAR_2019_merged.txt", sep = "\t", na = "", quote = F)
# 
# 
# ## Peak total
# df_PeakCount <- data.table(table(POSTAR_2019_merged$RBP_name, POSTAR_2019_merged$Cell_line, POSTAR_2019_merged$Method))
# df_PeakCount <- subset(df_PeakCount, N != 0)
# colnames(df_PeakCount) <- c("RBP_name", "Cell_line",
#                             "Method", "Number_peaks")
# 
# 
# fwrite(df_PeakCount, file = "Peaks_counts_POSTAR_2019_merged.txt", sep = "\t", na = "", quote = F, row.names = F)



#### Counts per gene ####
## Peak per gene
# First give a empty name for gene_id for those not matching genes so they can be counted later
# POSTAR_2019_merged[is.na(POSTAR_2019_merged$gene_id), "gene_id"] <- ""
# 
# # Create table with numbers of peaks per gene. This MUST ran on Abacus
# # df_PeakGene <- data.table(table(POSTAR_2019_merged$RBP_name, POSTAR_2019_merged$Cell_line, POSTAR_2019_merged$Method, POSTAR_2019_merged$gene_id))
# df_PeakGene <- subset(df_PeakGene, N != 0)
# colnames(df_PeakGene) <- c("RBP_name", "Cell_line",
#                             "Method", "gene_id", "Number_peaks")
# 
# 
# fwrite(df_PeakGene, file = "PeaksGene_counts_POSTAR_2019_merged.txt", sep = "\t", na = "", quote = F, row.names = F)



df_PeakGene <- fread(file = "PeaksGene_counts_POSTAR_2019_merged.txt", sep = "\t", na = "", quote = F, stringsAsFactors = F)
df_PeakGene <- df_PeakGene[!is.na(df_PeakGene$gene_id), ]
df_PeakGene[is.na(df_PeakGene$gene_id), "gene_id"] <- "NotOverlapped"




## Obtain absolute copy number per gene from Mukherjee et al 2016 (Ohler lab)
neel <- read.table("K:/Datasets/Mukherjee_2016_TranscriptomeCopyNumber_GSE84722_ST3.txt", header=T, sep="\t")  

# Remove version annotation (dot suffixes)
neel$Gene <- gsub("\\..*","",neel$Gene)

#neel <- neel[,c(1, grep("Locus", colnames(neel)))]
#neel[2:4] <- normalizeQuantiles(log10(neel[2:4]))

## Replace Inf for NA
#neel[sapply(neel, is.infinite)] <- NA

## Replace NaN for NA
#neel[sapply(neel, is.nan)] <- NA


# Remove only infinite values
#neel <- subset(neel, apply(neel[2:4], 1, function(x) {
#  ifelse(sum(is.na(x)) < 3, T, F)}))

#neel$Copies <- rowMeans(neel[2:4], na.rm = T)



## Add copy number info to data frame with peaks per gene number
df_PeakGene <- merge(df_PeakGene, neel[,c("Gene", "Copies")], by.x = "gene_id", by.y = "Gene", all.x = T)

## Substitute NA copy numbers per mean copy number of all genes
#df_PeakGene[is.na(df_PeakGene$Copies), "Copies"] <- median(neel$Copies, na.rm = T)

## Copies times number peaks
df_PeakGene$Number_peaksXcopies <- df_PeakGene$Number_peaks * 10**df_PeakGene$Copies




data_summary <- function(data = df_PeakGene, varname, groupnames){
  require(plyr)
  summary_func <- function(x, col){
    c(sum = sum(x[[col]], na.rm=TRUE))
  }
  data_sum<-ddply(data, groupnames, .fun=summary_func,
                  varname)
  data_sum <- rename(data_sum, c("sum" = "Number_peaks"))
  return(data_sum)
}

all_summary <- data_summary(df_PeakGene, varname="Number_peaksXcopies", 
                            groupnames=c("RBP_name",
                                         "Cell_line",
                                         "Method"))

fwrite(all_summary, file = "PeaksGene_counts_POSTAR_2019.txt", sep = "\t", na = "", quote = F, row.names = F)



### Average number of sites per trancript for reviewer #3
data_summary2 <- function(data = df_PeakGene, varname, groupnames){
  require(plyr)
  summary_func <- function(x, col){
    c(mean = mean(x[[col]], na.rm=TRUE))
  }
  data_sum<-ddply(data, groupnames, .fun=summary_func,
                  varname)
  data_sum <- rename(data_sum, c("mean" = "Number_peaks"))
  return(data_sum)
}

all_summary_averageNumberPeak <- data_summary2(df_PeakGene, varname="Number_peaks", 
                            groupnames=c("RBP_name",
                                         "Cell_line",
                                         "Method"))

fwrite(all_summary_averageNumberPeak, file = "PeaksMeanGene_counts_POSTAR_2019.txt", sep = "\t", na = "", quote = F, row.names = F)
