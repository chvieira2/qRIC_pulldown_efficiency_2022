setwd("K:/Results/20210903_Phospho_qRIC_revisionPaper/output")

library(ggplot2)
library(reshape2)
library(RColorBrewer)
library(gridExtra)
library(dplyr)
library(limma)
library(ggrepel)
library(gplots)
library(data.table)
library(GGally)





Requantify = T
ReproducibilityDifference = 2
DeltaAnalysisThresholdEfficiency <- -1





#### proteinGroups table load and preparation ####

PG <- fread("../txt_RQ/proteinGroups.txt", sep = "\t", stringsAsFactors = FALSE, data.table = F, check.names = T)
PG_noRQ <- fread("../txt_noRQ/proteinGroups.txt", sep = "\t", stringsAsFactors = FALSE, data.table = F, check.names = T)
if (!Requantify) {
  PG <- PG_noRQ
  
}
rm(Requantify)



PG$Gene.names <- sapply(strsplit(PG$Gene.names, ";"), "[", 1)
PG$Majority.protein.IDs <- sapply(strsplit(PG$Majority.protein.IDs, ";"), "[", 1)

PG_noRQ$Gene.names <- sapply(strsplit(PG_noRQ$Gene.names, ";"), "[", 1)
PG_noRQ$Majority.protein.IDs <- sapply(strsplit(PG_noRQ$Majority.protein.IDs, ";"), "[", 1)



#### Make sure both tables have same length
PG <- subset(PG, Majority.protein.IDs %in% PG_noRQ$Majority.protein.IDs)
PG_noRQ <- subset(PG_noRQ, Majority.protein.IDs %in% PG$Majority.protein.IDs)





intensities <- names(PG)[grepl("Intensity.[L|M|H].", names(PG))]
IBAQ <- names(PG)[grepl("iBAQ.[L|M|H].", names(PG))]
ratios <- names(PG)[grepl(".normalized.", names(PG))]

normalized <- F
if (!normalized) {
  ratios <- sub("normalized.","", ratios)  # Remove normalized from names
  
}
rm(normalized)



#### Define Requantified ratios and intensities ####

### Ratios
requantified.ratios <- c()
unscrupulous.ratios <- c()
for (i in ratios) {
  PG[paste0("Requantified.", i)] <- ifelse(is.na(PG_noRQ[i]) &
                                             !is.na(PG[i]), T, F)
  
  requantified.ratios <- c(requantified.ratios, paste0("Requantified.", i))
  unscrupulous.ratios <- c(unscrupulous.ratios, paste0("Unscrupulous.", i))
  rm(i)
}



### intensities
requantified.intensities <- c()
for (i in intensities) {
  PG[paste0("Requantified.", i)] <- ifelse(PG_noRQ[i] == 0 &
                                             PG[i] != 0, T, F)
  
  requantified.intensities <- c(requantified.intensities, paste0("Requantified.", i))
  rm(i)
}



#### Unscrupulous Requantification ####
# Select Unscrupulous Requantified ratios (ratios between two channels when both are requantified)
#Forward
PG$Unscrupulous.Ratio.H.L.Forward_II <- ifelse((PG$Requantified.Intensity.H.Forward_II & PG$Requantified.Intensity.L.Forward_II) & !is.na(PG$Ratio.H.L.Forward_II), T, F)

PG$Unscrupulous.Ratio.H.L.Forward_III <- ifelse((PG$Requantified.Intensity.H.Forward_III & PG$Requantified.Intensity.L.Forward_III) & !is.na(PG$Ratio.H.L.Forward_III), T, F)

PG$Unscrupulous.Ratio.H.L.Forward_IV <- ifelse((PG$Requantified.Intensity.H.Forward_IV & PG$Requantified.Intensity.L.Forward_IV) & !is.na(PG$Ratio.H.L.Forward_IV), T, F)




PG$Unscrupulous.Ratio.H.M.Forward_II <- ifelse((PG$Requantified.Intensity.H.Forward_II & PG$Requantified.Intensity.M.Forward_II) & !is.na(PG$Ratio.H.M.Forward_II), T, F)

PG$Unscrupulous.Ratio.H.M.Forward_III <- ifelse((PG$Requantified.Intensity.H.Forward_III & PG$Requantified.Intensity.M.Forward_III) & !is.na(PG$Ratio.H.M.Forward_III), T, F)

PG$Unscrupulous.Ratio.H.M.Forward_IV <- ifelse((PG$Requantified.Intensity.H.Forward_IV & PG$Requantified.Intensity.M.Forward_IV) & !is.na(PG$Ratio.H.M.Forward_IV), T, F)



PG$Unscrupulous.Ratio.M.L.Forward_II <- ifelse((PG$Requantified.Intensity.M.Forward_II & PG$Requantified.Intensity.L.Forward_II) & !is.na(PG$Ratio.M.L.Forward_II), T, F)

PG$Unscrupulous.Ratio.M.L.Forward_III <- ifelse((PG$Requantified.Intensity.M.Forward_III & PG$Requantified.Intensity.L.Forward_III) & !is.na(PG$Ratio.M.L.Forward_III), T, F)

PG$Unscrupulous.Ratio.M.L.Forward_IV <- ifelse((PG$Requantified.Intensity.M.Forward_IV & PG$Requantified.Intensity.L.Forward_IV) & !is.na(PG$Ratio.M.L.Forward_IV), T, F)







#Reverse
PG$Unscrupulous.Ratio.H.L.Reverse_II <- ifelse((PG$Requantified.Intensity.H.Reverse_II & PG$Requantified.Intensity.L.Reverse_II) & !is.na(PG$Ratio.H.L.Reverse_II), T, F)

PG$Unscrupulous.Ratio.H.L.Reverse_III <- ifelse((PG$Requantified.Intensity.H.Reverse_III & PG$Requantified.Intensity.L.Reverse_III) & !is.na(PG$Ratio.H.L.Reverse_III), T, F)

PG$Unscrupulous.Ratio.H.L.Reverse_IV <- ifelse((PG$Requantified.Intensity.H.Reverse_IV & PG$Requantified.Intensity.L.Reverse_IV) & !is.na(PG$Ratio.H.L.Reverse_IV), T, F)



PG$Unscrupulous.Ratio.H.M.Reverse_II <- ifelse((PG$Requantified.Intensity.H.Reverse_II & PG$Requantified.Intensity.M.Reverse_II) & !is.na(PG$Ratio.H.M.Reverse_II), T, F)

PG$Unscrupulous.Ratio.H.M.Reverse_III <- ifelse((PG$Requantified.Intensity.H.Reverse_III & PG$Requantified.Intensity.M.Reverse_III) & !is.na(PG$Ratio.H.M.Reverse_III), T, F)

PG$Unscrupulous.Ratio.H.M.Reverse_IV <- ifelse((PG$Requantified.Intensity.H.Reverse_IV & PG$Requantified.Intensity.M.Reverse_IV) & !is.na(PG$Ratio.H.M.Reverse_IV), T, F)



PG$Unscrupulous.Ratio.M.L.Reverse_II <- ifelse((PG$Requantified.Intensity.M.Reverse_II & PG$Requantified.Intensity.L.Reverse_II) & !is.na(PG$Ratio.M.L.Reverse_II), T, F)

PG$Unscrupulous.Ratio.M.L.Reverse_III <- ifelse((PG$Requantified.Intensity.M.Reverse_III & PG$Requantified.Intensity.L.Reverse_III) & !is.na(PG$Ratio.M.L.Reverse_III), T, F)

PG$Unscrupulous.Ratio.M.L.Reverse_IV <- ifelse((PG$Requantified.Intensity.M.Reverse_IV & PG$Requantified.Intensity.L.Reverse_IV) & !is.na(PG$Ratio.M.L.Reverse_IV), T, F)








# Delete Unscrupulous ratios
for (i in 1:length(ratios)) {
  PG[,ratios[i]] <- as.numeric(ifelse(PG[,unscrupulous.ratios[i]], NA, PG[,ratios[i]]))
  
  rm(i)
}



rm(PG_noRQ)


#### Filters ####


#filter out contaminants, reverse and only identified by site
PG <- subset(PG, Reverse != "+")
PG <- subset(PG, Potential.contaminant != "+")
PG <- subset(PG, Only.identified.by.site != "+")



#transform intensities and ratios to Log2 (or 10 if you prefer)
PG[c(intensities, IBAQ, ratios)] = log2(PG[c(intensities, IBAQ, ratios)])
# change Inf values for na
is.na(PG[c(intensities, IBAQ, ratios)]) <- sapply(PG[c(intensities, IBAQ, ratios)], is.infinite)

is.na(PG[c(intensities, IBAQ, ratios)]) <- sapply(PG[c(intensities, IBAQ, ratios)], is.nan)






# how many proteins were identified (have intensities values not NA) in each L-H group pair
for (i in 1:length(intensities)) {
  cat(intensities[i])
  cat("\t")
  cat("\t")
  cat(nrow(PG[!is.na(PG[intensities[i]]),]))
  cat("\t")
  cat(signif(mean(PG[,intensities[i]], na.rm = T), 4))
  cat("\t")
  cat(sum(PG[,requantified.intensities[i]], na.rm = T))
  cat("\t")
  cat(signif(100*(sum(PG[,requantified.intensities[i]], na.rm = T)/nrow(PG[!is.na(PG[intensities[i]]),])),3))
  cat("\n")
  rm(i)
}


# how many proteins were identified (have ratios values not NA) in each L-H group pair
for (i in 1:length(ratios)) {
  #group
  cat(ratios[i])
  cat("\t")
  cat("\t")
  #quantified ratios
  cat(nrow(PG[!is.na(PG[ratios[i]]),]))
  cat("\t")
  #ratio mean
  cat(signif(mean(PG[,ratios[i]], na.rm = T), 3))
  cat("\t")
  #Requantified ratios
  cat(sum(!is.na(PG[ratios[i]]) & PG[, requantified.ratios[i]], na.rm = T))
  cat("\t")
  #Percentage of requantified ratios
  cat(100*(sum(!is.na(PG[ratios[i]]) & PG[, requantified.ratios[i]], na.rm = T)) / nrow(PG[!is.na(PG[ratios[i]]),]))
  cat("\n")
  rm(i)
}


# how many proteins were identified (have intensities values not NA) in each L-H group pair
for (i in 1:length(unscrupulous.ratios)) {
  cat(unscrupulous.ratios[i])
  cat("\t")
  cat("\t")
  cat(sum(PG[,unscrupulous.ratios[i]]))
  cat("\t")
  cat(100*(sum(PG[,unscrupulous.ratios[i]])/sum(!PG[,unscrupulous.ratios[i]])))
  cat("\n")
  rm(i)
}


#rm(requantified.intensities, requantified.ratios, unscrupulous.ratios)
#nrow(intersect(PG[!is.na(PG[intensities[8]]),], PG[!is.na(PG[intensities[9]]),]))




#### Defining known RBPs ####

## RBPbase Gebauer et al 2020
load("K:/Datasets/20210903_RBPbase_RIC_STUDIEs.Rda")
RIC_STUDIES <- subset(RIC_STUDIES, Organism == "Hs")

load("K:/Datasets/20210903_RBPbase_COMPILED_TABLE.Rda")
COMPILED_TABLE_Hs <- COMPILED_TABLE[["Hs"]]

COMPILED_TABLE_Hs <- subset(COMPILED_TABLE_Hs, select = c("UNIQUE", "ID","Description",
                                                          names(COMPILED_TABLE_Hs)[names(COMPILED_TABLE_Hs) %in%
                                                                                     RIC_STUDIES$RBPBASEID],
                                                          "any_Hs", "hits_Hs"))

# Defining RBPs as having been detected in at least 2 RIC-like experiments
mRBPs_HEK <- subset(COMPILED_TABLE_Hs, RBPBASE000000007.1)$UNIQUE # Baltz et all

mRBPs <- subset(COMPILED_TABLE_Hs, RBPBASE000000007.1 | # Baltz et all
                  RBPBASE000000008.1 | # Hs_Beckmann2015
                  RBPBASE000000009.1 | # Hs_Castello2012
                  RBPBASE000000010.1 | # Hs_Castello2016
                  RBPBASE000000013.1 | # Hs_Kramer2014
                  #RBPBASE000000032.1 | # Hs_Mullari2017 excluded because of lack of NoCl control
                  RBPBASE000000034.1 | # Hs_Perez-Perri2018-RIC
                  RBPBASE000000035.1 | # Hs_Perez-Perri2018-eRIC
                  RBPBASE000000067.1 | # Hs_Garcia-Moreno
                  RBPBASE000000011.1 | # Hs_Conrad2016_chr
                  RBPBASE000000038.1 | # Hs_Backlund2020-cyto
                  RBPBASE000000039.1 | # Hs_Backlund2020-cyto-ars
                  RBPBASE000000012.1 | # Hs_Conrad2016
                  RBPBASE000000040.1 | # Hs_Backlund2020-nuc
                  RBPBASE000000041.1   # Hs_Backlund2020-nuc-ars
)$UNIQUE

allRBPs <- subset(COMPILED_TABLE_Hs, hits_Hs >= 2)$UNIQUE
allRBPs <- union(allRBPs, mRBPs) # Add mRBPs that were only identified once


rm(COMPILED_TABLE, COMPILED_TABLE_Hs, RIC_STUDIES)

# Define known RBPs and mRBPs
PG$Known_RBP <- ifelse(PG$Gene.names %in% allRBPs, T, F)
PG$Known_mRBP <- ifelse(PG$Gene.names %in% mRBPs, T, F)
PG$Known_mRBP_HEK <- ifelse(PG$Gene.names %in% mRBPs_HEK, T, F)





#### Define PG_summ as working table ####

PG_summ <- subset(PG,
                  #Known_RBP,
                  select = c("id", "Fasta.headers","Majority.protein.IDs","Gene.names",
                             "Known_mRBP", "Known_RBP", "Known_mRBP_HEK",
                             "Mol..weight..kDa.",
                             "Sequence.coverage....","Sequence.length", 
                             "Peptides",
                             "Razor...unique.peptides",
                             "Unique.peptides",
                             ratios, intensities, IBAQ,
                             requantified.intensities, requantified.ratios,
                             unscrupulous.ratios))



#### Invert label swap ratios ####
PG_summ$Ratio.H.L.Reverse_II <- -PG_summ$Ratio.H.L.Reverse_II
PG_summ$Ratio.H.L.Reverse_III <- -PG_summ$Ratio.H.L.Reverse_III
PG_summ$Ratio.H.L.Reverse_IV <- -PG_summ$Ratio.H.L.Reverse_IV

PG_summ$Ratio.M.L.Reverse_II <- -PG_summ$Ratio.M.L.Reverse_II

PG_summ$Ratio.H.M.Reverse_II <- -PG_summ$Ratio.H.M.Reverse_II


#### Correct ratios to input fractions ####

# I originally tried to correct for losses in the Ethanol precipitation of the input samples (50% loss) by multiplying all ratios by 0.5 (sum log2(0.5)). I decided to not do it to reduce data manipulation and because I don't really know how much of input is lost in ethanol precipitation

# Replicate II
PG_summ$Ratio.H.L.Forward_II <- (PG_summ$Ratio.H.L.Forward_II)
PG_summ$Ratio.M.L.Forward_II <- (PG_summ$Ratio.M.L.Forward_II) + log2(2) #I used half of Pulldown for no-crosslink samples

PG_summ$Ratio.H.L.Reverse_II <- (PG_summ$Ratio.H.L.Reverse_II)
PG_summ$Ratio.H.M.Reverse_II <- (PG_summ$Ratio.H.M.Reverse_II) + log2(2) #I used half of Pulldown for no-crosslink samples


# Replicate III
PG_summ$Ratio.H.L.Forward_III <- (PG_summ$Ratio.H.L.Forward_III)
PG_summ$Ratio.H.L.Reverse_III <- (PG_summ$Ratio.H.L.Reverse_III)



# Replicate IV
PG_summ$Ratio.H.L.Forward_IV <- (PG_summ$Ratio.H.L.Forward_IV)
PG_summ$Ratio.H.L.Reverse_IV <- (PG_summ$Ratio.H.L.Reverse_IV)





#### Define RBPs from NoCL ####
PG_summ$qRIC_mRBP <- ifelse((PG_summ$Ratio.H.M.Forward_II >= 1 &
                               PG_summ$Ratio.M.L.Reverse_II >= 1) &
                              (!is.na(PG_summ$Ratio.H.M.Forward_II) |
                                 !is.na(PG_summ$Ratio.M.L.Reverse_II)), T, F)

qRIC_mRBP <- subset(PG_summ, qRIC_mRBP)$Gene.names




#### Ratios into Percentages ####
# Replicate II
PG_summ$Percentage.H.L.Forward_II <- (2^PG_summ$Ratio.H.L.Forward_II)
PG_summ$Percentage.M.L.Forward_II <- (2^PG_summ$Ratio.M.L.Forward_II)

PG_summ$Percentage.H.L.Reverse_II <- (2^PG_summ$Ratio.H.L.Reverse_II)
PG_summ$Percentage.H.M.Reverse_II <- (2^PG_summ$Ratio.H.M.Reverse_II)


# Replicate III
PG_summ$Percentage.H.L.Forward_III <- (2^PG_summ$Ratio.H.L.Forward_III)
PG_summ$Percentage.H.L.Reverse_III <- (2^PG_summ$Ratio.H.L.Reverse_III)



# Replicate IV
PG_summ$Percentage.H.L.Forward_IV <- (2^PG_summ$Ratio.H.L.Forward_IV)
PG_summ$Percentage.H.L.Reverse_IV <- (2^PG_summ$Ratio.H.L.Reverse_IV)


percentages <- c("Percentage.H.L.Forward_II",
                 "Percentage.H.L.Forward_III",
                 "Percentage.H.L.Forward_IV",
                 "Percentage.H.L.Reverse_II",
                 "Percentage.H.L.Reverse_III",
                 "Percentage.H.L.Reverse_IV",
                 "Percentage.M.L.Forward_II",
                 "Percentage.H.M.Reverse_II")












#### Calculations ####
#Calculate mean making sure protein is identified in at least >1 replicates
MyMeanWithMinValues <- function(df, columns, MinNumberReplicates = 2) {
  foo <- apply(df[columns],
               1, 
               function(x) ifelse(sum(!is.na(x)) > MinNumberReplicates - 1,
                                  mean(x, na.rm = T),
                                  NA))
  
  return(foo)
}


# Ratios
PG_summ$Ratio.H.L.Forward <- MyMeanWithMinValues(PG_summ, c("Ratio.H.L.Forward_II",
                                                            "Ratio.H.L.Forward_III",
                                                            "Ratio.H.L.Forward_IV"),
                                                 MinNumberReplicates = 1)

PG_summ$Ratio.H.L.Reverse <- MyMeanWithMinValues(PG_summ, c("Ratio.H.L.Reverse_II",
                                                            "Ratio.H.L.Reverse_III",
                                                            "Ratio.H.L.Reverse_IV"),
                                                 MinNumberReplicates = 1)

PG_summ$Ratio.H.L.Mean <- MyMeanWithMinValues(PG_summ, c("Ratio.H.L.Forward",
                                                         "Ratio.H.L.Reverse"),
                                              MinNumberReplicates = 2)



PG_summ$Ratio.NoCl.Mean <- MyMeanWithMinValues(PG_summ, c("Ratio.M.L.Forward_II",
                                                          "Ratio.H.M.Reverse_II"),
                                               MinNumberReplicates = 2)


PG_summ$Ratio.ClXNoCl.Mean <- MyMeanWithMinValues(PG_summ, c("Ratio.H.M.Forward_II",
                                                             "Ratio.M.L.Reverse_II"),
                                                  MinNumberReplicates = 2)





#Intensities
PG_summ$Intensity.Input.Mean <- MyMeanWithMinValues(PG_summ, c("Intensity.L.Forward_II",
                                                               "Intensity.L.Forward_III",
                                                               "Intensity.L.Forward_IV",
                                                               "Intensity.H.Reverse_II",
                                                               "Intensity.H.Reverse_III",
                                                               "Intensity.H.Reverse_IV"),
                                                    MinNumberReplicates = 1)






#IBAQ
PG_summ$iBAQ.Input.Mean <- MyMeanWithMinValues(PG_summ, c("iBAQ.L.Forward_II",
                                                          "iBAQ.L.Forward_III",
                                                          "iBAQ.L.Forward_IV",
                                                          "iBAQ.H.Reverse_II",
                                                          "iBAQ.H.Reverse_III",
                                                          "iBAQ.H.Reverse_IV"),
                                               MinNumberReplicates = 1)



# Percentages
PG_summ$Percentage.H.L.Forward <- 2**PG_summ$Ratio.H.L.Forward

PG_summ$Percentage.H.L.Reverse <- 2**PG_summ$Ratio.H.L.Reverse

PG_summ$Percentage_Mean <- 2**PG_summ$Ratio.H.L.Mean

PG_summ$Percentage_Mean_NoCl <- 2**PG_summ$Ratio.NoCl.Mean




#### Reproducibility ####
PG_summ$Reproducible.logical_II <- ifelse((PG_summ$Ratio.H.L.Forward_II - PG_summ$Ratio.H.L.Reverse_II < ReproducibilityDifference) & (PG_summ$Ratio.H.L.Reverse_II - PG_summ$Ratio.H.L.Forward_II < ReproducibilityDifference), T, F)

PG_summ$Reproducible.logical_III <- ifelse((PG_summ$Ratio.H.L.Forward_III - PG_summ$Ratio.H.L.Reverse_III < ReproducibilityDifference) & (PG_summ$Ratio.H.L.Reverse_III - PG_summ$Ratio.H.L.Forward_III < ReproducibilityDifference), T, F)

PG_summ$Reproducible.logical_IV <- ifelse((PG_summ$Ratio.H.L.Forward_IV - PG_summ$Ratio.H.L.Reverse_IV < ReproducibilityDifference) & (PG_summ$Ratio.H.L.Reverse_IV - PG_summ$Ratio.H.L.Forward_IV < ReproducibilityDifference), T, F)




#### Write table ####
write.table(PG_summ, file = "PG_summary.txt", sep = "\t", na = "", quote = F, row.names = F)

#### Writing table ratios paper ####
## Do not further filter this table. It must be exactly the same as PG_summary.txt but with selected columns only, except for inverted reverse ratios
# Table 1 paper
foo_df <- subset(PG_summ,
                 select = c("Fasta.headers",
                            "Majority.protein.IDs",
                            "Gene.names",
                            
                            "Known_mRBP",
                            "Known_mRBP_HEK",
                            
                            "Peptides",
                            "Razor...unique.peptides",
                            "Unique.peptides",
                            
                            "Ratio.H.L.Forward_II",
                            "Ratio.H.L.Forward_III",
                            "Ratio.H.L.Forward_IV",
                            "Ratio.H.L.Reverse_II",
                            "Ratio.H.L.Reverse_III",
                            "Ratio.H.L.Reverse_IV",
                            "Ratio.H.M.Forward_II",
                            "Ratio.H.M.Reverse_II",
                            "Ratio.M.L.Forward_II",
                            "Ratio.M.L.Reverse_II",
                            
                            "Ratio.H.L.Forward",
                            "Ratio.H.L.Reverse",
                            
                            "Ratio.H.L.Mean",
                            
                            "Percentage.H.L.Forward",
                            "Percentage.H.L.Reverse",
                            
                            "Percentage_Mean"))


#Invert back ratios
foo_df[,c("Ratio.H.L.Reverse_II",
          "Ratio.H.L.Reverse_III",
          "Ratio.H.L.Reverse_IV",
          "Ratio.H.M.Reverse_II",
          "Ratio.M.L.Reverse_II")] <- -foo_df[,c("Ratio.H.L.Reverse_II",
                                                 "Ratio.H.L.Reverse_III",
                                                 "Ratio.H.L.Reverse_IV",
                                                 "Ratio.H.M.Reverse_II",
                                                 "Ratio.M.L.Reverse_II")]



# Adjust colum names
names(foo_df) <- c("Fasta.headers",
                   "Uniprot.IDs",
                   "Gene.names",
                   
                   "mRBP.(Gebauer)",
                   "HEK.mRBP.(Baltz)",
                   
                   "Peptides",
                   "Razor.unique.peptides",
                   "Unique.peptides",
                   
                   "log2.Ratio.H.L.Forward_Exp1",
                   "log2.Ratio.H.L.Forward_Exp2",
                   "log2.Ratio.H.L.Forward_Exp3",
                   "log2.Ratio.H.L.Reverse_Exp1",
                   "log2.Ratio.H.L.Reverse_Exp2",
                   "log2.Ratio.H.L.Reverse_Exp3",
                   "log2.Ratio.H.M.Forward_Exp1",
                   "log2.Ratio.H.M.Reverse_Exp1",
                   "log2.Ratio.M.L.Forward_Exp1",
                   "log2.Ratio.M.L.Reverse_Exp1",
                   
                   "log2.Ratio.H.L.Forward",
                   "log2.Ratio.H.L.Reverse",
                   
                   "log2.Ratio.H.L.Mean",
                   
                   "Pulldown.efficiency.Forward(%)",
                   "Pulldown.efficiency.Reverse(%)",
                   
                   "Pulldown.efficiency(%)")


write.table(foo_df, file = "Table1_ProteinEfficiency.txt", sep = "\t", na = "", quote = F, row.names = F)


#### Scatter plot function ####
MyScatterPlot <- function(df, X, Y,
                          xmin = -6, xmax = 6,
                          AxisSep.x = 2,
                          ymin = -6, ymax = 6,
                          AxisSep.y = 2,
                          bin.size = 128,
                          Title = NULL, xTitle = NULL, yTitle = NULL,
                          InPercentage.x = F, sec_breaks.x = c(0.1,0.5,2,8,32),
                          scaled_breaks.x = log2(sec_breaks.x),
                          InPercentage.y = F, sec_breaks.y = sec_breaks.x,
                          scaled_breaks.y = log2(sec_breaks.y),
                          ShowRBPs = F, ShowNotReproducible = F,
                          GradientColorMin = rgb(.5,.5,.5), GradientColorMax = rgb(0,0,0),
                          alpha_points = 1, size_points = 2.5) {
  
  plot <- ggplot(data = df, aes(x = eval(parse(text = X)), y = eval(parse(text = Y)))) +
    
    geom_point(na.rm = T,
               size = size_points,
               alpha = alpha_points,
               colour = densCols(x = df[,X],
                                 y = df[,Y],
                                 colramp = colorRampPalette(c(GradientColorMin,GradientColorMax)), nbin = bin.size)) +
    
    ggtitle(Title) +
    xlab(xTitle) +
    ylab(yTitle) +
    coord_cartesian(xlim = c(xmin,xmax),ylim = c(ymin,ymax)) +
    theme_bw() +
    theme(panel.grid.minor = element_blank(),
          plot.title = element_text(lineheight=.8,
                                    #face="bold", 
                                    vjust=0.5,
                                    hjust = 0.5,
                                    size=30),
          axis.title.x = element_text(#face="bold",
            size=25,
            hjust = 0.5,
            vjust = 1.5),
          axis.text.x  = element_text(#face = "bold",
            color = "black",
            angle=0, 
            vjust=1,
            size=30),
          axis.title.y = element_text(#face="bold",
            size=25,
            hjust = 0.5,
            vjust = 1.5),
          axis.text.y  = element_text(#face = "bold",
            color = "black",
            angle=0, 
            vjust=.4,
            size=30),
          #panel.grid=element_blank(),
          aspect.ratio=1)
  
  if (InPercentage.x) {
    
    plot <- plot +
      scale_x_continuous(breaks = scaled_breaks.x, labels = sec_breaks.x)
    
  } else {
    
    plot <- plot +
      scale_x_continuous(breaks=seq(xmin,xmax, AxisSep.x))
    
  }
  
  
  if (InPercentage.y) {
    
    plot <- plot +
      scale_y_continuous(breaks = scaled_breaks.y, labels = sec_breaks.y)
    
  } else {
    
    plot <- plot +
      scale_y_continuous(breaks=seq(ymin,ymax, AxisSep.y))
    
  }
  
  
  
  if(ShowNotReproducible) {
    
    plot <- plot +
      
      geom_point(data = subset(df, !Reproducible.logical),
                 na.rm = T,
                 size = size_points,
                 colour = rgb(.8,0,.8))
  }
  
  if(ShowRBPs) {
    
    plot <- plot +
      
      geom_point(data = subset(df, Known_mRBP),
                 na.rm = T,
                 shape = 21,
                 fill = NA,
                 size = size_points,
                 colour = rgb(.8,0,0))
    
    # stat_density_2d(aes(colour = Known_mRBP, alpha=..level..), size = 5, na.rm = T, bins = 5, show.legend = F)
  }
  
  return(plot)
}

MyScatterPlotPaper <- function(df = PG_summ, X = "Ratio.H.L.Forward", Y = "Ratio.H.L.Reverse",
                               xmin = -5, xmax = 6,
                               AxisSep.x = 2,
                               ymin = 0, ymax = 4, 
                               AxisSep.y = 1,
                               bin.size = 128,
                               Title = NULL, xTitle = NULL, yTitle = NULL,
                               InPercentage.x = F, sec_breaks.x = c(0.125,0.5,2,8,32),
                               scaled_breaks.x = log2(sec_breaks.x),
                               InPercentage.y = F, sec_breaks.y = 2**seq(ymin,ymax,1),
                               scaled_breaks.y = log2(sec_breaks.y),
                               ShowRBPs = F, ShowNotReproducible = F,
                               GradientColorMin = rgb(.5,.5,.5), GradientColorMax = rgb(0,0,0),
                               alpha_points = 1, size_points = 2.5) {
  
  plot <- ggplot(data = df, aes(x = eval(parse(text = X)), y = eval(parse(text = Y)))) +
    
    geom_point(na.rm = T,
               size = size_points,
               alpha = alpha_points,
               colour = densCols(x = df[,X],
                                 y = df[,Y],
                                 colramp = colorRampPalette(c(GradientColorMin,GradientColorMax)), nbin = bin.size)) +
    
    ggtitle(Title) +
    xlab(xTitle) +
    ylab(yTitle) +
    coord_cartesian(xlim = c(xmin,xmax),ylim = c(ymin,ymax)) +
    theme_bw() +
    theme(panel.grid.major = element_line(colour = rgb(0,0,0,.1), linetype = "dashed"),
          plot.title = element_text( size = 14,
                                     hjust = 0.5,
                                     vjust = 0.4),
          axis.title.x = element_text(size=14,
                                      hjust = 0.5,
                                      vjust = 0.4),
          axis.title.y = element_text(size=14,
                                      hjust = 0.5,
                                      vjust = 1.5),
          axis.text.x  = element_text(color = "black",
                                      angle=0, 
                                      vjust=0.5,
                                      hjust = 0.5,
                                      size=12),
          axis.text.y  = element_text(color = "black",
                                      angle=0, 
                                      vjust=0.5,
                                      hjust = 0.5,
                                      size=12),
          #panel.grid=element_blank(),
          aspect.ratio=1)
  
  
  
  
  
  if (InPercentage.x) {
    
    plot <- plot +
      scale_x_continuous(breaks = scaled_breaks.x, labels = sec_breaks.x, sec.axis = sec_axis(~ ., breaks = log2(sec_breaks.x)))
    
  } else {
    
    plot <- plot +
      scale_x_continuous(breaks=seq(xmin,xmax, AxisSep.x))
    
  }
  
  
  if (InPercentage.y) {
    
    plot <- plot +
      scale_y_continuous(breaks = scaled_breaks.y, labels = sec_breaks.y, sec.axis = sec_axis(~ ., breaks = log2(sec_breaks.y)))
    
  } else {
    
    plot <- plot +
      scale_y_continuous(breaks=seq(ymin,ymax, AxisSep.y))
    
  }
  
  
  
  if(ShowNotReproducible) {
    
    plot <- plot +
      
      geom_point(data = subset(df, !Reproducible.logical),
                 na.rm = T,
                 size = 5,
                 colour = rgb(.8,0,.8))
  }
  
  if(ShowRBPs) {
    
    plot <- plot +
      
      stat_density_2d(aes(colour = qRIC_mRBP, alpha=..level..),
                      size = 5, na.rm = T, bins = 5, show.legend = F)
  }
  
  return(plot)
}


#### Selected proteins ####

selected_protein <- c("RBM20", "SF3B1", "UPF1", "ELAVL1")
#selected_protein <- unique(subset(PG_summ, grepl("EIF", PG_summ$Gene.names))$Majority.protein.IDs)
#selected_protein <- unique(subset(PG_summ, Intensity.Input.Mean > 31 & Ratio.H.L.Mean > 1)$Majority.protein.IDs)
#selected_protein <- c("blank")

foo <- subset(PG_summ, Gene.names %in% selected_protein)


ratios_CL <- MyScatterPlot(PG_summ, xmin = -7,
                           ymin = 23, ymax = 36,
                           "Ratio.H.L.Mean", "Intensity.Input.Mean",
                           ShowRBPs = F, InPercentage.x = T, InPercentage.y = F) +
  
  geom_point(data = foo,
             na.rm = T,
             size = 5,
             colour = rgb(0,.4,.8)) +
  
  
  geom_text_repel(data = foo,
                  aes(label = Gene.names),
                  size = 7,
                  colour = rgb(0,0,1),
                  na.rm = T)

#png(filename = paste0("Protein_ratios_", paste(selected_protein, collapse = "_"), ".png"), height = 500, width = 500)
png(filename = paste0("Protein_ratios_", paste(selected_protein, collapse = "_"), ".png"), height = 500, width = 500)
grid.arrange(ratios_CL, 
             ncol = 1)
dev.off()
rm(foo, ratios_CL)







#### boxplot intensities ####
MyBoxplot <- function(df = PG_summ, variables = intensities,
                      AxisName_y = ("Log2"), LabelNames = variables,
                      limits_y_min = round(min(df[variables], na.rm = T))-1,
                      limits_y_max = round(max(df[variables], na.rm = T))+1,
                      limits_breaks = round((limits_y_max+limits_y_min)/20)) {
  
  all <- melt(df,
              id.vars = c("Majority.protein.IDs"), measure.vars=variables)
  
  
  
  bplot <- ggplot(all, aes(variable,value)) +
    
    geom_violin(fill = rgb(0,0,0,.25), na.rm = T) +
    
    geom_boxplot(fill = rgb(0,0,0,0), width = 0.75, na.rm = T, notch = T) +
    
    ylab(AxisName_y) +
    coord_flip(ylim = c(limits_y_min, limits_y_max)) +
    scale_y_continuous(breaks = seq(limits_y_min, limits_y_max, limits_breaks)) +
    scale_x_discrete(labels = LabelNames) +
    
    theme_bw() +
    theme(panel.grid = element_blank(),
          plot.title = element_text(lineheight=.8,
                                    face="bold", 
                                    vjust=0.5,
                                    hjust = 0.5,
                                    size=30),
          axis.title.y = element_blank(),
          axis.text.x  = element_text(face = "bold",
                                      color = "black",
                                      angle=0, 
                                      vjust=1,
                                      size=20),
          axis.title.x = element_text(face="bold",
                                      size=25,
                                      hjust = 0.5,
                                      vjust = 1.5),
          axis.text.y  = element_text(face = "bold",
                                      color = "black",
                                      angle=0, 
                                      vjust=1,
                                      size=20))
  
  return(bplot)
  
}




boxplot <- MyBoxplot(df = PG_summ, intensities, AxisName_y = "Intensities (Log2)",
                     LabelNames = intensities)






png(paste("Protein_Boxplot.png", sep = ""), width = 1000, height = 1000, pointsize = 25)
grid.arrange(boxplot, nrow = 1)
dev.off()
rm(boxplot)






#### Forward X Reverse ratios ####

#CL
ratios_CL_II <- MyScatterPlot(PG_summ,
                              "Ratio.H.L.Forward_II", "Ratio.H.L.Reverse_II",
                              ShowRBPs = T, InPercentage.x = T, InPercentage.y = T) + geom_abline(slope = 1)

ratios_CL_III <- MyScatterPlot(PG_summ,
                               "Ratio.H.L.Forward_III", "Ratio.H.L.Reverse_III",
                               ShowRBPs = T, InPercentage.x = T, InPercentage.y = T) + geom_abline(slope = 1)

ratios_CL_IV <- MyScatterPlot(PG_summ,
                              "Ratio.H.L.Forward_IV", "Ratio.H.L.Reverse_IV",
                              ShowRBPs = T, InPercentage.x = T, InPercentage.y = T) + geom_abline(slope = 1)



ratios_CL <- MyScatterPlot(PG_summ,
                           "Ratio.H.L.Forward", "Ratio.H.L.Reverse",
                           ShowRBPs = T, InPercentage.x = T, InPercentage.y = T)




#NoCL
ratios_NoCL_II <- MyScatterPlot(PG_summ,
                                "Ratio.M.L.Forward_II", "Ratio.H.M.Reverse_II",
                                ShowRBPs = T, InPercentage.x = T, InPercentage.y = T) + geom_abline(slope = 1)






blank <- ggplot() + theme(panel.background = element_blank(),
                          panel.border = element_blank())



png(filename = "Protein_ratios.png", height = 1000, width = 2000)
grid.arrange(ratios_CL_II, ratios_CL_III, ratios_CL_IV, ratios_CL,
             ratios_NoCL_II, blank, blank, blank,
             ncol = 4)
dev.off()
rm(ratios_CL_II, ratios_CL_III, ratios_CL_IV, ratios_CL,
   ratios_NoCL_II, blank)





#### Figure paper ####
MyScatterPlotPaper <- function(df = PG_summ, X = "Ratio.H.L.Forward", Y = "Ratio.H.L.Reverse",
                          xmin = -5, xmax = 6,
                          AxisSep.x = 2,
                          ymin = xmin, ymax = xmax,
                          AxisSep.y = AxisSep.x,
                          bin.size = 128,
                          Title = NULL, xTitle = NULL, yTitle = NULL,
                          InPercentage.x = F, sec_breaks.x = c(0.125,0.5,2,8,32),
                          scaled_breaks.x = log2(sec_breaks.x),
                          InPercentage.y = F, sec_breaks.y = sec_breaks.x,
                          scaled_breaks.y = log2(sec_breaks.y),
                          ShowRBPs = F, ShowNotReproducible = F,
                          GradientColorMin = rgb(.5,.5,.5), GradientColorMax = rgb(0,0,0),
                          alpha_points = 1, size_points = 3) {
  
  plot <- ggplot(data = df, aes(x = eval(parse(text = X)), y = eval(parse(text = Y)))) +
    
    geom_point(na.rm = T,
               size = size_points,
               alpha = alpha_points,
               colour = densCols(x = df[,X],
                                 y = df[,Y],
                                 colramp = colorRampPalette(c(GradientColorMin,GradientColorMax)), nbin = bin.size)) +
    
    ggtitle(Title) +
    xlab(xTitle) +
    ylab(yTitle) +
    coord_cartesian(xlim = c(xmin,xmax),ylim = c(ymin,ymax)) +
    theme_bw() +
    theme(panel.grid.minor = element_blank(),
          plot.title = element_text(lineheight=.8,
                                    #face="bold", 
                                    vjust=0.5,
                                    hjust = 0.5,
                                    size=30),
          axis.title.x = element_text(#face="bold",
            size=25,
            hjust = 0.5,
            vjust = 1.5),
          axis.text.x  = element_text(#face = "bold",
            color = "black",
            angle=0, 
            vjust=1,
            size=30),
          axis.title.y = element_text(#face="bold",
            size=25,
            hjust = 0.5,
            vjust = 1.5),
          axis.text.y  = element_text(#face = "bold",
            color = "black",
            angle=0, 
            vjust=.4,
            size=30),
          #panel.grid=element_blank(),
          aspect.ratio=1) +
    
    
    geom_hline(yintercept = DeltaAnalysisThresholdEfficiency,linetype = "dashed", color = rgb(1,0,0)) +
    
    geom_vline(xintercept = DeltaAnalysisThresholdEfficiency, linetype = "dashed", color = rgb(1,0,0))
  
    
    
    
    
  if (InPercentage.x) {
    
    plot <- plot +
      scale_x_continuous(breaks = scaled_breaks.x, labels = sec_breaks.x, sec.axis = sec_axis(~ ., breaks = log2(sec_breaks.x)))
    
  } else {
    
    plot <- plot +
      scale_x_continuous(breaks=seq(xmin,xmax, AxisSep.x))
    
  }
  
  
  if (InPercentage.y) {
    
    plot <- plot +
      scale_y_continuous(breaks = scaled_breaks.y, labels = sec_breaks.y, sec.axis = sec_axis(~ ., breaks = log2(sec_breaks.x)))
    
  } else {
    
    plot <- plot +
      scale_y_continuous(breaks=seq(ymin,ymax, AxisSep.y))
    
  }
  
  
  
  if(ShowNotReproducible) {
    
    plot <- plot +
      
      geom_point(data = subset(df, !Reproducible.logical),
                 na.rm = T,
                 size = 5,
                 colour = rgb(.8,0,.8))
  }
  
  if(ShowRBPs) {
    
    plot <- plot +
      
      stat_density_2d(aes(colour = Known_mRBP, alpha=..level..),
                      size = 5, na.rm = T, bins = 5, show.legend = F)
  }
  
  return(plot)
}



selected_protein <- c("RBM20", "SF3B1", "UPF1", "ELAVL1")
foo <- subset(PG_summ, Gene.names %in% selected_protein)


ratios_CL <- MyScatterPlotPaper(PG_summ,
                           "Ratio.H.L.Forward", "Ratio.H.L.Reverse",
                           ShowRBPs = F, InPercentage.x = T, InPercentage.y = T) +
  
  geom_point(data = foo,
             na.rm = T,
             size = 3,
             colour = rgb(0,.2,.8)) +
  
  
  geom_text_repel(data = foo,
                  aes(label = Gene.names),
                  size = 7,
                  colour = rgb(0,.2,.8),
                  na.rm = T)


ratios_NoCL <- MyScatterPlotPaper(PG_summ,
                             "Ratio.M.L.Forward_II", "Ratio.H.M.Reverse_II",
                             ShowRBPs = F, InPercentage.x = T, InPercentage.y = T) +
  
  geom_point(data = foo,
             na.rm = T,
             size = 3,
             colour = rgb(0,.2,.8)) +
  
  
  geom_text_repel(data = foo,
                  aes(label = Gene.names),
                  size = 7,
                  colour = rgb(0,.2,.8),
                  na.rm = T)

pdf(file = "Protein_ratios.pdf", height = 5, width = 10, useDingbats=FALSE)
grid.arrange(ratios_CL, ratios_NoCL,
             ncol = 2)
dev.off()
rm(ratios_CL, ratios_NoCL)





#### Figure paper Xl and no Xl controls ####
selected_protein <- c("HNRNPD", "ACTN1")
foo <- subset(PG_summ, Gene.names %in% selected_protein)


ratios_CL_II <- MyScatterPlotPaper(PG_summ,
                                   "Ratio.H.L.Forward_II", "Ratio.H.L.Reverse_II",
                                   ShowRBPs = F, InPercentage.x = T, InPercentage.y = T) +
  
  geom_point(data = foo,
             na.rm = T,
             size = 3,
             colour = rgb(0,.2,.8)) +
  
  geom_text_repel(data = foo,
                  aes(label = Gene.names),
                  size = 7,
                  colour = rgb(0,.2,.8),
                  na.rm = T)


ratios_CL_III <- MyScatterPlotPaper(PG_summ,
                                   "Ratio.H.L.Forward_III", "Ratio.H.L.Reverse_III",
                                   ShowRBPs = F, InPercentage.x = T, InPercentage.y = T) +
  
  geom_point(data = foo,
             na.rm = T,
             size = 3,
             colour = rgb(0,.2,.8)) +
  
  
  geom_text_repel(data = foo,
                  aes(label = Gene.names),
                  size = 7,
                  colour = rgb(0,.2,.8),
                  na.rm = T)


ratios_CL_IV <- MyScatterPlotPaper(PG_summ,
                                   "Ratio.H.L.Forward_IV", "Ratio.H.L.Reverse_IV",
                                   ShowRBPs = F, InPercentage.x = T, InPercentage.y = T) +
  
  geom_point(data = foo,
             na.rm = T,
             size = 3,
             colour = rgb(0,.2,.8)) +
  
  
  geom_text_repel(data = foo,
                  aes(label = Gene.names),
                  size = 7,
                  colour = rgb(0,.2,.8),
                  na.rm = T)


ratios_NoCL_II <- MyScatterPlotPaper(PG_summ,
                                  "Ratio.M.L.Forward_II", "Ratio.H.M.Reverse_II",
                                  ShowRBPs = F, InPercentage.x = T, InPercentage.y = T,
                                  size_points = 3) +
  
  geom_point(data = foo,
             na.rm = T,
             size = 3,
             colour = rgb(0,.2,.8)) +
  
  
  geom_text_repel(data = foo,
                  aes(label = Gene.names),
                  size = 7,
                  colour = rgb(0,.2,.8),
                  na.rm = T)



pdf(file = "Protein_ratios_repII.pdf", height = 10, width = 10, useDingbats=FALSE)
grid.arrange(ratios_CL_II, ratios_NoCL_II,
             ratios_CL_III, ratios_CL_IV,
             ncol = 2)
dev.off()
rm(ratios_CL_II, ratios_NoCL_II,
   ratios_CL_III, ratios_CL_IV)




#### Figure paper Ratios RBPs only ####


ratios_CL_RBPs <- MyScatterPlotPaper(subset(PG_summ, Known_mRBP),
                                     "Ratio.H.L.Forward", "Ratio.H.L.Reverse",
                                     ShowRBPs = F, InPercentage.x = T, InPercentage.y = T)

ratios_CL_notRBPs <- MyScatterPlotPaper(subset(PG_summ, !Known_mRBP),
                                        "Ratio.H.L.Forward", "Ratio.H.L.Reverse",
                                        ShowRBPs = F, InPercentage.x = T, InPercentage.y = T)

pdf(file = "Protein_ratios_mRBPs.pdf", height = 5, width = 10, useDingbats = F)
grid.arrange(ratios_CL_RBPs, ratios_CL_notRBPs,
             nrow = 1)
dev.off()
rm(ratios_CL_RBPs, ratios_CL_notRBPs)



#### Figure paper reproducibility ####
## All proteins
pdf(file = "Protein_Pairs_ratios.pdf", height = 5, width = 5, useDingbats=FALSE)
ggpairs(PG_summ[c("Ratio.H.L.Forward_II",
                  "Ratio.H.L.Forward_III",
                  "Ratio.H.L.Forward_IV",
                  "Ratio.H.L.Reverse_II",
                  "Ratio.H.L.Reverse_III",
                  "Ratio.H.L.Reverse_IV",
                  "Ratio.M.L.Forward_II",
                  "Ratio.H.M.Reverse_II")], 
        mapping = ggplot2::aes(alpha = 0.2),
        lower = list(continuous = wrap("points", size=0.1)),
        diag = list(continuous = "blankDiag")) +
  
  coord_cartesian(xlim = c(-8,10),ylim = c(-8,10))
dev.off()

## mRBPs only
pdf(file = "Protein_Pairs_ratios_mRBPs.pdf", height = 5, width = 5, useDingbats=FALSE)
ggpairs(subset(PG_summ, Known_mRBP)[c("Ratio.H.L.Forward_II",
                                      "Ratio.H.L.Forward_III",
                                      "Ratio.H.L.Forward_IV",
                                      "Ratio.H.L.Reverse_II",
                                      "Ratio.H.L.Reverse_III",
                                      "Ratio.H.L.Reverse_IV",
                                      "Ratio.M.L.Forward_II",
                                      "Ratio.H.M.Reverse_II")], 
        mapping = ggplot2::aes(alpha = 0.2),
        lower = list(continuous = wrap("points", size=0.1)),
        diag = list(continuous = "blankDiag")) +
  
  coord_cartesian(xlim = c(-8,10),ylim = c(-8,10))
dev.off()


#### Figure paper Histogram Ratios mRBPs ####
MyHistPlotPaper <- function(df = PG_summ,
                            X1 = "Ratio.H.L.Mean",
                            X2 = "Ratio.NoCl.Mean",
                            xmin = -5, xmax = 6,
                            AxisSep.x = 2,
                            ymin = 0, ymax = 0.5,
                            AxisSep.y = 0.1,
                            Title = NULL, xTitle = "Pull-down efficiency (%)", yTitle = "Density",
                            InPercentage.x = T, sec_breaks.x = c(0.125,0.5,2,8,32),
                            scaled_breaks.x = log2(sec_breaks.x)) {
  
  plot <- ggplot(data = df) +
    
    geom_density(aes(y = ..density.., x = eval(parse(text = X1))), na.rm = T, color = rgb(1,0,0)) +
    
    geom_density(aes(y = ..density.., x = eval(parse(text = X2))), na.rm = T) +
    
    ggtitle(Title) +
    xlab(xTitle) +
    ylab(yTitle) +
    coord_cartesian(xlim = c(xmin,xmax),ylim = c(ymin,ymax)) +
    scale_y_continuous(breaks=seq(ymin,ymax, AxisSep.y)) +
    theme_bw() +
    theme(panel.grid.major = element_line(colour = rgb(0,0,0,.1), linetype = "dashed"),
          panel.grid.minor = element_blank(),
          plot.title = element_text(lineheight=.8,
                                    vjust=0.5,
                                    hjust = 0.5,
                                    size=30),
          axis.title.x = element_text(size=25,
                                      hjust = 0.5,
                                      vjust = 1.5),
          axis.text.x  = element_text(color = "black",
                                      angle=0, 
                                      vjust=1,
                                      size=30),
          axis.title.y = element_text(size=25,
                                      hjust = 0.5,
                                      vjust = 1.5),
          axis.text.y  = element_text(color = "black",
                                      angle=0, 
                                      vjust=.4,
                                      size=30),
          panel.grid=element_blank(),
          aspect.ratio=1)
  
  
  
  
  
  if (InPercentage.x) {
    
    plot <- plot +
      scale_x_continuous(breaks = scaled_breaks.x, labels = sec_breaks.x, sec.axis = sec_axis(~ ., breaks = log2(sec_breaks.x)))
    
  } else {
    
    plot <- plot +
      scale_x_continuous(breaks=seq(xmin,xmax, AxisSep.x))
    
  }
  
  return(plot)
}

ratios_All_proteins <- MyHistPlotPaper(PG_summ, Title = "All proteins")
ratios_mRBPs <- MyHistPlotPaper(subset(PG_summ, Known_mRBP), Title = "mRBPs only")
ratios_nonmRBPs <- MyHistPlotPaper(subset(PG_summ, !Known_mRBP), Title = "Non mRBPs")

pdf(file = "Protein_Hist_mRBPs.pdf", height = 5, width = 15, useDingbats = F)
grid.arrange(ratios_All_proteins, ratios_mRBPs, ratios_nonmRBPs,
             nrow = 1)
dev.off()
rm(ratios_All_proteins, ratios_mRBPs, ratios_nonmRBPs)



#### Calculated RBP identification FDR ####
MySpecificFDR <- function(df = PG_summ, X = X, reference = reference, cutoff_specific = cutoff_specific) {
  
  Positives <- subset(df, df[,X] >= cutoff_specific)
  Negatives <- subset(df, df[,X] < cutoff_specific)
  
  N_Positives <- nrow(Positives)
  N_Negatives <- nrow(Negatives)
  
  N_TruePositive <- sum(Positives$Gene.names %in% reference)
  N_TrueNegatives <- sum(!Negatives$Gene.names %in% reference)
  
  N_FalsePositive <- sum(!Positives$Gene.names %in% reference)
  
  #sensitivity (true positive rate)
  TPR <- N_TruePositive/N_Positives
  
  #specificity (true negative rate)
  TNR <- N_TrueNegatives/N_Negatives
  
  # False Positive Rate (1- specificity)
  FPR <- 1 - TNR
  
  # False Discovery Rate
  FDR <- N_FalsePositive/N_Positives
  
  df_curve_specific <- data.frame(cutoff = round(cutoff_specific,3),
                                  FPR = round(FPR,3),
                                  TPR_sensitivity = round(TPR,3),
                                  FDR = round(FDR,3),
                                  TNR_specificity = round(TNR,3),
                                  Selected = N_Positives)
  
  return(df_curve_specific)
}

MyFDRCurve <- function(df = PG_summ, X, reference = allRBPs, co_min = 0.01, co_max = 2, co_breaks = 25, cutoff_specific = (2**DeltaAnalysisThresholdEfficiency)) {
  
  co <- co_min
  df_curve <- MySpecificFDR(df = df, X = X, reference = reference, cutoff_specific = cutoff_specific)
  
  while(co < co_max){
    
    df_curve <- rbind(df_curve, MySpecificFDR(df = df, X = X, reference = reference, cutoff_specific = co))
    
    co = co + ((co_max/co_min)/co_breaks)*co_min
    
  }
  
  df_curve$cutoff <- round(df_curve$cutoff, 3)
  
  return(df_curve)
}

MyFDRplot <- function(df_curve, breaks_names = 3, Title = NULL, cutoff_specific = (2**DeltaAnalysisThresholdEfficiency), co_min = 0.01, co_max = 2, co_breaks = 25){
  
  df_curve_cutoff <- subset(df_curve, cutoff == cutoff_specific)[1,]
  
  df_plot <- MyScatterPlot(df = df_curve, X = "cutoff", Y = "FDR",
                           xmin = 0, xmax = co_max, Title = Title,
                           AxisSep.x = 0.5,
                           ymin = 0, ymax = 1,
                           AxisSep.y = 0.2,
                           InPercentage.x = F, ShowRBPs = F,
                           ShowNotReproducible = F, size_points = 2,
                           GradientColorMin = rgb(1,0,0,0), GradientColorMax = rgb(1,0,0,0)) +
    
    geom_hline(yintercept = df_curve_cutoff[, "FDR"], linetype = "dashed", color = rgb(.5,.5,.5)) +
    
    geom_vline(xintercept = df_curve_cutoff[, "cutoff"], linetype = "dashed", color = rgb(.5,.5,.5)) +
    
    geom_text_repel(data = df_curve[seq(1, nrow(df_curve), breaks_names),], aes(cutoff, FDR, label = Selected), size = 8) +
    
    geom_text_repel(data = df_curve_cutoff, aes(cutoff, FDR, label = Selected), size = 8, color = rgb(1,0,0))
  
  return(df_plot)
}




FalsePositive_Cl <- MyFDRplot(MyFDRCurve(X = "Percentage_Mean"), Title = "Cl")
FalsePositive_NoCl <- MyFDRplot(MyFDRCurve(X = "Percentage_Mean_NoCl"), Title = "NoCl")

MySpecificFDR(X = "Percentage_Mean", reference = allRBPs, cutoff_specific = (2**DeltaAnalysisThresholdEfficiency))$FDR
MySpecificFDR(X = "Percentage_Mean", reference = mRBPs, cutoff_specific = (2**DeltaAnalysisThresholdEfficiency))$FDR
MySpecificFDR(X = "Percentage_Mean", reference = mRBPs_HEK, cutoff_specific = (2**DeltaAnalysisThresholdEfficiency))$FDR
MySpecificFDR(X = "Percentage_Mean", reference = qRIC_mRBP, cutoff_specific = (2**DeltaAnalysisThresholdEfficiency))$FDR



png(filename = "FalsePositiveRate.png", height = 500, width = 1000)
grid.arrange(FalsePositive_Cl, FalsePositive_NoCl, 
             ncol = 2)
dev.off()
rm(FalsePositive_Cl, FalsePositive_NoCl)



#### ROC curves ####
library(pROC)
MyCreate_ROC_df <- function(df, X, reference = c("Known_RBP",
                                                 "Known_mRBP", 
                                                 #"qRIC_mRBP",
                                                 "Known_mRBP_HEK")) {
  
  df <- subset(df, !is.na(eval(parse(text = X))))
  
  #Main df_interactors with Nuc
  pROC_obj <- roc(df[,reference[1]],
                  df[,X], direction = "<")
  
  df_pROC <- data.frame(sensitivities = pROC_obj$sensitivities,
                        thresholds = pROC_obj$thresholds,
                        specificities = pROC_obj$specificities)
  df_pROC$FPR <- 1- df_pROC$specificities
  df_pROC$Factor <- as.factor(reference[1])
  
  df_pROC$Number_Selected <- apply(df_pROC["thresholds"], 1, function(x) sum(PG_summ[,X] >= x, na.rm = T))
  
  cat("Reference", "\t", "AUROC", "\n")
  cat(reference[1], "\t", pROC_obj$auc, "\n")
  
  if (length(reference) > 1) {
    
    df_pROC <- rbind(df_pROC, MyCreate_ROC_df(df = df, X = X, reference = reference[2:length(reference)])) 
  }
  
  return(df_pROC)
  
}


df_Cl <- MyCreate_ROC_df(df = PG_summ, X = "Percentage_Mean")
df_Cl$Group <- "Crosslinked"
df_NoCl <- MyCreate_ROC_df(df = PG_summ, X = "Percentage_Mean_NoCl")
df_NoCl$Group <- "Control"

df_all <- rbind(df_Cl, df_NoCl)
df_all$Group <- factor(df_all$Group, levels = c("Crosslinked", "Control"))


myROCplot <- function(df = df_all, showLegend = T, Title = NULL, breaks_names = NULL, my.cutoff = (2**DeltaAnalysisThresholdEfficiency)) {
  
  ROC_plot <- ggplot(data = df, aes(x = FPR, y = sensitivities, color = Factor)) +
    
    facet_grid(cols = vars(Group)) +
    
    geom_abline(linetype = 1, color = "grey") +
    
    geom_line(size = 1, show.legend = showLegend) +
    
    ggtitle(Title) +
    ylab("TPR (sensitivity)") +
    xlab("FPR (1-specificity)") +
    coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) +
    scale_y_continuous(breaks = seq(0, 1, .2)) +
    scale_x_continuous(breaks = seq(0, 1, .2)) +
    theme_bw() +
    theme(panel.grid.minor = element_blank(),
          plot.title = element_text(lineheight=.8,
                                    #face="bold", 
                                    vjust=0.5,
                                    hjust = 0.5,
                                    size=16),
          axis.text.x  = element_text(#face = "bold",
            color = "black",
            angle=0, 
            hjust=.5,
            size=12),
          axis.title.x = element_text(#face="bold",
            size=16,
            hjust = 0.5,
            vjust = 1.5),
          axis.title.y = element_text(#face="bold",
            size=16,
            hjust = .5,
            vjust = 0.5),
          axis.text.y  = element_text(#face = "bold",
            color = "black",
            angle=0, 
            vjust=.5,
            size=12),
          aspect.ratio = 1)
  
  if(!is.null(breaks_names)) {
    ROC_plot <- ROC_plot  +
      
      geom_text_repel(data = df[seq(1, nrow(df), breaks_names),],
                      aes(x = FPR, y = sensitivities, label = round(thresholds, 3)),
                      size = 8, show.legend = F, nudge_x = 0.05, nudge_y = -0.05)
  }
  
  if(!is.null(my.cutoff)) {
    
    df_crosslink <- subset(df, Group == "Crosslinked")
    df_control <- subset(df, Group == "Control")
    
    index_crosslink <- which(abs(df_crosslink$thresholds-my.cutoff)==
                               min(abs(df_crosslink$thresholds-my.cutoff)))
    index_control <- which(abs(df_control$thresholds-my.cutoff)==
                                                                                                      min(abs(df_control$thresholds-my.cutoff)))
    
    df_sub_crosslink <- df_crosslink[index_crosslink,]
    df_sub_control <- df_control[index_control,]
    
    ROC_plot <- ROC_plot  +
      
      geom_point(data = df_sub_crosslink,
                 aes(x = FPR, y = sensitivities),
                 size = 3, show.legend = F) +
      
      geom_text_repel(data = df_sub_crosslink,
                      aes(x = FPR, y = sensitivities, 
                          label = paste0("(", round(specificities,2), ",",
                                         round(sensitivities,2), ")")),
                      size = 3, show.legend = F,
                      nudge_x = 0.05, nudge_y = -0.05) +
      
      geom_point(data = df_sub_crosslink,
                 aes(x = FPR, y = sensitivities),
                 size = 3, show.legend = F) +
      
      geom_text_repel(data = df_sub_control,
                      aes(x = FPR, y = sensitivities, 
                          label = paste0("(", round(specificities,2), ",",
                                         round(sensitivities,2), ")")),
                      size = 3, show.legend = F,
                      nudge_x = 0.05, nudge_y = -0.05)
  }
  
  return(ROC_plot)
  
  
}


ROC_plot <- myROCplot(df_all)



png("ROC_curves.png")
grid.arrange(ROC_plot, 
             ncol = 1)
dev.off()

pdf("ROC_curves.pdf", useDingbats = F)
grid.arrange(ROC_plot, 
             ncol = 1)
dev.off()
rm(ROC_plot)





#### Cutoff numbers ####
PG_summ_filtered <- subset(PG_summ, !is.na(Ratio.H.L.Mean) & Gene.names %in% mRBPs)

Identified_RBPs <- unique(PG_summ_filtered$Gene.names)

# Fraction of HEK mRBPs
100*length(intersect(Identified_RBPs, mRBPs_HEK))/length(mRBPs_HEK)

# Fraction of all mRBPs
100*length(intersect(Identified_RBPs, mRBPs))/length(mRBPs)




# how many proteins were identified (have intensities values not NA) in each L-H group pair
{cat("Group", "\t", "\t", "\t", "Number", "Mean", "Requantified", "% RQ", "\n")
  for (i in 1:length(intensities)) {
    cat(intensities[i])
    cat("\t")
    cat("\t")
    cat(nrow(PG_summ_filtered[!is.na(PG_summ_filtered[intensities[i]]),]))
    cat("\t")
    cat(signif(mean(PG_summ_filtered[,intensities[i]], na.rm = T), 4))
    cat("\t")
    cat(sum(PG_summ_filtered[,requantified.intensities[i]], na.rm = T))
    cat("\t")
    cat(signif(100*(sum(PG_summ_filtered[,requantified.intensities[i]], na.rm = T)/nrow(PG_summ_filtered[!is.na(PG_summ_filtered[intensities[i]]),])),3))
    cat("\n")
    rm(i)
  }}


# how many proteins were identified (have ratios values not NA) in each L-H group pair
{cat("Group", "\t", "\t", "\t", "\t", "Number", "Mean", "Requantified", "% RQ", "\n")
  for (i in 1:length(ratios)) {
    #group
    cat(ratios[i])
    cat("\t")
    cat("\t")
    #quantified ratios
    cat(nrow(PG_summ_filtered[!is.na(PG_summ_filtered[ratios[i]]),]))
    cat("\t")
    #ratio mean
    cat(signif(mean(PG_summ_filtered[,ratios[i]], na.rm = T), 3))
    cat("\t")
    #Requantified ratios
    cat(sum(!is.na(PG_summ_filtered[ratios[i]]) & PG_summ_filtered[, requantified.ratios[i]], na.rm = T))
    cat("\t")
    #Percentage of requantified ratios
    cat(100*(sum(!is.na(PG_summ_filtered[ratios[i]]) & PG_summ_filtered[, requantified.ratios[i]], na.rm = T)) / nrow(PG_summ_filtered[!is.na(PG_summ_filtered[ratios[i]]),]))
    cat("\n")
    rm(i)
  }
}




# how many proteins were identified (have intensities values not NA) in each L-H group pair
{cat("Group", "\t", "\t", "\t", "\t", "\t", "Unscrupulous", "% Unscrupulous","\n")
  for (i in 1:length(unscrupulous.ratios)) {
    cat(unscrupulous.ratios[i])
    cat("\t")
    cat("\t")
    cat(sum(PG_summ_filtered[,unscrupulous.ratios[i]]))
    cat("\t")
    cat(100*(sum(PG_summ_filtered[,unscrupulous.ratios[i]])/sum(!PG_summ_filtered[,unscrupulous.ratios[i]])))
    cat("\n")
    rm(i)
  }
}

