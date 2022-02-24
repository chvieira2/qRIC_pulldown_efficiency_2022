setwd("K:/Results/20210903_Phospho_qRIC_revisionPaper/output")

library(ggplot2)
library(reshape2)
library(RColorBrewer)
library(gridExtra)
library(dplyr)
library(limma)
library(ggrepel)
library(VennDiagram)
library(gplots)
library(data.table)
library(GGally)



Requantify = T
ReproducibilityDifference = 3
DeltaAnalysisThresholdEfficiency <- -1






#### PTM table load and preparation ####

PTMTable <- fread("../txt_RQ/Phospho (STY)Sites.txt", sep = "\t", stringsAsFactors = FALSE, data.table = F, check.names = T)
PTMTable_noRQ <- fread("../txt_noRQ/Phospho (STY)Sites.txt", sep = "\t", stringsAsFactors = FALSE, data.table = F, check.names = T)
if (!Requantify) {
  PTMTable <- PTMTable_noRQ
  
}
rm(Requantify)

#generate a new column concatenating Gene.name with PTM site information
#This code is ignoring the fact that some peptides are modified more than once......
PTMTable$Gene.names <- sapply(strsplit(PTMTable$Gene.names, ";"), "[", 1)
PTMTable$Protein <- sapply(strsplit(PTMTable$Protein, ";"), "[", 1)


PTMTable_noRQ$Gene.names <- sapply(strsplit(PTMTable_noRQ$Gene.names, ";"), "[", 1)
PTMTable_noRQ$Protein <- sapply(strsplit(PTMTable_noRQ$Protein, ";"), "[", 1)


#Generate a column with the gene name (or UniprotID when there is no gene name) plus modified site
PTMTable$Site.ID <- ifelse(is.na(PTMTable$Gene.names),
                           paste0(PTMTable$Protein, "_", PTMTable$Amino.acid, PTMTable$Position),
                           paste0(PTMTable$Gene.names, "_", PTMTable$Amino.acid, PTMTable$Position))

PTMTable_noRQ$Site.ID <- ifelse(is.na(PTMTable_noRQ$Gene.names),
                                paste0(PTMTable_noRQ$Protein, "_", PTMTable_noRQ$Amino.acid, PTMTable_noRQ$Position),
                                paste0(PTMTable_noRQ$Gene.names, "_", PTMTable_noRQ$Amino.acid, PTMTable_noRQ$Position))


#### Make sure both tables have same length
PTMTable <- subset(PTMTable, Site.ID %in% PTMTable_noRQ$Site.ID)
PTMTable_noRQ <- subset(PTMTable_noRQ, Site.ID %in% PTMTable$Site.ID)



intensities <- c("Intensity.L.Forward_II",
                 "Intensity.M.Forward_II",
                 "Intensity.H.Forward_II",
                 "Intensity.L.Forward_III",
                 "Intensity.H.Forward_III",
                 "Intensity.L.Forward_IV",
                 "Intensity.H.Forward_IV",
                 
                 "Intensity.L.Reverse_II",
                 "Intensity.M.Reverse_II",
                 "Intensity.H.Reverse_II",
                 "Intensity.L.Reverse_III",
                 "Intensity.H.Reverse_III",
                 "Intensity.L.Reverse_IV",
                 "Intensity.H.Reverse_IV")


ratios <- c("Ratio.H.L.Forward_II",
            "Ratio.H.L.Forward_III",
            "Ratio.H.L.Forward_IV",
            "Ratio.H.L.Reverse_II",
            "Ratio.H.L.Reverse_III",
            "Ratio.H.L.Reverse_IV",
            
            "Ratio.H.M.Forward_II",
            "Ratio.H.M.Reverse_II",
            
            "Ratio.M.L.Forward_II",
            "Ratio.M.L.Reverse_II")


PTMTable[,ratios] <- PTMTable[,c("Ratio.H.L.Forward_II___1",
                                 "Ratio.H.L.Forward_III___1",
                                 "Ratio.H.L.Forward_IV___1",
                                 "Ratio.H.L.Reverse_II___1",
                                 "Ratio.H.L.Reverse_III___1",
                                 "Ratio.H.L.Reverse_IV___1",
                                 
                                 "Ratio.H.M.Forward_II___1",
                                 "Ratio.H.M.Reverse_II___1",
                                 
                                 "Ratio.M.L.Forward_II___1",
                                 "Ratio.M.L.Reverse_II___1")]

PTMTable_noRQ[,ratios] <- PTMTable_noRQ[,c("Ratio.H.L.Forward_II___1",
                                           "Ratio.H.L.Forward_III___1",
                                           "Ratio.H.L.Forward_IV___1",
                                           "Ratio.H.L.Reverse_II___1",
                                           "Ratio.H.L.Reverse_III___1",
                                           "Ratio.H.L.Reverse_IV___1",
                                           
                                           "Ratio.H.M.Forward_II___1",
                                           "Ratio.H.M.Reverse_II___1",
                                           
                                           "Ratio.M.L.Forward_II___1",
                                           "Ratio.M.L.Reverse_II___1")]







ratios_unmod_pep <- c("Ratio.H.L.unmod..pep..Forward_II",
                      "Ratio.H.L.unmod..pep..Forward_III",
                      "Ratio.H.L.unmod..pep..Forward_IV",
                      "Ratio.H.L.unmod..pep..Reverse_II",
                      "Ratio.H.L.unmod..pep..Reverse_III",
                      "Ratio.H.L.unmod..pep..Reverse_IV",
                      
                      "Ratio.H.M.unmod..pep..Forward_II",
                      "Ratio.H.M.unmod..pep..Reverse_II",
                      
                      "Ratio.M.L.unmod..pep..Forward_II",
                      "Ratio.M.L.unmod..pep..Reverse_II")




#### Define Requantified ratios and intensities ####

### Ratios
requantified.ratios <- c()
unscrupulous.ratios <- c()
for (i in ratios) {
  PTMTable[paste0("Requantified.", i)] <- ifelse(is.na(PTMTable_noRQ[i]) &
                                                   !is.na(PTMTable[i]), T, F)
  
  requantified.ratios <- c(requantified.ratios, paste0("Requantified.", i))
  unscrupulous.ratios <- c(unscrupulous.ratios, paste0("Unscrupulous.", i))
  rm(i)
}



### intensities
requantified.intensities <- c()
for (i in intensities) {
  PTMTable[paste0("Requantified.", i)] <- ifelse(PTMTable_noRQ[i] == 0 &
                                                   PTMTable[i] != 0, T, F)
  
  requantified.intensities <- c(requantified.intensities, paste0("Requantified.", i))
  rm(i)
}









#### Unscrupulous Requantification ####
# Select Unscrupulous Requantified ratios (ratios between two channels when both are requantified)
#Forward
PTMTable$Unscrupulous.Ratio.H.L.Forward_II <- ifelse((PTMTable$Requantified.Intensity.H.Forward_II & PTMTable$Requantified.Intensity.L.Forward_II) & !is.na(PTMTable$Ratio.H.L.Forward_II), T, F)

PTMTable$Unscrupulous.Ratio.H.L.Forward_III <- ifelse((PTMTable$Requantified.Intensity.H.Forward_III & PTMTable$Requantified.Intensity.L.Forward_III) & !is.na(PTMTable$Ratio.H.L.Forward_III), T, F)

PTMTable$Unscrupulous.Ratio.H.L.Forward_IV <- ifelse((PTMTable$Requantified.Intensity.H.Forward_IV & PTMTable$Requantified.Intensity.L.Forward_IV) & !is.na(PTMTable$Ratio.H.L.Forward_IV), T, F)




PTMTable$Unscrupulous.Ratio.H.M.Forward_II <- ifelse((PTMTable$Requantified.Intensity.H.Forward_II & PTMTable$Requantified.Intensity.M.Forward_II) & !is.na(PTMTable$Ratio.H.M.Forward_II), T, F)




PTMTable$Unscrupulous.Ratio.M.L.Forward_II <- ifelse((PTMTable$Requantified.Intensity.M.Forward_II & PTMTable$Requantified.Intensity.L.Forward_II) & !is.na(PTMTable$Ratio.M.L.Forward_II), T, F)







#Reverse
PTMTable$Unscrupulous.Ratio.H.L.Reverse_II <- ifelse((PTMTable$Requantified.Intensity.H.Reverse_II & PTMTable$Requantified.Intensity.L.Reverse_II) & !is.na(PTMTable$Ratio.H.L.Reverse_II), T, F)

PTMTable$Unscrupulous.Ratio.H.L.Reverse_III <- ifelse((PTMTable$Requantified.Intensity.H.Reverse_III & PTMTable$Requantified.Intensity.L.Reverse_III) & !is.na(PTMTable$Ratio.H.L.Reverse_III), T, F)

PTMTable$Unscrupulous.Ratio.H.L.Reverse_IV <- ifelse((PTMTable$Requantified.Intensity.H.Reverse_IV & PTMTable$Requantified.Intensity.L.Reverse_IV) & !is.na(PTMTable$Ratio.H.L.Reverse_IV), T, F)




PTMTable$Unscrupulous.Ratio.H.M.Reverse_II <- ifelse((PTMTable$Requantified.Intensity.H.Reverse_II & PTMTable$Requantified.Intensity.M.Reverse_II) & !is.na(PTMTable$Ratio.H.M.Reverse_II), T, F)



PTMTable$Unscrupulous.Ratio.M.L.Reverse_II <- ifelse((PTMTable$Requantified.Intensity.M.Reverse_II & PTMTable$Requantified.Intensity.L.Reverse_II) & !is.na(PTMTable$Ratio.M.L.Reverse_II), T, F)








# Delet Unscrupulous ratios
for (i in 1:length(ratios)) {
  PTMTable[,ratios[i]] <- as.numeric(ifelse(PTMTable[,unscrupulous.ratios[i]], NA, PTMTable[,ratios[i]]))
  
  rm(i)
}







rm(PTMTable_noRQ)


#### Filters ####
PTMTable <- subset(PTMTable, PTMTable$Reverse != "+")
PTMTable <- subset(PTMTable, PTMTable$Potential.contaminant != "+")


#Localization probability threshold used in the field is usually 0.75 (class 1 peptides)
PTMTable <- subset(PTMTable, PTMTable$Localization.prob >= 0.75)

#Generate a column with the peptide sequence Site.ID of the given modification
PTMTable$Peptide.Site.ID <- gsub( " *\\(.*?\\) *", "", PTMTable$Phospho..STY..Probabilities)




# Line needed only if any column has all values = NA
#PTMTable[c(intensities, ratios, ratios_unmod_pep)] = apply(PTMTable[c(intensities, ratios, ratios_unmod_pep)], 1, function(x) as.numeric(x))

#transform intensities and ratios to Log2 (or 10 if you prefer)
PTMTable[c(intensities, ratios, ratios_unmod_pep)] = log2(PTMTable[c(intensities, ratios, ratios_unmod_pep)])
# change Inf values for na
is.na(PTMTable[c(intensities, ratios, ratios_unmod_pep)]) <- sapply(PTMTable[c(intensities, ratios, ratios_unmod_pep)], is.infinite)

is.na(PTMTable[c(intensities, ratios, ratios_unmod_pep)]) <- sapply(PTMTable[c(intensities, ratios, ratios_unmod_pep)], is.nan)






# how many proteins were identified (have intensities values not NA) in each L-H group pair
for (i in 1:length(intensities)) {
  cat(intensities[i])
  cat("\t")
  cat("\t")
  cat(nrow(PTMTable[!is.na(PTMTable[intensities[i]]),]))
  cat("\t")
  cat(signif(mean(PTMTable[,intensities[i]], na.rm = T), 4))
  cat("\t")
  cat(sum(PTMTable[,requantified.intensities[i]], na.rm = T))
  cat("\t")
  cat(signif(100*(sum(PTMTable[,requantified.intensities[i]], na.rm = T)/nrow(PTMTable[!is.na(PTMTable[intensities[i]]),])),3))
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
  cat(nrow(PTMTable[!is.na(PTMTable[ratios[i]]),]))
  cat("\t")
  #ratio mean
  cat(signif(mean(PTMTable[,ratios[i]], na.rm = T), 3))
  cat("\t")
  #Requantified ratios
  cat(sum(!is.na(PTMTable[ratios[i]]) & PTMTable[, requantified.ratios[i]], na.rm = T))
  cat("\t")
  #Percentage of requantified ratios
  cat(100*(sum(!is.na(PTMTable[ratios[i]]) & PTMTable[, requantified.ratios[i]], na.rm = T)) / nrow(PTMTable[!is.na(PTMTable[ratios[i]]),]))
  cat("\n")
  rm(i)
}


### Unmod peptides
for (i in 1:length(ratios_unmod_pep)) {
  #group
  cat(ratios_unmod_pep[i])
  cat("\t")
  cat("\t")
  #quantified ratios_unmod_pep
  cat(nrow(PTMTable[!is.na(PTMTable[ratios_unmod_pep[i]]),]))
  cat("\t")
  #ratio mean
  cat(signif(mean(PTMTable[,ratios_unmod_pep[i]], na.rm = T), 3))
  cat("\n")
  rm(i)
}


# how many proteins were identified (have intensities values not NA) in each L-H group pair
for (i in 1:length(unscrupulous.ratios)) {
  cat(unscrupulous.ratios[i])
  cat("\t")
  cat("\t")
  cat(sum(PTMTable[,unscrupulous.ratios[i]]))
  cat("\t")
  cat(100*(sum(PTMTable[,unscrupulous.ratios[i]])/sum(!PTMTable[,unscrupulous.ratios[i]])))
  cat("\n")
  rm(i)
}



#rm(requantified.intensities, requantified.ratios, unscrupulous.ratios)

#nrow(intersect(PTMTable[!is.na(PTMTable[intensities[4]]),], PTMTable[!is.na(PTMTable[intensities[3]]),]))







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
mRBPs_HEK <- subset(COMPILED_TABLE_Hs, RBPBASE000000007.1)$UNIQUE # Baltz et al

mRBPs <- subset(COMPILED_TABLE_Hs, RBPBASE000000007.1 | # Baltz et al
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
PTMTable$Known_RBP <- ifelse(PTMTable$Gene.names %in% allRBPs, T, F)
PTMTable$Known_mRBP <- ifelse(PTMTable$Gene.names %in% mRBPs, T, F)
PTMTable$Known_mRBP_HEK <- ifelse(PTMTable$Gene.names %in% mRBPs_HEK, T, F)








#### Defining PTM working table table ####

PTM_summ <- subset(PTMTable,
                   select = c("id","Protein.group.IDs", "Protein", "Gene.names", "Site.ID",
                              "Peptide.Site.ID", "Known_RBP", "Known_mRBP", "Known_mRBP_HEK",
                              "Amino.acid", "Position",
                              "Localization.prob",
                              "Score.diff", "PEP", "Score", "Delta.score",
                              "Score.for.localization",
                              intensities, ratios, ratios_unmod_pep,
                              requantified.intensities, requantified.ratios,
                              unscrupulous.ratios))








#### Invert label swap ratios ####
PTM_summ$Ratio.H.L.Reverse_II <- -PTM_summ$Ratio.H.L.Reverse_II
PTM_summ$Ratio.H.L.Reverse_III <- -PTM_summ$Ratio.H.L.Reverse_III
PTM_summ$Ratio.H.L.Reverse_IV <- -PTM_summ$Ratio.H.L.Reverse_IV

PTM_summ$Ratio.M.L.Reverse_II <- -PTM_summ$Ratio.M.L.Reverse_II

PTM_summ$Ratio.H.M.Reverse_II <- -PTM_summ$Ratio.H.M.Reverse_II






PTM_summ$Ratio.H.L.unmod..pep..Reverse_II <- -PTM_summ$Ratio.H.L.unmod..pep..Reverse_II
PTM_summ$Ratio.H.L.unmod..pep..Reverse_III <- -PTM_summ$Ratio.H.L.unmod..pep..Reverse_III
PTM_summ$Ratio.H.L.unmod..pep..Reverse_IV <- -PTM_summ$Ratio.H.L.unmod..pep..Reverse_IV

PTM_summ$Ratio.M.L.unmod..pep..Reverse_II <- -PTM_summ$Ratio.M.L.unmod..pep..Reverse_II

PTM_summ$Ratio.H.M.unmod..pep..Reverse_II <- -PTM_summ$Ratio.H.M.unmod..pep..Reverse_II











#### Correct ratios to input fractions ####

# I originally tried to correct for losses in the Ethanol precipitation of the input samples (50% loss) by multiplying all ratios by 1/2 (sum log2(0.5)). I decided to not do it to reduce data manipulation.

# Replicate II
PTM_summ$Ratio.H.L.Forward_II <- (PTM_summ$Ratio.H.L.Forward_II)
PTM_summ$Ratio.M.L.Forward_II <- (PTM_summ$Ratio.M.L.Forward_II) + log2(2) #I used half of Pulldown for no-crosslink samples

PTM_summ$Ratio.H.L.Reverse_II <- (PTM_summ$Ratio.H.L.Reverse_II)
PTM_summ$Ratio.H.M.Reverse_II <- (PTM_summ$Ratio.H.M.Reverse_II) + log2(2) #I used half of Pulldown for no-crosslink samples


# Replicate III
PTM_summ$Ratio.H.L.Forward_III <- (PTM_summ$Ratio.H.L.Forward_III)

PTM_summ$Ratio.H.L.Reverse_III <- (PTM_summ$Ratio.H.L.Reverse_III)



# Replicate IV
PTM_summ$Ratio.H.L.Forward_IV <- (PTM_summ$Ratio.H.L.Forward_IV)

PTM_summ$Ratio.H.L.Reverse_IV <- (PTM_summ$Ratio.H.L.Reverse_IV)




### Unmodified peptides
# Replicate II
PTM_summ$Ratio.H.L.unmod..pep..Forward_II <- (PTM_summ$Ratio.H.L.unmod..pep..Forward_II)
PTM_summ$Ratio.M.L.unmod..pep..Forward_II <- (PTM_summ$Ratio.M.L.unmod..pep..Forward_II) + log2(2) #I used half of Pulldown for no-crosslink samples

PTM_summ$Ratio.H.L.unmod..pep..Reverse_II <- (PTM_summ$Ratio.H.L.unmod..pep..Reverse_II)
PTM_summ$Ratio.H.M.unmod..pep..Reverse_II <- (PTM_summ$Ratio.H.M.unmod..pep..Reverse_II) + log2(2) #I used half of Pulldown for no-crosslink samples


# Replicate III
PTM_summ$Ratio.H.L.unmod..pep..Forward_III <- (PTM_summ$Ratio.H.L.unmod..pep..Forward_III)

PTM_summ$Ratio.H.L.unmod..pep..Reverse_III <- (PTM_summ$Ratio.H.L.unmod..pep..Reverse_III)



# Replicate IV
PTM_summ$Ratio.H.L.unmod..pep..Forward_IV <- (PTM_summ$Ratio.H.L.unmod..pep..Forward_IV)

PTM_summ$Ratio.H.L.unmod..pep..Reverse_IV <- (PTM_summ$Ratio.H.L.unmod..pep..Reverse_IV)


#### Define RBPs from NoCL ####
PTM_summ$qRIC_mRBP <- ifelse((PTM_summ$Ratio.H.M.Forward_II >= 1 &
                                PTM_summ$Ratio.M.L.Reverse_II >= 1) &
                               (!is.na(PTM_summ$Ratio.H.M.Forward_II) |
                                  !is.na(PTM_summ$Ratio.M.L.Reverse_II)), T, F)




#### Ratios into Percentages ####
# Replicate II
PTM_summ$Percentage.H.L.Forward_II <- (2^PTM_summ$Ratio.H.L.Forward_II)
PTM_summ$Percentage.M.L.Forward_II <- (2^PTM_summ$Ratio.M.L.Forward_II)

PTM_summ$Percentage.H.L.Reverse_II <- (2^PTM_summ$Ratio.H.L.Reverse_II)
PTM_summ$Percentage.H.M.Reverse_II <- (2^PTM_summ$Ratio.H.M.Reverse_II)


# Replicate III
PTM_summ$Percentage.H.L.Forward_III <- (2^PTM_summ$Ratio.H.L.Forward_III)

PTM_summ$Percentage.H.L.Reverse_III <- (2^PTM_summ$Ratio.H.L.Reverse_III)



# Replicate IV
PTM_summ$Percentage.H.L.Forward_IV <- (2^PTM_summ$Ratio.H.L.Forward_IV)

PTM_summ$Percentage.H.L.Reverse_IV <- (2^PTM_summ$Ratio.H.L.Reverse_IV)


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
PTM_summ$Ratio.H.L.Forward <- MyMeanWithMinValues(PTM_summ, c("Ratio.H.L.Forward_II",
                                                              "Ratio.H.L.Forward_III",
                                                              "Ratio.H.L.Forward_IV"),
                                                  MinNumberReplicates = 1)

PTM_summ$Ratio.H.L.Reverse <- MyMeanWithMinValues(PTM_summ, c("Ratio.H.L.Reverse_II",
                                                              "Ratio.H.L.Reverse_III",
                                                              "Ratio.H.L.Reverse_IV"),
                                                  MinNumberReplicates = 1)

PTM_summ$Ratio.H.L.Mean <- MyMeanWithMinValues(PTM_summ, c("Ratio.H.L.Forward",
                                                           "Ratio.H.L.Reverse"),
                                               MinNumberReplicates = 2)



PTM_summ$Ratio.NoCl.Mean <- MyMeanWithMinValues(PTM_summ, c("Ratio.M.L.Forward_II",
                                                          "Ratio.H.M.Reverse_II"),
                                               MinNumberReplicates = 2)


PTM_summ$Ratio.ClXNoCl.Mean <- MyMeanWithMinValues(PTM_summ, c("Ratio.H.M.Forward_II",
                                                             "Ratio.M.L.Reverse_II"),
                                                  MinNumberReplicates = 2)






## Unmod peptides
PTM_summ$Ratio.H.L.unmod..pep..Forward <- MyMeanWithMinValues(PTM_summ, c("Ratio.H.L.unmod..pep..Forward_II",
                                                              "Ratio.H.L.unmod..pep..Forward_III",
                                                              "Ratio.H.L.unmod..pep..Forward_IV"),
                                                  MinNumberReplicates = 1)

PTM_summ$Ratio.H.L.unmod..pep..Reverse <- MyMeanWithMinValues(PTM_summ, c("Ratio.H.L.unmod..pep..Reverse_II",
                                                              "Ratio.H.L.unmod..pep..Reverse_III",
                                                              "Ratio.H.L.unmod..pep..Reverse_IV"),
                                                  MinNumberReplicates = 1)

PTM_summ$Ratio.H.L.unmod..pep..Mean <- MyMeanWithMinValues(PTM_summ, c("Ratio.H.L.unmod..pep..Forward",
                                                           "Ratio.H.L.unmod..pep..Reverse"),
                                               MinNumberReplicates = 2)










#Intensities
PTM_summ$Intensity.Input.Mean <- MyMeanWithMinValues(PTM_summ, c("Intensity.L.Forward_II",
                                                                 "Intensity.L.Forward_III",
                                                                 "Intensity.L.Forward_IV",
                                                                 "Intensity.H.Reverse_II",
                                                                 "Intensity.H.Reverse_III",
                                                                 "Intensity.H.Reverse_IV"),
                                                     MinNumberReplicates = 2)









# Percentages
PTM_summ$Percentage.H.L.Forward <- 2**PTM_summ$Ratio.H.L.Forward

PTM_summ$Percentage.H.L.Reverse <- 2**PTM_summ$Ratio.H.L.Reverse

PTM_summ$Percentage.M.L.Forward <- 2**PTM_summ$Ratio.M.L.Forward
PTM_summ$Percentage.H.M.Reverse <- 2**PTM_summ$Ratio.H.M.Reverse



PTM_summ$Percentage_Mean <- 2**PTM_summ$Ratio.H.L.Mean

PTM_summ$Percentage_Mean_NoCl <- 2**PTM_summ$Ratio.NoCl.Mean




#### Reproducibility ####
PTM_summ$Reproducible.logical_II <- ifelse((PTM_summ$Ratio.H.L.Forward_II - PTM_summ$Ratio.H.L.Reverse_II < ReproducibilityDifference) & (PTM_summ$Ratio.H.L.Reverse_II - PTM_summ$Ratio.H.L.Forward_II < ReproducibilityDifference), T, F)

PTM_summ$Reproducible.logical_III <- ifelse((PTM_summ$Ratio.H.L.Forward_III - PTM_summ$Ratio.H.L.Reverse_III < ReproducibilityDifference) & (PTM_summ$Ratio.H.L.Reverse_III - PTM_summ$Ratio.H.L.Forward_III < ReproducibilityDifference), T, F)

PTM_summ$Reproducible.logical_IV <- ifelse((PTM_summ$Ratio.H.L.Forward_IV - PTM_summ$Ratio.H.L.Reverse_IV < ReproducibilityDifference) & (PTM_summ$Ratio.H.L.Reverse_IV - PTM_summ$Ratio.H.L.Forward_IV < ReproducibilityDifference), T, F)




#### Write table ####
write.table(PTM_summ, file = "PTM_summary.txt", sep = "\t", na = "", quote = F, row.names = F)

#### Writing table ratios paper ####
## Do not further filter this table. It must be exactly the same as PTM_summary.txt but with selected columns only, except for inverted reverse ratios
# Table 2 paper
foo_df <- subset(PTM_summ,
                 select = c("Protein",
                            "Gene.names",
                            "Site.ID",
                            
                            "Amino.acid",
                            "Position",
                            "Localization.prob",
                            "PEP", 
                            "Score", 
                            "Delta.score",
                            
                            "Known_mRBP",
                            "Known_mRBP_HEK",
                            
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
names(foo_df) <- c("Uniprot.IDs",
                   "Gene.names",
                   "Site.ID",
                   "Amino.acid",
                   "Position",
                   "Localization.prob",
                   "PEP", 
                   "Score", 
                   "Delta.score",
                   
                   "mRBP.(Gebauer)",
                   "HEK.mRBP.(Baltz)",
                   
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


write.table(foo_df, file = "Table2_PhosphoEfficiency.txt", sep = "\t", na = "", quote = F, row.names = F)
rm(foo_df)






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

selected_protein <- c(#"LARP1", 
  "RBM20", "SERBP1", 
  "SF3B1", "UPF1", "ELAVL1")
#selected_protein <- unique(subset(PTM_summ, grepl("EIF", PTM_summ$Gene.names))$Protein)
#selected_protein <- unique(subset(PTM_summ, Intensity.Input.Mean > 31 & Ratio.H.L.Mean > 1)$Protein)
#selected_protein <- c("blank")

foo <- subset(PTM_summ, Gene.names %in% selected_protein)


ratios_CL <- MyScatterPlot(PTM_summ, xmin = -7,
                           ymin = 20, ymax = 32,
                           "Ratio.H.L.Mean", "Intensity.Input.Mean",
                           ShowRBPs = F, InPercentage.x = T, InPercentage.y = F) +
  
  geom_point(data = foo,
             na.rm = T,
             size = 3,
             colour = rgb(0,.4,.8)) +
  
  
  geom_text_repel(data = foo,
                  aes(label = Site.ID),
                  size = 7,
                  colour = rgb(0,0,1),
                  na.rm = T) +
  
  geom_point(data = foo, aes(color = Gene.names),
             na.rm = T,
             size = 5)

#png(filename = paste0("PTM_ratios_", paste(selected_protein, collapse = "_"), ".png"), height = 500, width = 500)
png(filename = paste0("PTM_ratios_", paste(selected_protein, collapse = "_"), ".png"), height = 500, width = 500)
grid.arrange(ratios_CL, 
             ncol = 1)
dev.off()
rm(foo, ratios_CL)






#### boxplot intensities ####
MyBoxplot <- function(df, variables, id.vars,
                      AxisName_y = ("Log2"), LabelNames = variables,
                      limits_y_min = round(min(df[variables], na.rm = T))-1,
                      limits_y_max = round(max(df[variables], na.rm = T))+1,
                      limits_breaks = round((limits_y_max+limits_y_min)/20)) {
  
  all <- melt(df,
              id.vars = id.vars, measure.vars=variables)
  
  
  
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
                                      vjust=.4,
                                      size=20))
  
  return(bplot)
  
}




boxplot <- MyBoxplot(PTM_summ, intensities, id.vars = "Site.ID",
                     AxisName_y = "Intensities (Log2)",
                     LabelNames = intensities)






png(paste("PTM_Boxplot.png", sep = ""), width = 1000, height = 1000, pointsize = 25)
grid.arrange(boxplot, nrow = 1)
dev.off()
rm(boxplot)









#### Forward X Reverse ratios ####
#CL
ratios_CL_II <- MyScatterPlot(PTM_summ,
                              "Ratio.H.L.Forward_II", "Ratio.H.L.Reverse_II",
                              ShowRBPs = T, InPercentage.x = T, InPercentage.y = T) + geom_abline(slope = 1)

ratios_CL_III <- MyScatterPlot(PTM_summ,
                               "Ratio.H.L.Forward_III", "Ratio.H.L.Reverse_III",
                               ShowRBPs = T, InPercentage.x = T, InPercentage.y = T) + geom_abline(slope = 1)

ratios_CL_IV <- MyScatterPlot(PTM_summ,
                              "Ratio.H.L.Forward_IV", "Ratio.H.L.Reverse_IV",
                              ShowRBPs = T, InPercentage.x = T, InPercentage.y = T) + geom_abline(slope = 1)



ratios_CL <- MyScatterPlot(PTM_summ,
                           "Ratio.H.L.Forward", "Ratio.H.L.Reverse",
                           ShowRBPs = T, InPercentage.x = T, InPercentage.y = T)




#NoCL
ratios_NoCL_II <- MyScatterPlot(PTM_summ,
                                "Ratio.M.L.Forward_II", "Ratio.H.M.Reverse_II",
                                ShowRBPs = T, InPercentage.x = T, InPercentage.y = T) + geom_abline(slope = 1)





blank <- ggplot() + theme(panel.background = element_blank(),
                          panel.border = element_blank())



png(filename = "PTM_ratios.png", height = 1000, width = 2000)
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



selected_pSite <- c("SF3B1_T313", "UPF1_S1107")
foo <- subset(PTM_summ, Site.ID %in% selected_pSite)

ratios_CL <- MyScatterPlotPaper(PTM_summ,
                           "Ratio.H.L.Forward", "Ratio.H.L.Reverse",
                           InPercentage.x = T, InPercentage.y = T,
                           GradientColorMin = rgb(.9,.45,0), 
                           GradientColorMax = rgb(.45,.225,0)) +
  
  geom_point(data = foo,
             na.rm = T,
             size = 5,
             colour = rgb(0,.2,.8)) +
  
  
  geom_text_repel(data = foo,
                  aes(label = Site.ID),
                  size = 7,
                  colour = rgb(0,.2,.8),
                  na.rm = T)


ratios_NoCL <- MyScatterPlotPaper(PTM_summ,
                             "Ratio.M.L.Forward_II", "Ratio.H.M.Reverse_II",
                             InPercentage.x = T, InPercentage.y = T,
                             GradientColorMin = rgb(.9,.45,0), 
                             GradientColorMax = rgb(.45,.225,0)) +
  
  geom_point(data = foo,
             na.rm = T,
             size = 5,
             colour = rgb(0,.2,.8)) +
  
  
  geom_text_repel(data = foo,
                  aes(label = Site.ID),
                  size = 7,
                  colour = rgb(0,.2,.8),
                  na.rm = T)

pdf(file = "PTM_ratios.pdf", height = 5, width = 10, useDingbats = F)
grid.arrange(ratios_CL, ratios_NoCL,
             ncol = 2)
dev.off()
rm(foo, ratios_CL, ratios_NoCL)






#### Figure paper Ratios RBPs only ####
ratios_CL_RBPs <- MyScatterPlotPaper(subset(PTM_summ, Known_mRBP),
                                "Ratio.H.L.Forward", "Ratio.H.L.Reverse",
                                ShowRBPs = F, InPercentage.x = T, InPercentage.y = T,
                                GradientColorMin = rgb(.9,.45,0), 
                                GradientColorMax = rgb(.45,.225,0))

ratios_CL_notRBPs <- MyScatterPlotPaper(subset(PTM_summ, !Known_mRBP),
                                   "Ratio.H.L.Forward", "Ratio.H.L.Reverse",
                                   ShowRBPs = F, InPercentage.x = T, InPercentage.y = T,
                                   GradientColorMin = rgb(.9,.45,0), 
                                   GradientColorMax = rgb(.45,.225,0))

pdf(file = "PTM_ratios_mRBPs.pdf", height = 5, width = 10, useDingbats = F)
grid.arrange(ratios_CL_RBPs, ratios_CL_notRBPs,
             nrow = 1)
dev.off()
rm(ratios_CL_RBPs, ratios_CL_notRBPs)




#### Figure paper reproducibility ####
## All proteins
pdf(file = "PTM_Pairs_ratios.pdf", height = 5, width = 5, useDingbats=FALSE)
ggpairs(PTM_summ[c("Ratio.H.L.Forward_II",
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
pdf(file = "PTM_Pairs_ratios_mRBPs.pdf", height = 5, width = 5, useDingbats=FALSE)
ggpairs(subset(PTM_summ, Known_mRBP)[c("Ratio.H.L.Forward_II",
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

ratios_All_proteins <- MyHistPlotPaper(PTM_summ, Title = "All proteins")
ratios_mRBPs <- MyHistPlotPaper(subset(PTM_summ, Known_mRBP), Title = "mRBPs only")
ratios_nonmRBPs <- MyHistPlotPaper(subset(PTM_summ, !Known_mRBP), Title = "Non mRBPs")

pdf(file = "PTM_Hist_mRBPs.pdf", height = 5, width = 15, useDingbats = F)
grid.arrange(ratios_All_proteins, ratios_mRBPs, ratios_nonmRBPs,
             nrow = 1)
dev.off()
rm(ratios_All_proteins, ratios_mRBPs, ratios_nonmRBPs)



#### Cutoff numbers ####
PTM_summ_filtered <- subset(PTM_summ, !is.na(Ratio.H.L.Mean) & Gene.names %in% mRBPs)

Identified_RBPs <- unique(PTM_summ_filtered$Gene.names)

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
    cat(nrow(PTM_summ_filtered[!is.na(PTM_summ_filtered[intensities[i]]),]))
    cat("\t")
    cat(signif(mean(PTM_summ_filtered[,intensities[i]], na.rm = T), 4))
    cat("\t")
    cat(sum(PTM_summ_filtered[,requantified.intensities[i]], na.rm = T))
    cat("\t")
    cat(signif(100*(sum(PTM_summ_filtered[,requantified.intensities[i]], na.rm = T)/nrow(PTM_summ_filtered[!is.na(PTM_summ_filtered[intensities[i]]),])),3))
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
    cat(nrow(PTM_summ_filtered[!is.na(PTM_summ_filtered[ratios[i]]),]))
    cat("\t")
    #ratio mean
    cat(signif(mean(PTM_summ_filtered[,ratios[i]], na.rm = T), 3))
    cat("\t")
    #Requantified ratios
    cat(sum(!is.na(PTM_summ_filtered[ratios[i]]) & PTM_summ_filtered[, requantified.ratios[i]], na.rm = T))
    cat("\t")
    #Percentage of requantified ratios
    cat(100*(sum(!is.na(PTM_summ_filtered[ratios[i]]) & PTM_summ_filtered[, requantified.ratios[i]], na.rm = T)) / nrow(PTM_summ_filtered[!is.na(PTM_summ_filtered[ratios[i]]),]))
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
    cat(sum(PTM_summ_filtered[,unscrupulous.ratios[i]]))
    cat("\t")
    cat(100*(sum(PTM_summ_filtered[,unscrupulous.ratios[i]])/sum(!PTM_summ_filtered[,unscrupulous.ratios[i]])))
    cat("\n")
    rm(i)
  }
}



