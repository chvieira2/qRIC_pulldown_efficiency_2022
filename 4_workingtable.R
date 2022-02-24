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
library(stringr)
library(GGally)

#foo <- subset(working_table_PTM_filtered, Gene.names_PTM == "SF3B1")
#2**subset(working_table_PTM, Gene.names_PTM == "SF3B1")$Ratio.H.L.Mean_PTM[c(2,4)]

delete.irreproducible <- T
DeltaAnalysisThresholdEfficiency <- -1
p.val.method = "one_sample_T_test" # Unpaired_T_test, one_sample_T_test, eBayes

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






#### Load PG and PTM tables ####

PTM_summ <- read.csv("PTM_summary.txt", sep = "\t", stringsAsFactors = FALSE)

Peptides_summ <- read.csv("Peptides_summary.txt", sep = "\t", stringsAsFactors = FALSE)

PG_summ <- read.csv("PG_summary.txt", sep = "\t", stringsAsFactors = FALSE)







#### Reproducibility ####
PG_summ_raw <- PG_summ
PTM_summ_raw <- PTM_summ
if (delete.irreproducible) {
  PG_summ[,c("Ratio.H.L.Forward_II")] <- ifelse(!PG_summ$Reproducible.logical_II, NA, PG_summ[,c("Ratio.H.L.Forward_II")])
  PG_summ[,c("Ratio.H.L.Reverse_II")] <- ifelse(!PG_summ$Reproducible.logical_II, NA, PG_summ[,c("Ratio.H.L.Reverse_II")])
  
  PG_summ[,c("Ratio.H.L.Forward_III")] <- ifelse(!PG_summ$Reproducible.logical_III, NA, PG_summ[,c("Ratio.H.L.Forward_III")])
  PG_summ[,c("Ratio.H.L.Reverse_III")] <- ifelse(!PG_summ$Reproducible.logical_III, NA, PG_summ[,c("Ratio.H.L.Reverse_III")])
  
  PG_summ[,c("Ratio.H.L.Forward_IV")] <- ifelse(!PG_summ$Reproducible.logical_IV, NA, PG_summ[,c("Ratio.H.L.Forward_IV")])
  PG_summ[,c("Ratio.H.L.Reverse_IV")] <- ifelse(!PG_summ$Reproducible.logical_IV, NA, PG_summ[,c("Ratio.H.L.Reverse_IV")])
}


if (delete.irreproducible) {
  Peptides_summ[,c("Ratio.H.L.Forward_II")] <- ifelse(!Peptides_summ$Reproducible.logical_II, NA, Peptides_summ[,c("Ratio.H.L.Forward_II")])
  Peptides_summ[,c("Ratio.H.L.Reverse_II")] <- ifelse(!Peptides_summ$Reproducible.logical_II, NA, Peptides_summ[,c("Ratio.H.L.Reverse_II")])
  
  Peptides_summ[,c("Ratio.H.L.Forward_III")] <- ifelse(!Peptides_summ$Reproducible.logical_III, NA, Peptides_summ[,c("Ratio.H.L.Forward_III")])
  Peptides_summ[,c("Ratio.H.L.Reverse_III")] <- ifelse(!Peptides_summ$Reproducible.logical_III, NA, Peptides_summ[,c("Ratio.H.L.Reverse_III")])
  
  Peptides_summ[,c("Ratio.H.L.Forward_IV")] <- ifelse(!Peptides_summ$Reproducible.logical_IV, NA, Peptides_summ[,c("Ratio.H.L.Forward_IV")])
  Peptides_summ[,c("Ratio.H.L.Reverse_IV")] <- ifelse(!Peptides_summ$Reproducible.logical_IV, NA, Peptides_summ[,c("Ratio.H.L.Reverse_IV")])
}


if (delete.irreproducible) {
  PTM_summ[,c("Ratio.H.L.Forward_II")] <- ifelse(!PTM_summ$Reproducible.logical_II, NA, PTM_summ[,c("Ratio.H.L.Forward_II")])
  PTM_summ[,c("Ratio.H.L.Reverse_II")] <- ifelse(!PTM_summ$Reproducible.logical_II, NA, PTM_summ[,c("Ratio.H.L.Reverse_II")])
  
  PTM_summ[,c("Ratio.H.L.Forward_III")] <- ifelse(!PTM_summ$Reproducible.logical_III, NA, PTM_summ[,c("Ratio.H.L.Forward_III")])
  PTM_summ[,c("Ratio.H.L.Reverse_III")] <- ifelse(!PTM_summ$Reproducible.logical_III, NA, PTM_summ[,c("Ratio.H.L.Reverse_III")])
  
  PTM_summ[,c("Ratio.H.L.Forward_IV")] <- ifelse(!PTM_summ$Reproducible.logical_IV, NA, PTM_summ[,c("Ratio.H.L.Forward_IV")])
  PTM_summ[,c("Ratio.H.L.Reverse_IV")] <- ifelse(!PTM_summ$Reproducible.logical_IV, NA, PTM_summ[,c("Ratio.H.L.Reverse_IV")])
}




#### Numbers removed by irreproducibility threshold ####
#Calculate mean making sure protein is identified in at least >1 replicates
MyMeanWithMinValues <- function(df, columns, MinNumberReplicates = 2) {
  foo <- apply(df[columns],
               1, 
               function(x) ifelse(sum(!is.na(x)) > MinNumberReplicates - 1,
                                  mean(x, na.rm = T),
                                  NA))
  
  return(foo)
}

# Proteins
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


# PTM
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



## Proteins
# Reproducible
PG_summ_quant <- subset(PG_summ, !is.na(Ratio.H.L.Mean))
PG_summ_HighEff <- subset(PG_summ_quant, (Ratio.H.L.Forward >= DeltaAnalysisThresholdEfficiency &
                                                Ratio.H.L.Reverse >= DeltaAnalysisThresholdEfficiency))

# All
PG_summ_raw_quant <- subset(PG_summ_raw, !is.na(Ratio.H.L.Mean))
PG_summ_raw_HighEff <- subset(PG_summ_raw_quant, (Ratio.H.L.Forward >= DeltaAnalysisThresholdEfficiency &
                                                Ratio.H.L.Reverse >= DeltaAnalysisThresholdEfficiency))


## PTMs
# Reproducible
PTM_summ_quant <- subset(PTM_summ, !is.na(Ratio.H.L.Mean))
PTM_summ_HighEff <- subset(PTM_summ_quant, (Ratio.H.L.Forward >= DeltaAnalysisThresholdEfficiency &
                                               Ratio.H.L.Reverse >= DeltaAnalysisThresholdEfficiency))

# All
PTM_summ_raw_quant <- subset(PTM_summ_raw, !is.na(Ratio.H.L.Mean))
PTM_summ_raw_HighEff <- subset(PTM_summ_raw_quant, (Ratio.H.L.Forward >= DeltaAnalysisThresholdEfficiency &
                                               Ratio.H.L.Reverse >= DeltaAnalysisThresholdEfficiency))


# Total proteins QUANTIFIED
cat(" Quantified", "\n",
    "Reproducible", "\t", nrow(PG_summ_quant), "\n",
    "All", "\t", "\t", nrow(PG_summ_raw_quant))

# Total proteins QUANTIFIED above threshold
cat(" Quantified above threshold", "\n",
    "Reproducible", "\t", nrow(PG_summ_HighEff), "\n",
    "All", "\t", "\t", nrow(PG_summ_raw_HighEff))


# Total PTMs QUANTIFIED
cat(" Quantified", "\n",
    "Reproducible", "\t", nrow(PTM_summ_quant), "\n",
    "All", "\t", "\t", nrow(PTM_summ_raw_quant))

# Total PTMs QUANTIFIED above threshold
cat(" Quantified above threshold", "\n",
    "Reproducible", "\t", nrow(PTM_summ_HighEff), "\n",
    "All", "\t", "\t", nrow(PTM_summ_raw_HighEff))

rm(PG_summ_raw, PG_summ_quant, PG_summ_HighEff, PG_summ_raw_HighEff, PG_summ_raw_quant,
   PTM_summ_HighEff, PTM_summ_quant, PTM_summ_raw, PTM_summ_raw_HighEff, PTM_summ_raw_quant)



#### Working table ####

# Add tag for later manipulation of tables
names(PG_summ) <- paste0(names(PG_summ), "_PG")
names(Peptides_summ) <- paste0(names(Peptides_summ), "_PEP")
names(PTM_summ) <- paste0(names(PTM_summ), "_PTM")

## Generate working_table containing all information necessary for plots with Phosphopeptides and Proteins

#One alternative is to merge by Uniprot.ID
working_table_PTM <- merge(PTM_summ, PG_summ, by.x = "Protein_PTM", by.y = "Majority.protein.IDs_PG")
working_table_PEP <- merge(Peptides_summ, PG_summ, by.x = "Proteins_PEP", by.y = "Majority.protein.IDs_PG")






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

ratios_unmod_pep <- c("Ratio.H.L.unmod..pep..Forward_II",
                      "Ratio.H.L.unmod..pep..Forward_III",
                      "Ratio.H.L.unmod..pep..Forward_IV",
                      "Ratio.H.L.unmod..pep..Reverse_II",
                      "Ratio.H.L.unmod..pep..Reverse_III",
                      "Ratio.H.L.unmod..pep..Reverse_IV",
                      "Ratio.H.L.unmod..pep..Mean",
                      
                      
                      "Ratio.H.M.unmod..pep..Forward_II",
                      "Ratio.H.M.unmod..pep..Reverse_II",
                      
                      "Ratio.M.L.unmod..pep..Forward_II",
                      "Ratio.M.L.unmod..pep..Reverse_II")


intensities_PG <- paste0(intensities, "_PG")
intensities_PEP <- paste0(intensities, "_PEP")
intensities_PTM <- paste0(intensities, "_PTM")

ratios_PG <- paste0(ratios, "_PG")
ratios_PEP <- paste0(ratios, "_PEP")
ratios_PTM <- paste0(ratios, "_PTM")

ratios_unmod_pep_PTM <- paste0(ratios_unmod_pep, "_PTM")



rm(ratios, intensities, ratios_unmod_pep)




#### Simplify tables ####
working_table_PEP <- subset(working_table_PEP, select = c("Proteins_PEP",
                                                          "Gene.names_PEP",
                                                          "Sequence_PEP",
                                                          "Sequence.length_PG",
                                                          "Score_PEP",
                                                          "Start.position_PEP",
                                                          "End.position_PEP",
                                                          "Known_RBP_PG",
                                                          "Known_mRBP_PG",
                                                          "Unique..Groups._PEP",
                                                          "Unique..Proteins._PEP",
                                                          "Mol..weight..kDa._PG",
                                                          "Sequence.coverage...._PG",
                                                          "Sequence.length_PG",
                                                          ratios_PEP, ratios_PG))







working_table_PTM <- subset(working_table_PTM, select = c("Protein_PTM",
                                                          "Gene.names_PTM",
                                                          "Peptide.Site.ID_PTM",
                                                          "Site.ID_PTM",
                                                          "Amino.acid_PTM",
                                                          "Position_PTM",
                                                          "Known_RBP_PG",
                                                          "Known_mRBP_PG",
                                                          "Position_PTM",
                                                          "Mol..weight..kDa._PG",
                                                          "Sequence.coverage...._PG",
                                                          "Sequence.length_PG",
                                                          ratios_PTM, ratios_PG, ratios_unmod_pep_PTM))








#### Mean pulldown efficiency calculation ####
## I'm repeating here the calculation made in the other scripts for each PG, PTM and PEP tables. I repeat it because non-reproducible values might have been removed 

## Reverse ratios were already inverted before!!

#Calculate mean making sure protein is identified in at least >1 replicates
MyMeanWithMinValues <- function(df, columns, MinNumberReplicates = 1) {
  foo <- apply(df[columns],
               1, 
               function(x) ifelse(sum(!is.na(x)) > MinNumberReplicates - 1,
                                  mean(x, na.rm = T),
                                  NA))
  return(foo)
}



#PEP PEP
working_table_PEP$Ratio.H.L.Forward_PEP <- MyMeanWithMinValues(working_table_PEP, c("Ratio.H.L.Forward_II_PEP",
                                                            "Ratio.H.L.Forward_III_PEP",
                                                            "Ratio.H.L.Forward_IV_PEP"),
                                                 MinNumberReplicates = 1)

working_table_PEP$Ratio.H.L.Reverse_PEP <- MyMeanWithMinValues(working_table_PEP, c("Ratio.H.L.Reverse_II_PEP",
                                                            "Ratio.H.L.Reverse_III_PEP",
                                                            "Ratio.H.L.Reverse_IV_PEP"),
                                                 MinNumberReplicates = 1)

working_table_PEP$Ratio.H.L.Mean_PEP <- MyMeanWithMinValues(working_table_PEP, c("Ratio.H.L.Forward_PEP",
                                                         "Ratio.H.L.Reverse_PEP"),
                                              MinNumberReplicates = 2)



working_table_PEP$Ratio.NoCl.Mean_PEP <- MyMeanWithMinValues(working_table_PEP, c("Ratio.M.L.Forward_II_PEP",
                                                                                  "Ratio.H.M.Reverse_II_PEP"),
                                                             MinNumberReplicates = 2)


working_table_PEP$Ratio.ClXNoCl.Mean_PEP <- MyMeanWithMinValues(working_table_PEP, c("Ratio.H.M.Forward_II_PEP",
                                                             "Ratio.M.L.Reverse_II_PEP"),
                                                  MinNumberReplicates = 2)

#PEP PG
working_table_PEP$Ratio.H.L.Forward_PG <- MyMeanWithMinValues(working_table_PEP, c("Ratio.H.L.Forward_II_PG",
                                                                                    "Ratio.H.L.Forward_III_PG",
                                                                                    "Ratio.H.L.Forward_IV_PG"),
                                                               MinNumberReplicates = 1)

working_table_PEP$Ratio.H.L.Reverse_PG <- MyMeanWithMinValues(working_table_PEP, c("Ratio.H.L.Reverse_II_PG",
                                                                                    "Ratio.H.L.Reverse_III_PG",
                                                                                    "Ratio.H.L.Reverse_IV_PG"),
                                                               MinNumberReplicates = 1)

working_table_PEP$Ratio.H.L.Mean_PG <- MyMeanWithMinValues(working_table_PEP, c("Ratio.H.L.Forward_PG",
                                                                                 "Ratio.H.L.Reverse_PG"),
                                                            MinNumberReplicates = 2)



working_table_PEP$Ratio.NoCl.Mean_PG <- MyMeanWithMinValues(working_table_PEP, c("Ratio.M.L.Forward_II_PG",
                                                                                  "Ratio.H.M.Reverse_II_PG"),
                                                             MinNumberReplicates = 2)


working_table_PEP$Ratio.ClXNoCl.Mean_PG <- MyMeanWithMinValues(working_table_PEP, c("Ratio.H.M.Forward_II_PG",
                                                                                     "Ratio.M.L.Reverse_II_PG"),
                                                                MinNumberReplicates = 2)






#PTM PTM
working_table_PTM$Ratio.H.L.Forward_PTM <- MyMeanWithMinValues(working_table_PTM, c("Ratio.H.L.Forward_II_PTM",
                                                                                    "Ratio.H.L.Forward_III_PTM",
                                                                                    "Ratio.H.L.Forward_IV_PTM"),
                                                               MinNumberReplicates = 1)

working_table_PTM$Ratio.H.L.Reverse_PTM <- MyMeanWithMinValues(working_table_PTM, c("Ratio.H.L.Reverse_II_PTM",
                                                                                    "Ratio.H.L.Reverse_III_PTM",
                                                                                    "Ratio.H.L.Reverse_IV_PTM"),
                                                               MinNumberReplicates = 1)

working_table_PTM$Ratio.H.L.Mean_PTM <- MyMeanWithMinValues(working_table_PTM, c("Ratio.H.L.Forward_PTM",
                                                                                 "Ratio.H.L.Reverse_PTM"),
                                                            MinNumberReplicates = 2)



working_table_PTM$Ratio.NoCl.Mean_PTM <- MyMeanWithMinValues(working_table_PTM, c("Ratio.M.L.Forward_II_PTM",
                                                                                  "Ratio.H.M.Reverse_II_PTM"),
                                                             MinNumberReplicates = 2)


working_table_PTM$Ratio.ClXNoCl.Mean_PTM <- MyMeanWithMinValues(working_table_PTM, c("Ratio.H.M.Forward_II_PTM",
                                                                                     "Ratio.M.L.Reverse_II_PTM"),
                                                                MinNumberReplicates = 2)




#PTM PG
working_table_PTM$Ratio.H.L.Forward_PG <- MyMeanWithMinValues(working_table_PTM, c("Ratio.H.L.Forward_II_PG",
                                                                                   "Ratio.H.L.Forward_III_PG",
                                                                                   "Ratio.H.L.Forward_IV_PG"),
                                                              MinNumberReplicates = 1)

working_table_PTM$Ratio.H.L.Reverse_PG <- MyMeanWithMinValues(working_table_PTM, c("Ratio.H.L.Reverse_II_PG",
                                                                                   "Ratio.H.L.Reverse_III_PG",
                                                                                   "Ratio.H.L.Reverse_IV_PG"),
                                                              MinNumberReplicates = 1)

working_table_PTM$Ratio.H.L.Mean_PG <- MyMeanWithMinValues(working_table_PTM, c("Ratio.H.L.Forward_PG",
                                                                                "Ratio.H.L.Reverse_PG"),
                                                           MinNumberReplicates = 2)



working_table_PTM$Ratio.NoCl.Mean_PG <- MyMeanWithMinValues(working_table_PTM, c("Ratio.M.L.Forward_II_PG",
                                                                                 "Ratio.H.M.Reverse_II_PG"),
                                                            MinNumberReplicates = 2)


working_table_PTM$Ratio.ClXNoCl.Mean_PG <- MyMeanWithMinValues(working_table_PTM, c("Ratio.H.M.Forward_II_PG",
                                                                                    "Ratio.M.L.Reverse_II_PG"),
                                                               MinNumberReplicates = 2)









#### Delta calculations and Statistics ####
## Reverse ratios were already inverted!!

#PEP
# Relative funSILAC
working_table_PEP$delta_Forward_II <- working_table_PEP$Ratio.H.L.Forward_II_PEP - working_table_PEP$Ratio.H.L.Forward_II_PG

working_table_PEP$delta_Forward_III <- working_table_PEP$Ratio.H.L.Forward_III_PEP - working_table_PEP$Ratio.H.L.Forward_III_PG

working_table_PEP$delta_Forward_IV <- working_table_PEP$Ratio.H.L.Forward_IV_PEP - working_table_PEP$Ratio.H.L.Forward_IV_PG



working_table_PEP$delta_Reverse_II <- working_table_PEP$Ratio.H.L.Reverse_II_PEP - working_table_PEP$Ratio.H.L.Reverse_II_PG

working_table_PEP$delta_Reverse_III <- working_table_PEP$Ratio.H.L.Reverse_III_PEP - working_table_PEP$Ratio.H.L.Reverse_III_PG

working_table_PEP$delta_Reverse_IV <- working_table_PEP$Ratio.H.L.Reverse_IV_PEP - working_table_PEP$Ratio.H.L.Reverse_IV_PG




# Doing the calculation like this makes sure that the final value has been calculated from at least one value in For and one in Rev experiments
working_table_PEP$delta_Forward <- MyMeanWithMinValues(working_table_PEP, c("delta_Forward_II",
                                                                                        "delta_Forward_III",
                                                                                        "delta_Forward_IV"),
                                                                   MinNumberReplicates = 1)


working_table_PEP$delta_Reverse <- MyMeanWithMinValues(working_table_PEP, c("delta_Reverse_II",
                                                                                        "delta_Reverse_III",
                                                                                        "delta_Reverse_IV"),
                                                                   MinNumberReplicates = 1)

working_table_PEP$delta_Mean <- MyMeanWithMinValues(working_table_PEP, c("delta_Forward",
                                                                         "delta_Reverse"),
                                                    MinNumberReplicates = 2)


working_table_PEP$delta_Mean_all <- 
  MyMeanWithMinValues(working_table_PEP, c("delta_Forward_II",
                                           "delta_Forward_III",
                                           "delta_Forward_IV",
                                           "delta_Reverse_II",
                                           "delta_Reverse_III",
                                           "delta_Reverse_IV"),
                      MinNumberReplicates = 2)










#PTM
# Relative funSILAC
working_table_PTM$delta_Forward_II <- working_table_PTM$Ratio.H.L.Forward_II_PTM - working_table_PTM$Ratio.H.L.Forward_II_PG

working_table_PTM$delta_Forward_III <- working_table_PTM$Ratio.H.L.Forward_III_PTM - working_table_PTM$Ratio.H.L.Forward_III_PG

working_table_PTM$delta_Forward_IV <- working_table_PTM$Ratio.H.L.Forward_IV_PTM - working_table_PTM$Ratio.H.L.Forward_IV_PG



working_table_PTM$delta_Reverse_II <- working_table_PTM$Ratio.H.L.Reverse_II_PTM - working_table_PTM$Ratio.H.L.Reverse_II_PG

working_table_PTM$delta_Reverse_III <- working_table_PTM$Ratio.H.L.Reverse_III_PTM - working_table_PTM$Ratio.H.L.Reverse_III_PG

working_table_PTM$delta_Reverse_IV <- working_table_PTM$Ratio.H.L.Reverse_IV_PTM - working_table_PTM$Ratio.H.L.Reverse_IV_PG



# Doing the calculation like this makes sure that the final value has been calculated from at least one value in For and one in Rev experiments
working_table_PTM$delta_Forward <- MyMeanWithMinValues(working_table_PTM, c("delta_Forward_II",
                                                                                        "delta_Forward_III",
                                                                                        "delta_Forward_IV"),
                                                                   MinNumberReplicates = 1)


working_table_PTM$delta_Reverse <- MyMeanWithMinValues(working_table_PTM, c("delta_Reverse_II",
                                                                                        "delta_Reverse_III",
                                                                                        "delta_Reverse_IV"),
                                                                   MinNumberReplicates = 1)

working_table_PTM$delta_Mean <- MyMeanWithMinValues(working_table_PTM, c("delta_Forward",
                                                                                     "delta_Reverse"),
                                                                MinNumberReplicates = 2)


working_table_PTM$delta_Mean_all <-
  MyMeanWithMinValues(working_table_PTM, c("delta_Forward_II",
                                           "delta_Forward_III",
                                           "delta_Forward_IV",
                                           "delta_Reverse_II",
                                           "delta_Reverse_III",
                                           "delta_Reverse_IV"), 
                      MinNumberReplicates = 2)



# Relative funSILAC Unmod peptides
working_table_PTM$delta_unmod_Forward_II <- working_table_PTM$Ratio.H.L.unmod..pep..Forward_II_PTM - working_table_PTM$Ratio.H.L.Forward_II_PG

working_table_PTM$delta_unmod_Forward_III <- working_table_PTM$Ratio.H.L.unmod..pep..Forward_III_PTM - working_table_PTM$Ratio.H.L.Forward_III_PG

working_table_PTM$delta_unmod_Forward_IV <- working_table_PTM$Ratio.H.L.unmod..pep..Forward_IV_PTM - working_table_PTM$Ratio.H.L.Forward_IV_PG



working_table_PTM$delta_unmod_Reverse_II <- working_table_PTM$Ratio.H.L.unmod..pep..Reverse_II_PTM - working_table_PTM$Ratio.H.L.Reverse_II_PG

working_table_PTM$delta_unmod_Reverse_III <- working_table_PTM$Ratio.H.L.unmod..pep..Reverse_III_PTM - working_table_PTM$Ratio.H.L.Reverse_III_PG

working_table_PTM$delta_unmod_Reverse_IV <- working_table_PTM$Ratio.H.L.unmod..pep..Reverse_IV_PTM - working_table_PTM$Ratio.H.L.Reverse_IV_PG





# Doing the calculation like this makes sure that the final value has been calculated from at least one value in For and one in Rev experiments
working_table_PTM$delta_unmod_Forward <- MyMeanWithMinValues(working_table_PTM, c("delta_unmod_Forward_II",
                                                                                        "delta_unmod_Forward_III",
                                                                                        "delta_unmod_Forward_IV"),
                                                                   MinNumberReplicates = 1)


working_table_PTM$delta_unmod_Reverse <- MyMeanWithMinValues(working_table_PTM, c("delta_unmod_Reverse_II",
                                                                                        "delta_unmod_Reverse_III",
                                                                                        "delta_unmod_Reverse_IV"),
                                                                   MinNumberReplicates = 1)

working_table_PTM$delta_unmod_Mean <- MyMeanWithMinValues(working_table_PTM, c("delta_unmod_Forward",
                                                                                     "delta_unmod_Reverse"),
                                                                MinNumberReplicates = 2)










#### Delta efficiency p-values ####
### eBayes t-test
if (p.val.method == "eBayes") {
  MyModeratedTTest <- function(df = working_table_PTM, 
                               NumbTreatReps = 6, NumbCtlReps = 6,
                               Sample.type = "PTM",
                               Reps.value = c(ratios_PTM[grepl("H.L.", ratios_PTM)],
                                              ratios_PG[grepl("H.L.", ratios_PG)]),
                               Moderated.logic = F){
    
    # Moderated t-test (by Fritz)
    #Smyth, G. K. (2004). Linear Models and Empirical Bayes Methods for Assessing Differential Expression in Microarray Experiments. Statistical Applications in Genetics and Molecular Biology, 3(1), 1-25. doi:10.2202/1544-6115.1027 
    samples = factor(c(rep(Sample.type, NumbTreatReps), rep("PG", NumbCtlReps)))
    design = model.matrix(~ samples)
    fit = lmFit(as.matrix(df[Reps.value]), design)
    fit_ebayes = eBayes(fit,trend=Moderated.logic)
    # trend=T: intensity-trend allowed for the prior variance. Default is trend=F, i.e. prior variance is constant.
    
    #Check trend in data if in doubt about trend-variable:
    # plot(apply(df[Reps.value],1,mean,na.rm=T), apply(df[Reps.value],1,sd,na.rm=T))
    # plotSA(fit_ebayes)
    
    #Interpret level of moderation:
    # mean(fit_ebayes$df.prior/(fit_ebayes$df.prior+fit_ebayes$df.residual))
    #0.2 = different; 0.5 = balance; 0.99 = similar (interpretation from Smyth 2004 paper)
    
    #Volcanoplot
    # plot(fit_ebayes$coefficients[,2],
    #      -log10(p.adjust(fit_ebayes$p.value[,2], "BH")),
    #      col=(p.adjust(fit_ebayes$p.value[,2],"BH")>0.05)*4+4,pch=20,
    #      cex=1.1-log10(length(na.omit(fit_ebayes$p.value[,2])))/10,
    #     xlab="log2(fold-change)",ylab="-log10(p-value)",
    #     main=paste("Volcano (MTC=BH colouring cut-off=0.05)\n",colnames(fit_ebayes$coefficients)[2]))
    
    foo <- topTable(fit_ebayes,coef=2,number=nrow(df))
    foo$Neglog10pval_Cl_Mean <- -log10(foo$P.Value)
    foo$Neglog10padj_Cl_Mean <- -log10(foo$adj.P.Val)
    foo <- subset(foo, select = c(#"logFC",
                                  "P.Value", "Neglog10pval_Cl_Mean", "adj.P.Val", "Neglog10padj_Cl_Mean"))
    
    names(foo) <- c(#"delta_Mean",
                    "pval_Cl_Mean", "Neglog10pval_Cl_Mean", "padj_Cl_Mean", "Neglog10padj_Cl_Mean")
    return(foo)
    
  }
  
  working_table_PEP <- merge(working_table_PEP, MyModeratedTTest(working_table_PEP, Sample.type = "PEP", Reps.value = c(ratios_PEP[grepl("H.L.", ratios_PEP)], ratios_PG[grepl("H.L.", ratios_PG)])),by = "row.names")
  
  plot(working_table_PEP$delta_Mean, working_table_PEP$Neglog10pval_Cl_Mean)
  
  working_table_PTM <- merge(working_table_PTM, MyModeratedTTest(working_table_PTM, Sample.type = "PTM"), by = "row.names")
  
  plot(working_table_PTM$delta_Mean, working_table_PTM$Neglog10pval_Cl_Mean)
  
  working_table_PTM$delta_unmod_Mean <- MyModeratedTTest(working_table_PTM, Sample.type = "PTM",
                                                         Reps.value = c(ratios_unmod_pep_PTM[grepl("H.L.", ratios_unmod_pep_PTM) & !(ratios_unmod_pep_PTM %in% "Ratio.H.L.unmod..pep..Mean_PTM")],
                                                                        ratios_PG[grepl("H.L.", ratios_PG)]))[,"delta_Mean"]
  
}

### One Sample T-test
if (p.val.method == "one_sample_T_test"){
  MyOneSampleTTest <- function(df, columns, MinNumberReplicates = 2) {
    foo <- apply(df[columns],
                 1, 
                 function(x) ifelse(sum(!is.na(x)) > MinNumberReplicates - 1,
                                    t.test(x, alternative = "two.sided", mu = 0,
                                           na.rm = T)$p.value,
                                    NA))
    return(foo)
  }
  
  #PEP
  working_table_PEP$pval_Cl_Mean <- MyOneSampleTTest(working_table_PEP,
                                                     c("delta_Forward_II",
                                                       "delta_Forward_III",
                                                       "delta_Forward_IV",
                                                       "delta_Reverse_II",
                                                       "delta_Reverse_III",
                                                       "delta_Reverse_IV"))
  
  
  working_table_PEP$Neglog10pval_Cl_Mean <- -log10(working_table_PEP$pval_Cl_Mean)
  working_table_PEP$padj_Cl_Mean <- p.adjust(working_table_PEP$pval_Cl_Mean, method="fdr")
  working_table_PEP$Neglog10padj_Cl_Mean <- -log10(working_table_PEP$padj_Cl_Mean)
  
  plot(working_table_PEP$delta_Mean, working_table_PEP$Neglog10pval_Cl_Mean)
  
  
  #PTM
  working_table_PTM$pval_Cl_Mean <- MyOneSampleTTest(working_table_PTM, c("delta_Forward_II",
                                                                          "delta_Forward_III",
                                                                          "delta_Forward_IV",
                                                                          "delta_Reverse_II",
                                                                          "delta_Reverse_III",
                                                                          "delta_Reverse_IV"))
  
  
  working_table_PTM$Neglog10pval_Cl_Mean <- -log10(working_table_PTM$pval_Cl_Mean)
  working_table_PTM$padj_Cl_Mean <- p.adjust(working_table_PTM$pval_Cl_Mean, method="fdr")
  working_table_PTM$Neglog10padj_Cl_Mean <- -log10(working_table_PTM$padj_Cl_Mean)
  
  plot(working_table_PTM$delta_Mean, working_table_PTM$Neglog10pval_Cl_Mean)
}

### Unpaired T-test
if (p.val.method == "Unpaired_T_test"){
  MyUnpairedTTest <- function(df = working_table_PTM, 
                              Sample.type = "PTM",
                              columns = c(eval(parse(text = paste0("ratios_", Sample.type)))[grepl("H.L.", eval(parse(text = paste0("ratios_", Sample.type))))],
                                          ratios_PG[grepl("H.L.", ratios_PG)]),
                              MinNumberReplicates = 2) {
    foo <- apply(df[columns],
                 1, 
                 function(x) ifelse(sum(!is.na(x[1:(length(columns)/2)])) > MinNumberReplicates - 1 &
                                      sum(!is.na(x[(1+(length(columns)/2)):length(columns)])) > MinNumberReplicates - 1,
                                    t.test(x = x[1:(length(columns)/2)],
                                           y = x[(1+(length(columns)/2)):length(columns)],
                                           alternative = "two.sided", mu = 0, paired = F,
                                           na.rm = T)$p.value,
                                    NA))
    return(foo)
  }
  
  
  
  #PEP
  working_table_PEP$pval_Cl_Mean <- MyUnpairedTTest(df = working_table_PEP, Sample.type = "PEP")
  
  working_table_PEP$Neglog10pval_Cl_Mean <- -log10(working_table_PEP$pval_Cl_Mean)
  working_table_PEP$padj_Cl_Mean <- p.adjust(working_table_PEP$pval_Cl_Mean, method="fdr")
  working_table_PEP$Neglog10padj_Cl_Mean <- -log10(working_table_PEP$padj_Cl_Mean)
  
  plot(working_table_PEP$delta_Mean, working_table_PEP$Neglog10pval_Cl_Mean)
  
  
  #PTM
  working_table_PTM$pval_Cl_Mean <- MyUnpairedTTest(df = working_table_PTM, Sample.type = "PTM")
  
  working_table_PTM$Neglog10pval_Cl_Mean <- -log10(working_table_PTM$pval_Cl_Mean)
  working_table_PTM$padj_Cl_Mean <- p.adjust(working_table_PTM$pval_Cl_Mean, method="fdr")
  working_table_PTM$Neglog10padj_Cl_Mean <- -log10(working_table_PTM$padj_Cl_Mean)
  
  plot(working_table_PTM$delta_Mean, working_table_PTM$Neglog10pval_Cl_Mean)
  
}


#### Defining significant changes ####
#### Selection criteria: delta efficiency higher than 2 fold in mean
working_table_PEP$PEP_increase_Efficiency_Cl_Mean <-
  ifelse(working_table_PEP$delta_Mean >= log2(2), #&
           # ((working_table_PEP$Ratio.H.L.Forward_PG >= DeltaAnalysisThresholdEfficiency &
           #     working_table_PEP$Ratio.H.L.Reverse_PG >= DeltaAnalysisThresholdEfficiency) | 
           #    (working_table_PEP$Ratio.H.L.Forward_PEP >= DeltaAnalysisThresholdEfficiency &
           #       working_table_PEP$Ratio.H.L.Reverse_PEP >= DeltaAnalysisThresholdEfficiency)),
         T, F)
working_table_PEP["PEP_increase_Efficiency_Cl_Mean"] <- apply(working_table_PEP["PEP_increase_Efficiency_Cl_Mean"], 1, function(x) ifelse(is.na(x), FALSE, x))


working_table_PEP$PEP_decrease_Efficiency_Cl_Mean <-
  ifelse(working_table_PEP$delta_Mean <= -log2(2), #&
           # ((working_table_PEP$Ratio.H.L.Forward_PG >= DeltaAnalysisThresholdEfficiency &
           #     working_table_PEP$Ratio.H.L.Reverse_PG >= DeltaAnalysisThresholdEfficiency) | 
           #    (working_table_PEP$Ratio.H.L.Forward_PEP >= DeltaAnalysisThresholdEfficiency &
           #       working_table_PEP$Ratio.H.L.Reverse_PEP >= DeltaAnalysisThresholdEfficiency)),
         T, F)
working_table_PEP["PEP_decrease_Efficiency_Cl_Mean"] <- apply(working_table_PEP["PEP_decrease_Efficiency_Cl_Mean"], 1, function(x) ifelse(is.na(x), FALSE, x))






#PTM
#### Selection criteria: delta efficiency higher than 2 fold in mean
working_table_PTM$PTM_increase_Efficiency_Cl_Mean <-
  ifelse(working_table_PTM$delta_Mean >= log2(2), #&
           # ((working_table_PTM$Ratio.H.L.Forward_PG >= DeltaAnalysisThresholdEfficiency &
           #     working_table_PTM$Ratio.H.L.Reverse_PG >= DeltaAnalysisThresholdEfficiency) | 
           #    (working_table_PTM$Ratio.H.L.Forward_PTM >= DeltaAnalysisThresholdEfficiency &
           #       working_table_PTM$Ratio.H.L.Reverse_PTM >= DeltaAnalysisThresholdEfficiency)),
         T, F)
working_table_PTM["PTM_increase_Efficiency_Cl_Mean"] <- apply(working_table_PTM["PTM_increase_Efficiency_Cl_Mean"], 1, function(x) ifelse(is.na(x), FALSE, x))


working_table_PTM$PTM_decrease_Efficiency_Cl_Mean <-
  ifelse(working_table_PTM$delta_Mean <= -log2(2),# &
           # ((working_table_PTM$Ratio.H.L.Forward_PG >= DeltaAnalysisThresholdEfficiency &
           #     working_table_PTM$Ratio.H.L.Reverse_PG >= DeltaAnalysisThresholdEfficiency) | 
           #    (working_table_PTM$Ratio.H.L.Forward_PTM >= DeltaAnalysisThresholdEfficiency &
           #       working_table_PTM$Ratio.H.L.Reverse_PTM >= DeltaAnalysisThresholdEfficiency)),
         T, F)
working_table_PTM["PTM_decrease_Efficiency_Cl_Mean"] <- apply(working_table_PTM["PTM_decrease_Efficiency_Cl_Mean"], 1, function(x) ifelse(is.na(x), FALSE, x))





#### Write table ####
fwrite(working_table_PTM, file = "working_table_PTM.txt", sep = "\t", na = "", quote = F, row.names = F)
fwrite(working_table_PEP, file = "working_table_PEP.txt", sep = "\t", na = "", quote = F, row.names = F)












#### Writing table ratios paper ####
# Table 3 paper
foo_df <- subset(working_table_PTM, select = c("Protein_PTM",
                                      "Gene.names_PTM",
                                      "Site.ID_PTM",
                                      "Known_mRBP_PG",
                                      
                                      "Ratio.H.L.Mean_PTM",
                                      "Ratio.H.L.Mean_PG",
                                      
                                      "delta_Forward_II",
                                      "delta_Forward_III",
                                      "delta_Forward_IV",
                                      "delta_Reverse_II",
                                      "delta_Reverse_III",
                                      "delta_Reverse_IV",
                                      "delta_Forward",
                                      "delta_Reverse",
                                      "delta_Mean",
                                      "pval_Cl_Mean"))

## Calculate Percentages
foo_df[,c("Ratio.H.L.Mean_PTM",
          "Ratio.H.L.Mean_PG")] <- 2**foo_df[,c("Ratio.H.L.Mean_PTM",
                                                "Ratio.H.L.Mean_PG")]


# Adjust colum names
names(foo_df) <- c("Uniprot.IDs",
                   "Gene.names",
                   "Site.ID",
                   "Known.mRBP.Gebauer",
                   
                   "Phospho.Pulldown.efficiency(%)",
                   "Protein.Pulldown.efficiency(%)",
                   
                   "delta.efficiency.Forward_Exp1",
                   "delta.efficiency.Forward_Exp2",
                   "delta.efficiency.Forward_Exp3",
                   "delta.efficiency.Reverse_Exp1",
                   "delta.efficiency.Reverse_Exp2",
                   "delta.efficiency.Reverse_Exp3",
                   
                   "delta.efficiency.Forward",
                   "delta.efficiency.Reverse",
                   
                   "delta.efficiency.Mean",
                   "OneSampleTTest.p.value")


write.table(foo_df, file = "Table3_PhosphoDeltaEfficiency.txt", sep = "\t", na = "", quote = F, row.names = F)



#### Venn diagrams ####
PG_summ_filtered <- subset(PG_summ,
                           (Ratio.H.L.Forward_PG >= DeltaAnalysisThresholdEfficiency &
                              Ratio.H.L.Reverse_PG >= DeltaAnalysisThresholdEfficiency) &
                             !is.na(Ratio.H.L.Mean_PG))


#colect gene names
foo <- PG_summ
foo$Gene.names <- ifelse(is.na(foo$Gene.names) | foo$Gene.names == "", foo$Majority.protein.IDs, foo$Gene.names)

All_identified <- unique(foo$Gene.names)
All_quantified <- unique(subset(foo, !is.na(Ratio.H.L.Mean_PG))$Gene.names)
rm(foo)


High_eff_mRBPs <- unique(PG_summ_filtered$Gene.names)


MyVennDiagram <- function(A, B, C = NULL) {
  require(eulerr)
  
  
  if (is.null(C)) {
    A_B = length(intersect(A, B))
    len_A = length(A) - A_B
    len_B = length(B) - A_B
    
    venn_plot <- euler(c(A = len_A, B = len_B, "A&B" = A_B))
  } else {
    
    A_B_noC = length(intersect(A, B)) - length(intersect(A,intersect(B, C)))
    A_C_noB = length(intersect(A, C)) - length(intersect(A,intersect(B, C)))
    B_C_noA = length(intersect(B, C)) - length(intersect(A,intersect(B, C)))
    A_B_C = length(intersect(A,intersect(B, C)))
    
    len_A = length(A) - A_B_C - A_B_noC - A_C_noB
    len_B = length(B) - A_B_C - A_B_noC - B_C_noA
    len_C = length(C) - A_B_C - A_C_noB - B_C_noA
    
    venn_plot <- euler(c(A = len_A, B = len_B, C = len_C,
                         "A&B" = A_B_noC, "A&C" = A_C_noB, "B&C" = B_C_noA,
                         "A&B&C" = A_B_C))
  }
  
  return(venn_plot)
  
}

pdf(file = "Venn_plot_qRIC.pdf", height = 5, width = 5, useDingbats = F)
plot(MyVennDiagram(A = High_eff_mRBPs,
                   B = mRBPs_HEK,
                   C = mRBPs),
     legend = TRUE, quantities= TRUE)
dev.off()

pdf(file = "Venn_plot_Allquant.pdf", height = 5, width = 5, useDingbats = F)
plot(MyVennDiagram(A = All_quantified,
                   B = mRBPs_HEK,
                   C = mRBPs),
     legend = TRUE, quantities= TRUE)
dev.off()

# Percent of all mRBPs
100*length(High_eff_mRBPs)/length(mRBPs)

# Percent of all HEK mRBPs
100*length(High_eff_mRBPs)/length(mRBPs_HEK)

rm(PG_summ_filtered, High_eff_mRBPs, All_identified, All_quantified)




#### Numbers ####
PG_summ_all <- read.csv("PG_summary.txt", sep = "\t", stringsAsFactors = FALSE)
PG_summ_quant <- subset(PG_summ_all, !is.na(Ratio.H.L.Mean))
PG_summ_all_HighEff <- subset(PG_summ_all, 
                              (Ratio.H.L.Forward >= DeltaAnalysisThresholdEfficiency &
                                 Ratio.H.L.Reverse >= DeltaAnalysisThresholdEfficiency))

PG_summ_all_mRBPs <- subset(PG_summ_all, !is.na(Ratio.H.L.Mean) & Known_mRBP)
PG_summ_all_HEK <- subset(PG_summ_all, !is.na(Ratio.H.L.Mean) & Known_mRBP_HEK)


PTM_summ_all <- read.csv("PTM_summary.txt", sep = "\t", stringsAsFactors = FALSE)
PTM_summ_quant <- subset(PTM_summ_all, !is.na(Ratio.H.L.Mean))
PTM_summ_all_HighEff <- subset(PTM_summ_all,
                               (Ratio.H.L.Forward >= DeltaAnalysisThresholdEfficiency &
                                  Ratio.H.L.Reverse >= DeltaAnalysisThresholdEfficiency))

PTM_summ_all_mRBPs <- subset(PTM_summ_all, !is.na(Ratio.H.L.Mean) & Known_mRBP)
PTM_summ_all_HEK <- subset(PTM_summ_all, !is.na(Ratio.H.L.Mean) & Known_mRBP_HEK)


# Total proteins IDENTIFIED
nrow(PG_summ_all)

# Total proteins QUANTIFIED
nrow(PG_summ_quant)

# Total proteins QUANTIFIED in mRBPs
nrow(PG_summ_all_mRBPs)

# Total proteins QUANTIFIED above threshold
nrow(PG_summ_all_HighEff)

# HEK-specific mRBPs quantified in qRIC
mRBPs_HEK_quantified <- unique(mRBPs_HEK[(mRBPs_HEK %in% PG_summ_quant$Gene.names)])
length(mRBPs_HEK_quantified)

# mRBPs quantified in qRIC
mRBPs_quantified <- unique(mRBPs[(mRBPs %in% PG_summ_quant$Gene.names)])
length(mRBPs_quantified)



# Total p-sites IDENTIFIED
nrow(PTM_summ_all)

# Total p-sites QUANTIFIED
nrow(PTM_summ_quant)

# Total p-sites QUANTIFIED in mRBPs
nrow(PTM_summ_all_mRBPs)

# Total p-sites QUANTIFIED above threshold
nrow(PTM_summ_all_HighEff)

# Number of proteins from p-sites QUANTIFIED in mRBPs
length(unique(PTM_summ_all_mRBPs$Protein))


# Number of phosphosites in delta efficiency analysis
working_table_PTM_quant <- subset(working_table_PTM, Known_mRBP_PG & !is.na(delta_Mean))
nrow(working_table_PTM_quant)
# Number of proteins from phosphosites in delta efficiency analysis
length(unique(working_table_PTM_quant$Protein_PTM))

## High efficiency
working_table_PTM_HighEff <- 
  subset(working_table_PTM_quant, Known_mRBP_PG & 
           ((Ratio.H.L.Forward_PG >= DeltaAnalysisThresholdEfficiency &
               Ratio.H.L.Reverse_PG >= DeltaAnalysisThresholdEfficiency) | 
              (Ratio.H.L.Forward_PTM >= DeltaAnalysisThresholdEfficiency &
                 Ratio.H.L.Reverse_PTM >= DeltaAnalysisThresholdEfficiency)))

nrow(working_table_PTM_HighEff)

# Number of proteins from phosphosites in delta efficiency analysis
length(unique(working_table_PTM_HighEff$Protein_PTM))



# Number of pSite per protein
# 1 pSite
sum(table(working_table_PTM_quant$Protein_PTM) == 1)
# 2 pSite
sum(table(working_table_PTM_quant$Protein_PTM) == 2)
# 3 pSite
sum(table(working_table_PTM_quant$Protein_PTM) == 3)
# 4 or more pSite
sum(table(working_table_PTM_quant$Protein_PTM) > 3)




#Number of phosphosites quantified in total
table(PTM_summ_quant$Amino.acid)



# Number of phosphosites in RBPs
all_pSites_inproteins <- working_table_PTM_quant[, c("Site.ID_PTM", "Amino.acid_PTM")]
all_pSites_inproteins$Factor <- as.factor("inProteins")
table(all_pSites_inproteins$Amino.acid_PTM)


#Number of phosphosites changing
table(!is.na(subset(working_table_PTM_quant, PTM_increase_Efficiency_Cl_Mean | PTM_decrease_Efficiency_Cl_Mean)$delta_Mean))
changing_pSites <- subset(working_table_PTM_quant, PTM_increase_Efficiency_Cl_Mean | PTM_decrease_Efficiency_Cl_Mean)[, c("Site.ID_PTM", "Amino.acid_PTM", "PTM_increase_Efficiency_Cl_Mean", "PTM_decrease_Efficiency_Cl_Mean")]
changing_pSites$Factor <- as.factor("Changing")

table(changing_pSites$Amino.acid_PTM)
sum(changing_pSites$PTM_increase_Efficiency_Cl_Mean)
sum(changing_pSites$PTM_decrease_Efficiency_Cl_Mean)







#### Do RBPs below threshold are more abundant (background binding)
# PD efficiency > threshold HEK mRBPs
df_mRBPs_HEK_high <- subset(PG_summ_all, Gene.names %in% mRBPs_HEK  & 
                                 !is.na(Ratio.H.L.Mean) &
                                 (Ratio.H.L.Forward > DeltaAnalysisThresholdEfficiency &
                                    Ratio.H.L.Reverse > DeltaAnalysisThresholdEfficiency))

length(unique(df_mRBPs_HEK_high$Majority.protein.IDs))

# PD efficiency < threhsold HEK mRBPs
df_mRBPs_HEK_below0p5 <- subset(PG_summ_all, Gene.names %in% mRBPs_HEK  & 
                                  !is.na(Ratio.H.L.Mean) &
                                  (Ratio.H.L.Forward < DeltaAnalysisThresholdEfficiency & 
                                     Ratio.H.L.Reverse < DeltaAnalysisThresholdEfficiency))

length(unique(df_mRBPs_HEK_below0p5$Majority.protein.IDs))


# Mean log2 iBAQ values HEK mRBPs above 1% and below 0.5%
mean(df_mRBPs_HEK_high$iBAQ.Input.Mean)
sd(df_mRBPs_HEK_high$iBAQ.Input.Mean)
mean(df_mRBPs_HEK_below0p5$iBAQ.Input.Mean)
sd(df_mRBPs_HEK_below0p5$iBAQ.Input.Mean)
ggplot() +
  geom_density(data = df_mRBPs_HEK_high, aes(x = iBAQ.Input.Mean),
               color = rgb(0,0,0,1), fill = rgb(0,0,0,.25)) +
  geom_density(data = df_mRBPs_HEK_below0p5, aes(x = iBAQ.Input.Mean),
               color = rgb(1,0,0,1), fill = rgb(1,0,0,.25))

t.test(df_mRBPs_HEK_high$iBAQ.Input.Mean, df_mRBPs_HEK_below0p5$iBAQ.Input.Mean)$p.value

#fold change
2**(mean(df_mRBPs_HEK_below0p5$iBAQ.Input.Mean)-mean(df_mRBPs_HEK_high$iBAQ.Input.Mean))

