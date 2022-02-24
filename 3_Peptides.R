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






#### Peptides table load and preparation ####

Peptides <- fread("../txt_RQ/peptides.txt", sep = "\t", stringsAsFactors = FALSE, data.table = F, check.names = T)
Peptides_noRQ <- fread("../txt_noRQ/peptides.txt", sep = "\t", stringsAsFactors = FALSE, data.table = F, check.names = T)
if (!Requantify) {
  Peptides <- Peptides_noRQ
  
}
rm(Requantify)


#Peptides$Gene.names <- sapply(strsplit(Peptides$Gene.names, ";"), "[", 1)
Peptides$Proteins <- sapply(strsplit(Peptides$Proteins, ";"), "[", 1)

#Peptides_noRQ$Gene.names <- sapply(strsplit(Peptides_noRQ$Gene.names, ";"), "[", 1)
Peptides_noRQ$Proteins <- sapply(strsplit(Peptides_noRQ$Proteins, ";"), "[", 1)



#### Make sure both tables have same length
Peptides <- subset(Peptides, Sequence %in% Peptides_noRQ$Sequence)
Peptides_noRQ <- subset(Peptides_noRQ, Sequence %in% Peptides$Sequence)





intensities <- names(Peptides)[grepl("Intensity.[L|M|H].", names(Peptides))]
ratios <- names(Peptides)[grepl(".normalized.", names(Peptides))]

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
  Peptides[paste0("Requantified.", i)] <- ifelse(is.na(Peptides_noRQ[i]) &
                                                   !is.na(Peptides[i]), T, F)
  
  requantified.ratios <- c(requantified.ratios, paste0("Requantified.", i))
  unscrupulous.ratios <- c(unscrupulous.ratios, paste0("Unscrupulous.", i))
  rm(i)
}



### intensities
requantified.intensities <- c()
for (i in intensities) {
  Peptides[paste0("Requantified.", i)] <- ifelse(Peptides_noRQ[i] == 0 &
                                                   Peptides[i] != 0, T, F)
  
  requantified.intensities <- c(requantified.intensities, paste0("Requantified.", i))
  rm(i)
}









#### Unscrupulous Requantification ####
# Select Unscrupulous Requantified ratios (ratios between two channels when both are requantified)
#Forward
Peptides$Unscrupulous.Ratio.H.L.Forward_II <- ifelse((Peptides$Requantified.Intensity.H.Forward_II & Peptides$Requantified.Intensity.L.Forward_II) & !is.na(Peptides$Ratio.H.L.Forward_II), T, F)

Peptides$Unscrupulous.Ratio.H.L.Forward_III <- ifelse((Peptides$Requantified.Intensity.H.Forward_III & Peptides$Requantified.Intensity.L.Forward_III) & !is.na(Peptides$Ratio.H.L.Forward_III), T, F)

Peptides$Unscrupulous.Ratio.H.L.Forward_IV <- ifelse((Peptides$Requantified.Intensity.H.Forward_IV & Peptides$Requantified.Intensity.L.Forward_IV) & !is.na(Peptides$Ratio.H.L.Forward_IV), T, F)




Peptides$Unscrupulous.Ratio.H.M.Forward_II <- ifelse((Peptides$Requantified.Intensity.H.Forward_II & Peptides$Requantified.Intensity.M.Forward_II) & !is.na(Peptides$Ratio.H.M.Forward_II), T, F)

Peptides$Unscrupulous.Ratio.H.M.Forward_III <- ifelse((Peptides$Requantified.Intensity.H.Forward_III & Peptides$Requantified.Intensity.M.Forward_III) & !is.na(Peptides$Ratio.H.M.Forward_III), T, F)

Peptides$Unscrupulous.Ratio.H.M.Forward_IV <- ifelse((Peptides$Requantified.Intensity.H.Forward_IV & Peptides$Requantified.Intensity.M.Forward_IV) & !is.na(Peptides$Ratio.H.M.Forward_IV), T, F)



Peptides$Unscrupulous.Ratio.M.L.Forward_II <- ifelse((Peptides$Requantified.Intensity.M.Forward_II & Peptides$Requantified.Intensity.L.Forward_II) & !is.na(Peptides$Ratio.M.L.Forward_II), T, F)

Peptides$Unscrupulous.Ratio.M.L.Forward_III <- ifelse((Peptides$Requantified.Intensity.M.Forward_III & Peptides$Requantified.Intensity.L.Forward_III) & !is.na(Peptides$Ratio.M.L.Forward_III), T, F)

Peptides$Unscrupulous.Ratio.M.L.Forward_IV <- ifelse((Peptides$Requantified.Intensity.M.Forward_IV & Peptides$Requantified.Intensity.L.Forward_IV) & !is.na(Peptides$Ratio.M.L.Forward_IV), T, F)







#Reverse
Peptides$Unscrupulous.Ratio.H.L.Reverse_II <- ifelse((Peptides$Requantified.Intensity.H.Reverse_II & Peptides$Requantified.Intensity.L.Reverse_II) & !is.na(Peptides$Ratio.H.L.Reverse_II), T, F)

Peptides$Unscrupulous.Ratio.H.L.Reverse_III <- ifelse((Peptides$Requantified.Intensity.H.Reverse_III & Peptides$Requantified.Intensity.L.Reverse_III) & !is.na(Peptides$Ratio.H.L.Reverse_III), T, F)

Peptides$Unscrupulous.Ratio.H.L.Reverse_IV <- ifelse((Peptides$Requantified.Intensity.H.Reverse_IV & Peptides$Requantified.Intensity.L.Reverse_IV) & !is.na(Peptides$Ratio.H.L.Reverse_IV), T, F)



Peptides$Unscrupulous.Ratio.H.M.Reverse_II <- ifelse((Peptides$Requantified.Intensity.H.Reverse_II & Peptides$Requantified.Intensity.M.Reverse_II) & !is.na(Peptides$Ratio.H.M.Reverse_II), T, F)

Peptides$Unscrupulous.Ratio.H.M.Reverse_III <- ifelse((Peptides$Requantified.Intensity.H.Reverse_III & Peptides$Requantified.Intensity.M.Reverse_III) & !is.na(Peptides$Ratio.H.M.Reverse_III), T, F)

Peptides$Unscrupulous.Ratio.H.M.Reverse_IV <- ifelse((Peptides$Requantified.Intensity.H.Reverse_IV & Peptides$Requantified.Intensity.M.Reverse_IV) & !is.na(Peptides$Ratio.H.M.Reverse_IV), T, F)



Peptides$Unscrupulous.Ratio.M.L.Reverse_II <- ifelse((Peptides$Requantified.Intensity.M.Reverse_II & Peptides$Requantified.Intensity.L.Reverse_II) & !is.na(Peptides$Ratio.M.L.Reverse_II), T, F)

Peptides$Unscrupulous.Ratio.M.L.Reverse_III <- ifelse((Peptides$Requantified.Intensity.M.Reverse_III & Peptides$Requantified.Intensity.L.Reverse_III) & !is.na(Peptides$Ratio.M.L.Reverse_III), T, F)

Peptides$Unscrupulous.Ratio.M.L.Reverse_IV <- ifelse((Peptides$Requantified.Intensity.M.Reverse_IV & Peptides$Requantified.Intensity.L.Reverse_IV) & !is.na(Peptides$Ratio.M.L.Reverse_IV), T, F)








# Delete Unscrupulous ratios
for (i in 1:length(ratios)) {
  Peptides[,ratios[i]] <- as.numeric(ifelse(Peptides[,unscrupulous.ratios[i]], NA, Peptides[,ratios[i]]))
  
  rm(i)
}



rm(Peptides_noRQ)


#### Filters ####

Peptides <- subset(Peptides, Reverse != "+")
Peptides <- subset(Peptides, Potential.contaminant != "+")



#transform intensities and ratios to Log2 (or 10 if you prefer)
Peptides[c(intensities, ratios)] = log2(Peptides[c(intensities,ratios)])
# change Inf values for na
is.na(Peptides[c(intensities, ratios)]) <- sapply(Peptides[c(intensities, ratios)], is.infinite)
is.na(Peptides[c(intensities, ratios)]) <- sapply(Peptides[c(intensities, ratios)], is.nan)



# how many proteins were identified (have intensities values not NA) in each L-H group pair
for (i in 1:length(intensities)) {
  cat(intensities[i])
  cat("\t")
  cat("\t")
  cat(nrow(Peptides[!is.na(Peptides[intensities[i]]),]))
  cat("\t")
  cat(signif(mean(Peptides[,intensities[i]], na.rm = T), 4))
  cat("\t")
  cat(sum(Peptides[,requantified.intensities[i]], na.rm = T))
  cat("\t")
  cat(signif(100*(sum(Peptides[,requantified.intensities[i]], na.rm = T)/nrow(Peptides[!is.na(Peptides[intensities[i]]),])),3))
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
  cat(nrow(Peptides[!is.na(Peptides[ratios[i]]),]))
  cat("\t")
  #ratio mean
  cat(signif(mean(Peptides[,ratios[i]], na.rm = T), 3))
  cat("\t")
  #Requantified ratios
  cat(sum(!is.na(Peptides[ratios[i]]) & Peptides[, requantified.ratios[i]], na.rm = T))
  cat("\t")
  #Percentage of requantified ratios
  cat(100*(sum(!is.na(Peptides[ratios[i]]) & Peptides[, requantified.ratios[i]], na.rm = T)) / nrow(Peptides[!is.na(Peptides[ratios[i]]),]))
  cat("\n")
  rm(i)
}


# how many proteins were identified (have intensities values not NA) in each L-H group pair
for (i in 1:length(unscrupulous.ratios)) {
  cat(unscrupulous.ratios[i])
  cat("\t")
  cat("\t")
  cat(sum(Peptides[,unscrupulous.ratios[i]]))
  cat("\t")
  cat(100*(sum(Peptides[,unscrupulous.ratios[i]])/sum(!Peptides[,unscrupulous.ratios[i]])))
  cat("\n")
  rm(i)
}



#rm(requantified.intensities, requantified.ratios, unscrupulous.ratios)







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
Peptides$Known_RBP <- ifelse(Peptides$Gene.names %in% allRBPs, T, F)
Peptides$Known_mRBP <- ifelse(Peptides$Gene.names %in% mRBPs, T, F)
Peptides$Known_mRBP_HEK <- ifelse(Peptides$Gene.names %in% mRBPs_HEK, T, F)














#### Define Peptides_summ as working table ####

Peptides_summ <- subset(Peptides,
                        select = c("id", "Protein.group.IDs", "Proteins", "Gene.names",
                                   "Sequence",
                                   "Start.position", "End.position",
                                   "Length", "PEP", "Score",
                                   "Missed.cleavages", "Known_RBP", "Known_mRBP",
                                   "Unique..Groups.", "Unique..Proteins.",
                                   intensities, ratios,
                                   requantified.intensities, requantified.ratios,
                                   unscrupulous.ratios))






Peptides_summ$Unique..Groups. <- ifelse(Peptides_summ$Unique..Groups. == "yes", T, F)
Peptides_summ$Unique..Proteins. <- ifelse(Peptides_summ$Unique..Proteins. == "yes", T, F)










#### Invert label swap ratios ####
Peptides_summ$Ratio.H.L.Reverse_II <- -Peptides_summ$Ratio.H.L.Reverse_II
Peptides_summ$Ratio.H.L.Reverse_III <- -Peptides_summ$Ratio.H.L.Reverse_III
Peptides_summ$Ratio.H.L.Reverse_IV <- -Peptides_summ$Ratio.H.L.Reverse_IV

Peptides_summ$Ratio.M.L.Reverse_II <- -Peptides_summ$Ratio.M.L.Reverse_II

Peptides_summ$Ratio.H.M.Reverse_II <- -Peptides_summ$Ratio.H.M.Reverse_II











#### Correct ratios to input fractions ####

# I originally tried to correct for losses in the Ethanol precipitation of the input samples (50% loss) by multiplying all ratios by 1/2 (sum log2(0.5)). I decided to not do it to reduce data manipulation.

# Replicate II
Peptides_summ$Ratio.H.L.Forward_II <- (Peptides_summ$Ratio.H.L.Forward_II)
Peptides_summ$Ratio.M.L.Forward_II <- (Peptides_summ$Ratio.M.L.Forward_II) + log2(2) #I used half of Pulldown for no-crosslink samples

Peptides_summ$Ratio.H.L.Reverse_II <- (Peptides_summ$Ratio.H.L.Reverse_II)
Peptides_summ$Ratio.H.M.Reverse_II <- (Peptides_summ$Ratio.H.M.Reverse_II) + log2(2) #I used half of Pulldown for no-crosslink samples


# Replicate III
Peptides_summ$Ratio.H.L.Forward_III <- (Peptides_summ$Ratio.H.L.Forward_III)

Peptides_summ$Ratio.H.L.Reverse_III <- (Peptides_summ$Ratio.H.L.Reverse_III)



# Replicate IV
Peptides_summ$Ratio.H.L.Forward_IV <- (Peptides_summ$Ratio.H.L.Forward_IV)

Peptides_summ$Ratio.H.L.Reverse_IV <- (Peptides_summ$Ratio.H.L.Reverse_IV)





#### Define RBPs from NoCL ####
Peptides_summ$qRIC_mRBP <- ifelse((Peptides_summ$Ratio.H.M.Forward_II >= 1 &
                                     Peptides_summ$Ratio.M.L.Reverse_II >= 1) &
                                    (!is.na(Peptides_summ$Ratio.H.M.Forward_II) |
                                       !is.na(Peptides_summ$Ratio.M.L.Reverse_II)), T, F)




#### Ratios into Percentages ####
# Replicate II
Peptides_summ$Percentage.H.L.Forward_II <- (2^Peptides_summ$Ratio.H.L.Forward_II)
Peptides_summ$Percentage.M.L.Forward_II <- (2^Peptides_summ$Ratio.M.L.Forward_II)

Peptides_summ$Percentage.H.L.Reverse_II <- (2^Peptides_summ$Ratio.H.L.Reverse_II)
Peptides_summ$Percentage.H.M.Reverse_II <- (2^Peptides_summ$Ratio.H.M.Reverse_II)


# Replicate III
Peptides_summ$Percentage.H.L.Forward_III <- (2^Peptides_summ$Ratio.H.L.Forward_III)

Peptides_summ$Percentage.H.L.Reverse_III <- (2^Peptides_summ$Ratio.H.L.Reverse_III)



# Replicate IV
Peptides_summ$Percentage.H.L.Forward_IV <- (2^Peptides_summ$Ratio.H.L.Forward_IV)

Peptides_summ$Percentage.H.L.Reverse_IV <- (2^Peptides_summ$Ratio.H.L.Reverse_IV)


percentages <- c("Percentage.H.L.Forward_II",
                 "Percentage.H.L.Forward_III",
                 "Percentage.H.L.Forward_IV",
                 "Percentage.H.L.Reverse_II",
                 "Percentage.H.L.Reverse_III",
                 "Percentage.H.L.Reverse_IV",
                 "Percentage.M.L.Forward_II",
                 "Percentage.H.M.Reverse_II")







#### Reproducibility ####
Peptides_summ$Reproducible.logical_II <- ifelse((Peptides_summ$Ratio.H.L.Forward_II - Peptides_summ$Ratio.H.L.Reverse_II < ReproducibilityDifference) & (Peptides_summ$Ratio.H.L.Reverse_II - Peptides_summ$Ratio.H.L.Forward_II < ReproducibilityDifference), T, F)

Peptides_summ$Reproducible.logical_III <- ifelse((Peptides_summ$Ratio.H.L.Forward_III - Peptides_summ$Ratio.H.L.Reverse_III < ReproducibilityDifference) & (Peptides_summ$Ratio.H.L.Reverse_III - Peptides_summ$Ratio.H.L.Forward_III < ReproducibilityDifference), T, F)

Peptides_summ$Reproducible.logical_IV <- ifelse((Peptides_summ$Ratio.H.L.Forward_IV - Peptides_summ$Ratio.H.L.Reverse_IV < ReproducibilityDifference) & (Peptides_summ$Ratio.H.L.Reverse_IV - Peptides_summ$Ratio.H.L.Forward_IV < ReproducibilityDifference), T, F)






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
Peptides_summ$Ratio.H.L.Forward <- MyMeanWithMinValues(Peptides_summ, c("Ratio.H.L.Forward_II",
                                                              "Ratio.H.L.Forward_III",
                                                              "Ratio.H.L.Forward_IV"),
                                                  MinNumberReplicates = 1)

Peptides_summ$Ratio.H.L.Reverse <- MyMeanWithMinValues(Peptides_summ, c("Ratio.H.L.Reverse_II",
                                                              "Ratio.H.L.Reverse_III",
                                                              "Ratio.H.L.Reverse_IV"),
                                                  MinNumberReplicates = 1)

Peptides_summ$Ratio.H.L.Mean <- MyMeanWithMinValues(Peptides_summ, c("Ratio.H.L.Forward",
                                                           "Ratio.H.L.Reverse"),
                                               MinNumberReplicates = 2)



Peptides_summ$Ratio.NoCl.Mean <- MyMeanWithMinValues(Peptides_summ, c("Ratio.M.L.Forward_II",
                                                            "Ratio.H.M.Reverse_II"),
                                                MinNumberReplicates = 2)


Peptides_summ$Ratio.ClXNoCl.Mean <- MyMeanWithMinValues(Peptides_summ, c("Ratio.H.M.Forward_II",
                                                               "Ratio.M.L.Reverse_II"),
                                                   MinNumberReplicates = 2)







#Intensities
Peptides_summ$Intensity.Input.Mean <- MyMeanWithMinValues(Peptides_summ, c("Intensity.L.Forward_II",
                                                                           "Intensity.L.Forward_III",
                                                                           "Intensity.L.Forward_IV",
                                                                           "Intensity.H.Reverse_II",
                                                                           "Intensity.H.Reverse_III",
                                                                           "Intensity.H.Reverse_IV"),
                                                          MinNumberReplicates = 1)






# Percentages
Peptides_summ$Percentage.H.L.Forward <- 2**Peptides_summ$Ratio.H.L.Forward

Peptides_summ$Percentage.H.L.Reverse <- 2**Peptides_summ$Ratio.H.L.Reverse

Peptides_summ$Percentage_Mean <- 2**Peptides_summ$Ratio.H.L.Mean

Peptides_summ$Percentage_Mean_NoCl <- 2**Peptides_summ$Ratio.NoCl.Mean




#### Write table ####
write.table(Peptides_summ, file = "Peptides_summary.txt", sep = "\t", na = "", quote = F, row.names = F)

#### boxplot intensities ####
MyBoxplot <- function(df, variables,
                      AxisName_y = ("Log2"), LabelNames = variables,
                      limits_y_min = round(min(df[variables], na.rm = T))-1,
                      limits_y_max = round(max(df[variables], na.rm = T))+1,
                      limits_breaks = round((limits_y_max+limits_y_min)/20)) {
  
  all <- melt(df,
              id.vars = c("Sequence"), measure.vars=variables)
  
  
  
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




boxplot <- MyBoxplot(Peptides_summ, intensities, AxisName_y = "Intensities (Log2)",
                     LabelNames = intensities)






png(paste("PEP_Boxplot.png", sep = ""), width = 1000, height = 1000, pointsize = 25)
grid.arrange(boxplot, nrow = 1)
dev.off()
rm(boxplot)









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
  "RBM20",# "SERBP1", 
  "SF3B1", "UPF1", "ELAVL1")
#selected_protein <- unique(subset(Peptides_summ, grepl("EIF", Peptides_summ$Gene.names))$Proteins)
#selected_protein <- unique(subset(Peptides_summ, Intensity.Input.Mean > 31 & Ratio.H.L.Mean > 1)$Proteins)
#selected_protein <- c("blank")

foo <- subset(Peptides_summ, Gene.names %in% selected_protein)


ratios_CL <- MyScatterPlot(Peptides_summ, xmin = -7,
                           ymin = 20, ymax = 32,
                           "Ratio.H.L.Mean", "Intensity.Input.Mean",
                           ShowRBPs = F, InPercentage.x = T, InPercentage.y = F) +
  
  geom_point(data = foo, aes(color = Gene.names),
             na.rm = T,
             size = 3)

#png(filename = paste0("PEP_ratios_", paste(selected_protein, collapse = "_"), ".png"), height = 500, width = 500)
png(filename = paste0("PEP_ratios_", paste(selected_protein, collapse = "_"), ".png"), height = 500, width = 500)
grid.arrange(ratios_CL, 
             ncol = 1)
dev.off()
rm(foo, ratios_CL)






#### Figure paper Histogram Ratios mRBPs ####
MyHistPlotPaper <- function(df = Peptides_summ,
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

ratios_All_proteins <- MyHistPlotPaper(Peptides_summ, Title = "All proteins")
ratios_mRBPs <- MyHistPlotPaper(subset(Peptides_summ, Known_mRBP), Title = "mRBPs only")
ratios_nonmRBPs <- MyHistPlotPaper(subset(Peptides_summ, !Known_mRBP), Title = "Non mRBPs")

pdf(file = "PEP_Hist_mRBPs.pdf", height = 5, width = 15, useDingbats = F)
grid.arrange(ratios_All_proteins, ratios_mRBPs, ratios_nonmRBPs,
             nrow = 1)
dev.off()
rm(ratios_All_proteins, ratios_mRBPs, ratios_nonmRBPs)
