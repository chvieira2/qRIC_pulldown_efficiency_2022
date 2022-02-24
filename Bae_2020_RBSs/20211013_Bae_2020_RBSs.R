setwd("K:/Results/20210903_Phospho_qRIC_revisionPaper/Bae_2020_RBSs")

library(ggplot2)
library(reshape2)
library(RColorBrewer)
library(gridExtra)
library(dplyr)
library(limma)
library(ggrepel)
library(gplots)
library(data.table)

DeltaAnalysisThresholdEfficiency <- -1
LOG2_logical = F

#### Load Bae2020 files ####
RBSs <- fread("K:/Datasets/Bae_2020_NatStrucMolBio_RBS-ID_supp2_RBSs.txt", sep = "\t", stringsAsFactors = FALSE, data.table = F, check.names = T)
RBSs$Type <- as.factor("RBS")


RBSs_candidate <- fread("K:/Datasets/Bae_2020_NatStrucMolBio_RBS-ID_supp2_candidateRBSs.txt", sep = "\t", stringsAsFactors = FALSE, data.table = F, check.names = T)
# RBSs candidate has many NA rows that must be removed
RBSs_candidate <- RBSs_candidate[RBSs_candidate$RBS != "", ]

RBSs_candidate$Type <- as.factor("Candidate")

#### Add both qualities of data. In paper, candidates are not considered
#RBSs <- rbind(RBSs, RBSs_candidate) # I decided to not include it as it changes very little the final number of sites in my candidates that match to a protein with x-link site
rm(RBSs_candidate)




## Turn localization probabilities into numeric and only keep the highest
RBSs$Localizationscore_Total.RNA.RBS.Rep1 <- as.numeric(RBSs$Localizationscore_Total.RNA.RBS.Rep1)
RBSs$Localizationscore_Total.RNA.RBS.Rep2 <- as.numeric(RBSs$Localizationscore_Total.RNA.RBS.Rep2)
RBSs$Localizationscore_mRNA.RBS.Rep1 <- as.numeric(RBSs$Localizationscore_mRNA.RBS.Rep1)
RBSs$Localizationscore_mRNA.RBS.Rep2 <- as.numeric(RBSs$Localizationscore_mRNA.RBS.Rep2)

# Select maximal localization score
RBSs$Max.Localization.Prob <- as.numeric(apply(RBSs, 1, function(x) max(x[c("Localizationscore_Total.RNA.RBS.Rep1",
                                                             "Localizationscore_Total.RNA.RBS.Rep2",
                                                             "Localizationscore_mRNA.RBS.Rep1",
                                                             "Localizationscore_mRNA.RBS.Rep2")], na.rm = T)))


## Expand rows from RBSs ID separated by ";"
dt <- data.table(RBSs)
dt <- dt[ , list( ID = unlist( strsplit( RBS , ";" ) ) ) , by = eval(names(RBSs)[1:length(names(RBSs))]) ]

RBSs <- as.data.frame(dt)
rm(dt)


## ID into UNIPROT, gene name, Site.ID, amino acid position and amino acid type
#UNIPROT
RBSs$UNIPROT.ID <- sapply(RBSs$ID, function(x) unlist(strsplit(as.character(x),'\\|'))[2])
#Gene name
RBSs$Gene.name <- sapply(RBSs$ID, function(x) sub("^(?:.*\\|){2}(.*?)_.*", "\\1", x))
#Site.ID
RBSs$Site.ID <- sapply(RBSs$ID, function(x) sub("^(?:.*\\.)", "\\1", x))
#Amino acid type
RBSs$AA.type <- sapply(RBSs$Site.ID, function(x) unlist(strsplit(x, "[0-9]"))[1])
#Amino acid position
RBSs$AA.position <- sapply(RBSs$Site.ID, function(x) as.numeric(unlist(strsplit(x, "[A-Z]"))[2]))


#Adjust to pattern in my other tables
RBSs$Site.ID <- ifelse(!is.na(RBSs$Gene.name),
                       paste(RBSs[,"Gene.name"], RBSs[,"Site.ID"], sep = "_"),
                       paste(RBSs[,"UNIPROT.ID"], RBSs[,"Site.ID"], sep = "_"))



## Transform intensity into log2
if (LOG2_logical) {
  RBSs[,c("Intensity_Total.RNA.RBS.Rep1",
          "Intensity_Total.RNA.RBS.Rep2",
          "Intensity_mRNA.RBS.Rep1",
          "Intensity_mRNA.RBS.Rep2")] <- log2(RBSs[,c("Intensity_Total.RNA.RBS.Rep1",
                                                      "Intensity_Total.RNA.RBS.Rep2",
                                                      "Intensity_mRNA.RBS.Rep1",
                                                      "Intensity_mRNA.RBS.Rep2")])
  is.na(RBSs[c("Intensity_Total.RNA.RBS.Rep1",
               "Intensity_Total.RNA.RBS.Rep2",
               "Intensity_mRNA.RBS.Rep1",
               "Intensity_mRNA.RBS.Rep2")]) <- sapply(RBSs[c("Intensity_Total.RNA.RBS.Rep1",
                                                             "Intensity_Total.RNA.RBS.Rep2",
                                                             "Intensity_mRNA.RBS.Rep1",
                                                             "Intensity_mRNA.RBS.Rep2")], is.infinite)
}



## Select only columns of relevance
RBSs <- subset(RBSs, select = c("UNIPROT.ID", "Gene.name", "Site.ID", "AA.type", "AA.position",
                                "Max.Localization.Prob",
                                "Type",
                                "Intensity_Total.RNA.RBS.Rep1",
                                "Intensity_Total.RNA.RBS.Rep2",
                                "Intensity_mRNA.RBS.Rep1",
                                "Intensity_mRNA.RBS.Rep2",
                                "RNA.binding",
                                "Domain",
                                "RBSPTM_STYphospho",
                                "RBSPTM_Kacet",
                                "RBSPTM_Rmeth"))



write.table(RBSs, file = "RBSs.txt", sep = "\t", quote = F, row.names = F)



#### Min distance function ####
MinDist <- function(df = working_table_PTM, Site.ID) {
  if(Site.ID %in% df$Site.ID_PTM){
    # Filter candidate table for only that pSite
    x <- df[df$Site.ID_PTM == Site.ID,]
    
    # Filter RBSs for crosslinking sites in that protein
    if (x$Protein[1] %in% RBSs$UNIPROT.ID) {
      RBSs_filt <- RBSs[RBSs$UNIPROT.ID == x$Protein,]
      
      return(min(abs(x$Position_PTM - RBSs_filt$AA.position)))
      
      
    } else{return(NA)}
    
    
  } else{return(NA)}
  
}



#### qRIC files ####
working_table_PTM <- read.csv("../output/working_table_PTM.txt", sep = "\t", stringsAsFactors = FALSE)
working_table_PTM <- subset(working_table_PTM, !is.na(delta_Mean) &
                              Known_mRBP_PG)


working_table_PTM$Min.Dist.Xsite <- apply(working_table_PTM["Site.ID_PTM"], 1, function(x) MinDist(df = working_table_PTM, Site.ID = x))

working_table_PTM <- subset(working_table_PTM, !is.na(Min.Dist.Xsite))


if(LOG2_logical) {
  working_table_PTM$Sequence.length_PG <- log2(working_table_PTM$Sequence.length_PG)
  working_table_PTM$Min.Dist.Xsite <- log2(working_table_PTM$Min.Dist.Xsite)
}


## Create groups for plotting
working_table_PTM$Group <- factor(ifelse(working_table_PTM$PTM_increase_Efficiency_Cl_Mean, "Increase",
                                         ifelse(working_table_PTM$PTM_decrease_Efficiency_Cl_Mean, "Decrease",
                                                ifelse(!working_table_PTM$PTM_increase_Efficiency_Cl_Mean & 
                                                         !working_table_PTM$PTM_decrease_Efficiency_Cl_Mean, "Neutral",
                                                       NA))), levels = c("Decrease", "Neutral", "Increase"))





### Create filtered tables according to protein pull-down efficiency
working_table_PTM_high <- subset(working_table_PTM, Ratio.H.L.Forward_PG >= DeltaAnalysisThresholdEfficiency & 
                                                            Ratio.H.L.Reverse_PG >= DeltaAnalysisThresholdEfficiency)


working_table_PTM_low <- subset(working_table_PTM, Ratio.H.L.Forward_PG <= DeltaAnalysisThresholdEfficiency & 
                                                            Ratio.H.L.Reverse_PG <= DeltaAnalysisThresholdEfficiency)



### P.values ####
if (LOG2_logical) {
  ### Convert Inf values to NA for p.value calculation
  # All proteins
  working_table_PTM$Min.Dist.Xsite_noInf <- working_table_PTM$Min.Dist.Xsite
  is.na(working_table_PTM["Min.Dist.Xsite_noInf"]) <- sapply(working_table_PTM["Min.Dist.Xsite_noInf"], is.infinite)
  
  res.aov <- aov(Min.Dist.Xsite_noInf ~ Group, data = working_table_PTM)
  Tukey_test <- data.frame(TukeyHSD(res.aov, conf.level=0.95)$Group)
  
  
  # high efficiency proteins
  working_table_PTM_high$Min.Dist.Xsite_noInf <- working_table_PTM_high$Min.Dist.Xsite
  is.na(working_table_PTM_high["Min.Dist.Xsite_noInf"]) <- sapply(working_table_PTM_high["Min.Dist.Xsite_noInf"], is.infinite)
  
  res.aov <- aov(Min.Dist.Xsite_noInf ~ Group, data = working_table_PTM_high)
  Tukey_test_high <- data.frame(TukeyHSD(res.aov, conf.level=0.95)$Group)
  
  
  # low efficiency proteins
  working_table_PTM_low$Min.Dist.Xsite_noInf <- working_table_PTM_low$Min.Dist.Xsite
  is.na(working_table_PTM_low["Min.Dist.Xsite_noInf"]) <- sapply(working_table_PTM_low["Min.Dist.Xsite_noInf"], is.infinite)
  
  res.aov <- aov(Min.Dist.Xsite_noInf ~ Group, data = working_table_PTM_low)
  Tukey_test_low <- data.frame(TukeyHSD(res.aov, conf.level=0.95)$Group)
  
  
  
  
  ### Protein length
  # All proteins
  working_table_PTM_unique <- unique(working_table_PTM[,c("Protein_PTM", "Sequence.length_PG", "Group")])
  res.aov_length <- aov(Sequence.length_PG ~ Group, data = working_table_PTM_unique)
  Tukey_test_length <- data.frame(TukeyHSD(res.aov_length, conf.level=0.95)$Group)
  
  # high efficiency proteins
  working_table_PTM_high_unique <- unique(working_table_PTM_high[,c("Protein_PTM", "Sequence.length_PG", "Group")])
  res.aov_length <- aov(Sequence.length_PG ~ Group, data = working_table_PTM_high_unique)
  Tukey_test_length_high <- data.frame(TukeyHSD(res.aov_length, conf.level=0.95)$Group)
  
  # low efficiency proteins
  working_table_PTM_low_unique <- unique(working_table_PTM_low[,c("Protein_PTM", "Sequence.length_PG", "Group")])
  res.aov_length <- aov(Sequence.length_PG ~ Group, data = working_table_PTM_low_unique)
  Tukey_test_length_low <- data.frame(TukeyHSD(res.aov_length, conf.level=0.95)$Group)
} else {
  
  wiltest_Dist <- setNames(data.frame(matrix(ncol = 2, nrow = 0), stringsAsFactors = F), c("Comparison", "p.value"))
  
  wiltest_Dist <- rbind(wiltest_Dist,
                        list("DecreasexNeutral", 
                          signif(wilcox.test(subset(working_table_PTM, Group == "Decrease")[,"Min.Dist.Xsite"],
                                             subset(working_table_PTM, Group == "Neutral")[,"Min.Dist.Xsite"])$p.value, digits = 2)),
                        
                        list("DecreasexIncrease", 
                          signif(wilcox.test(subset(working_table_PTM, Group == "Decrease")[,"Min.Dist.Xsite"],
                                             subset(working_table_PTM, Group == "Increase")[,"Min.Dist.Xsite"])$p.value, digits = 2)),
                        
                        list("NeutralxIncrease", 
                          signif(wilcox.test(subset(working_table_PTM, Group == "Neutral")[,"Min.Dist.Xsite"],
                                             subset(working_table_PTM, Group == "Increase")[,"Min.Dist.Xsite"])$p.value, digits = 2)))
  
  
  
  
  wiltest_Length <- setNames(data.frame(matrix(ncol = 2, nrow = 0), stringsAsFactors = F), c("Comparison", "p.value"))
  
  wiltest_Length <- rbind(wiltest_Length,
                          list("DecreasexNeutral", 
                            signif(wilcox.test(subset(working_table_PTM, Group == "Decrease")[,"Sequence.length_PG"],
                                               subset(working_table_PTM, Group == "Neutral")[,"Sequence.length_PG"])$p.value, digits = 2)),
                          
                          list("DecreasexIncrease", 
                            signif(wilcox.test(subset(working_table_PTM, Group == "Decrease")[,"Sequence.length_PG"],
                                               subset(working_table_PTM, Group == "Increase")[,"Sequence.length_PG"])$p.value, digits = 2)),
                          
                          list("NeutralxIncrease", 
                            signif(wilcox.test(subset(working_table_PTM, Group == "Neutral")[,"Sequence.length_PG"],
                                               subset(working_table_PTM, Group == "Increase")[,"Sequence.length_PG"])$p.value, digits = 2)))
  
  
  ## Unique table just for plotting later protein length
  working_table_PTM_unique <- unique(working_table_PTM[,c("Protein_PTM", "Sequence.length_PG", "Group")])
  working_table_PTM_high_unique <- unique(working_table_PTM_high[,c("Protein_PTM", "Sequence.length_PG", "Group")])
  working_table_PTM_low_unique <- unique(working_table_PTM_low[,c("Protein_PTM", "Sequence.length_PG", "Group")])
  
  
  # Adjust col names
  colnames(wiltest_Dist) <- c("Comparison", "p.value")
  colnames(wiltest_Length) <- c("Comparison", "p.value")
}



#### Boxplot distances ####
MyBoxplot <- function(df = working_table_PTM, Test = NULL,
                      Y,
                      AxisName_y = Y, Factor = "Group",
                      ymin = 0,
                      ymax = 1500,
                      AxisSep.y = round((abs(ymax)+abs(ymin))/4)) {
  
  df <- melt(df,
              id.vars = c(Factor), measure.vars=Y)
  
  df <- subset(df, !is.na(eval(parse(text = Factor))) & !is.na(value))
  
  bplot <- ggplot(df, aes(eval(parse(text = Factor)), value)) +
    
    geom_jitter(position=position_jitter(0.2), na.rm = T, size = 3) +
    
    geom_boxplot(fill = rgb(0,0,0,0), width = 0.75, na.rm = T, notch = F, outlier.shape = NA)  +
    
    xlab(NULL) +
    ylab(AxisName_y) +
    coord_cartesian(ylim = c(ymin, ymax)) +
    scale_y_continuous(breaks = seq(ymin, ymax, AxisSep.y)) +
    
    theme_bw() +
    theme(plot.title = element_text(size = 35,
                                    hjust = 0.5,
                                    vjust = 0.4),
          axis.title.x = element_text(size=30,
                                      hjust = 0.5,
                                      vjust = 0.4),
          axis.text.x  = element_text(color = "black",
                                      angle=45, 
                                      #vjust=0,
                                      hjust = 1,
                                      size=25),
          axis.title.y = element_text(size=30,
                                      hjust = 0.5,
                                      vjust = 1.5),
          axis.text.y  = element_text(color = "black",
                                      angle=0, 
                                      vjust=0.5,
                                      hjust = 0.5,
                                      size=25),
          panel.grid=element_blank(),
          aspect.ratio=1)
  
  if (!is.null(Test)) {
    
    bplot <- bplot + annotate("text", x = 1, y = ymax*0.8,
                              label = paste0("N=", table(df$Group)[1])) +

      annotate("text", x = 2, y = ymax*0.8,
               label = paste0("N=", table(df$Group)[2])) +
      
      annotate("text", x = 3, y = ymax*0.8,
               label = paste0("N=", table(df$Group)[3]))
    
    
    if (LOG2_logical) {
      bplot <- bplot + annotate("text", x = 1.5, y = ymax*0.9,
                 label = paste0(row.names(Test)[1], "\n", "p-value=",  signif(Test$p.adj,2)[1])) +
        
        annotate("text", x = 2, y = ymax,
                 label = paste0(row.names(Test)[2], "\n", "p-value=",  signif(Test$p.adj,2)[2])) +
        
        annotate("text", x = 2.5, y = ymax*0.9,
                 label = paste0(row.names(Test)[3], "\n", "p-value=",  signif(Test$p.adj,2)[3]))
    } else {
      
      
      
      bplot <- bplot +
        annotate("text", x = 1.5, y = ymax*0.9,
                 label = paste0(Test[1,1], "\n", "p-value=",  Test[1,2])) +
        
        annotate("text", x = 2, y = ymax,
                 label = paste0(Test[2,1], "\n", "p-value=",  Test[2,2])) +
        
        annotate("text", x = 2.5, y = ymax*0.9,
                 label = paste0(Test[3,1], "\n", "p-value=",  Test[3,2]))
      
    }
  }
    
    
  
  return(bplot)
  
}





boxplot <- MyBoxplot(df = working_table_PTM, 
                     Test = if(LOG2_logical) Tukey_test else wiltest_Dist,
                     Y = "Min.Dist.Xsite",
                     AxisName_y = "Linear AAs distance\nto nearest RBS")
boxplot_high <- MyBoxplot(df = working_table_PTM_high, 
                     Y = "Min.Dist.Xsite", AxisName_y = "Linear AAs distance\nto nearest RBS (log2)")
boxplot_low <- MyBoxplot(df = working_table_PTM_low, 
                     Y = "Min.Dist.Xsite", AxisName_y = "Linear AAs distance\nto nearest RBS (log2)")





boxplot_length <- MyBoxplot(df = working_table_PTM_unique, 
                            Test = if(LOG2_logical) Tukey_test else wiltest_Length,
                            ymax = 2000,
                            Y = "Sequence.length_PG",
                            AxisName_y = "Protein AAs length")
boxplot_length_high <- MyBoxplot(df = working_table_PTM_high_unique, 
                                 Y = "Sequence.length_PG", AxisName_y = "Protein AAs length (log2)")
boxplot_length_low <- MyBoxplot(df = working_table_PTM_low_unique,
                                Y = "Sequence.length_PG", AxisName_y = "Protein AAs length (log2)")






png(paste("Boxplot_DistRBS.png", sep = ""), width = 1000, height = 500, pointsize = 25)
grid.arrange(boxplot, boxplot_length, 
             #boxplot_high, boxplot_length_high, 
             #boxplot_low, boxplot_length_low, 
             nrow = 1)
dev.off()


pdf(paste("Boxplot_DistRBS.pdf", sep = ""), width = 10, height = 15, useDingbats = F)
grid.arrange(boxplot, boxplot_length, 
             boxplot_high, boxplot_length_high, 
             boxplot_low, boxplot_length_low,
             nrow = 3)
dev.off()
rm(boxplot, boxplot_length, 
   boxplot_high, boxplot_length_high, 
   boxplot_low, boxplot_length_low)

rm(res.aov, res.aov_length,
   Tukey_test, Tukey_test_high, Tukey_test_length_high, Tukey_test_length_low, Tukey_test_low, Tukey_test_length)




#### Scatter Plot ####
MyScatterPlot <- function(df = working_table_PTM_high, X, Y,
                          xmin = round(min(df[,X], na.rm = T)),
                          xmax = round(max(df[,X], na.rm = T)),
                          AxisSep.x = round((abs(xmax)+abs(xmin))/20),
                          ymin = 0,
                          ymax = round(max(df[,Y], na.rm = T)),
                          AxisSep.y = round((abs(ymax)+abs(ymin))/4),
                          DotSize = 4,
                          Title = NULL,
                          xTitle = "Delta efficiency (log2)",
                          yTitle = "Linear AAs distance\nto nearest RBS (log2)",
                          bin.size = 128,
                          GradientColorMin = rgb(0,0,0), GradientColorMax = rgb(0,0,0)) {
  if(AxisSep.x <= 0) AxisSep.x <- 1
  if(AxisSep.y <= 0) AxisSep.y <- 1
  
  plot <- ggplot(data = df, aes(x = eval(parse(text = X)), y = eval(parse(text = Y)))) +
    
    ggtitle(Title) +
    xlab(xTitle) +
    ylab(yTitle) +
    coord_cartesian(xlim = c(xmin,xmax),ylim = c(ymin,ymax)) +
    scale_x_continuous(breaks=seq(xmin,xmax, AxisSep.x)) +
    scale_y_continuous(breaks=seq(ymin,ymax, AxisSep.y)) +
    theme_bw() +
    theme(plot.title = element_text(size = 35,
                                    hjust = 0.5,
                                    vjust = 0.4),
          axis.title.x = element_text(size=30,
                                      hjust = 0.5,
                                      vjust = 0.4),
          axis.text.x  = element_text(color = "black",
                                      angle=0, 
                                      vjust=0.5,
                                      hjust = 0.5,
                                      size=25),
          axis.title.y = element_text(size=30,
                                      hjust = 0.5,
                                      vjust = 1.5),
          axis.text.y  = element_text(color = "black",
                                      angle=0, 
                                      vjust=0.5,
                                      hjust = 0.5,
                                      size=25),
          panel.grid=element_blank(),
          aspect.ratio=1)
  
  foo1 <- subset(df, !PTM_increase_Efficiency_Cl_Mean & !PTM_decrease_Efficiency_Cl_Mean)
  plot <- plot +  
    geom_point(data = foo1,
               na.rm = T,
               size = DotSize,
               colour = densCols(x = foo1[,X],
                                 y = foo1[,Y],
                                 colramp = colorRampPalette(c(GradientColorMin, GradientColorMax)), nbin = bin.size)) +
    
    
    geom_point(data = subset(df, PTM_increase_Efficiency_Cl_Mean | PTM_decrease_Efficiency_Cl_Mean),
               na.rm = T,
               size = DotSize,
               colour = "blue")
  
  return(plot)
}



plot <- MyScatterPlot(df = working_table_PTM, X = "delta_Mean", Y = "Min.Dist.Xsite",
                      xmin = -3, xmax = 3,
                      ymax = 12)

plot_length <- MyScatterPlot(df = working_table_PTM, X = "delta_Mean", Y = "Sequence.length_PG",
                             xmin = -3, xmax = 3,
                             ymax = 12, yTitle = "Protein length (log2)")



plot_low <- MyScatterPlot(df = working_table_PTM_low, X = "delta_Mean", Y = "Min.Dist.Xsite",
                          xmin = -3, xmax = 3,
                          ymax = 12, Title = "Low efficiency proteins")

plot_length_low <- MyScatterPlot(df = working_table_PTM_low, X = "delta_Mean", Y = "Sequence.length_PG",
                                 xmin = -3, xmax = 3,
                                 ymax = 12, yTitle = "Protein length (log2)", Title = "Low efficiency proteins")



png("Scatter_GroupXDistRBS.png", width = 1000, height = 1000, pointsize = 25)
grid.arrange(plot, plot_length,
             plot_low, plot_length_low,
             ncol = 2)
dev.off()
rm(plot, plot_length, plot_low, plot_length_low)





#### Unquantified sites ####
PTM_summ <- fread("../output/PTM_summary.txt", sep = "\t", stringsAsFactors = FALSE, data.table = F, check.names = T)
PTM_summ <- subset(PTM_summ, Known_mRBP, select = c("Protein", "Gene.names", "Site.ID", "Position",
                                                    "Ratio.H.L.Forward",
                                                    "Ratio.H.L.Reverse",
                                                    "Ratio.H.L.Mean",
                                                    
                                                    "Intensity.L.Forward_II",
                                                    "Intensity.L.Forward_III",
                                                    "Intensity.L.Forward_IV",
                                                    "Intensity.H.Forward_II",
                                                    "Intensity.H.Forward_III",
                                                    "Intensity.H.Forward_IV",
                                                    
                                                    "Intensity.L.Reverse_II",
                                                    "Intensity.L.Reverse_III",
                                                    "Intensity.L.Reverse_IV",
                                                    "Intensity.H.Reverse_II",
                                                    "Intensity.H.Reverse_III",
                                                    "Intensity.H.Reverse_IV"))
# add sequence length info
PTM_summ <- merge(PTM_summ, working_table_PTM[,c("Protein_PTM","Sequence.length_PG")], 
                  by.x = "Protein", by.y = "Protein_PTM", all.x = T)

PTM_summ <- unique(subset(PTM_summ, !is.na(Sequence.length_PG)))


# Calculate distance to nearest cross-linking sites

MinDist2 <- function(df = PTM_summ, Site.ID) {
  #Site.ID = PTM_summ$Site.ID[73]
  if(Site.ID %in% df$Site.ID){
    # Filter candidate table for only that pSite
    x <- df[df$Site.ID == Site.ID,]
    
    # Filter RBSs for crosslinking sites in that protein
    if (x$Protein[1] %in% RBSs$UNIPROT.ID) {
      RBSs_filt <- RBSs[RBSs$UNIPROT.ID == x$Protein,]
      
      return(min(abs(x$Position - RBSs_filt$AA.position)))
      
      
    } else{return(NA)}
    
    
  } else{return(NA)}
  
}

PTM_summ$Min.Dist.Xsite <- log2(apply(PTM_summ["Site.ID"], 1, function(x) MinDist2(df = PTM_summ, Site.ID = x)))
print(c(nrow(PTM_summ), table(is.na(PTM_summ$Min.Dist.Xsite))))

PTM_summ <- subset(PTM_summ, !is.na(Min.Dist.Xsite))




### Create groups for plotting ####
PTM_summ$Group <- factor(ifelse(is.na(PTM_summ$Ratio.H.L.Mean),
                         "Not_quant",
                         
                         ifelse((PTM_summ$Ratio.H.L.Forward >= DeltaAnalysisThresholdEfficiency &
                                   PTM_summ$Ratio.H.L.Reverse >= DeltaAnalysisThresholdEfficiency),
                                "High",
                                "Low")), levels = c("Low", "High", "Not_quant"))

print(c(nrow(PTM_summ), table(PTM_summ$Group)))



### P.values ####

#Convert Inf values to NA for p.value calculation
PTM_summ$Min.Dist.Xsite_noInf <- PTM_summ$Min.Dist.Xsite
is.na(PTM_summ["Min.Dist.Xsite_noInf"]) <- sapply(PTM_summ["Min.Dist.Xsite_noInf"], is.infinite)

res.aov <- aov(Min.Dist.Xsite_noInf ~ Group, data = PTM_summ)
Tukey_test <- data.frame(TukeyHSD(res.aov, conf.level=0.95)$Group)



### Protein length
#Convert Inf values to NA for p.value calculation
PTM_summ_unique <- unique(PTM_summ[,c("Protein", "Sequence.length_PG", "Group")])
res.aov_length <- aov(Sequence.length_PG ~ Group, data = PTM_summ_unique)
Tukey_test_length <- data.frame(TukeyHSD(res.aov_length, conf.level=0.95)$Group)






## Boxplot
boxplot <- MyBoxplot(df = PTM_summ, Tukey = Tukey_test,
                     Y = "Min.Dist.Xsite", AxisName_y = "Linear AAs distance\nto nearest RBS (log2)", ymax = 15) 

boxplot_length <- MyBoxplot(df = PTM_summ_unique, Tukey = Tukey_test_length,
                            Y = "Sequence.length_PG", AxisName_y = "Protein AAs length (log2)", ymax = 15)



png(paste("Boxplot_DistRBS_pSites.png", sep = ""), width = 1000, height = 500, pointsize = 25)
grid.arrange(boxplot, boxplot_length, nrow = 1)
dev.off()
rm(boxplot, boxplot_length)












#### Input exclusive psites ####
MyExperimentLogic <- function(df = PTM_summ, Pulldown_col, Input_col) {
  
  vec <- as.factor(ifelse(is.na(df[,Pulldown_col]) & 
                  !is.na(df[,Input_col]),
                "Input only",
                
                ifelse(!is.na(df[,Pulldown_col]) & 
                         is.na(df[,Input_col]),
                       "Pulldown only",
                       
                       ifelse(is.na(df[,Pulldown_col]) & 
                                is.na(df[,Input_col]),
                              "Both",
                              
                              NA))))
  
  return(vec)
}



# Experiment II
PTM_summ$Experiment_Forward_II <- MyExperimentLogic(Pulldown_col = "Intensity.H.Forward_II",
                                                         Input_col = "Intensity.L.Forward_II")

PTM_summ$Experiment_Reverse_II <- MyExperimentLogic(Pulldown_col = "Intensity.L.Reverse_II",
                                                         Input_col = "Intensity.H.Reverse_II")

# Experiment III
PTM_summ$Experiment_Forward_III <- MyExperimentLogic(Pulldown_col = "Intensity.H.Forward_III",
                                                         Input_col = "Intensity.L.Forward_III")

PTM_summ$Experiment_Reverse_III <- MyExperimentLogic(Pulldown_col = "Intensity.L.Reverse_III",
                                                         Input_col = "Intensity.H.Reverse_III")

# Experiment IV
PTM_summ$Experiment_Forward_IV <- MyExperimentLogic(Pulldown_col = "Intensity.H.Forward_IV",
                                                         Input_col = "Intensity.L.Forward_IV")

PTM_summ$Experiment_Reverse_IV <- MyExperimentLogic(Pulldown_col = "Intensity.L.Reverse_IV",
                                                         Input_col = "Intensity.H.Reverse_IV")







# All Experiments
# Define Input only as those p-sites identified in at least only experiment as Input only and not as Both or Pulldown only in any other experiment


PTM_summ$Experiment_Forward <- apply(PTM_summ, 1,
               function(x) {
             # First select those in Both PD and input
            as.factor(ifelse("Both" %in% x[c("Experiment_Forward_II",
                                   "Experiment_Forward_III",
                                   "Experiment_Forward_IV")],
                                                                 
                   "Both",
                   
                   # Second ifelse for Input only
                   ifelse((sum(x[c("Experiment_Forward_II",
                                   "Experiment_Forward_III",
                                   "Experiment_Forward_IV")] == "Input only") >= 1) &
                            
                            (sum(x[c("Experiment_Forward_II",
                                     "Experiment_Forward_III",
                                     "Experiment_Forward_IV")] == "Pulldown only") < 1),
                          
                          "Input only",
                          
                          # Third ifelse for Pulldown only
                          ifelse((sum(x[c("Experiment_Forward_II",
                                          "Experiment_Forward_III",
                                          "Experiment_Forward_IV")] == "Pulldown only") >= 1) &
                                   
                                   (sum(x[c("Experiment_Forward_II",
                                            "Experiment_Forward_III",
                                            "Experiment_Forward_IV")] == "Input only") < 1),
                                 
                                 "Pulldown only",
                                 
                                 "Mixed"))))
                                          })

table(PTM_summ$Experiment_Forward)

PTM_summ$Experiment_Reverse <- apply(PTM_summ, 1,
                                          function(x) {
                                            # First select those in Both
                                            as.factor(ifelse("Both" %in% x[c("Experiment_Reverse_II",
                                                                   "Experiment_Reverse_III",
                                                                   "Experiment_Reverse_IV")],
                                                   
                                                   "Both",
                                                   
                                                   # Second ifelse for Input only
                                                   ifelse((sum(x[c("Experiment_Reverse_II",
                                                                   "Experiment_Reverse_III",
                                                                   "Experiment_Reverse_IV")] == "Input only") >= 1) &
                                                            
                                                            (sum(x[c("Experiment_Reverse_II",
                                                                     "Experiment_Reverse_III",
                                                                     "Experiment_Reverse_IV")] == "Pulldown only") < 1),
                                                          
                                                          "Input only",
                                                          
                                                          # Third ifelse for Pulldown only
                                                          ifelse((sum(x[c("Experiment_Reverse_II",
                                                                          "Experiment_Reverse_III",
                                                                          "Experiment_Reverse_IV")] == "Pulldown only") >= 1) &
                                                                   
                                                                   (sum(x[c("Experiment_Reverse_II",
                                                                            "Experiment_Reverse_III",
                                                                            "Experiment_Reverse_IV")] == "Input only") < 1),
                                                                 
                                                                 "Pulldown only",
                                                                 
                                                                 "Mixed"))))
                                          })

table(PTM_summ$Experiment_Reverse)




PTM_summ$Experiment <- as.factor(ifelse(PTM_summ$Experiment_Forward == "Input only" &
                                                       PTM_summ$Experiment_Reverse == "Input only",
                                                     "Input only",
                                                     
                                                     ifelse(PTM_summ$Experiment_Forward == "Pulldown only" &
                                                              PTM_summ$Experiment_Reverse == "Pulldown only",
                                                            "Pulldown only",
                                                            
                                                            ifelse(PTM_summ$Experiment_Forward == "Both" &
                                                                     PTM_summ$Experiment_Reverse == "Both",
                                                                   "Both",
                                                                   
                                                                   ifelse(!is.na(PTM_summ$Experiment_Forward) &
                                                                            !is.na(PTM_summ$Experiment_Reverse),
                                                                          "Both",NA)))))

table(PTM_summ$Experiment)




#### Boxplot Input only ####
boxplot_Forward <- MyBoxplot(PTM_summ, Y = "Min.Dist.Xsite", AxisName_y = "Linear AAs distance\nto nearest RBS (log2)", ymax = 15, Factor = "Experiment_Forward")

boxplot_Reverse <- MyBoxplot(PTM_summ, Y = "Min.Dist.Xsite", AxisName_y = "Linear AAs distance\nto nearest RBS (log2)", ymax = 15, Factor = "Experiment_Reverse")

boxplot <- MyBoxplot(PTM_summ, Y = "Min.Dist.Xsite", AxisName_y = "Linear AAs distance\nto nearest RBS (log2)", ymax = 15, Factor = "Experiment")






png(paste("Boxplot_InputOnly.png", sep = ""), width = 1500, height = 500, pointsize = 25)
grid.arrange(boxplot_Forward, boxplot_Reverse, boxplot, nrow = 1)
dev.off()
rm(boxplot)
