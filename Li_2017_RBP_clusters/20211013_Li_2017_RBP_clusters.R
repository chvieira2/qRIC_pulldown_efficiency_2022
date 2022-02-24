setwd("K:/Results/20210903_Phospho_qRIC_revisionPaper/Li_2017_RBP_clusters")

library(ggplot2)
library(reshape2)
library(RColorBrewer)
library(gridExtra)
library(dplyr)
library(limma)
library(ggrepel)
library(gplots)
library(data.table)




#### Load proteins ####
PG_summ <- read.csv("../output/PG_summary.txt", sep = "\t", stringsAsFactors = FALSE)


#PG_summ$Intensity_Input_Mean <- rowMeans(PG_summ[,c("Intensity.L.Forward", "Intensity.H.Reverse")], na.rm = F)

#PG_summ <- subset(PG_summ, (Ratio.H.L.Forward > 0 & Ratio.H.L.Reverse > 0) & Reproducible.logical)


PG_summ <- subset(PG_summ, !is.na(Ratio.H.L.Mean) & !is.na(iBAQ.Input.Mean),
                  select = c("Majority.protein.IDs","Gene.names",
                                                              "Ratio.H.L.Mean", "iBAQ.Input.Mean"))




#### Plot clusters of similar binders as defined in Li et al 2017 Genome Biol. ####
cluster_2 <- c("AGO4", "CPSF2", "ELAVL1", "LIN28A", "PUM2", "QKI", "RBPMS",
               "TNRC6A", "TNRC6B", "TNRC6C")
cluster_3 <- c("EWSR1", "FUS", "LIN28B", "TAF15")
cluster_4 <- c("HNRNPA1", "HNRNPA2B1", "HNRNPF", "HNRNPH1",
               "HNRNPM", "HNRNPU")
cluster_5 <- c("CPSF6", "NUDT21")
cluster_6 <- c("CSTF2T", "CSTF2")
cluster_7 <- c("C17orf85", "RTCB")
cluster_9 <- c("AGO1", "AGO2", "AGO3", "AGO4", "TNRC6A", "TNRC6C")
cluster_12 <- c("FBL", "NOP56", "NOP58")
cluster_13 <- c("ALKBH5", "CPSF1", "CPSF2", "CPSF3", "CPSF4")
cluster_15 <- c("CPSF6", "CPSF7")
cluster_16 <- c("EWSR1", "FXR1", "FXR2")
cluster_18 <- c("IGF2BP1", "IGF2BP2", "IGF2BP3")



PG_summ$Li_2017_Clusters <- ifelse(PG_summ$Gene.names %in% cluster_2, "2",
                                   ifelse(PG_summ$Gene.names %in% cluster_3, "3",
                                          ifelse(PG_summ$Gene.names %in% cluster_4, "4",
                                                 ifelse(PG_summ$Gene.names %in% cluster_5, "5",
                                                        ifelse(PG_summ$Gene.names %in% cluster_6, "6",
                                                               ifelse(PG_summ$Gene.names %in% cluster_7, "7",
                                                                      ifelse(PG_summ$Gene.names %in% cluster_9, "9",
                                                                             ifelse(PG_summ$Gene.names %in% cluster_12, "12",
                                                                                    ifelse(PG_summ$Gene.names %in% cluster_13, "13",
                                                                                           ifelse(PG_summ$Gene.names %in% cluster_15, "15",
                                                                                                  ifelse(PG_summ$Gene.names %in% cluster_16, "16",
                                                                                                         ifelse(PG_summ$Gene.names %in% cluster_18, "18",
                                                                                                                NA))))))))))))



# Remove clusters with less than 3 elements

for (i in as.character(c(2:7,9,12,13,15,16,18))) {
  
  if (nrow(subset(PG_summ, Li_2017_Clusters == i)) < 3) {
    
    PG_summ[PG_summ$Li_2017_Clusters == i & !is.na(PG_summ$Li_2017_Clusters), "Li_2017_Clusters"] <- NA
    
  }
}


#### Plots ####

MyScatterPlot <- function(df, X, Y,
                          xmin = -5, xmax = 6,
                          ymin = -1, ymax = 5,
                          AxisSep.x = 2,
                          AxisSep.y = 2,
                          dot.size = 2,
                          Title = NULL, xTitle = NULL, yTitle = NULL,
                          InPercentage.x = T, sec_breaks.x = c(0.1,0.5,1,2,4,8,16,32,64), scaled_breaks.x = log2(sec_breaks.x),
                          InPercentage.y = F, sec_breaks.y = sec_breaks.x, scaled_breaks.y = log2(sec_breaks.y),
                          GradientColorMin = rgb(.5,.5,.5), GradientColorMax = rgb(0,0,0), bin.size = 128) {
  
  plot <- ggplot(data = df, aes(x = eval(parse(text = X)), y = eval(parse(text = Y)))) +
    
    geom_point(na.rm = T,
               size = dot.size,
               colour = densCols(x = df[,X],
                                 y = df[,Y],
                                 colramp = colorRampPalette(c(GradientColorMin,GradientColorMax)), nbin = bin.size)) +
    
    ggtitle(Title) +
    xlab(xTitle) +
    ylab(yTitle) +
    coord_cartesian(xlim = c(xmin,xmax),ylim = c(ymin,ymax)) +
    theme_bw() +
    theme(plot.title = element_text(#face="bold",
      size = 35,
      hjust = 0.5,
      vjust = 0.4),
      axis.title.x = element_text(#face="bold",
        size=30,
        hjust = 0.5,
        vjust = 0.4),
      axis.text.x  = element_text(#face = "bold",
        color = "black",
        angle=0, 
        vjust=0.5,
        hjust = 0.5,
        size=25),
      axis.title.y = element_text(#face="bold",
        size=30,
        hjust = 0.5,
        vjust = 1.5),
      axis.text.y  = element_text(#face = "bold",
        color = "black",
        angle=0, 
        vjust=0.5,
        hjust = 0.5,
        size=25),
      panel.grid=element_blank(),
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
  
  return(plot)
}




main_plot <- MyScatterPlot(PG_summ,
                           "Ratio.H.L.Mean", "iBAQ.Input.Mean",
                           ymin = 17, ymax = 34)





plot_genes_in_Li_clusters <- main_plot +
  
  geom_point(data = subset(PG_summ, !is.na(Li_2017_Clusters)),
             aes(color = Li_2017_Clusters),
             na.rm = T,
             size = 3,
             show.legend = F) +
  
  geom_point(data = subset(PG_summ, Gene.names == "WDR33"),
             shape = 17,
             na.rm = T,
             size = 3,
             show.legend = F)




with_names <- plot_genes_in_Li_clusters +
  
  geom_text_repel(data = subset(PG_summ, !is.na(Li_2017_Clusters)),
                  aes(label = Gene.names,
                      color = Li_2017_Clusters),
                  na.rm = T,
                  #fontface = "bold",
                  size = 4,
                  show.legend = T, max.overlaps = 100)



png("Li_2017_clusters.png", width = 1500, height = 500, pointsize = 25)
grid.arrange(main_plot, plot_genes_in_Li_clusters, with_names, nrow = 1)
dev.off()




pdf("Li_2017_clusters.pdf", width = 5, height = 5, useDingbats = F)
grid.arrange(with_names, nrow = 1)
dev.off()

rm(main_plot, plot_genes_in_Li_clusters, with_names)
