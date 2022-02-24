setwd("K:/Results/20210903_Phospho_qRIC_revisionPaper/DYRK3_interactors")

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


#### Load Working table ####

working_table_PTM <- read.csv("../output/working_table_PTM.txt", sep = "\t", stringsAsFactors = FALSE)

working_table_PTM_filtered <- 
  subset(working_table_PTM, ((Ratio.H.L.Forward_PG > DeltaAnalysisThresholdEfficiency & 
                                Ratio.H.L.Reverse_PG > DeltaAnalysisThresholdEfficiency) |
                               (Ratio.H.L.Forward_PTM > DeltaAnalysisThresholdEfficiency & 
                                  Ratio.H.L.Reverse_PTM > DeltaAnalysisThresholdEfficiency)) &
           !is.na(delta_Mean) & Known_mRBP_PG)


# working_table_PTM_filtered <- 
#   subset(working_table_PTM, !is.na(delta_Mean) & Known_mRBP_PG)




#### Selected proteins ####
# CORUM complexes 
DYRK3_interactors <- fread("K:/Datasets/2018_Rai_DYRK3_interactors.txt", sep = "\t", stringsAsFactors = FALSE, data.table = F, check.names = T, skip = 19)


DYRK3_interactors_Uniprot <- DYRK3_interactors$Protein_ID  
DYRK3_interactors <- DYRK3_interactors$Protein_symbol



#### Plots ####
MyScatterPlot <- function(df = working_table_PTM_filtered, X = "delta_Mean",
                          Y = "Ratio.H.L.Mean_PG",
                          proteins_interest = Spliceosome,
                          xmin = -5, xmax = 5,
                          AxisSep.x = 2,
                          ymin = -5, ymax = 6,
                          AxisSep.y = 1,
                          DotSize = 3,
                          Title = NULL, xTitle = "Delta efficiency (log2)", yTitle = "Protein pulldown efficiency (%)",
                          InPercentage.x = F, sec_breaks.x = c(0.1,0.5,1,2,4,8,16,32,64), scaled_breaks.x = log2(sec_breaks.x),
                          InPercentage.y = T, sec_breaks.y = sec_breaks.x, scaled_breaks.y = log2(sec_breaks.y),
                          bin.size = 128,
                          GradientColorMin = rgb(.9,.45,0,.8), GradientColorMax = rgb(.3,.15,0,.8)) {
  
  df <- subset(df, !is.na(eval(parse(text = X))) & !is.na(eval(parse(text = Y))))
  
  foo_background <- subset(df, !(Gene.names_PTM %in% proteins_interest))
  
  foo_interest <- subset(df, Gene.names_PTM %in% proteins_interest)
  foo_interest_diff <- subset(foo_interest, PTM_increase_Efficiency_Cl_Mean | PTM_decrease_Efficiency_Cl_Mean)

  
  
  
{  plot <- ggplot(data = df, aes(x = eval(parse(text = X)), y = eval(parse(text = Y)))) +
    
    ggtitle(Title) +
    xlab(xTitle) +
    ylab(yTitle) +
    coord_cartesian(xlim = c(xmin,xmax),ylim = c(ymin,ymax)) +
    theme_bw() +
    theme(plot.title = element_text(size = 14,
                                    hjust = 0.5,
                                    vjust = 0.4),
          axis.title.x = element_text(size=14,
                                      hjust = 0.5,
                                      vjust = 0.4),
          axis.text.x  = element_text(color = "black",
                                      angle=0, 
                                      vjust=0.5,
                                      hjust = 0.5,
                                      size=12),
          axis.title.y = element_text(size=14,
                                      hjust = 0.5,
                                      vjust = 1.5),
          axis.text.y  = element_text(color = "black",
                                      angle=0, 
                                      vjust=0.5,
                                      hjust = 0.5,
                                      size=12),
          panel.grid.major = element_line(colour = rgb(0,0,0,.1), linetype = "dashed"),
          panel.grid.minor = element_blank(),
          aspect.ratio=1)
  }
  
    
  plot <- plot +
    
  geom_point(na.rm = T,
             size = DotSize,
             colour = densCols(x = df[,X],
                               y = df[,Y],
                               colramp = colorRampPalette(c(GradientColorMin, GradientColorMax)), nbin = bin.size)) +
    
    geom_point(data = foo_interest,
               na.rm = T,
               size = DotSize,
               colour = rgb(.4,0,.8))
  
  
    if(nrow(foo_interest_diff) <= 40) {
      plot <- plot +
        geom_text_repel(data = foo_interest_diff,
                        aes(label = Site.ID_PTM),
                        size = 3,
                        colour = rgb(0,0,0),
                        na.rm = T)
    }
    
  
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


MyCumulativePlot <- function(df = working_table_PTM_filtered, X = "delta_Mean",
                             proteins_interest = Spliceosome,
                          xmin = -3, xmax = 3,
                          AxisSep.x = 1,
                          ymin = 0, ymax = 1,
                          AxisSep.y = .2,
                          LineSize = 3,
                          Title = NULL, xTitle = "Delta efficiency (log2)", yTitle = "Cumulative probability",
                          InPercentage.x = F, sec_breaks.x = c(0.1,0.5,1,2,4,8,16,32,64), scaled_breaks.x = log2(sec_breaks.x),
                          InPercentage.y = F, sec_breaks.y = sec_breaks.x, scaled_breaks.y = log2(sec_breaks.y)) {
  
  df <- subset(df, !is.na(eval(parse(text = X))))
  
  foo_background <- subset(df, !(Gene.names_PTM %in% proteins_interest))
  
  foo_interest <- subset(df, Gene.names_PTM %in% proteins_interest)
  
  
  {  plot <- ggplot(foo_background,
                    aes(x = eval(parse(text = X)))) +
      geom_vline(xintercept = 0) +
      geom_hline(yintercept = 0.5) +
      
      coord_cartesian(xlim = c(xmin, xmax), ylim = c(ymin, ymax)) +
      ylab(yTitle) +
      xlab(xTitle) +
      labs(title = Title) +
      theme_bw() +
      theme(plot.title = element_text(size = 14,
                                      hjust = 0.5,
                                      vjust = 0.4),
            axis.title.x = element_text(size=14,
                                        hjust = 0.5,
                                        vjust = 0.4),
            axis.text.x  = element_text(color = "black",
                                        angle=0, 
                                        vjust=0.5,
                                        hjust = 0.5,
                                        size=12),
            axis.title.y = element_text(size=14,
                                        hjust = 0.5,
                                        vjust = 1.5),
            axis.text.y  = element_text(color = "black",
                                        angle=0, 
                                        vjust=0.5,
                                        hjust = 0.5,
                                        size=12),
            panel.grid.major = element_line(colour = rgb(0,0,0,.1), linetype = "dashed"),
            panel.grid.minor = element_blank(),
            aspect.ratio=1)
  }
  
  
  
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
  
  
  
{  plot <- plot +
    
    # Background
    annotate("text",
             x = xmin+2,
             y = ymax, 
             label = paste0("Background (", nrow(foo_background), ")"),
             colour = rgb(.9,.45,0,.8),
             size = LineSize) +
  
  
  geom_step(stat="ecdf",                         
              size = LineSize, 
              alpha= 1, 
              color = rgb(.9,.45,0,.8),
              na.rm = T) +
    
    # Interest
    annotate("text",
             x = xmin+2,
             y = ymax-0.1, 
             label = paste0("P-sites of interest (", nrow(foo_interest), ")"),
             colour = rgb(.4,0,.8),
             size = LineSize) +
    
    geom_step(data = foo_interest,            
              stat="ecdf",
              size = LineSize, 
              alpha= 1,
              color = rgb(.4,0,.8),
              na.rm = T) +
    
    # p.value
    annotate("text",
             x = 2,
             y = 0.5, 
             label = paste0("Wilcox test\n","p <= ",
                            signif(wilcox.test(foo_background[,X],
                                               foo_interest[,X])$p.value, 3)),
             colour = rgb(.4,0,.8),
             size = LineSize)
}
  return(plot)
  
}





# Xtree like plot
Xtree_plot_DYRK3_interactors <- MyScatterPlot(proteins_interest = DYRK3_interactors, Title = "DYRK3_interactors")

## Cumulative plot
CumPlot_delta_DYRK3_interactors <- MyCumulativePlot(proteins_interest = DYRK3_interactors, Title = "DYRK3_interactors")


png("Delta_funSILAC_DYRK3_interactors.png")
grid.arrange(Xtree_plot_DYRK3_interactors, CumPlot_delta_DYRK3_interactors,
             ncol = 2)
dev.off()



pdf("Delta_funSILAC_DYRK3_interactors.pdf", useDingbats = F)
grid.arrange(Xtree_plot_DYRK3_interactors, CumPlot_delta_DYRK3_interactors,
             ncol = 2)
dev.off()
rm(Xtree_plot_DYRK3_interactors, CumPlot_delta_DYRK3_interactors)