setwd("K:/Results/20210903_Phospho_qRIC_revisionPaper/RNABindingRegions/correlation")

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


name_meths <- c("PARCLIP_PARalyzer", "PARCLIP_Piranha") 
#"eCLIP", "iCLIP_Piranha", "iCLIP_CIMS")#, "HITSCLIP_Piranha", "HITSCLIP_CIMS")



#### Load number of peaks table ####
Peaks_POSTAR <- fread("../PeaksMeanGene_counts_POSTAR_2019.txt", sep = "\t", stringsAsFactors = FALSE, data.table = F, check.names = T)
Peaks_POSTAR <- subset(Peaks_POSTAR, Method == "PAR-CLIP,PARalyzer" |
                         Method == "PAR-CLIP,Piranha_0.01")

Peaks_ENCODE <- fread("../Peaks_counts_ENCODE.txt", sep = "\t", stringsAsFactors = FALSE, data.table = F, check.names = T)

# Filter ENCODE data to contain only RBPs also analysed by PAR-CLIP
#foo <- intersect(Peaks_POSTAR$RBP_name, Peaks_ENCODE$RBP_name)
#Peaks_POSTAR <- subset(Peaks_POSTAR, RBP_name %in% foo)
Peaks_ENCODE <- subset(Peaks_ENCODE, RBP_name %in% Peaks_POSTAR$RBP_name)


Peaks_summary <- rbind(Peaks_POSTAR, Peaks_ENCODE)
rm(Peaks_POSTAR, Peaks_ENCODE)




Peaks_summary[Peaks_summary$Method == "PAR-CLIP,Piranha_0.01", "Method"] <- "PARCLIP_Piranha"
Peaks_summary[Peaks_summary$Method == "iCLIP,Piranha_0.01", "Method"] <- "iCLIP_Piranha"
Peaks_summary[Peaks_summary$Method == "HITS-CLIP,Piranha_0.01", "Method"] <- "HITSCLIP_Piranha"
Peaks_summary[Peaks_summary$Method == "PAR-CLIP,PARalyzer", "Method"] <- "PARCLIP_PARalyzer"
Peaks_summary[Peaks_summary$Method == "iCLIP,CIMS", "Method"] <- "iCLIP_CIMS"
Peaks_summary[Peaks_summary$Method == "HITS-CLIP,CIMS", "Method"] <- "HITSCLIP_CIMS"



# Only select those to be analysed
Peaks_summary <- subset(Peaks_summary, Method %in% name_meths)


## Log transfor number of peaks
Peaks_summary$Number_peaks <- log10(Peaks_summary$Number_peaks)


## Organize cell line names
# Give HEK293 and HEK293T the same name
Peaks_summary$Cell_line <- ifelse(Peaks_summary$Cell_line == "HEK293T", "HEK293", Peaks_summary$Cell_line)



#### Load pull-down efficiency ####
PG_summ <- read.csv("../../output/PG_summary.txt", sep = "\t", stringsAsFactors = FALSE)

PG_summ <- subset(PG_summ, select = c("Majority.protein.IDs","Gene.names", "Known_mRBP",
                                      "Ratio.H.L.Forward",
                                      "Ratio.H.L.Reverse",
                                      "Ratio.H.L.Mean",
                                      "iBAQ.Input.Mean"))



#### Main functions ####
MyScatterPlot <- function(xmin = -5, xmax = 6,
                          AxisSep.x = 1,
                          ymin = .1, ymax = .8,
                          AxisSep.y = .1,
                          InPercentage.x = T, sec_breaks.x = c(0.1,0.5,1,2,4,8,16,32,64), scaled_breaks.x = log2(sec_breaks.x),
                          InPercentage.y = F, sec_breaks.y = sec_breaks.x, scaled_breaks.y = log2(sec_breaks.y),
                          Title = NULL, xTitle = NULL, yTitle = NULL) {
  
  plot <- ggplot() +
    
    ggtitle(Title) +
    xlab(xTitle) +
    ylab(yTitle) +
    coord_cartesian(xlim = c(xmin,xmax),ylim = c(ymin,ymax)) +
    theme_bw() +
    theme(plot.title = element_text(#face="bold",
      size = 25,
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
      aspect.ratio = 1)
  
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





MyHistogramPlot <- function(df, X = "r",
                            xTitle = NULL, yTitle = NULL,
                            data.intercept = NULL) {
  
  foo <- subset(df, type == "randomized")
  if (data.intercept >= 0) {
    final.pvalue_randomized <- (nrow(subset(foo, eval(parse(text = X)) >= data.intercept)) + 1)/(nrow(foo)+1)
  } else {
    final.pvalue_randomized <- (nrow(subset(foo, eval(parse(text = X)) <= data.intercept)) + 1)/(nrow(foo)+1)
  }
  
  
  foo <- subset(df, type == "boot")
  if (data.intercept >= 0) {
    final.pvalue_boot <- (nrow(subset(foo, eval(parse(text = X)) >= data.intercept)) + 1)/(nrow(foo)+1)
  } else {
    final.pvalue_boot <- (nrow(subset(foo, eval(parse(text = X)) <= data.intercept)) + 1)/(nrow(foo)+1)
  }
  
  
  
  plot <- ggplot(data = df, aes(x = eval(parse(text = X)))) +
    
    geom_density(aes(fill = type),
                 color = rgb(0,0,0,0), alpha = .5) +
    
    geom_vline(xintercept = 0) +
    
    geom_vline(xintercept = data.intercept, linetype = "dashed") +
    
    annotate(geom = "text", x = 0, y = 3, 
             label = paste0(X, " = ", signif(data.intercept, 2),
                   "\n",
                   "p.value (bootstrap) = ", signif(final.pvalue_boot, 3),
                   "\n",
                   "p.value (randomized) = ", signif(final.pvalue_randomized, 3))) +
    coord_cartesian(xlim = c(-1,1)) +
    xlab(xTitle) +
    ylab(yTitle) +
    theme_bw() +
    theme(plot.title = element_text(#face="bold",
      size = 25,
      hjust = 0.5,
      vjust = 0.4),
      axis.title.x = element_text(#face="bold",
        size=30,
        hjust = 0.5,
        vjust = 0.4),
      axis.text.x  = element_text(#face = "bold", color = "black",
        angle=0, 
        vjust=0.5,
        hjust = 0.5,
        size=25),
      axis.title.y = element_text(#face="bold",
        size=30,
        hjust = 0.5,
        vjust = 1.5),
      axis.text.y  = element_text(#face = "bold", color = "black",
        angle=0, 
        vjust=0.5,
        hjust = 0.5,
        size=25),
      panel.grid=element_blank(),
      aspect.ratio = 1)
  
  return(plot)
}




myRandomAndBoot <- function(df = PG_summ, StableVar = "Ratio.H.L.Mean", VariableVar,
                            pointer = 10000) {
  set.seed(1)
  value_r_random <- c()
  value_p_random <- c()
  
  value_r_boot <- c()
  value_p_boot <- c()
  
  while (pointer > 0) {
    foo <- subset(df, select = c(StableVar, VariableVar))
    
    foo <- subset(foo, !is.na(foo[,VariableVar]))
    
    #Random
    foo_random <- foo
    foo_random[,VariableVar] <- sample(foo_random[,VariableVar], nrow(foo_random), replace = F)
    
    
    xcor <- cor(foo_random[,StableVar], foo_random[,VariableVar], method = "spearman")
    value_r_random <- c(value_r_random, xcor)
    
    #fit_random <- lm(foo_random[,StableVar] ~ foo_random[,VariableVar], data = foo_random)
    #value_r_random <- c(value_r_random, summary(fit_random)$adj.r.squared)
    #value_p_random <- c(value_p_random, -log10(summary(fit_random)$coef[2,4]))
    
    
    #Bootstrap
    foo_boot <- foo
    foo_boot <- foo_boot[sample(1:nrow(foo_boot), nrow(foo_boot), replace = T),]
    
    
    xcor <- cor(foo_boot[,StableVar], foo_boot[,VariableVar], method = "spearman")
    value_r_boot <- c(value_r_boot, xcor)
    
    #fit_boot <- lm(foo_boot[,StableVar] ~ foo_boot[,VariableVar], data = foo_boot)
    #value_r_boot <- c(value_r_boot, summary(fit_boot)$adj.r.squared)
    #value_p_boot <- c(value_p_boot, -log10(summary(fit_boot)$coef[2,4]))
    
    
    pointer <- pointer -1
  }
  
  
  return(rbind(data.frame(r = value_r_random, type = "randomized"),
               data.frame(r = value_r_boot, type = "boot")))
  
}









#### plots correlation ####
MyPlotCorr <- function(Filt_method = "PARCLIP_Piranha", Filt_cell_line = NULL) {
  #### Receives a filtering method and returns a list with the correlation values
  
  ## Filter method
  df_filtered <- subset(Peaks_summary, Method == Filt_method)
  
  ## Filter by cell line
  if(!is.null(Filt_cell_line)) df_filtered <- subset(df_filtered, Cell_line == Filt_cell_line)
  
  ## Merge with efficiency
  PG_summ_merged <- merge(PG_summ, df_filtered, 
                          by.x = "Gene.names", by.y = "RBP_name", all.x = F, all.y = F)
  
  names(PG_summ_merged)[names(PG_summ_merged) == "Number_peaks"] <- paste0(Filt_method,
                                                                           "_Number_peaks")
  
PG_summ_merged <- subset(PG_summ_merged, !is.na(Ratio.H.L.Mean))


  ## Calculate correlation

foo <- subset(PG_summ_merged, !is.na(eval(parse(text = paste0(Filt_method,
                                                              "_Number_peaks")))))

x.cor <- cor.test(foo$Ratio.H.L.Mean, foo[,paste0(Filt_method, "_Number_peaks")], method = "spearman")



## Plot
pca_rotation <- prcomp(cbind(PG_summ_merged$Ratio.H.L.Mean,
                             PG_summ_merged[,paste0(Filt_method,"_Number_peaks")]))$rotation
beta <- as.numeric(pca_rotation[2,1]/pca_rotation[1,1])
intercept <- mean(PG_summ_merged[,paste0(Filt_method,"_Number_peaks")], na.rm = T) - (beta*mean(PG_summ_merged$Ratio.H.L.Mean, na.rm = T))


  plot_scatter_Number_peaks <- MyScatterPlot(ymin = -1 + round(min(PG_summ_merged[,paste0(Filt_method, "_Number_peaks")]),0), 
                                             ymax = 1 + round(max(PG_summ_merged[,paste0(Filt_method, "_Number_peaks")]),0),
                                             AxisSep.y = 1,
                                             InPercentage.y = F, sec_breaks.y = 10^seq(2,6,1)) +
    
    geom_point(data = PG_summ_merged,
               aes(x = Ratio.H.L.Mean,
                   y = eval(parse(text = paste0(Filt_method,
                                                "_Number_peaks")))),
               na.rm = T, size = 4, show.legend = T) +
    
    geom_segment(data = PG_summ_merged,
                 aes(x = min(Ratio.H.L.Mean, na.rm = T), 
                     xend = max(Ratio.H.L.Mean, na.rm = T),
                     y = -(intercept*DeltaAnalysisThresholdEfficiency + beta),
                     yend = intercept + beta*max(Ratio.H.L.Mean, na.rm = T))) +
    
    # geom_smooth(data = PG_summ_merged_filtered,
    #             aes(x = Ratio.H.L.Mean,
    #                 y = eval(parse(text = paste0(Filt_method,
    #                                              "_Number_peaks")))),
    #             method='lm',formula=y~x, na.rm = T) +
    
    labs(title = paste(paste0(Filt_method,"_",Filt_cell_line), "\n",
                       "rho = ",signif(x.cor$estimate, 3),
                       "p = ",signif(x.cor$p.value, 2))) #+
    # 
    # geom_text_repel(data = PG_summ_merged,
    #                 aes(label = Gene.names,
    #                     x = Ratio.H.L.Mean,
    #                     y = eval(parse(text = paste0(Filt_method,
    #                                                  "_Number_peaks")))),
    #                 na.rm = T, size = 5, show.legend = T)
  
  
  
  return(plot_scatter_Number_peaks)
}

#plot_PARCLIP_PARalyzer_HEK293T <- MyPlotCorr("PARCLIP_PARalyzer", Filt_cell_line = "HEK293T")
plot_PARCLIP_PARalyzer_HEK293 <- MyPlotCorr("PARCLIP_PARalyzer", Filt_cell_line = "HEK293")
#plot_PARCLIP_Piranha_HEK293T <- MyPlotCorr("PARCLIP_Piranha", Filt_cell_line = "HEK293T")
plot_PARCLIP_Piranha_HEK293 <- MyPlotCorr("PARCLIP_Piranha", Filt_cell_line = "HEK293")
#plot_eCLIP_HepG2 <- MyPlotCorr("eCLIP", Filt_cell_line = "HepG2")
#plot_eCLIP_K562 <- MyPlotCorr("eCLIP", Filt_cell_line = "K562")



png("Correlation_numberPeaks.png", width = 1000, height = 500, pointsize = 25)
grid.arrange(#plot_PARCLIP_PARalyzer_HEK293T,
             plot_PARCLIP_PARalyzer_HEK293,
             #plot_PARCLIP_Piranha_HEK293T, 
             plot_PARCLIP_Piranha_HEK293,
             #plot_eCLIP_HepG2,
             #plot_eCLIP_K562,
             nrow = 1)
dev.off()



pdf("Correlation_numberPeaks.pdf", width = 10, height = 5, useDingbats = F)
grid.arrange(#plot_PARCLIP_PARalyzer_HEK293T,
             plot_PARCLIP_PARalyzer_HEK293,
             #plot_PARCLIP_Piranha_HEK293T, 
             plot_PARCLIP_Piranha_HEK293,
             #plot_eCLIP_HepG2,
             #plot_eCLIP_K562,
             nrow = 1)
dev.off()



rm(#plot_PARCLIP_PARalyzer_HEK293T, 
   plot_PARCLIP_PARalyzer_HEK293,
   #plot_PARCLIP_Piranha_HEK293T,
   plot_PARCLIP_Piranha_HEK293,
   plot_eCLIP_HepG2, plot_eCLIP_K562)

#### Correlation p-values ####
MyPlotCorrPvalue <- function(Filt_method = "PARCLIP_PARalyzer",
                             Filt_cell_line = "HEK293") {
  
  #### Receives a filtering method and returns a list with the correlation values
  
  ## Filter method
  df_filtered <- subset(Peaks_summary, Method == Filt_method)
  
  ## Filter by cell line
  if(!is.null(Filt_cell_line)) df_filtered <- subset(df_filtered, Cell_line == Filt_cell_line)
  
  ## Merge with efficiency
  PG_summ_merged <- merge(PG_summ, df_filtered, 
                          by.x = "Gene.names", by.y = "RBP_name", all.x = F, all.y = F)
  
  names(PG_summ_merged)[names(PG_summ_merged) == "Number_peaks"] <- paste0(Filt_method,
                                                                           "_Number_peaks")
  
  PG_summ_merged <- subset(PG_summ_merged, !is.na(Ratio.H.L.Mean))
  
  
  
  ## Calculate correlation
  foo <- subset(PG_summ_merged, !is.na(eval(parse(text = paste0(Filt_method,
                                                                "_Number_peaks")))))
  
  x.cor <- cor(foo$Ratio.H.L.Mean, foo[,paste0(Filt_method, "_Number_peaks")], method = "spearman")
  
  
  ## Make plot
  plot_r <- MyHistogramPlot(df = myRandomAndBoot(df = PG_summ_merged,
                                                 VariableVar = paste0(Filt_method,
                                                                      "_Number_peaks")),
                            X = "r",
                            data.intercept = x.cor)
  
  return(plot_r + ggtitle(Filt_cell_line))
}



#plot_PARCLIP_PARalyzer_HEK293T <- MyPlotCorrPvalue("PARCLIP_PARalyzer", Filt_cell_line = "HEK293T")
plot_PARCLIP_PARalyzer_HEK293 <- MyPlotCorrPvalue("PARCLIP_PARalyzer", Filt_cell_line = "HEK293")
#plot_PARCLIP_Piranha_HEK293T <- MyPlotCorrPvalue("PARCLIP_Piranha", Filt_cell_line = "HEK293T")
plot_PARCLIP_Piranha_HEK293 <- MyPlotCorrPvalue("PARCLIP_Piranha", Filt_cell_line = "HEK293")
#plot_eCLIP_HepG2 <- MyPlotCorrPvalue("eCLIP", Filt_cell_line = "HepG2")
#plot_eCLIP_K562 <- MyPlotCorrPvalue("eCLIP", Filt_cell_line = "K562")



png("Correlation_numberPeaks_pvalue.png", width = 1000, height = 1000, pointsize = 25)
grid.arrange(#plot_PARCLIP_PARalyzer_HEK293T,
             plot_PARCLIP_PARalyzer_HEK293,
             #plot_PARCLIP_Piranha_HEK293T,
             plot_PARCLIP_Piranha_HEK293,
             #plot_eCLIP_HepG2,
             #plot_eCLIP_K562,
             nrow = 1)
dev.off()
rm(#plot_PARCLIP_PARalyzer_HEK293T,
   plot_PARCLIP_PARalyzer_HEK293,
   #plot_PARCLIP_Piranha_HEK293T,
   plot_PARCLIP_Piranha_HEK293,
   plot_eCLIP_HepG2, plot_eCLIP_K562)




