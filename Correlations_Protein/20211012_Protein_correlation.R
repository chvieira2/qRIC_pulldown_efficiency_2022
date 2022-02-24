setwd("K:/Results/20210903_Phospho_qRIC_revisionPaper/Correlations_Protein")

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



#### Genes in List ####
PG_summ <- read.csv("../output/PG_summary.txt", sep = "\t", stringsAsFactors = FALSE)

PG_summ <- subset(PG_summ, #Known_mRBP &
                  !is.na(Ratio.H.L.Mean),
                  select = c("Majority.protein.IDs","Gene.names",
                             "Ratio.H.L.Forward",
                             "Ratio.H.L.Reverse",
                             "Ratio.H.L.Mean",
                             "iBAQ.Input.Mean"))


PG_summ$iBAQ.Input.Mean <- log10(2**PG_summ$iBAQ.Input.Mean)



PG_summ_filtered <- subset(PG_summ, (Ratio.H.L.Forward > DeltaAnalysisThresholdEfficiency & 
                                       Ratio.H.L.Reverse > DeltaAnalysisThresholdEfficiency))



#### Plots ####
MyScatterPlot <- function(df, X, Y,
                          xmin = -5, xmax = 6, AxisSep.x = 1,
                          ymin = 0, ymax = 6, AxisSep.y = 1,
                          Title = NULL, xTitle = NULL, yTitle = NULL,
                          InPercentage.x = T, sec_breaks.x = c(0.1,0.5,1,2,4,8,16,32,64),
                          scaled_breaks.x = log2(sec_breaks.x),
                          InPercentage.y = F, sec_breaks.y = sec_breaks.x,
                          scaled_breaks.y = log2(sec_breaks.y),
                          GradientColorMin = rgb(.5,.5,.5), GradientColorMax = rgb(0,0,0),
                          bin.size = 128) {
  
  plot <- ggplot(data = df, aes(x = eval(parse(text = X)), y = eval(parse(text = Y)))) +
    
    geom_point(na.rm = T,
               size = 3,
               colour = densCols(x = df[,X],
                                 y = df[,Y],
                                 colramp = colorRampPalette(c(GradientColorMin,GradientColorMax)), nbin = bin.size)) +
    
    ggtitle(Title) +
    xlab(xTitle) +
    ylab(yTitle) +
    coord_cartesian(xlim = c(xmin,xmax),ylim = c(ymin,ymax)) +
    theme_bw() +
    theme(plot.title = element_text(size = 25,
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
          axis.title.y = element_text(face="bold",
                                      size=30,
                                      hjust = 0.5,
                                      vjust = 1.5),
          axis.text.y  = element_text(color = "black",
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





MyHistogramPlot <- function(df, X = "r",
                            Title = NULL, xTitle = NULL, yTitle = NULL,
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
    
    ggtitle(paste0(X, " = ", signif(data.intercept, 2),
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
                            pointer = 1000) {
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
    
    
    xcor <- cor(foo_random[,StableVar], foo_random[,VariableVar], method = "pearson")
    value_r_random <- c(value_r_random, xcor)
    
    #fit_random <- lm(foo_random[,StableVar] ~ foo_random[,VariableVar], data = foo_random)
    #value_r_random <- c(value_r_random, summary(fit_random)$adj.r.squared)
    #value_p_random <- c(value_p_random, -log10(summary(fit_random)$coef[2,4]))
    
    
    #Bootstrap
    foo_boot <- foo
    foo_boot <- foo_boot[sample(1:nrow(foo_boot), nrow(foo_boot), replace = T),]
    
    
    xcor <- cor(foo_boot[,StableVar], foo_boot[,VariableVar], method = "pearson")
    value_r_boot <- c(value_r_boot, xcor)
    
    #fit_boot <- lm(foo_boot[,StableVar] ~ foo_boot[,VariableVar], data = foo_boot)
    #value_r_boot <- c(value_r_boot, summary(fit_boot)$adj.r.squared)
    #value_p_boot <- c(value_p_boot, -log10(summary(fit_boot)$coef[2,4]))
    
    
    pointer <- pointer -1
  }
  
  
  return(rbind(data.frame(r = value_r_random, type = "randomized"),
               data.frame(r = value_r_boot, type = "boot")))
  
}





#### Input abundance ####
foo_abund <- subset(PG_summ_filtered, !is.na(iBAQ.Input.Mean))

x.cor_abund <- cor.test(foo_abund$Ratio.H.L.Mean, foo_abund$iBAQ.Input.Mean, method = "spearman")

pca_rotation <- prcomp(cbind(foo_abund$Ratio.H.L.Mean, foo_abund$iBAQ.Input.Mean))$rotation
beta <- as.numeric(pca_rotation[2,1]/pca_rotation[1,1])
intercept <- mean(foo_abund$iBAQ.Input.Mean, na.rm = T) - (beta*mean(foo_abund$Ratio.H.L.Mean, na.rm = T))




Abund_plot <- MyScatterPlot(PG_summ,
                            "Ratio.H.L.Mean", "iBAQ.Input.Mean",
                            AxisSep.y = 1, InPercentage.x = T,
                            ymin = 5, ymax = 10) +
  
  geom_segment(aes(x = DeltaAnalysisThresholdEfficiency, xend = max(Ratio.H.L.Mean, na.rm = T),
                   y = -(intercept*DeltaAnalysisThresholdEfficiency + beta), yend = intercept + beta*max(Ratio.H.L.Mean, na.rm = T))) +
  
  labs(title = paste("R = ",signif(x.cor_abund$estimate, 3),
                     "p = ",signif(x.cor_abund$p.value, 2)))

#labs(title = paste("R = ",signif(sqrt(summary(fit_abund)$adj.r.squared), 3),
#                   "Adj R2 = ",signif(summary(fit_abund)$adj.r.squared, 3), "\n",
#                  "P =",signif(summary(fit_abund)$coef[2,4], 3)))




plot_r <- MyHistogramPlot(df = myRandomAndBoot(df = foo_abund, VariableVar = "iBAQ.Input.Mean"), X = "r",
                          data.intercept = x.cor_abund$estimate)






png("Correlations_Abundance.png", width = 1000, height = 500)
grid.arrange(Abund_plot, plot_r,
             nrow = 1)
dev.off()

rm(fit_abund, plot_r)







### Figure paper
pdf("Correlations_Abundance.pdf", width = 5, height = 5, useDingbats = F)
grid.arrange(Abund_plot,
             nrow = 1)
dev.off()

rm(Abund_plot)

