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


DeltaAnalysisThresholdEfficiency <- -1


#### Load Working table ####

working_table_PTM <- read.csv("working_table_PTM.txt", sep = "\t", stringsAsFactors = FALSE)
working_table_PEP <- read.csv("working_table_PEP.txt", sep = "\t", stringsAsFactors = FALSE)

working_table_PEP <- subset(working_table_PEP, Proteins_PEP %in% working_table_PTM$Protein_PTM)

#### Filter working table ####
working_table_PEP_high <- subset(working_table_PEP, #Known_mRBP_PG & 
                                   ((Ratio.H.L.Forward_PG > DeltaAnalysisThresholdEfficiency &
                                       Ratio.H.L.Reverse_PG > DeltaAnalysisThresholdEfficiency) | 
                                      (Ratio.H.L.Forward_PEP > DeltaAnalysisThresholdEfficiency &
                                         Ratio.H.L.Reverse_PEP > DeltaAnalysisThresholdEfficiency)) & 
                              !is.na(delta_Mean))


working_table_PTM_high <- subset(working_table_PTM, #Known_mRBP_PG & 
                                   ((Ratio.H.L.Forward_PG > DeltaAnalysisThresholdEfficiency &
                                       Ratio.H.L.Reverse_PG > DeltaAnalysisThresholdEfficiency) | 
                                      (Ratio.H.L.Forward_PTM > DeltaAnalysisThresholdEfficiency &
                                         Ratio.H.L.Reverse_PTM > DeltaAnalysisThresholdEfficiency)) & 
                                   !is.na(delta_Mean))


working_table_PTM <- subset(working_table_PTM, Known_mRBP_PG)
working_table_PEP <- subset(working_table_PEP, Known_mRBP_PG)


#### Write table ####
Candidate_sites <- subset(working_table_PTM_high, select = c("Protein_PTM",
                                                        "Gene.names_PTM",
                                                        "Peptide.Site.ID_PTM",
                                                        "Site.ID_PTM",
                                                        "Amino.acid_PTM",
                                                        "Position_PTM",
                                                        "Known_RBP_PG",
                                                        "Known_mRBP_PG",
                                                        
                                                        "Mol..weight..kDa._PG",
                                                        "Sequence.length_PG",
                                                        
                                                        "Ratio.H.L.Forward_II_PTM",
                                                        "Ratio.H.L.Forward_III_PTM",
                                                        "Ratio.H.L.Forward_IV_PTM",
                                                        "Ratio.H.L.Reverse_II_PTM",
                                                        "Ratio.H.L.Reverse_III_PTM",
                                                        "Ratio.H.L.Reverse_IV_PTM",
                                                        
                                                        "Ratio.H.L.Mean_PTM",
                                                        
                                                        "Ratio.H.L.Forward_II_PG",
                                                        "Ratio.H.L.Forward_III_PG",
                                                        "Ratio.H.L.Forward_IV_PG",
                                                        "Ratio.H.L.Reverse_II_PG",
                                                        "Ratio.H.L.Reverse_III_PG",
                                                        "Ratio.H.L.Reverse_IV_PG",
                                                        
                                                        "Ratio.H.L.Mean_PG",
                                                        
                                                        
                                                        "delta_Forward_II",
                                                        "delta_Forward_III",
                                                        "delta_Forward_IV",
                                                        "delta_Reverse_II",
                                                        "delta_Reverse_III",
                                                        "delta_Reverse_IV",
                                                        
                                                        "delta_Mean",
                                                        
                                                        "PTM_increase_Efficiency_Cl_Mean",
                                                        "PTM_decrease_Efficiency_Cl_Mean"))


names(Candidate_sites) <- c("Protein",
                            "Gene.names",
                            "Peptide.Site.ID",
                            "Site.ID",
                            "Amino.acid",
                            "Position",
                            "Known_RBP",
                            "Known_mRBP",
                            
                            "Mol.weight.kDa",
                            "Sequence.length",
                            
                            "Eff.log2.Phospho.Forward_II",
                            "Eff.log2.Phospho.Forward_III",
                            "Eff.log2.Phospho.Forward_IV",
                            "Eff.log2.Phospho.Reverse_II",
                            "Eff.log2.Phospho.Reverse_III",
                            "Eff.log2.Phospho.Reverse_IV",
                            
                            "Eff.log2.Phospho.Mean",
                            
                            "Eff.log2.Protein.Forward_II",
                            "Eff.log2.Protein.Forward_III",
                            "Eff.log2.Protein.Forward_IV",
                            "Eff.log2.Protein.Reverse_II",
                            "Eff.log2.Protein.Reverse_III",
                            "Eff.log2.Protein.Reverse_IV",
                            
                            "Eff.log2.Protein.Mean",
                            
                            
                            "delta_Eff_Forward_II",
                            "delta_Eff_Forward_III",
                            "delta_Eff_Forward_IV",
                            "delta_Eff_Reverse_II",
                            "delta_Eff_Reverse_III",
                            "delta_Eff_Reverse_IV",
                            
                            "delta_Eff_Mean",
                            
                            "PTM.Increase.Eff_FC2",
                            "PTM.Decrease.Eff_FC2")




fwrite(Candidate_sites, file = "Candidate_sites.txt", sep = "\t", na = "", quote = F, row.names = F)
rm(Candidate_sites)






#### Plots ####
MyScatterPlot <- function(df, X, Y,
                          xmin = -1, xmax = 5,
                          AxisSep.x = 2,
                          ymin = -1, ymax = 5,
                          AxisSep.y = 2,
                          DotSize = 2,
                          Title = NULL, xTitle = NULL, yTitle = NULL,
                          InPercentage.x = F, sec_breaks.x = c(0.1,0.5,1,2,4,8,16,32), scaled_breaks.x = log2(sec_breaks.x),
                          InPercentage.y = F, sec_breaks.y = sec_breaks.x, scaled_breaks.y = log2(sec_breaks.y),
                          DotShape = "Amino.acid_PTM",
                          bin.size = 128,
                          GradientColorMin = rgb(.9,.45,0,.8),
                          GradientColorMax = rgb(.3,.15,0,.8)) {
  
  plot <- ggplot(data = df, aes(x = eval(parse(text = X)), y = eval(parse(text = Y)))) +
    
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
  
  
  
  if(!is.null(DotShape)) {
    plot <- plot +
      geom_point(aes(shape = eval(parse(text = DotShape))),
                 na.rm = T,
                 size = DotSize,
                 colour = densCols(x = df[,X],
                                   y = df[,Y],
                                   colramp = colorRampPalette(c(GradientColorMin, GradientColorMax)), nbin = bin.size)) +
      
      theme(legend.title = element_blank())
  }
  return(plot)
}









#### Figure paper xTree and Vulcano ####
Maybe_interesting_pSites <- c(working_table_PTM[working_table_PTM$Gene.names_PTM %in% c(#"SERBP1",
  "RBM20",
  #"LARP1"#,
  "ELAVL1",
  "UPF1",
  "SF3B1"
),
"Site.ID_PTM"])

foo <- subset(working_table_PTM, !(PTM_increase_Efficiency_Cl_Mean | PTM_decrease_Efficiency_Cl_Mean))
foo1 <- subset(working_table_PTM, PTM_increase_Efficiency_Cl_Mean | PTM_decrease_Efficiency_Cl_Mean)
foo2 <- subset(working_table_PTM, Site.ID_PTM %in% Maybe_interesting_pSites)

# Vulcano plot
Vulcano_plot <- MyScatterPlot(foo,
                              X = "delta_Mean", Y = "Neglog10pval_Cl_Mean",
                              xmin = -5, xmax = 5,
                              ymin = 0, ymax = 4, AxisSep.y = 1) +
  
  
  geom_segment(aes(x = -log2(2), y = 1, xend = -1, yend = 10), size = .5, color = rgb(0,.4,.8), linetype = 2) +
  geom_segment(aes(x = log2(2), y = 1, xend = 1, yend = 10), size = .5, color = rgb(0,.4,.8), linetype = 2) +
  
  geom_segment(aes(x = -10, y = 1, xend = -1, yend = 1), size = .5, color = rgb(0,.4,.8), linetype = 2) +
  geom_segment(aes(x = 1, y = 1, xend = 10, yend = 1), size = .5, color = rgb(0,.4,.8), linetype = 2) +
  
  
  geom_point(data = foo1,
             aes(shape = Amino.acid_PTM),
             na.rm = T,
             size = 2,
             colour = rgb(0,.4,.8)) +
  
  
  geom_point(data = foo2,
             aes(shape = Amino.acid_PTM),
             na.rm = T,
             size = 2,
             stroke = 1,
             #shape = 21,
             colour = rgb(0,0,0)) +
  
  
  geom_text_repel(data = foo2,
                  aes(label = Site.ID_PTM),
                  size = 4,
                  colour = rgb(0,0,0),
                  na.rm = T)




# Xtree like plot
Xtree_plot <- MyScatterPlot(foo,
                            X = "delta_Mean", Y = "Ratio.H.L.Mean_PG",
                            xmin = -5, xmax = 5,
                            ymin = -5, ymax = 6, AxisSep.y = 1, InPercentage.y = T) +
  
  
  geom_segment(aes(x = -1, y = -10, xend = -1, yend = 10), size = .5, color = rgb(0,.4,.8), linetype = 2) +
  
  geom_segment(aes(x = 1, y = -10, xend = 1, yend = 10), size = .5, color = rgb(0,.4,.8), linetype = 2) +
  
  
  geom_point(data = foo1,
             aes(shape = Amino.acid_PTM),
             na.rm = T,
             size = 2,
             colour = rgb(0,.4,.8)) +
  
  
  geom_point(data = foo2,
             aes(shape = Amino.acid_PTM),
             na.rm = T,
             size = 2,
             stroke = 1,
             #shape = 21,
             colour = rgb(0,0,0)) +
  
  
  geom_text_repel(data = foo2,
                  aes(label = Site.ID_PTM),
                  size = 4,
                  colour = rgb(0,0,0),
                  na.rm = T)



pdf("Xtree_Vulcano_DeltaEff.pdf", width = 10, height = 5, useDingbats = F)
grid.arrange(Xtree_plot, Vulcano_plot,
             nrow = 1)
dev.off()
rm(Xtree_plot, Vulcano_plot, foo, foo1, foo2)




#### Vulcano plot per number samples ####
## Removing irreproducible values in working table has the side effect of only letting forward and reverse pairs in the table, so delta efficiencies are measured in 2, 4 or 6 samples

foo <- subset(working_table_PTM, select = c("delta_Mean","Neglog10pval_Cl_Mean","delta_Forward_II",
                                            "delta_Forward_III",
                                            "delta_Forward_IV",
                                            "delta_Reverse_II",
                                            "delta_Reverse_III",
                                            "delta_Reverse_IV"))

foo$N_delta <- apply(foo, 1,
             function(x) factor(sum(!is.na(x[c("delta_Forward_II",
                                        "delta_Forward_III",
                                        "delta_Forward_IV",
                                        "delta_Reverse_II",
                                        "delta_Reverse_III",
                                        "delta_Reverse_IV")])),
                                   levels = c(2,3,4,5,6)))

foo_2 <- subset(foo, !is.na(delta_Mean) & !is.na(Neglog10pval_Cl_Mean))


MyVulcanoNvalues <- function(df = foo,
                             X = "delta_Mean_all", Y = "Neglog10pval_Cl_Mean",
                          xmin = -5, xmax = 5,
                          AxisSep.x = 2,
                          ymin = 0, ymax = 4,
                          AxisSep.y = 1,
                          DotSize = 2,
                          Title = NULL, xTitle = X, yTitle = Y) {
  
  plot <- ggplot(data = df, aes(x = eval(parse(text = X)), y = eval(parse(text = Y)))) +
    geom_point(aes(color = N_delta), alpha = 0.75,
               na.rm = T,
               size = DotSize) +
    
    ggtitle(Title) +
    xlab(xTitle) +
    ylab(yTitle) +
    coord_cartesian(xlim = c(xmin,xmax),ylim = c(ymin,ymax))+
    scale_x_continuous(breaks=seq(xmin,xmax, AxisSep.x)) +
    scale_y_continuous(breaks=seq(ymin,ymax, AxisSep.y)) +
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
          panel.grid=element_blank(),
          aspect.ratio=1)
  
  
  return(plot)
}




Vulcano_plot <- MyVulcanoNvalues()


png("Vulcano.png", width = 500, height = 500, pointsize = 25)
grid.arrange(Vulcano_plot,
             nrow = 1)
dev.off()

pdf("Vulcano.pdf", width = 5, height = 5, useDingbats = F)
grid.arrange(Vulcano_plot,
             nrow = 1)
dev.off()
rm(Vulcano_plot)

#### Figure paper reproducibility ####
pdf(file = "Delta_Pairs.pdf", height = 5, width = 5, useDingbats=FALSE)
ggpairs(working_table_PTM[c("delta_Forward_II",
                            "delta_Forward_III",
                            "delta_Forward_IV",
                            "delta_Reverse_II",
                            "delta_Reverse_III",
                            "delta_Reverse_IV")], 
        mapping = ggplot2::aes(alpha = 0.2),
        lower = list(continuous = wrap("points", size=0.1)),
        diag = list(continuous = "blankDiag")) +
  
  coord_cartesian(xlim = c(-5,5),ylim = c(-5,5))
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



#### Pep x pSites cummulative Paper ####
foo_PTM_high <- subset(working_table_PTM, !is.na(delta_Mean) & Ratio.H.L.Mean_PG >= DeltaAnalysisThresholdEfficiency)
foo_PEP_high <- subset(working_table_PEP, Proteins_PEP %in% foo_PTM_high$Protein_PTM &
                         !is.na(delta_Mean & 
                                  Ratio.H.L.Mean_PG >= DeltaAnalysisThresholdEfficiency))

foo_PTM_low <- subset(working_table_PTM, !is.na(delta_Mean) & Ratio.H.L.Mean_PG <= DeltaAnalysisThresholdEfficiency)
foo_PEP_low <- subset(working_table_PEP, Proteins_PEP %in% foo_PTM_low$Protein_PTM &
                         !is.na(delta_Mean & 
                                  Ratio.H.L.Mean_PG <= DeltaAnalysisThresholdEfficiency))


CumPlot_Cl_Mean_PTM <- ggplot(foo_PEP_high,
                              aes(x = delta_Mean)) +
  #geom_vline(xintercept = 0) +
  geom_hline(yintercept = 0.5, linetype = 2) +
  
  coord_cartesian(xlim = c(-3, 3), ylim = c(0, 1)) +
  scale_x_continuous(breaks=seq(-3, 3, 1)) +
  scale_y_continuous(breaks=seq(0, 1, 0.25)) +
  ylab(NULL) +
  xlab(NULL) +
  labs(title = NULL) +
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
    aspect.ratio = 1) +
  
  ## Above threshold
  geom_step(data = foo_PEP_high,            
            stat="ecdf",                         # Peptides
            size = 3, 
            alpha= 1, 
            color = rgb(.5,.5,.5),
            na.rm = T) +
  
  geom_step(data = foo_PTM_high,            # PTM
            stat="ecdf",
            size = 3, 
            alpha= 1,
            color = rgb(.9,.45,0,.8),
            na.rm = T) +
  
  annotate("text",
           x = 2,
           y = 0.5, 
           label = paste0("High efficiency\n","p <= ",
                          signif(wilcox.test(foo_PEP_high$delta_Mean,
                                             foo_PTM_high$delta_Mean)$p.value, 3)),
           colour = rgb(.9,.45,0,.8),
           size = 3) +
  
  
  
  ## Below threshold
  geom_step(data = foo_PEP_low,            
            stat="ecdf",                         # Peptides
            size = 3, 
            alpha= 1, 
            color = rgb(0,0,0),
            na.rm = T) +
  
  geom_step(data = foo_PTM_low,            # PTM
            stat="ecdf",
            size = 3, 
            alpha= 1,
            color = rgb(.6,.3,0,.8),
            na.rm = T) +
  
  annotate("text",
           x = 2,
           y = 0.4, 
           label = paste0("Low efficiency\n","p <= ",
                          signif(wilcox.test(foo_PEP_low$delta_Mean,
                                             foo_PTM_low$delta_Mean)$p.value, 3)),
           colour = rgb(.9,.45,0,.8),
           size = 3)


pdf("PTM_CumPlot_Mean.pdf", width = 5, height = 5, useDingbats = F)
grid.arrange(CumPlot_Cl_Mean_PTM,
             nrow = 1)
dev.off()
rm(CumPlot_Cl_Mean_PTM, foo_PEP_high, foo_PEP_low, foo_PTM_high, foo_PTM_low)













#### Unmodified peptides ####
Maybe_interesting_pSites <- c(working_table_PTM[working_table_PTM$Gene.names_PTM %in% c(#"SERBP1",
  "RBM20",
  #"LARP1",
  "ELAVL1",
  "UPF1",
  "SF3B1"),
  "Site.ID_PTM"])


foo <- subset(working_table_PTM, !(PTM_increase_Efficiency_Cl_Mean | PTM_decrease_Efficiency_Cl_Mean))
foo1 <- subset(working_table_PTM, PTM_increase_Efficiency_Cl_Mean | PTM_decrease_Efficiency_Cl_Mean)
foo2 <- subset(working_table_PTM, Site.ID_PTM %in% Maybe_interesting_pSites)


# against unmodified peptides
# Delta efficiency
Unmod_delta_plot <- MyScatterPlot(foo,
                                  X = "delta_Mean", Y = "delta_unmod_Mean",
                                  xmin = -5, xmax = 5,
                                  ymin = -5, ymax = 5,
                                  GradientColorMin = rgb(.9,.45,0,.8),
                                  GradientColorMax = rgb(.3,.15,0,.8),
                                  DotSize = 3,
                                  xTitle = "delta_Mean", yTitle = "delta_unmod_Mean") +
  
  theme(panel.grid.major = element_line(colour = rgb(0,0,0,.1), linetype = "dashed")) +
  
  
  geom_point(data = foo1,
             na.rm = T,
             size = 3,
             colour = rgb(0,.4,.8))


# geom_point(data = foo2,
#            na.rm = T,
#            size = 3,
#            stroke = 1,
#            shape = 21,
#            colour = rgb(0,0,0)) +
#   
#   
#   geom_text_repel(data = foo2,
#                   aes(label = Site.ID_PTM),
#                   size = 4,
#                   colour = rgb(0,0,0),
#                   na.rm = T)



# Xtree Delta efficiency
Unmod_xtree_plot <- MyScatterPlot(foo,
                                  X = "delta_unmod_Mean", Y = "Ratio.H.L.Mean_PG",
                                  xmin = -5, xmax = 5,
                                  ymin = -5, ymax = 6, AxisSep.y = 1, InPercentage.y = T,
                                  GradientColorMin = rgb(.9,.45,0,.8),
                                  GradientColorMax = rgb(.3,.15,0,.8),
                                  DotSize = 3,
                                  xTitle = "delta_unmod_Mean", yTitle = "Ratio.H.L.Mean_PG") +
  
  theme(panel.grid.major = element_line(colour = rgb(0,0,0,.1), linetype = "dashed")) +
  
  
  geom_point(data = foo1,
             na.rm = T,
             size = 3,
             colour = rgb(0,.4,.8)) +
  
  
  geom_point(data = foo2,
             na.rm = T,
             size = 3,
             stroke = 1,
             shape = 21,
             colour = rgb(0,0,0)) +
  
  
  geom_text_repel(data = foo2,
                  aes(label = Site.ID_PTM),
                  size = 4,
                  colour = rgb(0,0,0),
                  na.rm = T)


sum(!is.na(working_table_PTM$Ratio.H.L.Mean_PTM))
sum(!is.na(working_table_PTM$Ratio.H.L.unmod..pep..Mean_PTM))


pdf("UnmodPeptides_qRIC.pdf", width = 10, height = 5, useDingbats = F)
grid.arrange(Unmod_delta_plot, Unmod_xtree_plot,
             nrow = 1)
dev.off()
rm(Unmod_delta_plot, Unmod_xtree_plot)





#### Figure for paper Unique and scores ####
all <- melt(working_table_PEP,
            id.vars = c("Sequence_PEP"), measure.vars = "delta_Mean")
all$Factor <- as.factor("All")

Unique <- melt(subset(working_table_PEP, Unique..Proteins._PEP & !is.na(delta_Mean)),
               id.vars = c("Sequence_PEP"), measure.vars = "delta_Mean")
Unique$Factor <- as.factor("Unique")

Non_Unique <- melt(subset(working_table_PEP, !Unique..Proteins._PEP & !is.na(delta_Mean)),
                   id.vars = c("Sequence_PEP"), measure.vars = "delta_Mean")
Non_Unique$Factor <- as.factor("Non-unique")

Low <- melt(subset(working_table_PEP, Score_PEP <= 100 & !is.na(delta_Mean)),
            id.vars = c("Sequence_PEP"), measure.vars = "delta_Mean")
Low$Factor <- as.factor("Low")


Medium_low <- melt(subset(working_table_PEP, Score_PEP > 100 & Score_PEP < 125 & !is.na(delta_Mean)),
                   id.vars = c("Sequence_PEP"), measure.vars = "delta_Mean")
Medium_low$Factor <- as.factor("Mid-Low")


Medium <- melt(subset(working_table_PEP, Score_PEP > 125 & Score_PEP < 155 & !is.na(delta_Mean)),
               id.vars = c("Sequence_PEP"), measure.vars = "delta_Mean")
Medium$Factor <- as.factor("Medium")


Medium_high <- melt(subset(working_table_PEP, Score_PEP > 155 & Score_PEP < 200 & !is.na(delta_Mean)),
                    id.vars = c("Sequence_PEP"), measure.vars = "delta_Mean")
Medium_high$Factor <- as.factor("Mid-High")


High <- melt(subset(working_table_PEP, Score_PEP >= 200 & !is.na(delta_Mean)),
             id.vars = c("Sequence_PEP"), measure.vars = "delta_Mean")
High$Factor <- as.factor("High")

pSites <- melt(subset(working_table_PTM, !is.na(delta_Mean)),
               id.vars = c("Site.ID_PTM"), measure.vars = "delta_Mean")
pSites$Factor <- as.factor("p-sites")



df <- rbind(#all[,2:4],
  Unique[,2:4],
  Non_Unique[,2:4],
  Low[,2:4],
  Medium_low[,2:4],
  Medium[,2:4],
  Medium_high[,2:4],
  High[,2:4],
  pSites[,2:4])
rm(all, Unique, Non_Unique, Low, Medium_low, Medium, Medium_high, High, pSites)


n_fun <- function(x){
  return(data.frame(y = 2,
                    label = length(x)))
}



bplot <- ggplot(data = df, mapping = aes(x = Factor, y = value)) +
  
  geom_hline(yintercept = 0, linetype = 1, color = "grey") +
  
  stat_summary(fun.data = n_fun, geom = "text", hjust = 0.5, na.rm = T) +
  
  geom_violin(aes(fill = Factor), color = rgb(0,0,0,0), na.rm = T, show.legend = F) +
  
  ylab(NULL) +
  xlab(NULL) +
  coord_cartesian(ylim = c(-2, 2)) +
  scale_y_continuous(breaks = seq(-2, 2, 1)) +
  scale_fill_manual(values = c(#rgb(0,0,0),
    brewer.pal(n = 8, name = "Dark2")[c(1,4)],
    brewer.pal(n = 8, name = "Greys")[4:8],
    rgb(.9,.45,0))) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        plot.title = element_text(lineheight=.8,
                                  vjust=0.5,
                                  hjust = 0.5,
                                  size=30),
        axis.text.x  = element_text(color = "black",
                                    angle=0, 
                                    hjust=.5,
                                    size=20),
        axis.title.x = element_text(size=25,
                                    hjust = 0.5,
                                    vjust = 1.5),
        axis.title.y = element_blank(),
        axis.text.y  = element_text(color = "black",
                                    angle=0, 
                                    vjust=.5,
                                    size=20),
        aspect.ratio = 1)




library(ggridges)

Ridgelines_plot <- ggplot(df, aes(x = value, y = Factor, fill = Factor)) +
  geom_density_ridges() +
  theme_ridges() + 
  theme(legend.position = "none") +
  scale_fill_manual(values = c(#rgb(0,0,0),
    brewer.pal(n = 8, name = "Dark2")[c(1,4)],
    brewer.pal(n = 8, name = "Greys")[4:8],
    rgb(.9,.45,0))) +
  ylab(NULL) +
  xlab("Delta efficiency (log2)") +
  coord_cartesian(xlim = c(-3, 3)) +
  scale_x_continuous(breaks = seq(-3, 3, 1))



pdf("Boxplot_PEP_PTM.pdf", width = 10, height = 5, useDingbats = F)
grid.arrange(bplot, Ridgelines_plot,
             nrow = 1)
dev.off()

rm(bplot, Ridgelines_plot, df)


#### Selected proteins and correlation ####

fit1 <- lm(delta_Forward ~ delta_Reverse, data = working_table_PTM)



foo <- subset(working_table_PTM, PTM_increase_Efficiency_Cl_Mean | PTM_decrease_Efficiency_Cl_Mean)
foo1 <- subset(foo, Gene.names_PTM %in% c("RBM20"))

#PTM
Mean_Cl_PTM_SILAC <- MyScatterPlot(working_table_PTM,
                                   "delta_Forward", "delta_Reverse",
                                   xmin = -7, xmax = 7, AxisSep.x = 2,
                                   ymin = -7, ymax = 7, AxisSep.y = 2, 
                                   GradientColorMin = rgb(.9,.45,0,.8),
                                   GradientColorMax = rgb(.3,.15,0,.8)) +
  
#  geom_smooth(data = working_table_PTM,
#              aes(x = delta_Forward,
#                  y = delta_Reverse),
#              method='lm',formula=y~x, na.rm = T) +
  
#  labs(title = paste("R = ",signif(sqrt(summary(fit1)$adj.r.squared), 3),
#                     "Adj R2 = ",signif(summary(fit1)$adj.r.squared, 3),
#                     "P =",signif(summary(fit1)$coef[2,4], 3))) +
  
  geom_point(data = foo,
             na.rm = T,
             size = 3,
             colour = rgb(0,.4,.8)) +
  
  
  geom_text_repel(data = foo1,
                  aes(label = Site.ID_PTM),
                  size = 6,
                  colour = rgb(0,0,1),
                  na.rm = T)


pdf("ForwardxReverse_Delta.pdf", width = 5, height = 5, useDingbats = F)
grid.arrange(Mean_Cl_PTM_SILAC, nrow = 1)
dev.off()
rm(foo, Mean_Cl_PTM_SILAC)
