setwd("K:/Results/20210903_Phospho_qRIC_revisionPaper/PFAM_RBDs")

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


#### Load proteins ####
PG_summ <- read.csv("../output/PG_summary.txt", sep = "\t", stringsAsFactors = FALSE)


#PG_summ <- subset(PG_summ, (Ratio.H.L.Forward > 0 & Ratio.H.L.Reverse > 0))

PG_summ <- subset(PG_summ,  #Known_mRBP &
                    !is.na(Ratio.H.L.Mean),
                  select = c("Majority.protein.IDs","Gene.names", "Sequence.length", "Ratio.H.L.Mean", "iBAQ.Input.Mean"))



#### Matched Pfam features ####

#Load PFAM Human proteins information
Pfam <- read.table("K:/Datasets/PFAM_human_domains.tsv", sep = "\t", stringsAsFactors = FALSE)
names(Pfam) = c("seq_id", "alignment_start", "alignment_end", "envelope_start", "envelope_end", "hmm_acc", "hmm_name", "type", "hmm_start", "hmm_end", "hmm_length" ,"bit_score", "E-value", "clan")

Pfam <- subset(Pfam, type == "Domain")



#### Label known RBDs

### Pfam features known to be RBD in Gerstberger
# Load Gerstberger 2014
gerstberger <- read.csv("K:/Datasets/Gerstberger_RBP_Table_2014.txt", sep = "\t", stringsAsFactors = FALSE)

##Generate list of all known RBDs from Gerstberger et al 2014
RBD.Gerstberger = gerstberger["domains.count."]

#split domain column using ","
RBD.Gerstberger <- data.frame(do.call('rbind', strsplit(as.character(RBD.Gerstberger$domains.count.),',', fixed=TRUE)), stringsAsFactors = FALSE)

#replace brackets with anything inside for "" (nothing) in all collumns, apply to all collumns and add it up as a data frame
RBD.Gerstberger = data.frame(sapply(RBD.Gerstberger, FUN =  gsub, pattern = "\\[.\\]", replacement = "", "[", 1), stringsAsFactors = FALSE)

#append all collumns from RBD.Gerstberger, keeping only the unique names
RBD.Gerstberger = sort(unique(append(RBD.Gerstberger$X1, c(RBD.Gerstberger$X2, RBD.Gerstberger$X3, RBD.Gerstberger$X4, RBD.Gerstberger$X5, RBD.Gerstberger$X6, RBD.Gerstberger$X7, RBD.Gerstberger$X8, RBD.Gerstberger$X9))))

#Remove empty RBD that was originally named "."
RBD.Gerstberger <- subset(RBD.Gerstberger, RBD.Gerstberger != ".")

#Mark Pfam features known to be RBDs based on Gerstberger
Pfam$RBD.logical <- ifelse(Pfam$hmm_name %in% RBD.Gerstberger, T, F) 



rm(gerstberger)



#Filter Pfam to contain only relevant proteins from PG_summ
Pfam <- subset(Pfam, seq_id %in% PG_summ$Majority.protein.IDs)

#Count number of domains and number of RBDs as defined in gerstberger
PG_summ[,"Number_Domains"] <- apply(PG_summ["Majority.protein.IDs"], 1, function(x) nrow(subset(Pfam, seq_id == x)))
PG_summ[,"Number_RBDs"] <- apply(PG_summ["Majority.protein.IDs"], 1, function(x) nrow(subset(Pfam, seq_id == x & RBD.logical)))




lvls = c("0", "1", "2", "3 or more")





PG_summ[, "Number_Domains_factors"] <- sapply(PG_summ[,"Number_Domains"],
                                              function(x)  factor(ifelse(x == 1, "1",
                                                                         ifelse(x == 2, "2",
                                                                                ifelse(x >= 3 , "3 or more", "0"))),
                                                                  levels = lvls))

PG_summ[, "Number_RBDs_factors"] <- sapply(PG_summ[,"Number_RBDs"],
                                           function(x)  factor(ifelse(x == 1, "1",
                                                                      ifelse(x == 2, "2",
                                                                             ifelse(x >= 3 , "3 or more", "0"))),
                                                               levels = lvls))





#### Cumulative Plots ####
MyCumPlot <- function(xmin = -5, xmax = 6, AxisSep.x = 1,
                      ymin = 0, ymax = 1, AxisSep.y = .2,
                      InPercentage.x = T, sec_breaks.x = c(0.1,0.5,1,2,4,8,16,32,64),
                      scaled_breaks.x = log2(sec_breaks.x),
                      Title = NULL, xlabel = NULL, ylabel = NULL) {
  
  plot <- ggplot() +
    
    geom_hline(yintercept = .5, linetype = 2, color = rgb(.25,.25,.25)) + 
    
    coord_cartesian(xlim = c(xmin, xmax)) +
    ylab(ylabel) +
    xlab(xlabel) +
    labs(title = Title) +
    scale_y_continuous(breaks=seq(ymin,ymax, AxisSep.y)) +
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
  
  if (InPercentage.x) {
    
    plot <- plot +
      scale_x_continuous(breaks = scaled_breaks.x, labels = sec_breaks.x)
    
  } else {
    
    plot <- plot +
      scale_x_continuous(breaks=seq(xmin,xmax, AxisSep.x))
    
  }
  
  
  return(plot)
}





#Number_Domains_factors
wiltest_Number_Domains_factors <- wilcox.test(subset(PG_summ, Number_Domains_factors == lvls[1])[,"Ratio.H.L.Mean"],
                                              subset(PG_summ, Number_Domains_factors == lvls[length(lvls)])[,"Ratio.H.L.Mean"])
plot_Cumu_Number_Domains <- MyCumPlot() +
  
  geom_step(data = PG_summ,
            aes(x = Ratio.H.L.Mean, color = Number_Domains_factors), stat="ecdf", size = 3, na.rm = T, show.legend = T) +
  
  labs(title = paste0("HighXLow, p-value <= ", signif(wiltest_Number_Domains_factors$p.value, digits = 3)))






#Number_RBDs_factors
wiltest_Number_RBDs_factors <- wilcox.test(subset(PG_summ, Number_RBDs_factors == lvls[1])[,"Ratio.H.L.Mean"],
                                           subset(PG_summ, Number_RBDs_factors == lvls[4])[,"Ratio.H.L.Mean"])
plot_Cumu_Number_RBDs <- MyCumPlot() +
  
  geom_step(data = PG_summ,
            aes(x = Ratio.H.L.Mean, color = Number_RBDs_factors), stat="ecdf", size = 3, na.rm = T, show.legend = F) +
  
  labs(title = paste0("HighXLow, p-value <= ", signif(wiltest_Number_RBDs_factors$p.value, digits = 3)))






png("Pfam_RBP_domains.png", width = 1000, height = 500, pointsize = 25)
grid.arrange(plot_Cumu_Number_Domains, plot_Cumu_Number_RBDs,
             nrow = 2)
dev.off()






#### Figure paper 

#colors <- brewer.pal(n = 8, name = "Set2")
colors <- c("black", brewer.pal(n = 8, name = "PuBu")[c(4,6,8)])
  
plot_Cumu_Number_RBDs <- MyCumPlot(AxisSep.y = .25, xmax = 6) +
  
  scale_color_manual(values=colors[1:4]) +
  
  geom_step(data = PG_summ,
            aes(x = Ratio.H.L.Mean, color = Number_RBDs_factors), stat="ecdf", size = 3, na.rm = T, show.legend = T)



pdf("Pfam_RBP_domains.pdf", width = 5, height = 5, useDingbats = F)
grid.arrange(plot_Cumu_Number_RBDs,
             nrow = 1)
dev.off()



png("Pfam_RBP_domains.png", width = 500, height = 500)
grid.arrange(plot_Cumu_Number_RBDs,
             nrow = 1)
dev.off()



rm(plot_Cumu_Number_Domains, plot_Cumu_Number_RBDs)

rm(scaffold_plot, wiltest_Number_Domains_factors, wiltest_Number_RBDs_factors, lvls)




#### Types of RBDs ####
# Give collection of RBDs per gene
PG_summ$RBDs <- as.character(lapply(PG_summ[, "Majority.protein.IDs"], function(x) paste(unique(subset(Pfam, RBD.logical & seq_id == x)$hmm_name), collapse = ";")))

PG_summ[PG_summ[ ,"RBDs"] == "" ,"RBDs"] <- "None"



sum(grepl("RRM_1", PG_summ$RBDs))

# exclude repeated domains in same protein
Pfam_RBD_count <- unique(subset(Pfam, RBD.logical & seq_id %in% PG_summ$Majority.protein.IDs,
                                select = c("seq_id", "hmm_name")))

# Must abundant RBDs
abundant_RBDs <- as.data.frame(sort(table(Pfam_RBD_count$hmm_name), decreasing = T), stringsAsFactors = F)



top_abundant_domains <- head(abundant_RBDs$Var1, 4)
top_abundant_domains <- c("RRM_1", "KH_1", "DEAD")
PG_summ[!PG_summ[ ,"RBDs"] %in% c(top_abundant_domains, "None"),"RBDs"] <- "Other"

PG_summ$RBDs <- factor(PG_summ$RBDs, levels = c("RRM_1",
                                                     "KH_1",
                                                     "DEAD",
                                                     #"zf-CCHC",
                                                     "Other",
                                                     "None"))

# Colors for plotting
gg_color_hue <- brewer.pal(n = 8, name = "Set2") #hcl(h = seq(15, 375, length = 10 + 1), l = 65, c = 100)[1:10]



MyCumPlotPaper <- function(df = PG_summ,
                      xmin = -5, xmax = 6, AxisSep.x = 1,
                      ymin = 0, ymax = 1, AxisSep.y = .2,
                      InPercentage.x = T, sec_breaks.x = c(0.1,0.5,1,2,4,8,16,32,64),
                      scaled_breaks.x = log2(sec_breaks.x),
                      Title = NULL, xlabel = NULL, ylabel = NULL) {
  
  plot <- ggplot() +
    
    geom_hline(yintercept = .5, linetype = 2, color = rgb(.25,.25,.25)) + 
    
    geom_step(data = df,
              aes(x = Ratio.H.L.Mean,
                  colour = RBDs),
              stat="ecdf", size = 3, na.rm = T) +
    
    scale_color_manual(values = gg_color_hue) +
    
    coord_cartesian(xlim = c(xmin, xmax)) +
    ylab(ylabel) +
    xlab(xlabel) +
    labs(title = Title) +
    scale_y_continuous(breaks=seq(ymin,ymax, AxisSep.y)) +
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
  
  if (InPercentage.x) {
    
    plot <- plot +
      scale_x_continuous(breaks = scaled_breaks.x, labels = sec_breaks.x)
    
  } else {
    
    plot <- plot +
      scale_x_continuous(breaks=seq(xmin,xmax, AxisSep.x))
    
  }
  
  
  return(plot)
}


plot_Cummu_RBDs <- MyCumPlotPaper(df = PG_summ, Title = "At least 1 RBD")


##  Select proteins with the given top_abundant_domain and without any other domains.
# I'm doing this for p-value calculation because otherwise subgroups won't be unique, as some proteins have more than 1 RBD

PG_summ_singleRBD <- subset(PG_summ, Number_RBDs == 1 | Number_RBDs == 0)



foo <- subset(Pfam_RBD_count, hmm_name == top_abundant_domains[1])

foo1 <- subset(PG_summ_singleRBD, Majority.protein.IDs %in% foo$seq_id)


foo <- subset(Pfam_RBD_count, hmm_name == top_abundant_domains[2])

foo2 <- subset(PG_summ_singleRBD, Majority.protein.IDs %in% foo$seq_id)


foo <- subset(Pfam_RBD_count, hmm_name == top_abundant_domains[3])

foo3 <- subset(PG_summ_singleRBD, Majority.protein.IDs %in% foo$seq_id)


# foo <- subset(Pfam_RBD_count, hmm_name == top_abundant_domains[4])
# 
# foo4 <- subset(PG_summ_singleRBD, Majority.protein.IDs %in% foo$seq_id)


foo5 <- subset(PG_summ_singleRBD, RBDs == "Other")


foo6 <- subset(PG_summ_singleRBD, RBDs == "None")




plot_Cummu_RBDs_pvalues <- MyCumPlotPaper(df = PG_summ_singleRBD, Title = "Single RBD") +
  
  ## RRM_1
  annotate("text", x = 3, y = .7, hjust = 0,
           label = paste0(top_abundant_domains[1], " X ", top_abundant_domains[2], " = ", 
                          signif(wilcox.test(foo1$Ratio.H.L.Mean,
                                             foo2$Ratio.H.L.Mean)$p.value,
                                 digits = 2)),
           size = 5) +
  
  annotate("text", x = 3, y = .65, hjust = 0,
           label = paste0(top_abundant_domains[1], " X ", top_abundant_domains[3], " = ", 
                          signif(wilcox.test(foo1$Ratio.H.L.Mean,
                                             foo3$Ratio.H.L.Mean)$p.value,
                                 digits = 2)),
           size = 5) +
  
  # annotate("text", x = 3, y = .6, hjust = 0,
  #          label = paste0(top_abundant_domains[1], " X ", top_abundant_domains[4], " = ", 
  #                         signif(wilcox.test(foo1$Ratio.H.L.Mean,
  #                                            foo4$Ratio.H.L.Mean)$p.value,
  #                                digits = 2)),
  #          size = 5) +
  
  annotate("text", x = 3, y = .55, hjust = 0,
           label = paste0(top_abundant_domains[1], " X ", "Other", " = ", 
                          signif(wilcox.test(foo1$Ratio.H.L.Mean,
                                             foo5$Ratio.H.L.Mean)$p.value,
                                 digits = 2)),
           size = 5) +

  annotate("text", x = 3, y = .5, hjust = 0,
           label = paste0(top_abundant_domains[1], " X ", "None", " = ",
                          signif(wilcox.test(foo1$Ratio.H.L.Mean,
                                             foo6$Ratio.H.L.Mean)$p.value,
                                 digits = 2)),
           size = 5) +
  
  
  ## KH_1
  annotate("text", x = 3, y = .45, hjust = 0,
           label = paste0(top_abundant_domains[2], " X ", top_abundant_domains[3], " = ", 
                          signif(wilcox.test(foo2$Ratio.H.L.Mean,
                                             foo3$Ratio.H.L.Mean)$p.value,
                                 digits = 2)),
           size = 5) +
  
  # annotate("text", x = 3, y = .4, hjust = 0,
  #          label = paste0(top_abundant_domains[2], " X ", top_abundant_domains[4], " = ", 
  #                         signif(wilcox.test(foo2$Ratio.H.L.Mean,
  #                                            foo4$Ratio.H.L.Mean)$p.value,
  #                                digits = 2)),
  #          size = 5) +
  
  annotate("text", x = 3, y = .35, hjust = 0,
           label = paste0(top_abundant_domains[2], " X ", "Other", " = ", 
                          signif(wilcox.test(foo2$Ratio.H.L.Mean,
                                             foo5$Ratio.H.L.Mean)$p.value,
                                 digits = 2)),
           size = 5) +

  annotate("text", x = 3, y = .3, hjust = 0,
           label = paste0(top_abundant_domains[2], " X ", "None", " = ",
                          signif(wilcox.test(foo2$Ratio.H.L.Mean,
                                             foo6$Ratio.H.L.Mean)$p.value,
                                 digits = 2)),
           size = 5) +
  
  
  ### Dead 
  # annotate("text", x = 3, y = .25, hjust = 0,
  #          label = paste0(top_abundant_domains[3], " X ", top_abundant_domains[4], " = ", 
  #                         signif(wilcox.test(foo3$Ratio.H.L.Mean,
  #                                            foo4$Ratio.H.L.Mean)$p.value,
  #                                digits = 2)),
  #          size = 5) +
  
  annotate("text", x = 3, y = .2, hjust = 0,
           label = paste0(top_abundant_domains[3], " X ", "Other", " = ", 
                          signif(wilcox.test(foo3$Ratio.H.L.Mean,
                                             foo5$Ratio.H.L.Mean)$p.value,
                                 digits = 2)),
           size = 5) +

  annotate("text", x = 3, y = .15, hjust = 0,
           label = paste0(top_abundant_domains[3], " X ", "None", " = ",
                          signif(wilcox.test(foo3$Ratio.H.L.Mean,
                                             foo6$Ratio.H.L.Mean)$p.value,
                                 digits = 2)),
           size = 5) +
  
  
  ## zf-CCHC
  # annotate("text", x = 3, y = .1, hjust = 0,
  #          label = paste0(top_abundant_domains[4], " X ", "Other", " = ", 
  #                         signif(wilcox.test(foo4$Ratio.H.L.Mean,
  #                                            foo5$Ratio.H.L.Mean)$p.value,
  #                                digits = 2)),
  #          size = 5) +
  # 
  # annotate("text", x = 3, y = .05, hjust = 0,
  #          label = paste0(top_abundant_domains[4], " X ", "None", " = ",
  #                         signif(wilcox.test(foo4$Ratio.H.L.Mean,
  #                                            foo6$Ratio.H.L.Mean)$p.value,
  #                                digits = 2)),
  #          size = 5) +

  ## Other
  annotate("text", x = 3, y = .0, hjust = 0,
           label = paste0("Other", " X ", "None", " = ",
                          signif(wilcox.test(foo5$Ratio.H.L.Mean,
                                             foo6$Ratio.H.L.Mean)$p.value,
                                 digits = 2)),
           size = 5)





png("Cumm_RBDs.png", width = 500, height = 1000, pointsize = 25)
grid.arrange(plot_Cummu_RBDs, plot_Cummu_RBDs_pvalues, nrow = 2)
dev.off()

pdf("Cumm_RBDs.pdf", width = 5, height = 10)
grid.arrange(plot_Cummu_RBDs, plot_Cummu_RBDs_pvalues, nrow = 2)
dev.off()

rm(plot_Cummu_RBDs_pvalues, plot_Cummu_RBDs, foo, foo1, foo2, foo3, foo4, foo5, foo6)



















#### With or without RBDs ####

PG_summ$Type.Domain <- ifelse((PG_summ$Number_Domains != 0) & (PG_summ$Number_RBDs != 0), "RBD",
                              ifelse((PG_summ$Number_Domains != 0), "Other",
                                     "None"))

plot_Cummu_all <- MyCumPlot(xmax = 5, AxisSep.y = .25) +
  
  geom_step(data = PG_summ,
            aes(x = Ratio.H.L.Mean, color = as.factor(Type.Domain)),
            stat="ecdf", size = 3, na.rm = T, show.legend = F) +
  
  scale_color_manual(values=c(rgb(.8,0,.8), rgb(0,.8,0), rgb(0,.8,.8))) +
  
  annotate("text", x = 0, y = 0, hjust = 0,
           label = paste0("p-value <= ", 
                          signif(wilcox.test(subset(PG_summ, Type.Domain == "None")$Ratio.H.L.Mean,
                                             subset(PG_summ, Type.Domain == "Other")$Ratio.H.L.Mean)$p.value,
                                 digits = 7)),
           #fontface = "bold",
           colour = rgb(0,0,0), size = 6) +
  
  annotate("text", x = 0, y = 0.2, hjust = 0,
           label = paste0("p-value <= ", 
                          signif(wilcox.test(subset(PG_summ, Type.Domain == "RBD")$Ratio.H.L.Mean,
                                             subset(PG_summ, Type.Domain == "Other")$Ratio.H.L.Mean)$p.value,
                                 digits = 7)),
           #fontface = "bold",
           colour = rgb(0,.8,0), size = 6) +
  
  annotate("text", x = 0, y = 0.1, hjust = 0,
           label = paste0("p-value <= ", 
                          signif(wilcox.test(subset(PG_summ, Type.Domain == "RBD")$Ratio.H.L.Mean,
                                             subset(PG_summ, Type.Domain == "None")$Ratio.H.L.Mean)$p.value,
                                 digits = 7)),
           #fontface = "bold", 
           colour = rgb(.8,0,.8), size = 6)




png("Pfam_CumPlot_RBDs.png", width = 500, height = 500, pointsize = 25)
grid.arrange(plot_Cummu_all,
             nrow = 2)
dev.off()





rm(Pfam, Pfam_RBD_count, PG_summ, abundant_RBDs, gg_color_hue, top_abundant_domains, RBD.Gerstberger)
