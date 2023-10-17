library(dplyr)
library(reshape)
library(ggplot2)
library(forcats)
library(utils)
library(reader)

# This script creates graph that summarizes for every specie for every motif the precent of motifs found in specific 
# position around the TSS and the mean score of all motifs in the specific postion
#

main_path = "C:\\Main\\Path\\"
plot_sufix <-".png"
plot_sufix <-".pdf"

output_sub_path = "corePlots"
output_dir =  paste(main_path,output_sub_path,sep="")
setwd(output_dir)

sub_path = "4CoreGraph\\"
lib_pref <- paste(main_path,sub_path,sep="")


###############################
plotgraph <- function(percent_fn, mean_fn, spname, f_type, motif, high_color, low_color,plot_sufix)
{ 
  if (motif == "GAGA") {
    break_l = c(5, 10, 15)
  } else {
    break_l = c(10, 30, 50)
  }
  if (substr(motif,1,3) == "DPE") {
    x_label <- "position relative to the A+1 of dInr"
  } else {
    x_label <- "position relative to TSS"
  }
  percent_fn = paste(lib_pref,spname,"_",f_type, "_",motif,"_percent.txt",sep="")
  percent <- read.table(percent_fn,header = TRUE)
    
  indOfUnderscore <-  unlist(gregexpr('_', colnames(percent)[1]))
  current_motif_title <- substr(colnames(percent)[1],1,indOfUnderscore-1)
  current_motif_title <- paste(current_motif_title,f_type)
  colnames(percent)[1] <- "percent"
  colnames(percent) <- gsub("X", "", colnames(percent))
  
  mean_fn = paste(lib_pref,spname,"_",f_type, "_",motif,"_meanScore.txt",sep="")
  Mean <- read.table(mean_fn,header=TRUE)
  colnames(Mean)[1] <-"MeanScore"
  colnames(Mean) <- gsub("X", "", colnames(Mean))
  #melting and renaming
  mper <- melt(percent, id=c("percent"))
  colnames(mper) <- c('pos','devStage','percent')
  mscr <- melt(Mean, id=c("MeanScore"))
  colnames(mscr) <- c('pos','devStage','MeanScore')
  
  #creating dataframe
  df <- merge(mper, mscr, by=c("pos","devStage"))
  df$devStage <- as.factor(df$devStage)
  df$pos <- as.factor(df$pos)
  #x_label <- "position relative to TSS"
  
  #prints all plots per devStage
  ggplot(df, aes(x = pos, y = percent, fill = MeanScore)) + 
    geom_bar(stat = "identity", width = 1) + 
    scale_fill_gradient(high = high_color, low = low_color ) + #consider fixed 0 to 1 scale: limits=c(0, 1))
    facet_grid(rows = vars(devStage)) +  #switch="both" to change facet position
    #scale_y_continuous(breaks=c(10, 30, 50)) + #manually specify y axis ticks: breaks=c(10, 30)
    scale_y_continuous(breaks=break_l) +
    xlab(x_label) +  ylab("% of transcripts containing the element") +
    theme_classic()   +
    theme(axis.ticks.x = element_blank(),
          axis.title.x = element_text(vjust = -0.9,size = 30),
          axis.title.y = element_text(margin = unit(c(0, 7, 0, 0), "mm"),size = 33),
          text = element_text(size=18),
          strip.text.y = element_text(size = 20),
          plot.title = element_text(hjust = 0.5, size = 40)) +
    ggtitle(motif) #specify name
  sfn = paste(motif,"_",spname,plot_sufix,sep="")
  ggsave(sfn,plot=last_plot(),height = 15, width = 11, dpi = 300) 
}




######################
#read percent and mean score
#
specie_l <- list("all")
f_type <- 'Core'
motif <- 'dInr'
#high_color <- "#073801"
#low_color <- "#bbf2bb"

for (spname in specie_l){
  ##############################################
  #   dInr
  ##############################################
  motif <- 'dInr'
  high_color <- "#073801"
  low_color <- "#bbf2bb"
  plotgraph(percent_fn, mean_fn, spname, f_type, motif, high_color, low_color,plot_sufix)
  
  ##############################################
  #   hInr
  ##############################################
  motif <- 'hInr'
  high_color <- "#073801"
  low_color <- "#bbf2bb"
  plotgraph(percent_fn, mean_fn, spname, f_type, motif, high_color, low_color,plot_sufix)
  
  ##############################################
  #   BBCABW
  ##############################################
  motif <- 'BBCABW'
  high_color <- "#073801"
  low_color <- "#bbf2bb"
  plotgraph(percent_fn, mean_fn, spname, f_type, motif, high_color, low_color,plot_sufix)
  
  
  ##############################################
  #   DPE
  ##############################################
  
  motif <- 'DPE'
  high_color <- "#2d0357"
  low_color <- "#d8c5eb"
  plotgraph(percent_fn, mean_fn, spname, f_type, motif, high_color, low_color,plot_sufix)
  
  ##############################################
  #   DPEhInr
  ##############################################
  
  motif <- 'DPEhInr'
  high_color <- "#2d0357"
  low_color <- "#d8c5eb"
  plotgraph(percent_fn, mean_fn, spname, f_type, motif, high_color, low_color,plot_sufix)
  
  ##############################################
  #   DPEdInr
  ##############################################
  
  motif <- 'DPEdInr'
  high_color <- "#2d0357"
  low_color <- "#d8c5eb"
  plotgraph(percent_fn, mean_fn, spname, f_type, motif, high_color, low_color,plot_sufix)
  
  ##############################################
  #   DPEBBCABW
  ##############################################
  
  motif <- 'DPEBBCABW'
  high_color <- "#2d0357"
  low_color <- "#d8c5eb"
  plotgraph(percent_fn, mean_fn, spname, f_type, motif, high_color, low_color,plot_sufix)
  
  
  
  ##########################################################
  #        TATA
  ##########################################################
  motif <- 'TATA'
  high_color <- "#132B43"
  low_color <- "#bed5ed"
  plotgraph(percent_fn, mean_fn, spname, f_type, motif, high_color, low_color,plot_sufix)
  
  ##########################################################
  #        PB
  ##########################################################
  motif <- 'PB'
  high_color <- "#DC0085"
  low_color <- #DB72C4"
  plotgraph(percent_fn, mean_fn, spname, f_type, motif, high_color, low_color,plot_sufix)
  
  ##########################################################
  #        GAGA
  ##########################################################
  motif <- 'GAGA'
  high_color <-  "#3182bd"
  low_color <-  "#deebf7"
  plotgraph(percent_fn, mean_fn, spname, f_type, motif, high_color, low_color,plot_sufix)
  
  ##########################################################
  #        Motif1
  ##########################################################
  motif <- 'Motif1'
  high_color <-  "#31a354"
  low_color <- "#e5f5e0"
  plotgraph(percent_fn, mean_fn, spname, f_type, motif, high_color, low_color,plot_sufix)
  
  ##########################################################
  #        dTCT
  ##########################################################
  motif <- 'dTCT'
  high_color <- "#753304"
  low_color <- "#edbe9d"
  plotgraph(percent_fn, mean_fn, spname, f_type, motif, high_color, low_color,plot_sufix)
  
}
   
