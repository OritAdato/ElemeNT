library(dplyr)
library(reshape)
library(ggplot2)
library(forcats)
library(utils)
library(reader)

# This script creates graph that summarizes for every specie for every motif the precent of motifs found in specific 
# position around the TSS and the mean score of all motifs in the specific postion
#

main_path = "C:\\Main\\Path"
output_sub_path = "nascentPlots"
output_dir =  paste(main_path,output_sub_path,sep="")
setwd(output_dir)

sub_path = "4nascentGraph\\"
lib_pref <- paste(main_path,sub_path,sep="")


plot_sufix <-".png"
#plot_sufix <-".pdf"

###############################
plotgraph <- function(percent_fn, mean_fn, spname, f_type, motif, high_color, low_color ,plot_sufix, g_width, g_height)
{ 
  if (motif == "GAGA") {
    break_l = c(5, 10, 15)
  } else {
    break_l = c(0, 10, 20, 30, 40, 50, 60, 70, 80)
  }
  if (motif == "DPE") {
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
  
  #prints all plots per devStage
  ggplot(df, aes(x = pos, y = percent, fill = MeanScore)) + 
    geom_bar(stat = "identity", width = 1) + 
    scale_fill_gradient(high = high_color, low = low_color ) + #consider fixed 0 to 1 scale: limits=c(0, 1))
    #facet_grid defines the title of the y right axis
    #facet_grid(rows = vars(devStage)) +  #switch="both" to change facet position
    # in scale_y_continuous the expand = c(0, 0), limits = c(0, NA) makes sure that the y axis starts at 0
    #scale_y_continuous(breaks=c(0, 10, 20, 30, 40, 50, 60, 70, 80), expand = c(0, 0), limits = c(0, NA)) + #manually specify y axis ticks: breaks=c(10, 20, 30, 40, 50, 60)
    scale_y_continuous(breaks=break_l, expand = c(0, 0), limits = c(0, NA)) +
    xlab(x_label) +  ylab("% of transcripts containing the element") +
    theme_classic()   +
    theme(axis.ticks.x = element_blank(),
          axis.title.x = element_text(vjust = -0.9,size = 33),
          axis.title.y = element_text(margin = unit(c(0, 7, 0, 0), "mm"),size = 33),
          text = element_text(size=26),
          legend.text = element_text(size = 15 ),
          strip.text.y = element_text(size = 26, hjust = 0.8),
          plot.title = element_text(hjust = 0.5, size = 60)) +
   # ggtitle(current_motif_title) #specify name
    
    ggtitle(motif) #specify name
    theme(plot.title = element_text(size = 80, face = "bold"))
  sfn = paste(motif,"_",spname,plot_sufix,sep="")
  ggsave(sfn,plot=last_plot(),height = g_height, width = g_width,dpi = 300) 
}




######################
#read percent and mean score
specie_l <- list("dm6")
# for separate stages
f_type <- 'nascentRNA'
g_height = 15
g_width =10
# for all stages together
f_type <- 'nascentRNAallDS'
g_height = 10
g_width =14

motif <- 'dInr'
high_color <- "#073801"
low_color <- "#bbf2bb"

for (spname in specie_l){
  ##############################################
  #   dInr
  ##############################################
  motif <- 'dInr'
  high_color <- "#073801"
  low_color <- "#bbf2bb"
  plotgraph(percent_fn, mean_fn, spname, f_type, motif, high_color, low_color, plot_sufix, g_width,  g_height)
  
  ##############################################
  #   DPE
  ##############################################
  
  motif <- 'DPE'
  high_color <- "#2d0357"
  low_color <- "#d8c5eb"
  plotgraph(percent_fn, mean_fn, spname, f_type, motif, high_color, low_color, plot_sufix, g_width,  g_height)
  ##########################################################
  #        TATA
  ##########################################################
  motif <- 'TATA'
  high_color <- "#132B43"
  low_color <- "#bed5ed"
  plotgraph(percent_fn, mean_fn, spname, f_type, motif, high_color, low_color, plot_sufix, g_width,  g_height)
  
  ##########################################################
  #        PB
  ##########################################################
  motif <- 'PB'
  high_color <- "#DC0085"
  low_color <- #DB72C4"
  plotgraph(percent_fn, mean_fn, spname, f_type, motif, high_color, low_color, plot_sufix, g_width,  g_height)
  
  ##########################################################
  #        dTCT
  ##########################################################
  motif <- 'dTCT'
  high_color <- "#753304"
  low_color <- "#edbe9d"
  plotgraph(percent_fn, mean_fn, spname, f_type, motif, high_color, low_color, plot_sufix, g_width,  g_height)
 
  ##########################################################
  #        Motif1
  ##########################################################
  motif <- 'Motif1'
  high_color <-  "#31a354"
  low_color <- "#e5f5e0"
  plotgraph(percent_fn, mean_fn, spname, f_type, motif, high_color, low_color, plot_sufix, g_width,  g_height)
  ##########################################################
  #        GAGA
  ##########################################################
  motif <- 'GAGA'
  high_color <-  "#3182bd"
  low_color <-  "#deebf7"
  #high_color <- "#000000"
  #low_color <- "#e8e8e8"
  plotgraph(percent_fn, mean_fn, spname, f_type, motif, high_color, low_color, plot_sufix, g_width,  g_height)
  
  
  
}
  
