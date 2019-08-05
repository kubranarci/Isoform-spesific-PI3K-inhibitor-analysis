#step 1 prioritization per HCC line
setwd("~/Desktop ")
data.path1 <- paste0(getwd(), "/100Randomization/random_nodeattributes/")
randomization.files <- list.files(data.path1, ".tsv$", recursive = T, full.names = T)
data.path2 <- paste0(getwd(), "/13-2-centrality_measures_")
centrality.files <- list.files(data.path2, ".txt$", recursive = T, full.names = T)
merge_final =list()
for (random in randomization.files)  {
  random_temp <- read.delim(random, header = TRUE, stringsAsFactors = F, sep = "\t")
  random_name <- gsub("_randomTerminals_nodeattributes.tsv", "", basename(random))
  random_temp <- data.frame("Node"=random_temp$Protein, random_temp[-1], stringsAsFactors = F ) 
  for (central in centrality.files){
    central_temp <- read.delim(central, header = TRUE, stringsAsFactors = F, sep = "\t")
    central_name <- gsub("_centrality.txt", "", basename(central))
    central_temp <- data.frame("Node"=central_temp$Node.Name, "Set"=central_name, central_temp[-1], stringsAsFactors = F )   
     if (random_name == central_name){
      merged_set <- merge(random_temp, central_temp, by = "Node") 
      merged_set <- (merged_set)[merged_set$FractionOfOptimalForestsContaining == 0.01,] 
      merged_set <- (merged_set)[merged_set$BetweennessCentrality < 0.001,] 
      merged_set <- (merged_set)[merged_set$Degree.Centrality > 0.001,] 
      merged_set <- (merged_set)[merged_set$Betweenness.Centrality > 0.001,] 
      merged_set <- (merged_set)[merged_set$EigenVector.Centrality > 0.001,]
      if (nrow(merged_set) > 5) {merged_set <- merged_set[1:5,]} 
      merge_final <- rbind(merge_final, merged_set, stringsAsFactors = F)   }}}
for(i in 1:length(merge_final$TerminalType)) {
  if(merge_final$TerminalType[i] == "") {
    merge_final$TerminalType[i]  <- "Steiner"}
  else 
    merge_final$TerminalType[i]  <- "Transcriptomic"}} 
data <- data.frame("Inhibitor Type" = merge_final$Set, Node = merge_final$Node, Random_Betweenness = merge_final$BetweennessCentrality, Degree = merge_final$Degree.Centrality, Eigen = merge_final$EigenVector.Centrality, "Terminal Type"= merge_final$TerminalType, HCC= "Mahlavu")
MV_data <- merge_final
library(tidyr)   
library(tidyverse)
#perform step 1 for Huh7 also and merge them in data3
data3 <- rbind(Huh7_data, MV_data, stringsAsFactors = F)
data3$Node = with(data3, reorder(Node,Random_Betweenness))
#visualize the results
  ggplot(data = data3) + 
    geom_point(aes(Random_Betweenness,Node,color= Inhibitor.Type, shape= Terminal.Type), size=3) +
    facet_grid(.~HCC) +
    xlab("Betweenness Centralities of Random Networks") + ylab("Gene Nodes") +
    theme(legend.title= element_blank(), text = element_text(size=14), axis.text.x = element_text( size=10, angle=45))