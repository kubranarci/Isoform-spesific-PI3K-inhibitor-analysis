#kubra narci kbrnrc@gmail.com

#reads a list of .bgo results from BINGO results, merges, annotates and filtrates by evidence codes. 

#usage
#bingo_functions.r working_directory GO_terms.txt goa_human.gaf output.txt
#Go_terms.txt lists  GO ids and descriptions
#goa_human.gaf file includes evidence codes

#necessary packages
library(stringr)
library(GO.db)

args = commandArgs(trailingOnly=TRUE)

setwd(args[1])
data.path <- paste0(getwd(), "/")
count.files <- list.files(data.path, ".bgo$", recursive = T, full.names = T)
snames <- gsub(".bgo", "", basename(count.files))
merged = data.frame("GOID"="")
n <- ncol(merged)
for (countf in count.files) {
  temp <- read.delim(countf, header = TRUE, stringsAsFactors = F)
  name <- gsub(".bgo", "", basename(countf))
  temp <- temp[temp$corr.p.value < 0.05, ]
  if (nrow(temp) > 0){  
    temp=data.frame("GOID"= temp$GO.ID, temp$corr.p.value, stringsAsFactors = F )
    temp$GOID <- str_pad(temp$GOID, 7, pad = "0")
    temp$GOID <- paste0("GO:", temp$GOID)
    merged=merge(merged,temp, by="GOID", all=TRUE)
    n=n+1
    colnames(merged)[n] <-  paste("", name, sep="")    }
  else {
    temp=data.frame("GOID"= "NA", 0, stringsAsFactors = F )
    merged=merge(merged,temp, by="GOID", all=TRUE)
    n=n+1
    colnames(merged)[n] <-  paste("", name, sep="")}}
GO_all <- read.delim(args[2], header = TRUE, stringsAsFactors = F)
converted <- merge(GO_all, merged, by= "GOID", all=FALSE )
#for filtered cluster tables
deflist <- select(GO.db, keys= as.character(converted$GOID), columns=c("TERM","ONTOLOGY"), keytype="GOID")
GO_evidence <- read.csv(file= args[3], header = FALSE, sep = "\t", skip = 31)
GO_evidence_2 <- data.frame("GOID" = GO_evidence$V5, "Evidence" = GO_evidence$V7 )
merged_GO <- merge(converted, GO_evidence_2, by= "GOID", all= FALSE)
filtered_GO <- dplyr::filter(merged_GO, grepl('EXP|IDA|IPI|IMP|IGI|IEP', Evidence))
uniqu_GO <- filtered_GO[!duplicated(filtered_GO$GOID), ]
write.table(uniqu_GO, file= args[4], sep = "\t" )
