#kubra narci kbrnrc@gmail.com

#read the file containing ensemble gene ids of housekeeping genes and the count matrix for Huh7 or Mahlavu treatments

#usage
#common_dispersion.r housekeeping.list.txt countmatrix.txt
#housekeeping.list.txt includes the list of house keeping genes with ENSG id
#countmatrix is a tab limited file including expression values for RNA-seq results for example from the same cell line for different tratments. Row names should be for ENSG. 
args = commandArgs(trailingOnly=TRUE)
ensembl_gene_id = read.table(file=args[1], sep="\t", header=T, stringsAsFactors = F)
count_matrix = read.table(file=args[2], sep="\t", header=T, stringsAsFactors = F)
#merge ensemble gene ids to the count matrix
merged=merge(ensembl_gene_id, count_matrix, by="ENSEMBL_GENE_ID", all=FALSE)
#install edgeR
source("http://bioconductor.org/biocLite.R")
biocLite("edgeR")
a #update all packages
y #yes
library(edgeR)
#create count table
gene <- merged[1]
rownames(merged) <- gene[,1]
merged <- merged[,-1]
number=ncol(merged)
mobDataGroups <- c( as.integer(runif(number, min = 1, max = (3))) )
#specific to edgeR keep data in DGElist format
d <- DGEList(counts=merged ,group=factor(mobDataGroups))
d <- calcNormFactors(d)
design <- model.matrix(~factor(mobDataGroups))
d <- estimateDisp(d,design)
#print BVC value
d$common.dispersion
