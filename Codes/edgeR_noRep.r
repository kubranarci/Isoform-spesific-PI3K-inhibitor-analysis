
#Kubra Narci kbrnrc@gmail.com
#Given a count taable including control vs sample expression values and filraation criterias for cpm, dispersion, FDR and P value Differentially Expressed Genes calculated and annotated and written into the output file.  

#usage
#edgeR_noRep.R CountTable.txt cpm dispersion  FDR_value_up FDR_value_down  P_value Output
#CountTable.txt will be in form of > 1. CONTROL 2. SAMPLE

args = commandArgs(trailingOnly=TRUE)
#install packages
source("http://bioconductor.org/biocLite.R")
biocLite("lattice")
biocLite("edgeR")
library(lattice)
library(edgeR)
biocLite("org.Hs.eg.db")
library(org.Hs.eg.db)
#create count table
countsTable=read.delim(file=args[1],header=TRUE, row.names = 1)
mobDataGroups <- c("Mt", "Mu")
#specific to edgeR keep data in DGElist format
d <- DGEList(counts=countsTable ,group=factor(mobDataGroups))
cpm= as.numeric(args[2])
#filering
d.full <- d # keep the old one in case we mess up
head(cpm(d))
apply(d$counts, 2, sum)
keep <- rowSums(cpm(d)>cpm) >= 2
d <- d[keep,]
dim(d)
d$samples$lib.size <- colSums(d$counts)
d$samples
#normalizing the data
d <- calcNormFactors(d)
#dispersion is handgiven
disp=as.numeric(args[3])
et <- exactTest(d, dispersion= disp)
pvaluesort <- topTags( et , n = nrow( et$table ) , sort.by = "p.value" )$table
pvaluesort $symbol <- mapIds(org.Hs.eg.db,
                             keys=row.names(pvaluesort ),
                             column="SYMBOL",
                             keytype="ENSEMBL",
                             multiVals="first")
pvaluesort $EntrezID <- mapIds(org.Hs.eg.db,
                              keys=row.names(pvaluesort ),
                               column="ENTREZID",
                               keytype="ENSEMBL",
                               multiVals="first")				 
pvaluesort $Unigene <- mapIds(org.Hs.eg.db,
                              keys=row.names(pvaluesort ),
                              column="UNIGENE",
                              keytype="ENSEMBL",
                              multiVals="first")				
#data export
FDR_up=args[4]
FDR_down=args[5]
P_value=args[6]
FDRsort= (pvaluesort)[pvaluesort$FDR <= P_value,] 
Psort= (FDRsort)[FDRsort$PValue <= P_value,] 
Upregule= (Psort  )[Psort $logFC>= FDR_up,] 
Downregule= (Psort  )[Psort $logFC<= FDR_down,] 
y=rbind(Upregule,Downregule)
write.csv( y, file=args[7], quote = F, row.names = F)
