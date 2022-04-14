Input the fasta file of Hydractinia symbiolongicarpus aa sequences into eggnog5 mapper at 

http://eggnog5.embl.de/

download and extract the output "EggNOGout.txt" into the working directory of the R.

#In R install packages required to run the analysis

if (!requireNamespace("BiocManager", quietly=TRUE))
+ install.packages("BiocManager")
BiocManager::install("clusterProfiler")
BiocManager::install("enrichplot")
BiocManager::install("dplyr")
install.packages("data.table")
install.packages("dplyr")
install.packages("tidyverse")

library(clusterProfiler)
library(enrichplot)
library(data.table)

## Get the necessary information from the eggNOG output file

# read the data
EggData <- read.delim("EggNOGout.txt", header = TRUE)
head(EggData)

                           -
         EC   KEGG_ko KEGG_Pathway KEGG_Module KEGG_Reaction KEGG_rclass                   BRITE           KEGG_TC CAZy BiGG_Reaction
1         - ko:K15193            -           -             -           -         ko00000,ko03021                 -    -             -
2         - ko:K15193            -           -             -           -         ko00000,ko03021                 -    -             -
3 2.3.1.225 ko:K20032            -           -             -           - ko00000,ko01000,ko04131 9.B.37.1,9.B.37.3    -             -
4         -         -            -           -             -           -                       -                 -    -             -
5         -         -            -           -             -           -                       -                 -    -             -
6         -         -            -           -             -           -                       -                 -    -             -
                         PFAMs
1                         SPT2
2                         SPT2
3 Ank_2,Ank_3,Ank_4,Ank_5,DHHC
4                        RVT_1
5                            -
6                          SET


#get columns 1 (query) and 12 (KO terms)
KEggDataP <- EggData[c(1,13)]
head(KEggDataP)

query KEGG_Pathway
1:  HyS0001.3            -
2:  HyS0001.4            -
3:  HyS0001.9            -
4: HyS0001.13            -
5: HyS0001.21            -
6: HyS0001.23            -

# clean up by removing the "ko:" in front of every KO term
KEggDataP$KEGG_ko <- gsub( "ko:", "", as.character(KEggDataP$KEGG_Pathway))

# expand, since some genes/proteins will have multiple assigned KO terms
KEggDataP <- data.table(KEggDataP)
KEggDataP <- KEggDataP[, list(KEGG_Pathway = unlist(strsplit(KEGG_Pathway , ","))), by = query]

# select the needed columns
KEggDataPFinal <- KEggDataP[,c(2,1)]
head(KEggDataPFinal)

 KEGG_Pathway      query
1:            -  HyS0001.3
2:            -  HyS0001.4
3:            -  HyS0001.9
4:            - HyS0001.13
5:            - HyS0001.21
6:            - HyS0001.23

HelDE1=read.table(file = "HelDE36DnLFC1.txt", header = TRUE, quote = "")
HelDE1=data.frame(HelDE1)

library(dplyr)
library(tidyverse)

HelDE1GO = HelDE1 %>%
dplyr::select(Gene,Log2FC,Pval)

HelDE1Sort <- HelDE1GO[order(-HelDE1GO$Log2FC), ]
str(HelDE1Sort)

'data.frame':	14 obs. of  3 variables:
 $ Gene  : chr  "HyS0407.3" "HyS0007.80" "HyS0002.525" "HyS0033.9" ...
 $ Log2FC: num  -2.18 -2.35 -2.51 -2.68 -3 ...
 $ Pval  : num  3.83e-09 1.84e-06 1.32e-09 7.09e-05 9.30e-06 3.32e-05 2.25e-05 2.43e-05 1.73e-05 1.72e-09 ...

geneList <- factor(as.integer(geneNames %in% HelDE1Sort$Gene))
str(geneList)

 Factor w/ 2 levels "0","1": 1 1 1 1 1 1 1 1 1 1 ...

names(geneList) <- HelDE1Sort$Gene
str(geneList)

Factor w/ 2 levels "0","1": 1 1 1 1 1 1 1 1 1 1 ...
 - attr(*, "names")= chr [1:11127] "HyS0407.3" "HyS0007.80" "HyS0002.525" "HyS0033.9" ...

enr_HelDE1 <- enricher(HelDE1Sort$Gene, TERM2GENE=KEggDataPFinal, pvalueCutoff = 0.5 , pAdjustMethod = "BH", qvalueCutoff = 0.5 , minGSSize = 10)

# write the results of the analysis
write.table(enr_HelDE1, file = "Hel36DnLfc1_KOEnr_results.csv", sep = " ") 
# write table
