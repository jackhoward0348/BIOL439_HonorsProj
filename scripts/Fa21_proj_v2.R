#In this script, I did the analysis for the mutant and the wild type all the way through, getting the clusters and running david. 
#This yielded a result that showed the same clusters on the same genes
#confused, Alenka informed me that it is not the clusters that differ between the treatments, it is the count of sRNAs in each cluster
#now the project will focus on doing a fold change, seeing the difference in the integer column for each cluster 

####Loading Packages
  ```{r}
library(dbplyr)
library(MethylIT)
library(MethylIT.utils)
library(dplyr)
library(ggplot2)
library(tidyr)
library(tidyverse)
library(reshape2)
library(biomaRt)
```

####Importing Genome
```{r}
AG = import("/Users/jackhoward/Desktop/Fa21_Project/gffannotation.Arabidopsis_thaliana.TAIR10.38.gtf")
#AG means arabidopsis genome
#araport11 is updated version of TAIR10
#gtf files are g ranges with metadata

AT_genes = AG[ AG$type == "gene", c( "gene_id", "gene_biotype" ) ]
#filtering for type == 'gene', filter() function is not usable for objects of class GRange
#within this we also have mitochondrial DNA which we do not need

AT_genes = AT_genes[ AT_genes$gene_biotype == "protein_coding", "gene_id" ]
#put gene_id after column to keep metadata
#get 27628 genes, that is what arabidopsis has

seqlevels(AT_genes, pruning.mode = "coarse") <- c("1", "2", "3", "4", "5")
#pruning: trimming off data we dont want, in our case, plastids and mitochondiral genes
#selecting chromosome 1-5 
#seqlevels knows to do this action on the first column

seqlevels(AT_genes)<- paste0("", 1:5)
#in the event that the file has chr1, chr2, we change to number

AT_genes = sortBySeqnameAndStart(AT_genes)
#sorts it 1-5, specific function to sort by chromosome and start site
```

##Wild Type
#sRNA counts
```{r}
sRNAcounts_wild <- read.table(file = "/Users/jackhoward/Desktop/Fa21_Project/GSM3933516_WT-1_counts.txt", sep = "\t", header = TRUE )
#txt file from the paper

dim(sRNAcounts_wild)

#factor(RNAcounts_wild$seqnames)

sRNAcount_wild_gr <- makeGRangesFromDataFrame(sRNAcounts_wild,
                                              keep.extra.columns = TRUE)

#sRNAcount_wild_gr$integer is how many sRNAs map into each cluster
```

#Overlaps between
```{r}
Overlaps_sRNA_wild <- findOverlaps(sRNAcount_wild_gr, AT_genes, ignore.strand = TRUE, maxgap = 0L) # maxgap is how far away can the cluster be - read the help on this function
#RNAcount_gr is query, #gene is the subject
#literature says 2kb away may still have power to influence
#we are just looking at direct overlap

Overlaps_sRNA_wild <-AT_genes[queryHits(Overlaps_sRNA_wild)]
#this step picks each of the genes that the query hits are on, we need all the information of genes

Overlaps_sRNA_wild <-unique(Overlaps_sRNA_wild)
#just want the unique genes
#8262 ranges

gene_id_wild <- Overlaps_sRNA_wild$gene_id
```

#Matching genes to descriptions
```{r}
listDatasets(useMart(biomart="plants_mart",host="plants.ensembl.org"))
#retrieves gene descriptions from an online database

mart <- useDataset(dataset = "athaliana_eg_gene", 
                   mart    = useMart("plants_mart",host = "plants.ensembl.org"))
#specifies that we want arabidopsis from the list of plants

head(listFilters(useDataset(dataset = "athaliana_eg_gene", 
                            mart    = useMart("plants_mart",
                                              host = "plants.ensembl.org"))), 100)
#tons of information about the arabidopsis genome that we can choose from
#visualizing data source

resultTable <- getBM(attributes = c("ensembl_gene_id", "external_gene_name",
                                    "description"),
                     filters = "ensembl_gene_id", values = gene_id_wild, mart = mart)
#other things TAIR gene name, TARI symbol, go_parent_term
#look into mart website

dim(resultTable)
names(resultTable)

write.csv(resultTable, file = "/Users/jackhoward/Desktop/Fa21_Project/wild_geneID_description.csv")
```

In between these steps, took the CSV file and isolated just the names of each of the genes, the names are analyzed in DAVID, an online software

#David analysis
```{r}
david_wild <- read.delim("/Users/jackhoward/Desktop/Fa21_Project/wild_davidresult.txt", header = TRUE)

colnames(david_wild) [4] <- "Percent_Genes_Involved"

#cleaning david_wild dataframe
ncols <- max(stringr::str_count(david_wild$Term, "~")) + 1
colmn <- paste("col", 1:ncols)

david_wild <- cbind(david_wild, reshape2::colsplit(david_wild$Term, "~", colmn))

names(david_wild)[names(david_wild) == "col 1"] <- "GO_Term"
names(david_wild)[names(david_wild) == "col 2"] <- "Function"
```

##Mutant
#sRNA counts
```{r}
sRNAcounts_mut <- read.table(file = "/Users/jackhoward/Desktop/Fa21_Project/GSM3933520_memory_gen1-1_counts.txt", sep = "\t", header = TRUE )

dim(sRNAcounts_mut)

sRNAcount_mut_gr <- makeGRangesFromDataFrame(sRNAcounts_mut,
                                             keep.extra.columns = TRUE)
```

#Overlaps between
```{r}
Overlaps_sRNA_mut <- findOverlaps(sRNAcount_mut_gr, AT_genes, ignore.strand = TRUE, maxgap = 0L) 

Overlaps_sRNA_mut <-AT_genes[queryHits(Overlaps_sRNA_mut)]

Overlaps_sRNA_mut <-unique(Overlaps_sRNA_mut)

gene_id_mut <- Overlaps_sRNA_mut$gene_id
```

#Matching genes to descriptions
```{r}
resultTable2 <- getBM(attributes = c("ensembl_gene_id", "external_gene_name",
                                     "description"),
                      filters = "ensembl_gene_id", values = gene_id_mut, mart = mart)

write.csv(resultTable, file = "/Users/jackhoward/Desktop/Fa21_Project/mut_geneID_description.csv")
```

In between these steps, took the CSV file and isolated just the names of each of the genes, the names are analyzed in DAVID, an online software

#David analysis
```{r}
david_mut <- read.delim("/Users/jackhoward/Desktop/Fa21_Project/mut_davidresult.txt", header = TRUE)

colnames(david_mut) [4] <- "Percent_Genes_Involved"

#cleaning david_wild dataframe
ncols <- max(stringr::str_count(david_mut$Term, "~")) + 1
colmn <- paste("col", 1:ncols)

david_mut <- cbind(david_mut, reshape2::colsplit(david_mut$Term, "~", colmn))

names(david_mut)[names(david_mut) == "col 1"] <- "GO_Term"
names(david_mut)[names(david_mut) == "col 2"] <- "Function"
```





