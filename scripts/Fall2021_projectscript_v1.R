library(dbplyr)
library(MethylIT)
library(MethylIT.utils)
library(dplyr)
library(ggplot2)
library(tidyr)
library(tidyverse)
library(reshape2)

########### TAIR 10 gene annotation ################
AG = import("/Users/jackhoward/Desktop/Fa21_Project/gffannotation.Arabidopsis_thaliana.TAIR10.38.gtf")
#AG means arabidopsis genome
#araport11 is updated version of TAIR10
#gtf files are g ranges with metadata
gene = AG[ AG$type == "gene", c( "gene_id", "gene_biotype" ) ]
#filtering for type == 'gene'
#within this we also have mitochondrial DNA which we do not need
gene = gene[ gene$gene_biotype == "protein_coding", "gene_id" ]
#put gene_id after column to keep metadata
#get 27628 genes, that is what arabidopsis has
seqlevels(gene, pruning.mode = "coarse") <- c("1", "2", "3", "4", "5")
#pruning: trimming off data we dont want, in our case, plastids and mitochondiral genes
#selecting chromosome 1-5 
#GRanges format, something to look up
#seqlevels knows to do this action on the first column
seqlevels(gene)<- paste0("", 1:5)
#in the event that the file has chr1, chr2, we change to number
gene = sortBySeqnameAndStart(gene)
#sorts it 1-5, specific function to sort by chromosome and start site


############ Make a gr object out of sRNA cluster counts ############
RNAcounts <- read.table(file = "/Users/jackhoward/Desktop/Fa21_Project/GSM3933520_memory_gen1-1_counts.txt", sep = "\t", header = TRUE )
dim(RNAcounts)

#factor(RNAcounts$seqnames)

RNAcount_gr <- makeGRangesFromDataFrame(RNAcounts,
                                        keep.extra.columns = TRUE)

#RNAcount_gr$integer is how many sRNAs map into each cluster

#GRanges, type of data objects

############ Find overlap between the gene ranges and sRNA cluster ranges ########
Overlap_gene_RNA <- findOverlaps(RNAcount_gr, gene, ignore.strand = TRUE, maxgap = 0L) # maxgap is how far away can the cluster be - read the help on this function
#RNAcount_gr is query, #gene is the subject
#literature says 2kb away may still have power to influence
#we are just looking at direct overlap

Overlap_gene_RNA <-gene[queryHits(Overlap_gene_RNA)]
#this step picks each of the genes that the query hits are on, we need all the information of genes

Overlap_gene_RNA <-unique(Overlap_gene_RNA)
#just want the unique genes
#8262 ranges

gene_id <- Overlap_gene_RNA$gene_id



#######have the genes of interest, need their descriptions#######

library(biomaRt)

listDatasets(useMart(biomart="plants_mart",host="plants.ensembl.org"))
#retrieves gene descriptions from an online database

mart <- useDataset(dataset = "athaliana_eg_gene", 
                   mart    = useMart("plants_mart",host = "plants.ensembl.org"))
#specifies that we want arabidopsis from the list of plants

head(listFilters(useDataset(dataset = "athaliana_eg_gene", 
                           mart    = useMart("plants_mart",
                                            host = "plants.ensembl.org"))), 100)
#tons of information about the arabidopsis genome that we can choose from

resultTable <- getBM(attributes = c("ensembl_gene_id", "external_gene_name",
                                    "description"),
                     filters = "ensembl_gene_id", values = gene_id, mart = mart)
#other things TAIR gene name, TARI symbol, go_parent_term
#look into mart website

dim(resultTable)
names(resultTable)

write.csv(resultTable, file = "/Users/jackhoward/Desktop/fall21project.csv")


#####Ran David for control, analyzing data#####

david <- read.delim("/Users/jackhoward/Desktop/Fa21_Project/fall21project_davidresult.txt", header = TRUE)

colnames(david) [4] <- "Percent_Genes_Involved"


#clean up term column

ncols <- max(stringr::str_count(david$Term, "~")) + 1
colmn <- paste("col", 1:ncols)

david <- cbind(david, reshape2::colsplit(david$Term, "~", colmn))

names(david)

names(david)[names(david) == "col 1"] <- "GO_Term"
names(david)[names(david) == "col 2"] <- "Function"

#select the most activated functions
primary <-
  david %>%
  filter(Percent_Genes_Involved > 0.5)

#plot of function and percent of genes involved

primary %>%
  ggplot(aes(x = Function, y= Percent_Genes_Involved)) +
  geom_col()
  




