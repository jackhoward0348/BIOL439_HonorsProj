####Loading Packages
library(dbplyr)
library(MethylIT)
library(MethylIT.utils)
library(dplyr)
library(ggplot2)
library(tidyr)
library(tidyverse)
library(reshape2)
library(biomaRt)
library(scales)

####Importing Genome
AG <-data.table::fread("https://www.arabidopsis.org/download_files/Genes/TAIR10_genome_release/TAIR10_gff3/TAIR10_GFF3_genes.gff")
#AG means arabidopsis genome
#araport11 is updated version of TAIR10
#gtf files are g ranges with metadata

AG <- AG[ ,c(1,3,4,5,9)]
colnames(AG) <- c("chromosome", "type", "start", "end", "metadata")

AT_genes <-
  AG %>%
  filter(type == 'gene')
#filtering for type == 'gene', getting rid of mRNA, exons, CDS, mitochondrial DNA 

AT_genes <- separate(AT_genes, metadata, c("gene_id", "biotype", "name"), ";")
AT_genes$gene_id <- gsub(pattern = ".*=", replacement = "", AT_genes$gene_id)
AT_genes$biotype <- gsub(pattern = ".*=", replacement = "", AT_genes$biotype)
AT_genes$name <- gsub(pattern = ".*=", replacement = "", AT_genes$name)
#regex used to separate and clean up columns

AT_genes <-
  AT_genes %>%
  filter(biotype == "protein_coding_gene")

AT_genes$chromosome <- gsub(pattern = "[Chr]", replacement = "", AT_genes$chromosome)

AT_genes <-
  AT_genes %>%
  filter(chromosome %in% c(1,2,3,4,5))
#getting rid of mitochondrial and plastid locations

AT_genes <- makeGRangesFromDataFrame(AT_genes, keep.extra.columns = TRUE)


####sRNA counts####
sRNAcounts_wild_1 <- data.table::fread("https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM3933nnn/GSM3933516/suppl/GSM3933516_WT-1_counts.txt.gz")
#txt file from the paper, copy link address of ftp download button and import from the online source

sRNAcounts_wild_2 <- data.table::fread("https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM3933nnn/GSM3933517/suppl/GSM3933517_WT-2_counts.txt.gz")

sRNAcounts_wild_3 <- data.table::fread("https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM3933nnn/GSM3933518/suppl/GSM3933518_WT-3_counts.txt.gz")

sRNAcounts_wild_4 <- data.table::fread("https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM3933nnn/GSM3933519/suppl/GSM3933519_WT-4_counts.txt.gz")

sRNAcounts_memory_1 <- data.table::fread("https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM3933nnn/GSM3933520/suppl/GSM3933520_memory_gen1-1_counts.txt.gz")

sRNAcounts_memory_2 <- data.table::fread("https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM3933nnn/GSM3933521/suppl/GSM3933521_memory_gen1-2_counts.txt.gz")

#combining and cleaning
sRNAcounts_comb <- merge(sRNAcounts_wild_1, sRNAcounts_wild_2, by = "Name")
sRNAcounts_comb <-
  sRNAcounts_comb %>%
  dplyr::select(Name, start.x, end.x, "WT-1", "WT-2")

sRNAcounts_comb <- merge(sRNAcounts_comb, sRNAcounts_wild_3, by = "Name")
sRNAcounts_comb <-
  sRNAcounts_comb %>%
  dplyr::select(Name, start.x, end.x, "WT-1", "WT-2", "WT-3")

sRNAcounts_comb <- merge(sRNAcounts_comb, sRNAcounts_wild_4, by = "Name")
sRNAcounts_comb <-
  sRNAcounts_comb %>%
  dplyr::select(Name, start.x, end.x, "WT-1", "WT-2", "WT-3", "WT-4")

sRNAcounts_comb <- merge(sRNAcounts_comb, sRNAcounts_memory_1, by = "Name")
sRNAcounts_comb <-
  sRNAcounts_comb %>%
  dplyr::select(Name, "WT-1", "WT-2", "WT-3", "WT-4", "memory_gen1-1", start.x, end.x)

sRNAcounts_comb <- merge(sRNAcounts_comb, sRNAcounts_memory_2, by = "Name")
sRNAcounts_comb <-
  sRNAcounts_comb %>%
  dplyr::select(Name, "WT-1", "WT-2", "WT-3", "WT-4", "memory_gen1-1", "memory_gen1-2", seqnames, start.x, end.x)

colnames(sRNAcounts_comb) <- c('Name', 'WT-1', 'WT-2', 'WT-3', 'WT-4', 'MEM-1', 'MEM-2', 'seqnames', 'start', 'end')

####determining clusters of interest
#find average sRNA cluster for each of the treatments and mutate into a new column

sRNAcounts_comb <-
  sRNAcounts_comb %>%
  mutate('WT_avg' = rowMeans(dplyr::select(sRNAcounts_comb, starts_with("WT")), na.rm = TRUE)) %>%
  mutate('MEM_avg' = rowMeans(dplyr::select(sRNAcounts_comb, starts_with("MEM")), na.rm = TRUE))

sRNAcounts_comb <-
  sRNAcounts_comb %>%
  mutate(foldchange = as.numeric(WT_avg/MEM_avg)) %>%
  mutate(the_rank = rank(foldchange))

sRNAcounts_higherWTexp <-
  sRNAcounts_comb %>%
  filter(foldchange > 4) %>%
  dplyr::select(Name, start, end, seqnames, foldchange)

sRNAcounts_higherMEMexp <-
  sRNAcounts_comb %>%
  filter(foldchange < 0.35) %>%
  dplyr::select(Name, start, end, seqnames, foldchange)
#thresholds are arbitrary

####make a GRange####
sRNAcounts_higherWTexp_gr <- makeGRangesFromDataFrame(sRNAcounts_higherWTexp, keep.extra.columns = TRUE)

Overlaps_sRNAWTexp <- findOverlaps(sRNAcounts_higherWTexp_gr, AT_genes, ignore.strand = TRUE, maxgap = 0L)
#RNAcount_gr is query, #gene is the subject
#literature says 2kb away may still have power to influence
#we are just looking at direct overlap

Overlaps_sRNAWTexp <-AT_genes[queryHits(Overlaps_sRNAWTexp)]
#this step picks each of the genes that the query hits are on, we need all the information of genes

Overlaps_sRNAWTexp <-unique(Overlaps_sRNAWTexp)
#just want the unique genes

geneid_sRNA_higherWTexp <- Overlaps_sRNAWTexp$gene_id
#list of genes for higher wild type expressions based on fold change

###memory high exp###
sRNAcounts_higherMEMexp_gr <- makeGRangesFromDataFrame(sRNAcounts_higherMEMexp, keep.extra.columns = TRUE)

Overlaps_sRNAMEMexp <- findOverlaps(sRNAcounts_higherMEMexp_gr, AT_genes, ignore.strand = TRUE, maxgap = 0L)
#RNAcount_gr is query, #gene is the subject
#literature says 2kb away may still have power to influence
#we are just looking at direct overlap

Overlaps_sRNAMEMexp <-AT_genes[queryHits(Overlaps_sRNAMEMexp)]
#this step picks each of the genes that the query hits are on, we need all the information of genes

Overlaps_sRNAMEMexp <-unique(Overlaps_sRNAMEMexp)
#just want the unique genes

geneid_sRNA_higherMEMexp <- Overlaps_sRNAMEMexp$gene_id

####gene descriptions
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

resultTable_higherWTexp <- getBM(attributes = c("ensembl_gene_id", "external_gene_name",
                                    "description", "name_1006", "definition_1006", "namespace_1003"),
                     filters = "ensembl_gene_id", values = geneid_sRNA_higherWTexp, mart = mart)


resultTable_higherMEMexp <- getBM(attributes = c("ensembl_gene_id", "external_gene_name",
                                                "description", "name_1006", "definition_1006", "namespace_1003"),
                                 filters = "ensembl_gene_id", values = geneid_sRNA_higherMEMexp, mart = mart)

###david
write.csv(geneid_sRNA_higherWTexp, file = "/Users/jackhoward/Desktop/geneid_sRNA_higherWTexp.csv")

write.csv(geneid_sRNA_higherMEMexp, file = "/Users/jackhoward/Desktop/geneid_sRNA_higherMEMexp.csv")










