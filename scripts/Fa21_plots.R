####plot of fold change
sRNAcounts_comb_plot3 <- separate(sRNAcounts_comb, Name, c("Cluster", "Number"))
#need to change the name of the first column so we can get a numeric value for each cluster
sRNAcounts_comb_plot3$Number <- as.numeric(sRNAcounts_comb_plot3$Number)
#change to numeric
sRNAcounts_comb_plot3$Chromosome <- as.numeric(sRNAcounts_comb_plot3$Chromosome)


sRNAcounts_comb_plot3 <-
  sRNAcounts_comb_plot3 %>%
  pivot_longer(cols = c(WT_avg, MEM_avg), names_to = "Treatment", values_to = "Avg_Count") %>%
  dplyr::select(Number, seqnames, Treatment, Avg_Count)

sRNAcounts_comb_plot %>%
  ggplot(aes(x = Number, y= Avg_Count)) +
  geom_point(aes(color= Treatment)) +
  ylim(0, 10000)

####plot of differences between treatments for each cluster
sRNAcounts_comb_plot2 <- separate(sRNAcounts_comb, Name, c("Cluster", "Number"))
sRNAcounts_comb_plot2$Number <- as.numeric(sRNAcounts_comb_plot2$Number)
names(sRNAcounts_comb_plot2)[names(sRNAcounts_comb_plot2) == 'seqnames'] <- 'Chromosome'






###violin plot
sRNAcounts_comb_plot2 %>%
  ggplot(aes(x = Chromosome, y= difference_ofavg)) +
  geom_violin() +
  scale_y_continuous(name="sRNA Count Difference", labels = comma) +
  xlab("Cluster") +
  ggtitle("Difference in sRNA Count per Cluster")


###for the results, show count of gene difference for each funtion
WT_result_plot <- resultTable_higherWTexp[ , c(1,4,6)]
colnames(WT_result_plot) <- c("ID", "funct", "class")
WT_result_plot <-
  WT_result_plot %>%
  filter(class == "molecular_function" | class == "biological_process")

WT_result_plot2 <-
  WT_result_plot %>%
  count(funct)

WT_result_plot %>%
  ggplot(aes())

# break it down by GO terms bc straight up function is too much


resultTable_higherWTexp2 <- unique(resultTable_higherWTexp$external_gene_name)
write.csv(resultTable_higherWTexp2, file = "/Users/jackhoward/Desktop/geneid_sRNA_higherWTexp.csv")

resultTable_higherMEMexp2 <- unique(resultTable_higherMEMexp$external_gene_name)
write.csv(resultTable_higherMEMexp2, file = "/Users/jackhoward/Desktop/geneid_sRNA_higherMEMexp.csv")


WT_princeton_results <- "https://go.princeton.edu/tmp/7362500//query.html"
WT_table <- web_page %>%
  read_html() %>%
  html_nodes(css = "table") %>%
  html_table(fill = TRUE)

WT_GOparent <- WT_table[[1]]

###fig 3
GOparent_terms$GOterm <- str_wrap(GOparent_terms$GOterm, width = 24)

GOparent_terms %>%
  ggplot(aes(x = GOterm, y = freq, fill = scale)) +
  geom_bar(stat = 'identity', position = position_dodge()) +
  theme(axis.text.x = element_text(angle = 55, hjust = 1), legend.title = element_blank()) +
  ylab("Frequency of Genes (%)") +
  xlab("Gene Ontology Parent Term") +
  ggtitle("Frequency of Genes Per GO Term in Significant Genes vs Genome")

#fig 2
sRNAcounts_comb_plot2 %>%
  ggplot(aes(x = Chromosome, fill = significance)) + 
  geom_bar() +
  coord_flip() +
  ylab("Number of Clusters") +
  ggtitle("Significance of Fold Change per Cluster")

#figure 1
sRNAcounts_comb_plot %>%
  ggplot(aes(x = Number, y= difference_ofavg)) +
  geom_point(aes(color= Chromosome)) +
  scale_y_continuous(name="sRNA Count Difference", labels = comma) +
  xlab("Cluster") +
  ggtitle("Difference in sRNA Count per Cluster") +
  ylim(-180000, 150000) +
  ylab("Difference in sRNA counts (WT - MEM)")


