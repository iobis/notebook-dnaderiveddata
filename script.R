library(dplyr)
library(robis)
library(ape)
library(msa)
library(dendextend)

# get data from OBIS

occ <- occurrence(hasextensions = "DNADerivedData", extensions = "DNADerivedData") %>%
  filter(!is.na(phylum))

table(occ$phylum)
occ %>%
  group_by(phylum) %>%
  summarize(species = length(unique(speciesid)), records = n()) %>%
  arrange(desc(species))

dna <- dna_records(occ, fields = c("id", "scientificName", "phylum", "class", "order", "family", "genus", "species", "eventDate"))

# random subset data

dna_subset <- dna[sample(1:nrow(dna), 150),]

# align with clustal

aligned <- msa::msa(dna_subset$DNA_sequence, type = "dna")

# convert to ape DNAbin

sequences <- t(sapply(strsplit(as.character(aligned), ""), tolower))
rownames(sequences) <- paste0(substr(dna_subset$id, 1, 5), " - ", dna_subset$phylum)
sequences <- ape::as.DNAbin(sequences)
sequences

# dist.dna and hclust

distances <- ape::dist.dna(sequences, model = "raw")
h_cluster <- hclust(distances, method = "average")

dend <- as.dendrogram(h_cluster) %>% set("labels_cex", 0.6)
dend <- dend %>% set("labels_col", as.numeric(as.factor(substr(labels(dend), 9, 100))))

plot(dend)
