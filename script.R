library(dplyr)
library(robis)
library(ape)
library(ips)
library(dendextend)
library(ggplot2)
library(ggtree)

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


# convert to ape DNAbin

sequences <- t(sapply(strsplit(as.character(dna_subset), ""), tolower))
rownames(sequences) <- paste0(substr(dna_subset$id, 1, 5), " - ", dna_subset$phylum)
sequences <- ape::as.DNAbin(sequences)
sequences

#align with mafft

ali_mafft=ips::mafft(sequences, method = "auto")


# dist.dna and hclust

distances <- ape::dist.dna(ali_mafft, model = "raw")
h_cluster <- hclust(distances, method = "average")

dend <- as.dendrogram(h_cluster) %>% set("labels_cex", 0.6)
dend <- dend %>% set("labels_col", as.numeric(as.factor(substr(labels(dend), 9, 100))))

plot(dend)

#For flexible phylogenetic tree visualization, can build on this new package:

tree<-njs(distances)
ggt<-ggtree(tree, cex = 0.8, aes(color=branch.length))+geom_tiplab(size=2)
ggt

#Also bootstrapping

#Bootstrapping
myBoots <- boot.phylo(tree, ali_mafft, function(e) njs(dist.dna(e, model = "TN93")))
bp2 <- tibble(node=1:Nnode(tree) + Ntip(tree), bootstrap = myBoots)
bptree <- full_join(tree, bp2, by="node")

ggtree(bptree) +
  geom_tiplab(size=2) +
  geom_nodelab(aes(label=bootstrap), size=2)
