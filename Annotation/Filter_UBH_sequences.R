library(dplyr)
library(phylotools)

uniprot_ubh <- 
  read.table("uniprot_vs_macro_UBH.txt", header = T) %>% 'colnames<-'(c("uni_id"))

head(uniprot_ubh)

uniprot_fasta <- 
  read.fasta("uniprot_sprot.fasta") 
nrow(uniprot_ubh) == nrow(uniprot_fasta)

colnames(uniprot_fasta)
head(uniprot_fasta$seq.name)

uniprot_fasta <- uniprot_fasta %>% 
  mutate(seq.name= gsub(" .*","",seq.name)) %>% 
  filter(seq.name %in% uniprot_ubh$uni_id)

#Check results
length(unique(uniprot_ubh$uni_id)) == nrow(uniprot_fasta)

#Save UBH Uniprot protein sequences
dat2fasta(uniprot_fasta, "uniprot_fasta_ubh.fasta")
