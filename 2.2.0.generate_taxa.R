fasta = seqinr::read.fasta(file = 'datasets/full_length.fasta', as.string = TRUE,
                           forceDNAtolower = FALSE, whole.header = TRUE)
fasta = data.frame(unlist(fasta))
colnames(fasta) = "sequence"
fasta$accession = row.names(fasta)
fasta = fasta[(duplicated(fasta$sequence==F)),]
input = fasta$sequence
for (database in c("rdp_train_set_18.fa.gz","silva_nr99_v138.1_train_set.fa.gz","pr2_version_4.14.0_SSU_dada2.fasta.gz")){
  tag = strsplit(database,"_")[[1]][1]
  taxa <- dada2::assignTaxonomy(input, paste0("raw_data/",database))
  taxa = data.frame(taxa)
  taxa$sequence = row.names(taxa)
  reads = read.csv("datasets/full_length_reads.csv")
  reads = dplyr::left_join(reads,taxa,by = "sequence")
  write.csv(reads,paste0("taxa/",tag,"_full_length_taxa.csv"),row.names=F)
}