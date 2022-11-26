reads = read.csv("datasets/oriented_reads.csv")
fasta = seqinr::read.fasta(file = 'datasets/full_length.fasta', as.string = TRUE,
                                   forceDNAtolower = FALSE, whole.header = TRUE)
fasta = data.frame(unlist(fasta))
colnames(fasta) = "sequence"
fasta$accession = row.names(fasta)
full_length = dplyr::left_join(fasta,reads[,c(1,3:5)],by = "accession")
write.csv(full_length,"datasets/full_length_reads.csv",row.names = F)



