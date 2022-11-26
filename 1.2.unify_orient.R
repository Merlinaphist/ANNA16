library(readr)
library(dplyr)
library(stringr)
library(Biostrings)
reads = read.csv("datasets/cleaned_reads.csv")
blast_results = read_tsv("blast/results.tsv")
colnames(blast_results)[1] = "accession"
reads = left_join(reads,blast_results, by = "accession")

reads = reads%>%
  filter(evalue<0.01)

rev_reads = reads[(reads$sstart>reads$send),1:5]
fwd_reads = reads[(reads$sstart<reads$send),1:5]

rc = c()
for (seq in rev_reads$sequence){
  rc = append(rc,as.character(reverseComplement(DNAString(seq))))
}

rev_reads$sequence = rc

reads = rbind(fwd_reads,rev_reads)

set.seed(114514)
reads = reads[sample(1:nrow(reads)),]

write.csv(reads,"datasets/oriented_reads.csv",row.names = F)

fa = character(2*nrow(reads))
fa[c(TRUE, FALSE)] = sprintf(">%s", reads$accession)
fa[c(FALSE, TRUE)] = reads$sequence
writeLines(fa,"datasets/oriented_reads.fasta")


