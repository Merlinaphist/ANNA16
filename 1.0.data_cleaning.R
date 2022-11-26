library(seqinr)
library(readr)
library(dplyr)
library(Biostrings)
library(DECIPHER)
library(ggplot2)
library(stringr)
reads = seqinr::read.fasta(file = 'raw_data/rrnDB-5.7_16S_rRNA.fasta', as.string = TRUE,
                           forceDNAtolower = FALSE, whole.header = TRUE)
reads = data.frame(unlist(reads))
reads= cbind(row.names(reads),reads)
row.names(reads) = seq(1,nrow(reads))
colnames(reads) = c('name',"sequence")
reads$accession = sapply(strsplit(reads$name,split = "|", fixed = T),`[`,2)
reads$name = sapply(strsplit(reads$name,split = "|", fixed = T),`[`,1)

reads$name = str_replace_all(reads$name,"'","")
reads$name = str_replace_all(reads$name,"\\]","")
reads$name = str_replace_all(reads$name,"\\[","")


meta = read_tsv("raw_data/rrnDB-5.7.tsv")
meta = select(meta,1:2,12)
colnames(meta) = c("accession","taxid","copy_number")
reads = left_join(reads,meta,by = "accession")
rm(meta)
reads = reads[(duplicated(reads$accession)==F),]
set.seed(114514)
reads = reads[sample(1:nrow(reads)),]

write.csv(reads,"datasets/cleaned_reads.csv",row.names = F)

fa = character(2*nrow(reads))
fa[c(TRUE, FALSE)] = sprintf(">%s", reads$accession)
fa[c(FALSE, TRUE)] = reads$sequence
writeLines(fa,"blast/cleaned_reads.fasta")
