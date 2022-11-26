cd blast
makeblastdb -dbtype nucl -in marker.fasta -out test
blastn -db test -query cleaned_reads.fasta -task blastn-short -outfmt "6 qseqid qlen qstart qend stitle slen sstart send pident bitscore evalue" -out output.tsv -max_hsps 1
echo -e "qseqid\tqlen\tqstart\tqend\tstitle\tslen\tsstart\tsend\tpident\tbitscore\tevalue" | cat - output.tsv > results.tsv