hyperex -p datasets/tmp_full_length -m 2 --forward-primer AGAGTTTGATCCTGGCTCAG --reverse-primer TACGGYTACCTTGTTACGACT datasets/oriented_reads.fasta 
hyperex -p datasets/tmp_V1-V2 -m 2 --forward-primer AGAGTTTGATCCTGGCTCAG --reverse-primer GCTGCCTCCCGTAGGAGT datasets/oriented_reads.fasta 
hyperex -p datasets/tmp_V1-V3 -m 2 --forward-primer AGAGTTTGATCCTGGCTCAG --reverse-primer ATTACCGCGGCTGCTGG datasets/oriented_reads.fasta 
hyperex -p datasets/tmp_V3-V4 -m 2 --forward-primer CCTACGGGAGGCAGCAG --reverse-primer GACTACHVGGGTATCTAATCC datasets/oriented_reads.fasta 
hyperex -p datasets/tmp_V4 -m 2 --forward-primer GTGCCAGCMGCCGCGGTAA --reverse-primer GGACTACHVGGGTWTCTAAT datasets/oriented_reads.fasta 
hyperex -p datasets/tmp_V4-V5 -m 2 --forward-primer GTGCCAGCMGCCGCGGTAA --reverse-primer CCGYCAATTYMTTTRAGTTT datasets/oriented_reads.fasta 
hyperex -p datasets/tmp_V6-V8 -m 2 --forward-primer GAATTGACGGGGGCCCGCACAAG --reverse-primer CGGTGTGTACAAGGCCCGGGAACG datasets/oriented_reads.fasta 
hyperex -p datasets/tmp_V7-V9 -m 2 --forward-primer CAACGAGCGCAACCCT --reverse-primer TACGGYTACCTTGTTACGACT datasets/oriented_reads.fasta 
rm datasets/*.gff
awk -F" " '{print $1}' datasets/tmp_full_length.fa > datasets/full_length.fasta
awk -F" " '{print $1}' datasets/tmp_V1-V2.fa > datasets/V1-V2.fasta
awk -F" " '{print $1}' datasets/tmp_V1-V3.fa > datasets/V1-V3.fasta
awk -F" " '{print $1}' datasets/tmp_V3-V4.fa > datasets/V3-V4.fasta
awk -F" " '{print $1}' datasets/tmp_V4.fa > datasets/V4.fasta
awk -F" " '{print $1}' datasets/tmp_V4-V5.fa > datasets/V4-V5.fasta
awk -F" " '{print $1}' datasets/tmp_V6-V8.fa > datasets/V6-V8.fasta
awk -F" " '{print $1}' datasets/tmp_V7-V9.fa > datasets/V7-V9.fasta
rm datasets/*.fa

hyperex -p datasets/tmp_V1_simu -m 2 --forward-primer AGAGTTTGATCATGGCTCA --reverse-primer TTACTCACCCGTBCGCC datasets/oriented_reads.fasta 
hyperex -p datasets/tmp_V2_simu -m 2 --forward-primer GGCGVACGGGTGAGTAA --reverse-primer GGRCCGTGTCTCAGTYCCARTG datasets/oriented_reads.fasta 
hyperex -p datasets/tmp_V3_simu -m 2 --forward-primer CAYTGGRACTGAGACACGGYCC --reverse-primer GCTGCTGGCACGDAGTTAGCC datasets/oriented_reads.fasta 
hyperex -p datasets/tmp_V4_simu -m 2 --forward-primer GGCTAACTHCGTGCCAGCAGC --reverse-primer CCTGYGMTMCCCACRCTTTCG datasets/oriented_reads.fasta 
hyperex -p datasets/tmp_V5_simu -m 2 --forward-primer CGAAAGYGTGGGKAKCRCAGG --reverse-primer GYCCCCGTCAATTCMTTTGAGT datasets/oriented_reads.fasta 
hyperex -p datasets/tmp_V6_simu -m 2 --forward-primer ACTCAAAKGAATTGACGGGGRC --reverse-primer AGCTGACGACARCCATGCASCAC datasets/oriented_reads.fasta 
hyperex -p datasets/tmp_V7_simu -m 2 --forward-primer GTGSTGCATGGYTGTCGTCAGCT --reverse-primer TTGACGTCRTCCCCRCCTTCC datasets/oriented_reads.fasta 
hyperex -p datasets/tmp_V8_simu -m 2 --forward-primer GGAAGGYGGGGAYGACGTCAA --reverse-primer GTGTGTGACGGGCGGTGTGTACA datasets/oriented_reads.fasta 
hyperex -p datasets/tmp_V9_simu -m 2 --forward-primer TGTACACACCGCCCGTCACACAC --reverse-primer TACGGTTACCTTGTTACGACTT datasets/oriented_reads.fasta 
rm datasets/*.gff
awk -F" " '{print $1}' datasets/tmp_V1_simu.fa > datasets/V1_simu.fasta
awk -F" " '{print $1}' datasets/tmp_V2_simu.fa > datasets/V2_simu.fasta
awk -F" " '{print $1}' datasets/tmp_V3_simu.fa > datasets/V3_simu.fasta
awk -F" " '{print $1}' datasets/tmp_V4_simu.fa > datasets/V4_simu.fasta
awk -F" " '{print $1}' datasets/tmp_V5_simu.fa > datasets/V5_simu.fasta
awk -F" " '{print $1}' datasets/tmp_V6_simu.fa > datasets/V6_simu.fasta
awk -F" " '{print $1}' datasets/tmp_V7_simu.fa > datasets/V7_simu.fasta
awk -F" " '{print $1}' datasets/tmp_V8_simu.fa > datasets/V8_simu.fasta
awk -F" " '{print $1}' datasets/tmp_V9_simu.fa > datasets/V9_simu.fasta
rm datasets/*.fa