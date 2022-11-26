# full_length
mafft datasets/full_length.fasta > tree/full_length_aln.fasta
FastTree -spr 4 -gamma -fastest -no2nd -constraintWeight 100 -nt tree/full_length_aln.fasta > tree/full_length.tre
# V1-V2
mafft datasets/V1-V2.fasta > tree/V1-V2_aln.fasta
FastTree -spr 4 -gamma -fastest -no2nd -constraintWeight 100 -nt tree/V1-V2_aln.fasta > tree/V1-V2.tre
# V1-V3
mafft datasets/V1-V3.fasta > tree/V1-V3_aln.fasta
FastTree -spr 4 -gamma -fastest -no2nd -constraintWeight 100 -nt tree/V1-V3_aln.fasta > tree/V1-V3.tre
# V3-V4
mafft datasets/V3-V4.fasta > tree/V3-V4_aln.fasta
FastTree -spr 4 -gamma -fastest -no2nd -constraintWeight 100 -nt tree/V3-V4_aln.fasta > tree/V3-V4.tre
# V4
mafft datasets/V4.fasta > tree/V4_aln.fasta
FastTree -spr 4 -gamma -fastest -no2nd -constraintWeight 100 -nt tree/V4_aln.fasta > tree/V4.tre
# V4-V5
mafft datasets/V4-V5.fasta > tree/V4-V5_aln.fasta
FastTree -spr 4 -gamma -fastest -no2nd -constraintWeight 100 -nt tree/V4-V5_aln.fasta > tree/V4-V5.tre
# V6-V8
mafft datasets/V6-V8.fasta > tree/V6-V8_aln.fasta
FastTree -spr 4 -gamma -fastest -no2nd -constraintWeight 100 -nt tree/V6-V8_aln.fasta > tree/V6-V8.tre
# V7-V9
mafft datasets/V7-V9.fasta > tree/V7-V9_aln.fasta
FastTree -spr 4 -gamma -fastest -no2nd -constraintWeight 100 -nt tree/V7-V9_aln.fasta > tree/V7-V9.tre