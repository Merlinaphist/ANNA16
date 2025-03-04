#!/bin/bash

# Function to reverse complement a DNA sequence
revcomp() {
    echo "$1" | tr "ACGTMRWSYKVHDBNacgtmrwsykvhdbn" "TGCAKYWSRMBDHVNtgcakysrwmbdhvn" | rev
}

# Function to display usage
usage() {
    echo "Usage: [-h] -i <input_file> -o <output_prefix> -s <start> -e <end> -t <tmp_dir>"
    echo "  -i  Input FASTA file"
    echo "  -o  Prefix of the output file"
    echo "  -s  Start of hypervariable region (options: V1, V3, V4, V6, V7)"
    echo "  -e  End of hypervariable region (options: V2, V3, V4, V5, V8, V9)"
    echo "  -t  Folder for storing temporary files"
    echo "  -c  Number of cores"
    exit 1
}

START_OPTIONS=("V1" "V3" "V4" "V6" "V7")
END_OPTIONS=("V2" "V3" "V4" "V5" "V8" "V9")

# Parse command-line arguments
while getopts ":hi:o:s:e:t:c:" opt; do
    case $opt in
        i) INPUT_FILE="$OPTARG" ;;
        o) OUTPUT_FILE="$OPTARG" ;;
        s) START="$OPTARG" ;;
        e) END="$OPTARG" ;;
        t) TMP_PATH="$OPTARG" ;;
        c) NUM_CORES="$OPTARG" ;;
        h) usage ;;
        *) usage ;;
    esac
done

# Check for required arguments
if [[ -z "$INPUT_FILE" || -z "$OUTPUT_FILE" || -z "$START" || -z "$END" ]]; then
    usage
fi

# Check if the region is valid
START_IN="FALSE"
for item in "${START_OPTIONS[@]}"
do
    if [ "$START" = "$item" ]; then
        START_IN="TRUE"
    fi
done


if [ "$START_IN" = "FALSE" ]; then
    echo "Error: Invalid starting region specified. Available options: ${START_OPTIONS[@]}"
    exit 1
fi

END_IN="FALSE"
for item in "${END_OPTIONS[@]}"
do
    if [ "$END" = "$item" ]; then
        END_IN="TRUE"
    fi
done


if [ "$END_IN" = "FALSE" ]; then
    echo "Error: Invalid ending region specified. Available options: ${END_OPTIONS[@]}"
    exit 1
fi

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

START_PRIMER_FWD=`cat ${SCRIPT_DIR}/primers/${START}_start_fwd.txt`
END_PRIMER_FWD=`cat ${SCRIPT_DIR}/primers/${END}_end_fwd.txt`
START_PRIMER_REV=`cat ${SCRIPT_DIR}/primers/${START}_start_rev.txt`
END_PRIMER_REV=`cat ${SCRIPT_DIR}/primers/${END}_end_rev.txt`

# Trim sequences using forward primers
cutadapt -g "$START_PRIMER_FWD"..."$END_PRIMER_FWD" \
         -e 0.2 -j $NUM_CORES \
         --untrimmed-output $TMP_PATH/tmp_fw_untrimmed.fasta -o $TMP_PATH/tmp_fw.fasta $INPUT_FILE

cutadapt -g "$END_PRIMER_REV"..."$START_PRIMER_REV" \
         -e 0.2 -j $NUM_CORES \
         --untrimmed-output ${OUTPUT_FILE}_unidentified.fasta -o $TMP_PATH/tmp_rc.fasta $TMP_PATH/tmp_fw_untrimmed.fasta

awk '
function revcomp(seq, rev, i, base, complement) {
    complement["A"] = "T"; complement["T"] = "A";
    complement["C"] = "G"; complement["G"] = "C";
    complement["M"] = "K"; complement["K"] = "M";
    complement["R"] = "Y"; complement["Y"] = "R";
    complement["W"] = "W"; complement["S"] = "S";
    complement["V"] = "B"; complement["B"] = "V";
    complement["H"] = "D"; complement["D"] = "H";
    complement["N"] = "N";

    rev = "";
    for (i = length(seq); i > 0; i--) {
        base = substr(seq, i, 1);
        rev = rev complement[toupper(base)];
    }
    return rev;
}

/^>/ { print; next }
{ print revcomp($0) }
' "$TMP_PATH/tmp_rc.fasta" > "$TMP_PATH/tmp_rc_rc.fasta"


cat $TMP_PATH/tmp_fw.fasta $TMP_PATH/tmp_rc_rc.fasta > ${OUTPUT_FILE}.fasta

rm $TMP_PATH/tmp_fw.fasta $TMP_PATH/tmp_rc_rc.fasta $TMP_PATH/tmp_rc.fasta