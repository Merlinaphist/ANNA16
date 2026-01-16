

## Table of Contents

- [User Guide](#user_guide)
    - [Installation](#installation-)
    - [Preprocessing](#preprocessing-)
    - [Usage](#usage-)
- [About ANNA16](#about)
- [Cite ANNA16](#citation)

# User Guide <a name="user_guide"></a>
The current page shows ANNA16 (Version 2.0). To install the older versions or learn about the updates in ANNA16, please visit [the wiki](https://github.com/Merlinaphist/ANNA16/wiki).

## Installation <a name="installation"></a>

We recommend create a separate environment to host ANNA16:

```bash
wget https://github.com/Merlinaphist/ANNA16/archive/refs/tags/v1.1.0.zip
unzip v1.1.0.zip
cd ANNA16-1.1.0
conda create -n anna16
conda activate anna16
conda install python=3.9 cutadapt=5.0
pip install scikit-learn==1.1.2 tensorflow[and-cuda]==2.17.0 numpy==1.26.4 pandas==2.0.2
pip install -e .
chmod +x "`pwd`/bin/extract_regions.sh"
chmod +x "`pwd`/bin/run_anna16.py"
echo "export PATH=$PATH:`pwd`/bin" >> ~/.bashrc
source ~/.bashrc
```

## Preprocessing <a name="preprocessing"></a>

The input to ANNA16 needs to be:

1. FASTA, FNA, or FA files of **DNA** sequences;
2. **Positive** strands sequences;
3. (Recommended) Trimmed by our `extract_regions.sh` script. This script calls [cutadapt](https://cutadapt.readthedocs.io/en/stable/) and the [Primer Table](#primers-) to identify the regions of interest.

```bash
extract_regions.sh -i <input_file> -o <output_prefix> \
                   -s <start> -e <end> -t <tmp_dir> -c <num_core>
```
`-i` - Input FASTA file

`-o` - Prefix of the output file

`-s` - Start of hypervariable region (options: V1, V3, V4, V6, V7)

`-e` - End of hypervariable region (options: V2, V3, V4, V5, V8, V9)

`-t` - Folder for storing temporary files

`-c` - Number of cores

An example command is:

```bash
mkdir tmp
extract_regions.sh -i raw_data/input.fasta \
                   -o intermediate_data/trimmed_full_length \
                   -s V1 -e V9 -t tmp -c 4
```

## Usage <a name="usage"></a>

### ANNA16 as a Python Library

```python
from anna16 import Preprocessing, CopyNumberPredictor
region = "V1-V2" #Options: [full_length, V1-V2, V1-V3, V3-V4, V4-V5, V4, V6-V8, V7-V9]
pp = Preprocessing()
seqs = pp.ReadFASTA(filename)
kmer_counts = pp.CountKmers(seqs)
model = CopyNumberPredictor(region=region) 
model.load(trimmed=True)
copy_number_pred = model.predict(kmer_counts)
```

### ANNA16 as a Command-Line Tool:

```bash
run_anna16.py -r <REGION> -t <TRIM> -i <INPUT_FILE(S)> -o <OUTPUT_FILE(S)>
```

**Required Parameters:**

`-r` - Region on 16S rRNA gene that the input sequences belong to. Options: `[full_length, V1-V2, V1-V3, V3-V4, V4-V5, V4, V6-V8, V7-V9]`

`-t` - Whether or not the primers of the amplicons have been trimmed off. 
Options: `[True, False]`

`-i` - A list of files of the input sequences.

**Optional Parameter:**

`-o` - The filename of the output. Output will be saved in `csv` format. If not provided, output file(s) will take the same name(s) with the input.

An example command is:

```bash
run_anna16.py -r full_length -t True -i input0.fasta input1.fasta -o pred0 pred1
```

# About ANNA16 <a name="about"></a>

ANNA16 is an end-to-end tool that predicts 16S rRNA gene copy number (GCN) from 16S rRNA gene sequence. The tool utilizes an ensembled architecture of Multi-layer Perceptron (MLP), Support Vector Machine (SVM), and Ridge Regression. This repository releases the model weights of ANNA16.

![Summary of ANNA16](assets/ANNA16_summary.png)

### Primer Table <a name="primers"></a>
The regions used in training were obtained from the following primers: 

| Region | Forward Primer Name | Forward Primer Sequence | Reverse Primer Name |Reverse Primer |
|-------:|--------------------:|------------------------:|-----------:|--------------:|
| Full Length | 27F | AGA GTT TGA TCC TGG CTC AG     | 1492R | TAC GGY TAC CTT GTT ACG ACT T    |
| V1-V2       | 27F | AGA GTT TGA TCC TGG CTC AG     | 338R | GCT GCC TCC CGT AGG AGT         |
| V1-V3       | 27F | AGA GTT TGA TCC TGG CTC AG     | 534R | ATT ACC GCG GCT GCT GG          |
| V3-V4       | 341F | CCT ACG GGA GGC AGC AG         | 785R | GAC TAC HVG GGT ATC TAA TCC     |
| V4          | 515F | GTG CCA GCM GCC GCG GTA A      | 806R | GGA CTA CHV GGG TWT CTA AT      |
| V4-V5       | 515F | GTG CCA GCM GCC GCG GTA A      | 926R | CCG YCA ATT YMT TTR AGT TT      |
| V6-V8       | 939F | GAA TTG ACG GGG GCC CGC ACA AG | 1378R | CGG TGT GTA CAA GGC CCG GGA ACG |
| V7-V9       | 1115F | CAA CGA GCG CAA CCC T          | 1492R | TAC GGY TAC CTT GTT ACG ACT T    |



# Cite ANNA16 <a name="citation"></a>

Miao, J., Chen, T., Misir, M., & Lin, Y. (2024). Deep Learning for Predicting 16S rRNA Gene Copy Number. *Scientific reports*, 14(14282). https://doi.org/10.1038/s41598-024-64658-5 
