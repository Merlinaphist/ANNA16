# About ANNA16
ANNA16 is an end-to-end tool that predicts 16S rRNA gene copy number (GCN) from 16S rRNA gene sequence. The tool utilizes an ensembled architecture of Multi-layer Perceptron (MLP), Support Vector Machine (SVM), and Ridge Regression. This repository releases the model weights of ANNA16.

![Summary of ANNA16](ANNA16_summary.png)

# Use ANNA16

## Data Preprocessing

+ Prepare FASTA, FNA, or FA files of **DNA** sequences that encode your target 16S rRNA genes.
+ The sequences must be **positive** strands.
+ We recommend trimming the sequences with the primers listed in the following table:

| Region | Forward Primer Name | Forward Primer Sequence | Reverse Primer Name |Reverse Primer |
|-------:|--------------------:|------------------------:|-----------:|--------------:|
| Full Length | 27F | AGA GTT TGA TCC TGG CTC AG     | 1492R | TAC GGY TAC CTT GTT ACG ACT     |
| V1-V2       | 27F | AGA GTT TGA TCC TGG CTC AG     | 338R | GCT GCC TCC CGT AGG AGT         |
| V1-V3       | 27F | AGA GTT TGA TCC TGG CTC AG     | 534R | ATT ACC GCG GCT GCT GG          |
| V3-V4       | 341F | CCT ACG GGA GGC AGC AG         | 785R | GAC TAC HVG GGT ATC TAA TCC     |
| V4          | 515F | GTG CCA GCM GCC GCG GTA A      | 806R | GGA CTA CHV GGG TWT CTA AT      |
| V4-V5       | 515F | GTG CCA GCM GCC GCG GTA A      | 926R | CCG YCA ATT YMT TTR AGT TT      |
| V6-V8       | 939F | GAA TTG ACG GGG GCC CGC ACA AG | 1378R | CGG TGT GTA CAA GGC CCG GGA ACG |
| V7-V9       | 1115F | CAA CGA GCG CAA CCC T          | 1492R | TAC GGY TAC CTT GTT ACG ACT     |

## Colab Notebook

ANNA16 can be run on Google Colab. Please visit https://colab.research.google.com/drive/1XwpTMCHSfTmzpHyKrmiD8aC8C_1nndUV#scrollTo=V_Fe5x8g8sCv.

## Local Server

Alternatively, ANNA16 can be installed on a local server.

**(1) Download ANNA16**

```bash
git clone https://github.com/Merlinaphist/ANNA16.git
chmod +x /path/to/ANNA16/run_anna16.py
export PATH=$PATH:/path/to/ANNA16/
```

**(2) Python Dependencies**

```bash
pickle
numpy
pandas==2.7.3
zipfile
scikit-learn==1.1.2
tensorflow==2.9.0
```

We recommend create an environment for ANNA16 to host the Python dependencies:

```bash
conda create --name anna16 python=3.8
conda activate anna16
pip install -r ANNA16/requirements.txt
```

**(3) Run ANNA16**

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
./run_anna16.py -r full_length -t True -i input0.fasta input1.fasta -o pred0 pred1
```

# Cite ANNA16:

Miao, J., Chen, T., Misir, M., & Lin, Y. (2022). Deep Learning for Predicting 16S rRNA Copy Number. *bioRxiv*, 2022.2011.2026.518038. https://doi.org/10.1101/2022.11.26.518038
