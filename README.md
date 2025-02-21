[Work in Progress]
# About ANNA16

- [Update Log](#updates-)
- [User Guide](#user_guide-)
    - [Colab Notebook](#colab-)
    - [Local Device](#local_device-)
        - [Installation](#installation-)
        - [Preprocessing](#preprocessing-)
        - [Usage](#usage-)
- [Cite ANNA16](#citation-)

ANNA16 is an end-to-end tool that predicts 16S rRNA gene copy number (GCN) from 16S rRNA gene sequence. The tool utilizes an ensembled architecture of Multi-layer Perceptron (MLP), Support Vector Machine (SVM), and Ridge Regression. This repository releases the model weights of ANNA16.

![Summary of ANNA16](assets/ANNA16_summary.png)


# Update Log <a name="updates"></a>

`v1.1`:
1. Change the primer sequence of `1492R` to the more common format: `TAC GGY TAC CTT GTT ACG ACT T` (one additional `T`at the end). Re-train the model weights for trimmed full-length and V7-V9.
2. Add `extract_regions.sh` for sequence preprocessing. This script extracts the 16S rRNA full-length or subregions from the input FASTA file and unifies the sequence orientation.
3. Update `tensorflow` version requirements to `2.17.0` to accommodate the development of new GPUs and CUDA framework.
4. `anna16` is now a Python package that can be imported


# User Guide <a name="user_guide"></a>

## Colab Notebook <a name="colab"></a>

ANNA16 can be run on Google Colab. Please visit https://colab.research.google.com/drive/1XwpTMCHSfTmzpHyKrmiD8aC8C_1nndUV#scrollTo=V_Fe5x8g8sCv

## Local Device <a name="local_device"></a>

Alternatively, ANNA16 can be installed on a local device.

### Installation <a name="installation"></a>

We recommend create a separate environment to host ANNA16:

```bash
conda env create -f anna16_env.yaml
conda activate anna16
git clone https://github.com/Merlinaphist/ANNA16.git
cd ANNA16
pip install -e .
chmod +x "`pwd`/bin/extract_regions.sh"
chmod +x "`pwd`/bin/run_anna16.py"
echo "export PATH=$PATH:`pwd`/bin" >> ~/.bashrc
source ~/.bashrc
```

### Preprocessing <a name="preprocessing"></a>

The input to ANNA16 needs to be:

1. FASTA, FNA, or FA files of **DNA** sequences;
2. **Positive** strands sequences;
3. Trimmed with the following primers:

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

We recommend using the `extract_regions.sh` script to preprocess the input files.

```bash
extract_regions.sh -i <input_file> -o <output_prefix> -s <start> -e <end> -t <tmp_dir>
```
`-i` - Input FASTA file

`-o` - Prefix of the output file

`-s` - Start of hypervariable region (options: V1, V3, V4, V6, V7)

`-e` - End of hypervariable region (options: V2, V3, V4, V5, V8, V9)

`-t` - Folder for storing temporary files

An example command is:

```bash
mkdir tmp
extract_regions.sh -i raw_data/input.fasta -o intermediate_data/trimmed_full_length -s V1 -e V9 -t tmp
```

### Usage <a name="usage"></a>

#### To use ANNA16 in a `py` or `ipynb` script:

```python
from anna16 import Preprocessing, CopyNumberPredictor
region = "V1-V2" #Options: [full_length, V1-V2, V1-V3, V3-V4, V4-V5, V4, V6-V8, V7-V9]
pp = Preprocessing()
seqs = pp.ReadFASTA(filename)
kmer_counts = pp.CountKmers(seqs)
model = CopyNumberPredictor(region=region) 
model.load(filename=f"model_files/trimmed/{region}.zip")
copy_number_pred = model.predict(kmer_counts)
```

#### To use ANNA16 as a command-line tool:

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

# Cite ANNA16 <a name="citation"></a>

Miao, J., Chen, T., Misir, M., & Lin, Y. (2024). Deep Learning for Predicting 16S rRNA Gene Copy Number. *Scientific reports*, 14(14282). https://doi.org/10.1038/s41598-024-64658-5 
