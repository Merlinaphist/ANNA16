#!/usr/bin/env python3
from anna16 import Preprocessing, get_model  # , predict_from_fasta
import pandas as pd
import numpy as np
import warnings
warnings.filterwarnings("ignore")
import argparse
import os
# from multiprocessing import Pool

parser = argparse.ArgumentParser()
parser.add_argument('-r', '--region', type=str,
                    choices=["full_length", "V1-V2", "V1-V3", "V3-V4",
                             "V4-V5", "V4", "V6-V8", "V7-V9"],
                    help='Target Region', required=True, dest='region')
parser.add_argument('-m', '--ml_type', type=str,
                    choices=["cuml", "sklearn"],
                    help='Machine Learning Framework Type', default="cuml", dest='ml_type')
parser.add_argument('-i', '--input', type=str, nargs='+',
                    help='Input FASTA file(s)', required=True, dest='input')
parser.add_argument('-o', '--output', type=str, nargs='*',
                    help='Output file name(s) without extension', dest='output')
# parser.add_argument('--streaming', help='Process the input file sequence by sequence', action="store_true")

args = parser.parse_args()

# Validate output count if provided
if args.output is not None and len(args.output) != len(args.input):
    parser.error(
        f"Number of -o/--output names ({len(args.output)}) must match "
        f"number of -i/--input files ({len(args.input)})."
    )

# Initialize the Model
model = get_model(ml_type=args.ml_type)
model.load(args.region)
pp = Preprocessing()

# Input FASTA files
uploaded_seqs = {}
for filename in args.input:
    if filename.split(".")[-1] not in ["fasta", "fna", "fa"]:
        raise ValueError(
            'Invalid file format. Expected formats are ["fasta","fna","fa"].'
        )
    uploaded_seqs[filename] = pp.ReadFASTA(filename)

# Run Predictions
task_names = list(uploaded_seqs.keys())
results = {}
for task_name in task_names:
    X = pp.CountKmers(uploaded_seqs[task_name]['sequence'])
    results[task_name] = model.predict(X)
    results[task_name] = pd.DataFrame(
        results[task_name],
        index=uploaded_seqs[task_name]['seqid'],
        columns=["predicted_copy_number"],
    )

# Save Prediction Results
for i, task_name in enumerate(task_names):
    if args.output is None:
        suffix = task_name.split(".")[-1]
        downloaded_filename = task_name[:-(len(suffix) + 1)] + ".csv"
    else:
        downloaded_filename = args.output[i] + ".csv"
    results[task_name].to_csv(downloaded_filename, index=True)