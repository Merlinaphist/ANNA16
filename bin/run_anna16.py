#!/usr/bin/env python3
from anna16 import Preprocessing, get_model #, predict_from_fasta
import pandas as pd
import numpy as np
import warnings
warnings.filterwarnings("ignore")
import argparse
import os
# from multiprocessing import Pool

parser = argparse.ArgumentParser()
parser.add_argument('-r', '--region', type=str, choices=["full_length", "V1-V2", "V1-V3", "V3-V4", "V4-V5", "V4", "V6-V8", "V7-V9"], help='Target Region', required=True)
parser.add_argument('-m', '--ml_type', type=str, choices=["cuml", "sklearn"], help='Machine Learning Framework Type', default="cuml")
parser.add_argument('-i', '--input', action='append', nargs='+', help='Input File', required=True)
parser.add_argument('-o', '--output', action='append', nargs='*', help='Output File')
# parser.add_argument('--streaming', help='Process the input file sequence by sequence', action="store_true")

args = parser.parse_args()

#Initialize the Model

model = get_model(ml_type=args.ml_type)
model.load(args.region)
pp = Preprocessing()

#Input FASTA files
uploaded_seqs = {}
for filename in args.i[0]:
    if filename.split(".")[-1] not in ["fasta","fna","fa"]:
        raise ValueError('Invalid file format. Expected formats are ["fasta","fna","fa"].')
    else:
        uploaded_seqs[filename] = pp.ReadFASTA(filename)

#Run Predictions
task_names = list(uploaded_seqs.keys())
results = {}
for task_name in task_names:
    X = pp.CountKmers(uploaded_seqs[task_name])
    results[task_name] = model.predict(X)
    results[task_name] = pd.DataFrame(results[task_name],
                                    index = uploaded_seqs[task_name].index,
                                    columns = ["predicted_copy_number"])

#Save Prediction Results
for i in range(len(task_names)):
    task_name = task_names[i]
    if args.o == None:
        suffix = filename.split(".")[-1]
        downloaded_filename = task_name.split("."+suffix)[0]+".csv"
    else:
        downloaded_filename = args.o[0][i]+".csv"
    output_file = results[[task_name]].to_csv(downloaded_filename, index=True)
