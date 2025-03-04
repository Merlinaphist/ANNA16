#!/usr/bin/env python3
from anna16 import Preprocessing, CopyNumberPredictor, predict_from_fasta
import pandas as pd
import numpy as np
import warnings
warnings.filterwarnings("ignore")
import argparse
import os
# from multiprocessing import Pool

parser = argparse.ArgumentParser()
parser.add_argument('-r', type=str, choices=["full_length", "V1-V2", "V1-V3", "V3-V4", "V4-V5", "V4", "V6-V8", "V7-V9"], help='Target Region', required=True)
parser.add_argument('-t', type=str, choices=["True", "False"], help='Primers Trimmed', required=True)
parser.add_argument('-i', action='append', nargs='+', help='Input File', required=True)
parser.add_argument('-o', action='append', nargs='*', help='Output File')

args = parser.parse_args()

#Initialize the Model
model = CopyNumberPredictor(region=args.r)
model.load(args.t == "True")
pp = Preprocessing()

for i in range(len(args.i[0])):
    input_name = args.i[0][i]
    if input_name.split(".")[-1] not in ["fasta","fna","fa"]:
        raise ValueError('Invalid file format. Expected formats are ["fasta","fna","fa"].')
    if args.o == None:
        suffix = input_name.split(".")[0]
        output_name = suffix+".csv"
    else:
        output_name = args.o[0][i]+".csv"
    predict_from_fasta(input_name, output_name, model)