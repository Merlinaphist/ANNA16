#!/usr/bin/env python3
from anna16 import Preprocessing, CopyNumberPredictor
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



def predict_from_fasta(args):
    input_name, output_name = args
    input_file = open(input_name,"r")
    output_file = open(output_name, 'w')
    output_file.write("index,predicted_copy_number\n")
    tmp_seq = ""
    for line in input_file:
        if line[0] != ">":
            tmp_seq+=line.split("\n")[0]
        else:
            if tmp_seq != "":
                X = pp.CountKmers([tmp_seq])
                pred = model.predict(X)[0]
                output_file.write(f"{seqid},{pred}\n")
            tmp_seq = ""
            seqid = line.split("\n")[0][1:]       
    input_file.close()
    output_file.close()


for i in range(args.i[0]):
    input_name = args.i[0][i]
    if input_name.split(".")[-1] not in ["fasta","fna","fa"]:
        raise ValueError('Invalid file format. Expected formats are ["fasta","fna","fa"].')
    if args.o == None:
        suffix = input_name.split(".")[0]
        output_name = suffix+".csv"
    else:
        output_name = args.o[0][i]+".csv"
    predict_from_fasta([input_name, output_name])


# #Input FASTA files
# uploaded_seqs = {}
# for filename in args.i[0]:
#     if filename.split(".")[-1] not in ["fasta","fna","fa"]:
#         raise ValueError('Invalid file format. Expected formats are ["fasta","fna","fa"].')
#     else:
#         uploaded_seqs[filename] = pp.ReadFASTA(filename)

# #Run Predictions
# task_names = list(uploaded_seqs.keys())
# results = {}
# for task_name in task_names:
#     X = pp.CountKmers(uploaded_seqs[task_name])
#     results[task_name] = model.predict(X)
#     results[task_name] = pd.DataFrame(results[task_name],
#                                       index = uploaded_seqs[task_name].index,
#                                       columns = ["predicted_copy_number"])

# #Save Prediction Results
# for i in range(len(task_names)):
#     task_name = task_names[i]
#     if args.o == None:
#         suffix = filename.split(".")[-1]
#         downloaded_filename = task_name.split("."+suffix)[0]+".csv"
#     else:
#         downloaded_filename = args.o[0][i]+".csv"
#     values = results[task_name]["predicted_copy_number"].values.tolist()
#     indices = results[task_name].index.tolist()
#     with open(downloaded_filename, 'w') as f:
#         f.write('index,predicted_copy_number\n')
#         for i in range(len(values)):
#             line = str(indices[i])+","+str(values[i])+"\n"
#             f.write(line)