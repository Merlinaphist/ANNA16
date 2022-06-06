from tkinter import *
from tkinter import filedialog

window = Tk()
window.title('16S Gene Copy Number Predictor')
window.geometry("300x300")
window.config(background = "white")

import pandas as pd
from tensorflow.keras.models import load_model
from sklearn.feature_extraction.text import CountVectorizer


class Predict():
    def __init__(self):
        self.input_filename = ''
        self.output_filename = ''
        self.seq = []
        self.id = []
        self.res = []
        self.response = None
        self.model = load_model("mlp.h5")
        self.vocal = []
        file_handle = open("vector_vocal.csv","r")
        for line in file_handle:
            self.vocal.append(line.split("\n")[0])
        file_handle.close()

    def browseFiles(self):
        self.input_filename = filedialog.askopenfilename(initialdir = "/",title = "Select a File",filetypes = (("FASTA files","*.fasta"),("all files","*.*")))
        file_handle = open(self.input_filename,"r")
        tmp_seq = ""
        for line in file_handle:
            if (line[0] == ">"):
                if tmp_seq != "":
                    self.seq.append(tmp_seq)
                self.id.append(line.split("\n")[0])
                tmp_seq = ""
            else:
                tmp_seq+=line.split("\n")[0]
        file_handle.close()
        global window
        self.response = Label(window,text = "Input Successful.")
        self.response.grid(column = 99,row = 70)

    def vectorize(self,X):
        vectorizer = CountVectorizer(vocabulary = self.vocal)
        kmer_str = self.generate_kmer_multiple(X, 6)
        return vectorizer.fit_transform(kmer_str).toarray()

    def run(self):
        global window
        self.response.destroy()
        self.response = Label(window,text = "Running Data Now...")
        self.response.grid(column = 99,row = 70)
        X = self.vectorize(self.seq)
        Y = self.model.predict(X)
        for i in range(0,len(Y),1):
            self.res.append([str(self.id[i])[1:],self.seq[i],str(Y[i])[1:-1]])
        self.res = pd.DataFrame(self.res)
        self.response.destroy()
        self.response = Label(window,text = "Prediction Finished.")
        self.response.grid(column = 99,row = 70)

    def save(self):
        self.output_filename = filedialog.asksaveasfilename(initialdir = "/",title = "SAVE",filetypes = (("CSV files","*.csv"),("all files","*.*")))
        self.res.to_csv(path_or_buf = self.output_filename, index = False,header = ["id","sequence","copy_number"])
        global window
        self.response.destroy()
        self.response = Label(window,text = "Save Successful.")
        self.response.grid(column = 99,row = 70)

    def generate_kmer_multiple(self,seqlist,k):
        kmer_list = []
        for seq in seqlist:
            kmer_list.append(self.generate_kmer_single(seq,k))
        return kmer_list
    
    def generate_kmer_single(self,seq,k):
        kmer = ""
        for i in range(0,len(seq)-k,1):
            kmer += seq[i:i+k]+" "
        return kmer[:-1]
																							
main = Predict()

button_input = Button(window,text = "Input",command = main.browseFiles)
button_input.grid(column = 99, row = 20)

button_save = Button(window,text = "Run",command = main.run)
button_save.grid(column = 99,row = 30)

button_save = Button(window,text = "Save",command = main.save)
button_save.grid(column = 99,row = 40)

button_exit = Button(window,text = "Exit",command = exit)
button_exit.grid(column = 99,row = 50)

version = Label(window,text = "Version: 6mer+CountVectorize+MLP")
version.grid(column = 99,row = 80)

window.mainloop()
