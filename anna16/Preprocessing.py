from sklearn.feature_extraction.text import CountVectorizer
import pandas as pd

class Preprocessing():
    def __init__(self, k_size=6):
        self.k_size = k_size
        kmers = self.ref_kmers("", self.k_size)
        self.vectorizer = CountVectorizer(vocabulary = kmers)
        self.seqs = []

    def ref_kmers(self, current_kmer, current_depth):
        if current_depth == 1:
            return [current_kmer+"a",current_kmer+"u",current_kmer+"c",current_kmer+"g"]
        else:
            ret = self.ref_kmers(current_kmer+"a",current_depth-1)
            for nt in ['u','c','g']:
                ret += self.ref_kmers(current_kmer+nt,current_depth-1)
            return ret

    def seq2kmer(self, seq, k):
        kmer = ""
        for i in range(0,len(seq)-k,1):
            kmer += seq[i:i+k]+" "
        return kmer[:-1]

    def CountKmers(self,seqs):
        if type(seqs) in [type([]),type(pd.core.series.Series([1]))]:
            kmer = pd.Series(seqs).apply(lambda x: self.seq2kmer(x, self.k_size))
            transformed_X = self.vectorizer.transform(kmer).toarray()
            return transformed_X
        else:
            raise ValueError("Invalid 'seqs' format. Expected formats are 'list' or 'pandas.core.series.Series'.")

    def ReadFASTA(self,filename,as_pd=True):
        if filename.split(".")[-1] not in ["fasta","fna","fa"]:
            raise ValueError('Invalid file format. Expected formats are ["fasta","fna","fa"].')
        file_handle = open(filename,"r")
        seqs = []
        seqid = []
        tmp_seq = ""
        for line in file_handle:
            if (line[0] == ">"):
                if tmp_seq != "":
                    seqs.append(tmp_seq)
                seqid.append(line.split("\n")[0][1:])
                tmp_seq = ""
            else:
                tmp_seq+=line.split("\n")[0]
        seqs.append(tmp_seq)
        file_handle.close()
        if as_pd:
            fasta = {}
            for i in range(len(seqs)):
                fasta[seqid[i]] = seqs[i]
            return pd.DataFrame(fasta,index=["sequence"]).transpose()["sequence"]
        else:
            return seqs, seqid