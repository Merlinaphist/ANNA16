#!/usr/bin/env python3
import pandas as pd
from Bio import SeqIO
from sklearn.feature_extraction.text import CountVectorizer
import warnings
import torch
import torch.nn as nn
import torch.optim as optim
from torch.utils.data import DataLoader
import numpy as np
import pandas as pd
import pickle, os
from multiprocessing import Pool, cpu_count


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

    def seq2kmer(self, seq):
        kmer = ""
        for i in range(0,len(seq)-self.k_size,1):
            kmer += seq[i:i+self.k_size]+" "
        return kmer[:-1]

    def CountKmers(self,seqs):
        if type(seqs) in [type([]),type(pd.core.series.Series([1]))]:
            with Pool(cpu_count()) as p:
                kmer = p.map(self.seq2kmer, seqs)
            # kmer = pd.Series(seqs).apply(lambda x: self.seq2kmer(x, self.k_size))
            transformed_X = self.vectorizer.transform(kmer).toarray()
            return transformed_X
        else:
            raise ValueError("Invalid 'seqs' format. Expected formats are 'list' or 'pandas.core.series.Series'.")

    def ReadFASTA(self,filename):
        if filename.split(".")[-1] not in ["fasta","fna","fa"]:
            raise ValueError('Invalid file format. Expected formats are ["fasta","fna","fa"].')
        records = list(SeqIO.parse(filename, "fasta"))
        data = []
        for record in records:
            name_parts = record.description
            data.append([name_parts, str(record.seq)])
        return pd.DataFrame(data, columns=["seqid", "sequence"])



class MLP(nn.Module):
    """PyTorch MLP model equivalent to the Keras Sequential model"""
    def __init__(self, input_dim, hidden_config, output_dim, device='cpu', batch_size=100):
        super(MLP, self).__init__()
        layers = []
        prev_dim = input_dim
        for dim, name in hidden_config:
            layers.extend(self.get_layer(prev_dim, dim, name))
            prev_dim = dim
        layers.extend([
            nn.Linear(prev_dim, output_dim),
            nn.Softplus()
        ])
        self.layers = nn.Sequential(*layers)
        self.device = device
        self.batch_size = batch_size
    
    def get_layer(self, prev_dim, dim, name):
        if name == "relu":
            return [nn.Linear(prev_dim, dim), nn.ReLU()]
        elif name == "gelu":
            return [nn.Linear(prev_dim, dim), nn.GELU()]
        elif name == "selu":
            return [nn.Linear(prev_dim, dim), nn.SELU()]
        elif name == "elu":
            return [nn.Linear(prev_dim, dim), nn.ELU()]
        elif name == "linear":
            return [nn.Linear(prev_dim, dim)]
        else:
            raise ValueError(f"Invalid activation function: {name}")
        
    def forward(self, x):
        return self.layers(x)

    def predict(self, x):
        self.eval()
        outputs = []
        dataloader = DataLoader(RRNDataset(x, device=self.device), batch_size=self.batch_size, shuffle=False)
        with torch.no_grad():
            for batch_X in dataloader:
                outputs.append(self.forward(batch_X))
        return torch.cat(outputs, dim=0).cpu().numpy()
    
    def save(self, filename):
        torch.save(self.state_dict(), filename)
    
    def load(self, filename, map_location="cpu"):
        self.load_state_dict(torch.load(filename, map_location=map_location))
        return self
    
class RMSELoss(nn.Module):
    def __init__(self):
        super().__init__()
        self.mse = nn.MSELoss()
    def forward(self,yhat,y):
        return torch.sqrt(self.mse(yhat,y))

class RRNDataset():
    def __init__(self, x, y=None, device="cpu"):
        self.x = torch.tensor(x, dtype=torch.float32)
        if y is not None:
            self.y = torch.tensor(y, dtype=torch.float32)
        else:
            self.y = None
        self.device = device
    def __len__(self):
        return len(self.x)
    def __getitem__(self, idx):
        if self.y is not None:
            return self.x[idx].to(self.device), self.y[idx].to(self.device)
        else:
            return self.x[idx].to(self.device)


default_mlp_params = {
    # "hidden_dims": [489, 926, 645, 929, 582, 82],
    # "activations": ["gelu", "relu", "relu", "elu", 'gelu', "linear"],
    'hidden_config':[[489, "gelu"], [926, "relu"], [645, "relu"], [929, "elu"], [582, "gelu"], [82, "linear"]],
    "learning_rate": 0.00016770313599,
    "epochs": 59,
    "batch_size": 100
}

default_ml_params = {
    "ridge_alpha": 49,
    "pca_n_components": 100,
    "svr_kernel": "rbf",
    "svr_C": 11,
    "svr_gamma": "scale"
}


class ANNA16Model():
    def __init__(self, pca, svr, ridge, device="cpu", mlp_params=default_mlp_params):
        self.hidden_config = mlp_params["hidden_config"]
        self.learning_rate = mlp_params["learning_rate"]
        self.epochs = mlp_params["epochs"]
        self.batch_size = mlp_params["batch_size"]
        self.device = device
        self.mlp = MLP(4096, self.hidden_config, 1, device=self.device, batch_size=self.batch_size)
        self.mlp.to(self.device)
        self.ml_models = {
            "pca": pca,
            "svr": svr,
            "ridge": ridge
        }
        self.loss_func = RMSELoss()
    
    def fit_mlp(self, X_train, Y_train, X_val=None, Y_val=None):
        """Train the PyTorch MLP model"""
        train_dataloader = DataLoader(RRNDataset(X_train, Y_train, device=self.device), batch_size=self.batch_size, shuffle=True)
        if X_val is not None and Y_val is not None:
            val_dataloader = DataLoader(RRNDataset(X_val, Y_val, device=self.device), batch_size=self.batch_size, shuffle=False)
        else:
            val_dataloader = None
        # Optimizer and loss
        optimizer = optim.Adam(self.mlp.parameters(), lr=self.learning_rate)
        # Training loop
        loss_history = {"train_loss": []} if val_dataloader is None else {"train_loss": [], "val_loss": []}
        for epoch in range(self.epochs):
            train_epoch_loss = 0.0
            self.mlp.train()
            for batch_X, batch_Y in train_dataloader:
                optimizer.zero_grad()
                outputs = self.mlp(batch_X)
                loss = self.loss_func(outputs, batch_Y)
                loss.backward()
                optimizer.step()
                train_epoch_loss += loss.item()
            loss_history["train_loss"].append(train_epoch_loss / len(train_dataloader))
            # Validation (optional, for monitoring)
            if val_dataloader is not None:
                self.mlp.eval()
                val_epoch_loss = 0.0
                with torch.no_grad():
                    for batch_X, batch_Y in val_dataloader:
                        val_outputs = self.mlp(batch_X)
                        val_loss = self.loss_func(val_outputs, batch_Y).item()
                        val_epoch_loss += val_loss
                loss_history["val_loss"].append(val_epoch_loss / len(val_dataloader))
        return loss_history
    
    def fit(self, X_train, Y_train, X_val=None, Y_val=None, verbose=True, return_history=False):
        self.echo(text="------Training Starts------", verbose=verbose)
        
        X_train_pca = self.ml_models["pca"].fit_transform(X_train)
        self.ml_models["svr"].fit(X_train_pca, Y_train)
        svr_pred = self.ml_models["svr"].predict(X_train_pca)
        
        self.echo(text="SVR done.", verbose=verbose)

        loss_history = self.fit_mlp(X_train, Y_train, X_val, Y_val)
        mlp_pred = self.mlp.predict(X_train)
        self.echo(text="MLP done.", verbose=verbose)
        
        new_X_train = np.concatenate([X_train_pca, mlp_pred, svr_pred.reshape(-1,1)], axis=1)
        
        self.ml_models["ridge"].fit(new_X_train, Y_train)
        self.echo(text="Meta-Model: Ridge done.", verbose=verbose)
        if return_history:
            return loss_history

    def predict(self, X_test):
        X_test_pca = self.ml_models["pca"].transform(X_test)
        mlp_pred = self.mlp.predict(X_test)
        svr_pred = self.ml_models["svr"].predict(X_test_pca)
        new_X_test = np.concatenate([X_test_pca, mlp_pred, svr_pred.reshape(-1,1)], axis=1)
        final_pred = self.ml_models["ridge"].predict(new_X_test)
        return final_pred
    
    def echo(self, text, verbose):
        if verbose not in [True, False]:
            raise ValueError('verbose must be True or False')
        if verbose:
            print(text)


class ANNA16ModelSKLEARN(ANNA16Model):
    def __init__(self, pca, svr, ridge, device="cpu", mlp_params=default_mlp_params):
        super().__init__(pca, svr, ridge, device=device, mlp_params=mlp_params)
    
    def save(self, file_prefix):
        self.mlp.save(f"{file_prefix}_mlp.pt")
        for model_name in self.ml_models:
            with open(f"{file_prefix}_{model_name}_skl.pkl", 'wb') as file:
                pickle.dump(self.ml_models[model_name], file)
    
    def load(self, region):
        current_dir = os.path.dirname(os.path.abspath(__file__))
        file_prefix=f"{current_dir}/model_files/{region}/{region}"
        self.mlp.load(f"{file_prefix}_mlp.pt", map_location=self.device)
        for model_name in self.ml_models:
            with open(f"{file_prefix}_{model_name}_skl.pkl", 'rb') as file:
                self.ml_models[model_name] = pickle.load(file)
        return self


class ANNA16ModelCUML(ANNA16Model):
    def __init__(self, pca, svr, ridge, device="cpu", mlp_params=default_mlp_params):
        super().__init__(pca, svr, ridge, device=device, mlp_params=mlp_params)

    def save(self, file_prefix, as_sklearn=True):
        if as_sklearn not in [True, False]:
            raise ValueError('as_sklearn must be True or False')
        self.mlp.save(f"{file_prefix}_mlp.pt")
        if not as_sklearn:
            for model_name in self.ml_models:
                with open(f"{file_prefix}_{model_name}_cuml.pkl", 'wb') as file:
                    pickle.dump(self.ml_models[model_name], file)
        else:
            for model_name in self.ml_models:
                with open(f"{file_prefix}_{model_name}_skl.pkl", 'wb') as file:
                    pickle.dump(self.ml_models[model_name].as_sklearn(), file)

    def load(self, region, from_sklearn=True):
        current_dir = os.path.dirname(os.path.abspath(__file__))
        file_prefix=f"{current_dir}/model_files/{region}/{region}"
        if from_sklearn not in [True, False]:
            raise ValueError('from_sklearn must be True or False')
        self.mlp.load(f"{file_prefix}_mlp.pt")
        if from_sklearn:
            for model_name in self.ml_models:
                with open(f"{file_prefix}_{model_name}_skl.pkl", 'rb') as file:
                    self.ml_models[model_name] = self.ml_models[model_name].from_sklearn(pickle.load(file))
        else:
            for model_name in self.ml_models:
                with open(f"{file_prefix}_{model_name}_cuml.pkl", 'rb') as file:
                    self.ml_models[model_name] = pickle.load(file)



def get_model(ml_type="cuml", device="cpu", mlp_params=default_mlp_params, ml_params=default_ml_params):
    if ml_type == "cuml":
        from cuml import Ridge, PCA
        from cuml.svm import SVR
        pca = PCA(n_components=ml_params["pca_n_components"])
        svr = SVR(kernel=ml_params["svr_kernel"], C=ml_params["svr_C"], gamma=ml_params["svr_gamma"])
        ridge = Ridge(alpha=ml_params["ridge_alpha"])
        return ANNA16ModelCUML(pca, svr, ridge, device=device, mlp_params=mlp_params)
    elif ml_type == "sklearn":
        from sklearn.linear_model import Ridge
        from sklearn.decomposition import PCA
        from sklearn.svm import SVR
        pca = PCA(n_components=ml_params["pca_n_components"])
        svr = SVR(kernel=ml_params["svr_kernel"], C=ml_params["svr_C"], gamma=ml_params["svr_gamma"])
        ridge = Ridge(alpha=ml_params["ridge_alpha"])
        return ANNA16ModelSKLEARN(pca, svr, ridge, device=device, mlp_params=mlp_params)
    else:
        raise ValueError(f"Invalid ML type: {ml_type}")
