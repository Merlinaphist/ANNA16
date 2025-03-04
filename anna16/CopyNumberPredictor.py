import os
import shutil
import pickle
import numpy as np
import pandas as pd
from zipfile import ZipFile
from sklearn.svm import SVR
from sklearn.decomposition import PCA
from sklearn.linear_model import Ridge
import tensorflow as tf
from tensorflow import cast,float32
from tensorflow.keras.layers import Dense
from tensorflow.keras.optimizers import Adam
from tensorflow.keras.models import Sequential, load_model

class CopyNumberPredictor():
    def __init__(self,region):
        self.region = region
        self.state = [59, 0.00016770313599, 0.11, 100, 489, 1, 926, 0, 645, 0, 929, 3, 582, 1, 82, 4]
        self.activation_indices={0:"relu",1:"gelu",2:"selu",3:"elu",4:"linear"}
        self.mlp = self.create_mlp()
        self.ridge = Ridge(alpha = 49)
        self.pca = PCA(n_components=100)
        self.svr = SVR(kernel='rbf',C=11,gamma='auto')
        
    def save(self,filename):
        if filename[-4:] != ".zip":
            raise ValueError('Invalid filename. Expect a zip file.')
            
        path = filename[:-4]
        prefix = path + "/" + self.region
            
        if not os.path.exists(path):
            os.makedirs(path)
        
        with open(prefix+'_pca.pkl','wb') as file:
            pickle.dump(self.pca,file)
        with open(prefix+'_ridge.pkl','wb') as file:
            pickle.dump(self.ridge,file)
        with open(prefix+'_svr.pkl','wb') as file:
            pickle.dump(self.svr,file)
        self.mlp.save(prefix+"_mlp.h5")
        
        shutil.make_archive(path, 'zip', path)
        shutil.rmtree(path)
            
    def load(self,trimmed=True):
        # if filename[-4:] != ".zip":
        #     raise ValueError('Invalid input file. Expect a zip file.')
            
        # path = filename[:-4]
        
        # if not os.path.exists(path):
        #     os.makedirs(path)
        trim_or_not={True:"trimmed", False:"untrimmed"}

        path = f"{__file__.split("CopyNumerPredictor.py")[0]}model_files/{trim_or_not[trimmed]}"
        filename = f"{path}/{self.region}.zip"

        with ZipFile(filename,'r') as zObject:
            zObject.extractall(path=path)
        
        prefix = path + "/" + self.region
        with open(prefix+'_pca.pkl', 'rb') as file:
            self.pca = pickle.load(file) 
        with open(prefix+'_ridge.pkl', 'rb') as file:
            self.ridge = pickle.load(file)
        with open(prefix+'_svr.pkl', 'rb') as file:
            self.svr = pickle.load(file)
        self.mlp = load_model(prefix + "_mlp.h5",
                              custom_objects={"root_mean_squared_error": self.root_mean_squared_error})
        shutil.rmtree(path)
        
    def fit(self,X_train,Y_train,verbose=True):
        self.echo(text = "------Training Starts------", verbose = verbose)
        X_train_pca = self.pca.fit_transform(X_train)
        self.fit_mlp(X_train,Y_train)
        mlp_pred = self.mlp.predict(X_train,verbose=0)
        self.echo(text = "Model 1: MLP done.", verbose = verbose)
        
        reshaped_Y_train = Y_train.values.reshape(Y_train.shape[0])
        new_X_train = pd.concat([pd.DataFrame(X_train_pca),pd.DataFrame(mlp_pred)], axis = 1)
        
        signals = ["Model 2: SVR done."]
        for model in [self.svr]:
            model.fit(X_train_pca,reshaped_Y_train)
            pred = model.predict(X_train_pca)
            new_X_train = pd.concat([new_X_train, pd.DataFrame(pred)], axis = 1)
            self.echo(text = signals[0], verbose = verbose)
            del signals[0]
        
        self.ridge.fit(new_X_train,reshaped_Y_train)
        self.echo(text = "Meta-Model: Ridge done.", verbose = verbose)
    
    def echo(self,text,verbose):
        if verbose not in [True,False]:
            raise ValueError('verbose must be True or False')
        if verbose:
            print(text)
        
    def predict(self, X_test):
        X_test_pca = self.pca.transform(X_test)
        mlp_pred = self.mlp.predict(X_test,verbose=0)
        new_X_test = pd.concat([pd.DataFrame(X_test_pca),pd.DataFrame(mlp_pred)], axis = 1)
        
        for model in [self.svr]:
            pred = model.predict(X_test_pca)
            new_X_test = pd.concat([new_X_test, pd.DataFrame(pred)], axis = 1)
            
        final_pred = self.ridge.predict(new_X_test)
        return final_pred
    
    def create_mlp(self):
        learning_rate = self.state[1]
        model = Sequential()
        for i in range(4,len(self.state),2):
            n_neurons = self.state[i]
            activation = self.activation_indices[self.state[i+1]]
            if n_neurons != 0:
                model.add(Dense(n_neurons,activation=activation))
        model.add(Dense(1, activation="linear"))
        model.compile(loss=self.root_mean_squared_error,optimizer=Adam(learning_rate))
        return model
    
    def root_mean_squared_error(self, y_true, y_pred):
        y_true = cast(y_true,float32)
        return tf.sqrt(tf.reduce_mean(tf.square(y_pred - y_true)))
    
    def fit_mlp(self,X_train,Y_train):
        epochs = self.state[0]
        validation_split = self.state[2]
        batch_size = self.state[3]
        self.mlp.fit(X_train,Y_train,validation_split=validation_split,
                     batch_size=batch_size,epochs=epochs,verbose=0)
