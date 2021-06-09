import os
import sys
import pickle
import numpy as np
import pandas as pd
import argparse
from sklearn import preprocessing
from sklearn import metrics
from sklearn.ensemble import RandomForestRegressor
from sklearn.metrics import mean_absolute_error
from extractFeatures import extractCT,get_resname_char


def load_ss2cs_model(nucleus, DIR_PATH):
  ''' load save model '''
  filename = DIR_PATH + '/model/RF_' + nucleus + '.sav'
  model = pickle.load(open(filename, 'rb'))
  return(model)

def main():
    # configure parser
    parser = argparse.ArgumentParser()
    parser.add_argument("-i","--input", help="path to read input CT secondary structure file", required = True)
    parser.add_argument("-o","--output", help="path to save output chemical shift file", required = True)
    parser.add_argument("-s","--ss2cs_path", help="path to SS2CS repo", required = True)

    # parse command line
    a = parser.parse_args()  

    # initialize    
    inFile = a.input
    outFile = a.output
    DIR_PATH = a.ss2cs_path   
    rna = "user"
    nuclei = ["C1'", "C2'", "C3'", "C4'", "C5'","C2","C5","C6","C8", "H1'", "H2'", "H3'","H4'", "H5'","H5''","H2","H5","H6","H8"]

    # featurization
    features = extractCT(DIR_PATH+"/"+inFile, rna)
    features.drop('i_resname_char', axis=1, inplace=True)

    # fit one hot encoder
    train_X = pd.read_csv(DIR_PATH+"/data/train_X_NEW.csv",sep=' ',header=0)
    train_X = train_X.drop(['id','length','resid'],axis = 1)
    enc = preprocessing.OneHotEncoder(sparse = False)
    enc.fit(train_X)
    
    # fit model for each nucleus type
    results = pd.DataFrame([])
    for nucleus in nuclei:
    # one hot encoding testing data
        features_resname = features.drop(['id', 'length', 'resid'],axis=1)
        features_info = features['length']
        features_resname_enc = pd.DataFrame(enc.transform(features_resname))
        features_enc = pd.concat([features_info, features_resname_enc],axis = 1)

        # model prediction
        model = load_ss2cs_model(nucleus, DIR_PATH)
        y_pred = model.predict(features_enc)

        # format prediction
        output_resname = features['i_resname'].apply(lambda x: get_resname_char(x))
        output_resid = features['resid']
        output_nucleus = pd.Series([nucleus]*len(features))
        output_cs = pd.Series(y_pred)
        output_error = pd.Series(["."]*len(features))
        result = pd.concat([output_resname, output_resid, output_nucleus, output_cs, output_error],axis=1)
        results = pd.concat([results, result],ignore_index=True)
        
    results.to_csv(DIR_PATH+"/"+outFile, sep=' ', header=None, index=False)

if __name__ == "__main__":
    main()