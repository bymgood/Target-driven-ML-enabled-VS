#!/usr/bin/env python
# coding: utf-8
# Author: Yuemin Bian

print('Module 5: Virtual screening')

import pandas as pd
import pickle
import argparse
import numpy as np
import matplotlib.pyplot as plt

from scipy import interp
from itertools import cycle
from sklearn.naive_bayes import BernoulliNB
from sklearn.metrics import roc_curve, auc
from sklearn.model_selection import StratifiedKFold
from sklearn.model_selection import train_test_split
from sklearn.cluster import KMeans

# # Parse input

parser = argparse.ArgumentParser()
parser.add_argument('-m','--model', help = 'Load the prepared machine learning model', required=True)
parser.add_argument('-t','--model_type', help = 'Specify the algorithm of model to be loaded - choose from MLP and RF', required=True)
parser.add_argument('-s','--screen_set', help='load the compound set to screen', required=True)
parser.add_argument('-f','--file_name', help="filename to save the result to", required=True)
args = parser.parse_args()

# # Make functions

def read_model(name):
    model = pickle.load(open(name, 'rb'))
    print (name + " loaded")
    return model

def read_library(name):
    screen_set = pd.read_csv(name, header=0)
    print('screen set loaded')
    print('screen set includes '+str(screen_set.shape[0])+' compounds')
    return screen_set

def screening(model, screen_set, name):
    print('In progress of screening...')
    comp_id = screen_set.iloc[:,0]
    smiles = screen_set.iloc[:,1]
    fps = screen_set.iloc[:,2:]
    test_set_probas_ = model.predict_proba(fps)
    fps = []
    outcome=pd.DataFrame()
    outcome['comp_id'] = comp_id
    outcome['smiles'] = smiles
    col_name = name + '_prediction_score'
    outcome[col_name] = test_set_probas_[:,1]
    return outcome

def write_out(outcome,name):
    file_name=str(name)
    outcome.to_csv(file_name+'.csv',index=False)
    print ("Done. The library file with scores is saved to the disk.")
    return

# # Use functions

model = read_model(args.model)
screen_set=read_library(args.screen_set)
outcome = screening(model, screen_set, args.model_type)
write_out(outcome,args.file_name)

