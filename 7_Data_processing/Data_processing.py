#!/usr/bin/env python
# coding: utf-8
## Author: Yuemin Bian

print('Module 7: Data processing')

import sys
import numpy as np
import pandas as pd
np.set_printoptions(threshold=sys.maxsize)
import argparse

# # Parse input

parser = argparse.ArgumentParser()
parser.add_argument('-r','--result_rf', help = 'The outcome rf file from the module 6', required=True)
parser.add_argument('-m','--result_mlp', help = 'The outcome mlp file from the module 6', required=True)
parser.add_argument('-f','--file_name', help='The output file name', required=True)
args = parser.parse_args()

# # Make functions

def read_screen_result(name):
    print ('Read '+ name +'...')
    outcome = pd.read_csv(name, header=0)
    return outcome

def add_rank(outcome,model_type):
    i=1
    rank=[]
    while i < len(outcome)+1:
        rank.append(i)
        i+=1
    outcome[model_type+'_rank']=rank
    return outcome

def merge_two(mlp,rf):
    print ('Merge two files...')
    RF_score=[]
    RF_rank=[]
    mlp = mlp[['comp_id','MLP_prediction_score','mlp_rank']]
    combined = rf.merge(mlp,on='comp_id')
    RF_score = combined['RF_prediction_score']
    RF_rank = combined['rf_rank']
    combined = combined.drop('RF_prediction_score', axis=1)
    combined['RF_prediction_score'] = RF_score
    combined = combined.drop('rf_rank', axis=1)
    combined['rf_rank'] = RF_rank
    return combined

def normalize_scores(combined):
    norm_rf = (combined['RF_prediction_score']-combined['RF_prediction_score'].min())/(combined['RF_prediction_score'].max()-combined['RF_prediction_score'].min())
    norm_mlp = (combined['MLP_prediction_score']-combined['MLP_prediction_score'].min())/(combined['MLP_prediction_score'].max()-combined['MLP_prediction_score'].min())
    combined['norm_mlp_score'] = norm_mlp
    combined['norm_RF_score'] = norm_rf
    return combined

def ensemble_score(normalized_combined):
    print ('Add ensemble ranking as well...')
    i=0
    ensemble=[]
    while i < len(normalized_combined):
        ensemble.append((normalized_combined.iloc[i,:]['rf_rank']+normalized_combined.iloc[i,:]['mlp_rank'])/2)
        i+=1
    normalized_combined['ensemble_rank']=ensemble
    normalized_combined = normalized_combined.sort_values(by=['ensemble_rank'],ascending=True)
    return normalized_combined

def top_1_percent_combined(normalized_combined_ensemble):
    print ('Select non-duplicated top 1% VS hits from MLP, RF, and ensemble...')
    length = int(len(normalized_combined_ensemble)/100)
    normalized_combined_ensemble_1 = normalized_combined_ensemble.iloc[:length]
    mlp_1 = normalized_combined_ensemble.sort_values(by=['mlp_rank'],ascending=True).iloc[:length]
    rf_1 = normalized_combined_ensemble.sort_values(by=['rf_rank'],ascending=True).iloc[:length]
    
    frames= [normalized_combined_ensemble_1, mlp_1, rf_1]
    percent_1_combined = pd.concat(frames)
    percent_1_combined = percent_1_combined.drop_duplicates(subset=['comp_id'])
    
    return percent_1_combined

def write_out(outcome,name):
    file_name=str(name)
    outcome.to_csv(file_name+'_top_1_percent_VS_hits'+'.csv',index=False)
    print ("Done. Top selected compounds are written to the disk")
    return

def write_out_all(outcome,name):
    file_name=str(name)
    outcome.to_csv(file_name+'_VS_all'+'.csv',index=False)
    print ("Done. All screened compounds are written to the disk")
    return

# # Use functions
# Load in outcome
mlp=read_screen_result(args.result_mlp)
rf=read_screen_result(args.result_rf)

mlp = add_rank(mlp,'mlp')
rf = add_rank(rf,'rf')
combined=merge_two(mlp,rf)
combined_ensemble=ensemble_score(combined)
outcome = top_1_percent_combined(combined_ensemble)

write_out(outcome,args.file_name)
write_out_all(combined_ensemble,args.file_name)
