#!/usr/bin/env python
# coding: utf-8
# Author: Yuemin Bian

print('Module 2: Compound collection')

#Python3 is required for the chembl api
#It is so sad that I have to jump from python2 to python3 because of this...

import argparse
from chembl_webresource_client.new_client import new_client
import pandas as pd
import numpy as np
#pd.set_option('display.max_columns', None)

# # Parse input

parser = argparse.ArgumentParser()
parser.add_argument('-i','--input', help = "The target list for compound collection. Or simply the generated file from the Module 1", required=True)
parser.add_argument('-f','--file_name', help = "A name to save out files...", required=True)
args = parser.parse_args()

# # Make functions

def uniprot_id_chembl_id(uniprot_id):
    print('Transfering the uniprot id to the chembl id...')

    target = new_client.target
    uniprot_id = uniprot_id
    res = target.filter(target_components__accession=uniprot_id)

    res = pd.DataFrame(res)
    return res

def protein_chembl_id(res):
    i=0
    while i < res.shape[0]:
        if res['target_type'][i] == 'SINGLE PROTEIN':
            target_chembl_id = res['target_chembl_id'][i]
        i=i+1
    return target_chembl_id

def get_compounds(retrieved_id):
    print('Collecting compounds for '+retrieved_id)
    activity = new_client.activity
    res = activity.filter(target_chembl_id=retrieved_id)
    res = pd.DataFrame(res)
    print (str(len(res))+' compounds colllected')
    return res

def analysis_to_collected_cmpds(combined):
    test_type=[]
    number=[]
    assay=pd.DataFrame()
    for i in pd.unique(combined['standard_type']):
        test_type.append(i)
        number.append(combined.loc[combined['standard_type'] == i].shape[0])
        #print (i, combined.loc[combined['standard_type'] == i].shape[0])
    assay['type']=test_type
    assay['count']=number
    assay = assay.sort_values('count', ascending=False, ignore_index=True)
    print(str(assay.shape[0])+' values types')
    print(assay.iloc[0,0] + ' has the most record(s) - ' + str(assay.iloc[0,1]) + ' records')
    return assay

def remove_dup(combined,assay):
    combined = combined.loc[combined['standard_type'] == assay.iloc[0,0]]
    combined_one = combined.drop_duplicates(subset=['molecule_chembl_id'])
    print('After removing duplicates, ' + str(combined_one.shape[0]) + ' compounds remains with ' + assay.iloc[0,0]) 
    return combined_one

def active_cmpds(training_data):
    training_data["standard_value"] = pd.to_numeric(training_data["standard_value"])
    if str(training_data["standard_type"].iloc[0]) == 'Inhibition':
        training_data_active = training_data.loc[training_data['standard_value'] >= 50]
        print (str(training_data_active.shape[0]) + " active compounds collected with standard value > 50%") 
    else:
        training_data_active = training_data.loc[training_data['standard_value'] <= 1000]
        print (str(training_data_active.shape[0]) + " active compounds collected with standard value < 1000 nM") 
    return training_data_active

def inactive_cmpds(training_data):
    training_data["standard_value"] = pd.to_numeric(training_data["standard_value"])
    if str(training_data["standard_type"].iloc[0]) == 'Inhibition':
        training_data_inactive = training_data.loc[training_data['standard_value'] < 50]
        print (str(training_data_inactive.shape[0]) + " inactive compounds collected with standard value < 50%")         
    else:
        training_data_inactive = training_data.loc[training_data['standard_value'] > 1000]
        print (str(training_data_inactive.shape[0]) + " inactive compounds collected with standard value > 1000 nM") 
    return training_data_inactive

def write_out(table,name):
    file_name=str(name)
    table.to_csv(file_name+'.csv',index=False)
    print ("Collected info is written to the disk.")
    return

# # Use functions

df=pd.read_csv(args.input)
combined=pd.DataFrame()
print (args.input + ' is recognized')

for uniprot_id_one in df['uniprot id']:
    res=uniprot_id_chembl_id(uniprot_id_one)
    if np.array(res) != []:
        chembl_id = protein_chembl_id(res)
        table=get_compounds(chembl_id)
        combined = pd.concat([combined,table],ignore_index=True)
    
write_out(combined,args.file_name)

assay=analysis_to_collected_cmpds(combined)
write_out(assay,'assay_info_' + args.file_name)

training_data=remove_dup(combined,assay)
actives=active_cmpds(training_data)
write_out(actives,'actives_' + args.file_name)
inactives=inactive_cmpds(training_data)
write_out(inactives,'inactives_' + args.file_name)
