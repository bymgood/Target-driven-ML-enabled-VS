#!/usr/bin/env python
# coding: utf-8
## Author: Yuemin Bian

print('Module 3: Fingerprint calculation')

import sys
import numpy as np
import pandas as pd
import argparse
from rdkit import Chem
from rdkit import DataStructs
from rdkit.Chem import AllChem
from rdkit.Chem import MACCSkeys
from rdkit.Chem.AtomPairs import Pairs
from rdkit.Chem.AtomPairs import Torsions
#np.set_printoptions(threshold=sys.maxsize)

# # Parse input

parser = argparse.ArgumentParser()
parser.add_argument('-i','--input', help = "The compound list for FP calculation. Or simply the output file from Module 2", required=True)
parser.add_argument('-a','--active_inactive', help = "Optional. Select from active and inactive. Is the input file for active or inactive molecules? Default=blank", required=False)
parser.add_argument('-t','--fp_type', help = "Optional. Select from Morgan, AtomPair, Topological, and MACCS. Default=Morgan", required=False)
parser.add_argument('-n','--number_of_bits', help = "Optional. The number of bits for Morgan, AtomPair and Topological FPs. Default=1024", required=False)
parser.add_argument('-f','--file_name', help = "A name to save out files...", required=True)
args = parser.parse_args()


# # Make functions

def read_cmpds(name):
    print("Read input compounds file")
    df=pd.read_csv(name)
    original=df.shape[0]
    df['canonical_smiles'].replace('', np.nan, inplace=True)
    df.dropna(subset=['canonical_smiles'], inplace=True)
    recognized=df.shape[0]
    diff=original-recognized
    print(str(recognized) + ' compounds recognized. ' + str(diff) + ' unrecognized SMILES are dropped.')
    df.reset_index(drop=True, inplace=True)
    df=df['canonical_smiles']
    return df

def get_morganfp(df,bits):
    print("Calcuating Morgan FP...")
    bit_fp=[]
    for i in df:
        m = Chem.MolFromSmiles(i)
        n = AllChem.GetMorganFingerprintAsBitVect(m,2,bits)
        fp=np.zeros((1,))
        DataStructs.ConvertToNumpyArray(n, fp)
        bit_fp.append(fp)

    cmpd_fp=pd.concat([pd.DataFrame(df),pd.DataFrame(bit_fp).astype(int)],axis=1)
    return cmpd_fp

def get_AtomPairfp(df,bits):
    print("Calcuating Atom Pair FP...")
    bit_fp=[]
    for i in df:
        m = Chem.MolFromSmiles(i)
        n = Pairs.GetHashedAtomPairFingerprint(m,bits)
        fp=np.zeros((1,))
        DataStructs.ConvertToNumpyArray(n, fp)
        bit_fp.append(fp)

    cmpd_fp=pd.concat([pd.DataFrame(df),pd.DataFrame(bit_fp).astype(int)],axis=1)
    return cmpd_fp

def get_TopologicalTorsionfp(df,bits):
    print("Calcuating Topological FP...")    
    bit_fp=[]
    for i in df:
        m = Chem.MolFromSmiles(i)
        n = Torsions.GetHashedTopologicalTorsionFingerprint(m,bits)
        fp = np.zeros((1,))
        DataStructs.ConvertToNumpyArray(n, fp)
        bit_fp.append(fp)

    cmpd_fp=pd.concat([pd.DataFrame(df),pd.DataFrame(bit_fp).astype(int)],axis=1)
    return cmpd_fp

def get_MACCS(df):
    print("Calcuating MACCS FP...")    
    bit_fp=[]
    for i in df:
        m = Chem.MolFromSmiles(i)
        n = MACCSkeys.GenMACCSKeys(m)
        fp = np.zeros((1,))
        DataStructs.ConvertToNumpyArray(n, fp)
        bit_fp.append(fp)

    cmpd_fp=pd.concat([pd.DataFrame(df),pd.DataFrame(bit_fp).astype(int)],axis=1)
    return cmpd_fp

def label_active(cmpd_and_fp):
    print("Compounds are labeled as active")
    label=[]
    count=0
    while count < cmpd_and_fp.shape[0]:
        label.append(int(1))
        count+=1
        
    cmpd_and_fp['label']=label    
    return cmpd_and_fp

def label_inactive(cmpd_and_fp):
    print("Compounds are labeled as inactive")
    label=[]
    count=0
    while count < cmpd_and_fp.shape[0]:
        label.append(int(0))
        count+=1
        
    cmpd_and_fp['label']=label    
    return cmpd_and_fp

def write_out(table,name):
    file_name=str(name)
    table.to_csv(file_name+'.csv',index=False)
    print ("Calculation finished. Files saved to the disk.")
    return

# # Use functions

df=read_cmpds(args.input)

if args.number_of_bits:
    bits = int(args.number_of_bits)
else:
    bits = 1024

if args.fp_type == 'Morgan':
    cmpd_fp = get_morganfp(df,bits)
elif args.fp_type == 'AtomPair':
    cmpd_fp = get_AtomPairfp(df,bits)
elif args.fp_type == 'Topological':
    cmpd_fp = get_TopologicalTorsionfp(df,bits)
elif args.fp_type == 'MACCS':
    cmpd_fp = get_MACCS(df)
else:
    cmpd_fp = get_morganfp(df,bits)

if args.active_inactive == "active":
    cmpd_and_fp = label_active(cmpd_fp)
elif args.active_inactive == "inactive":
    cmpd_and_fp = label_inactive(cmpd_fp)
else:
    cmpd_and_fp = cmpd_fp

write_out(cmpd_and_fp,args.file_name)
