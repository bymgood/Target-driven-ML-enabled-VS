#!/usr/bin/env python
# coding: utf-8
# Author: Yuemin Bian

print('Module of library preparation')

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
parser.add_argument('-i','--input', help = "The library file to be prepared", required=True)
parser.add_argument('-s','--smiles', help = "Which column is the SMILES column? Starting from 1", required=True)
parser.add_argument('-c','--cmpd_id', help = "Which column is the compound id column? Starting from 1", required=True)
parser.add_argument('-t','--fp_type', help = "Optional. Select from Morgan, AtomPair, Topological, and MACCS. Default=Morgan", required=False)
parser.add_argument('-n','--number_of_bits', help = "Optional. The number of bits for Morgan, AtomPair and Topological FPs. Default=1024", required=False)
parser.add_argument('-f','--file_name', help = "A name to save out files...", required=True)
args = parser.parse_args()


# # Make functions

def read_cmpds(name,smi_col, id_col):
    print("Read input library file")
    df=pd.read_csv(name)
    smiles=df.iloc[:,smi_col]
    cmpd_id=df.iloc[:,id_col]
    #smiles = pd.DataFrame(smiles)
    #cmpd_id = pd.DataFrame(cmpd_id)
    return smiles, cmpd_id

def get_morganfp(df,bits):
    print("Calcuating Morgan FP...")
    bit_fp=[]
    count=0
    for i in df:
        m = Chem.MolFromSmiles(i)
        n = AllChem.GetMorganFingerprintAsBitVect(m,2,bits)
        fp=np.zeros((1,))
        DataStructs.ConvertToNumpyArray(n, fp)
        bit_fp.append(fp)
        count = count+1
        print (count)
    cmpd_fp=pd.concat([pd.DataFrame(df),pd.DataFrame(bit_fp).astype(int)],axis=1)
    print('FP calcuated')
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

def write_out(table,name):
    print('In progress...')
    file_name=str(name)
    table.to_csv(file_name+'.csv',index=False)
    print ("Calculation finished. Files saved to the disk.")
    return

# # Use functions

df,cmpd_id = read_cmpds(args.input,(int(args.smiles)-1),(int(args.cmpd_id)-1))

if args.number_of_bits:
    bits = int(args.number_of_bits)
else:
    bits = 1024

if args.fp_type == 'Morgan':
    cmpd_fp = get_morganfp(df,bits)
    cmpd_fp = pd.concat([cmpd_id, cmpd_fp],axis=1)
    print('Preparing write out')
elif args.fp_type == 'AtomPair':
    cmpd_fp = get_AtomPairfp(df,bits)
    cmpd_fp = pd.concat([cmpd_id, cmpd_fp],axis=1)
    print('Preparing write out')
elif args.fp_type == 'Topological':
    cmpd_fp = get_TopologicalTorsionfp(df,bits)
    cmpd_fp = pd.concat([cmpd_id, cmpd_fp],axis=1)
    print('Preparing write out')
elif args.fp_type == 'MACCS':
    cmpd_fp = get_MACCS(df)
    cmpd_fp = pd.concat([cmpd_id, cmpd_fp],axis=1)
    print('Preparing write out')
else:
    cmpd_fp = get_morganfp(df,bits)
    cmpd_fp = pd.concat([cmpd_id, cmpd_fp],axis=1)
    print('Preparing write out')

write_out(cmpd_fp,args.file_name)
