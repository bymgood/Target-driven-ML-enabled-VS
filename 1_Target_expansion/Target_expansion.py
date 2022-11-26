#!/usr/bin/env python
# coding: utf-8
# Author: Yuemin Bian

print('Module 1: Target expansion')

import argparse
import pandas as pd
import numpy as np
from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML
from Bio import SeqIO
from Bio import ExPASy
from Bio import Align
from Bio import pairwise2
from Bio.SubsMat.MatrixInfo import blosum62

# # Parse input

parser = argparse.ArgumentParser()
parser.add_argument('-i','--input', help = "Entering the protein Uniprot ID to start", required=True)
parser.add_argument('-f','--file_name', help = "A name to save out files...", required=True)
parser.add_argument('-p','--percent', help = "Optional. A percentage identity of selecting similar proteins. e.g. 0.20", required=False)
args = parser.parse_args()

# # Make functions

def get_blast_record(uniprot_id,size):
    print('Get targets through blastp search')

    print('Connecting NCBI server...')

    #Run blast
    #fasta_string = open(args.input).read()
    result_handle = NCBIWWW.qblast(program="blastp",database= "nr",sequence= uniprot_id, ncbi_gi=True, hitlist_size=size, entrez_query="txid9606[ORGN]")
    name='raw'+"_blast"+".xml"
    name=str(name)
    with open(name, "w") as out_handle:
        out_handle.write(result_handle.read())
    result_handle.close()
    result_handle = open(name)

    print('Data extracted. Preparing the results...')
    
    blast_record = NCBIXML.read(result_handle)
    return blast_record

def get_dirty(blast_record):
    print('Blast result prepared')
    E_VALUE_THRESH = 0.04
    result=[]
    length=[]
    for alignment in blast_record.alignments:
        length.append(alignment.length)
    length=int(length[0])
    for alignment in blast_record.alignments:
         for hsp in alignment.hsps:    
            if hsp.expect < E_VALUE_THRESH:
                per_ident=float(hsp.identities)/float(length)
                one=(alignment.title, per_ident, hsp.bits)
                result.append(one)
    
    result=pd.DataFrame(result, columns=['description','percent identity','bit score'])
    #result.to_csv('blast_dirty.csv',index=False)
    return result

def blast_result(result):
    res=np.array(result)
    l=0
    gi=[]
    sc=[]
    uni=[]
    identity_score=float(0)
    for row in res:
        gi_one=str(row[0])
        sc_one=row[1]
        if ">sp|" in gi_one:
            if sc_one > identity_score:
                start = gi_one.find('>sp|') + 4
                end = start+6
                uniprot_id=gi_one[start:end]
                gi.append(gi_one)
                sc.append(sc_one)
                uni.append(uniprot_id)
                l=l+1
                
    gi=pd.DataFrame(gi, columns=['description'])
    sc=pd.DataFrame(sc, columns=['percent score'])
    uni=pd.DataFrame(uni, columns=['uniprot id'])
    
    df=pd.concat([gi, sc, uni], axis=1)
    df=df.sort_values(by=['percent score'], ascending = False)
    df=pd.DataFrame(np.array(df))
    df.columns=['description','percent score', 'uniprot id']
    df=df.assign(**dict.fromkeys(['source'], 'blastp'))
    return df

def write_blast_result(df,name):
    df.to_csv(name+'.csv',index=False)
    print('Blast result is written to the disk')
    return

# # Use functions

blast_record=get_blast_record(args.input,500)
result=get_dirty(blast_record)
df=blast_result(result)

if args.percent:
    df=df.loc[df['percent score'] >= float(args.percent)]

print (str(df.shape[0]) + ' records collected')
write_blast_result(df,args.file_name)

