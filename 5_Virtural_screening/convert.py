#!/usr/bin/env python
# coding: utf-8
# Author: Cong Liu, Yuemin Bian

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.Descriptors import ExactMolWt
from rdkit.Chem.Crippen import MolLogP
from rdkit.Chem.rdMolDescriptors import CalcNumHBA, CalcNumHBD, CalcTPSA
from rdkit.Chem.Descriptors import NumRotatableBonds
from tqdm import tqdm
import argparse

parser=argparse.ArgumentParser(description="Convert SDF file to CSV format")
parser.add_argument("--input",type=str)
parser.add_argument("--output",type=str,default="output.csv")

args=parser.parse_args()
sdfinput = args.input
csvoutput = args.output

suppl = AllChem.SDMolSupplier(sdfinput)
with open(csvoutput,"w") as fl:
    fl.write(
        "Molecule (RDKit Mol),Catalog ID,MW,MW (desalted),ClogP,HBD,HBA,TPSA,RotBonds\n"
    )
    for mol in tqdm(suppl):
        # compute molecular properties
        smi = Chem.MolToSmiles(mol)
        name = mol.GetProp("Catalog ID")
        MW1 = round(ExactMolWt(mol),3)
        MW2 = mol.GetProp("MW (desalted)")
        clogP = round(MolLogP(mol),3)
        HBD = CalcNumHBD(mol)
        HBA = CalcNumHBD(mol)
        TPSA = round(CalcTPSA(mol),1)
        RotBonds = NumRotatableBonds(mol)

        # write output.
        fl.write(
            f"{smi},{name},{MW1},{MW2},{clogP},{HBD},{HBA},{TPSA},{RotBonds}\n"
        )
