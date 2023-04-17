#!/usr/bin/env python
# coding: utf-8

import pandas as pd
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem
import mordred
from mordred import Calculator, descriptors
import sys

def embed(mol):
    #print (Chem.MolToSmiles(mol))
    mol_with_H = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol_with_H)
    AllChem.MMFFOptimizeMolecule(mol_with_H)
    return mol_with_H

if __name__ == '__main__':
    path="./data/"
    filename=sys.argv[1]
    mols = pd.read_csv(path+str(filename))
    mols['rdmol'] = mols['SMILES'].map(lambda x: Chem.MolFromSmiles(x))    # drop duplicates based on inchi
    mols = mols.drop_duplicates(subset="rdmol")
    mols['rdmol_optimized'] = mols['rdmol'].map(embed)
    calc = Calculator(descriptors) # create calculator for all mordred descriptors (can also specify subtype)
    df=calc.pandas(mols['rdmol_optimized'])
    df=df.applymap(lambda x: np.nan if type(x) in [mordred.error.Missing,mordred.error.Error] else x)
    df=df.dropna(axis=1)
    non_zero_std = df.std() != 0
    df = df [non_zero_std[non_zero_std].index]
    threshold=0.98
    df_corr = df.corr().abs()
    upper = df_corr.where(np.triu(np.ones(df_corr.shape), k=1).astype(np.bool))
    to_drop = [column for column in upper.columns if any(upper[column] > threshold)]
    df = df.drop(to_drop, axis=1)
    to_save=mols[["ID","SMILES"]].join(df)
    to_save.to_csv(path+"Des_Mordred.csv")