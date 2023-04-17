#!/usr/bin/env python
# coding: utf-8
import stk
import rdkit.Chem,rdkit.Chem.AllChem as rdkit
from rdkit import Chem
from rdkit.Chem import rdqueries
from rdkit.Chem import RWMol
from rdkit.Chem import Lipinski
from rdkit import RDLogger
import numpy as np
import pandas as pd
import sys

# Some predefined functions for calculating VDW and MW
def Get_molecular(molecular):
    m=molecular.replace("[*]", "").replace("()", "")
    m1= rdkit.MolFromSmiles(m)
    m1= rdkit.AddHs(m1)
    return m1

def Count_number(molecular,index):
    q = rdqueries.AtomNumEqualsQueryAtom(index)
    return len(molecular.GetAtomsMatchingQuery(q))

def Index_atom(atom):
    if atom=='H':
        index=1
    elif atom=='B':
        index=5        
    elif atom=='C':
        index=6
    elif atom=='N':
        index=7
    elif atom=='O':
        index=8
    elif atom=='F':
        index=9
    elif atom=='Si':
        index=14
    elif atom=='P':
        index=15
    elif atom=='S':
        index=16
    elif atom=='Cl':
        index=17
    elif atom=='Se':
        index=34
    elif atom=='Br':
        index=35
    elif atom=='As':
        index=33
    elif atom=='I':
        index=53
    return index

def Count_Aromatic(molecular):
    return Lipinski.NumAromaticRings(molecular)

def Count_Noaromatic(molecular):
    return Lipinski.NumAliphaticRings(molecular)

# Some predefined functions for calculating monomer length
#Note that we took iodine atoms as monomeric linkage sites because they do not appear in the backbone
def Group(smi):
    moment=smi.replace("I", "[3H]")
    moment=moment.replace("[*]", "I")
    moment=moment.replace("*", "I")
    return moment

def Line_polymer(monomer,number):
    bb=stk.BuildingBlock(monomer,[stk.IodoFactory()])
    polymer= stk.ConstructedMolecule(
    topology_graph=stk.polymer.Linear((bb, ), 'A', number),)
    return  polymer.to_rdkit_mol()

def polymer_Optimize(polymer):
    rdkit.SanitizeMol(polymer)
    rdkit.MMFFOptimizeMolecule(polymer)
    return polymer

def Get_Connect(polymer):
    Connect_index=0
    atom_Num=[]
    for atom in polymer.GetAtoms():
        Num=atom.GetAtomicNum()
        atom_Num.append(Num)
    for i in range (0,len(atom_Num)):
        if (atom_Num[i-1]==1 and atom_Num[i]>1 and i>0):
            Connect_index=i
    return Connect_index 

def Get_NoHydrogen_Connect(smiles):
    smiles=smiles.replace("[*]", "")
    smiles=smiles.replace("()", "")
    mol=rdkit.MolFromSmiles(smiles)
    Atoms=[]
    for atom in mol.GetAtoms():
        Atoms.append(atom.GetIdx())
    return (len (Atoms)+1)

def Cal_distance(monomer,index1,index2):
    x1,y1,z1=monomer.GetConformer().GetAtomPosition(index1)
    x2,y2,z2=monomer.GetConformer().GetAtomPosition(index2)
    distance=np.sqrt((x1-x2)**2+(y1-y2)**2+(z1-z2)**2)
    return distance

# Some predefined functions for making mainchain
def set_linker_flag(mol, reverse=False):
    """
    poly.set_linker_flag
    Args:
        mol: RDkit Mol object
        reverse: Reversing head and tail (boolean)
    Returns:
        boolean
    """

    flag = False
    mol.SetIntProp('head_idx', -1)
    mol.SetIntProp('tail_idx', -1)
    mol.SetIntProp('head_ne_idx', -1)
    mol.SetIntProp('tail_ne_idx', -1)

    for atom in mol.GetAtoms():
        atom.SetBoolProp('linker', False)
        atom.SetBoolProp('head', False)
        atom.SetBoolProp('tail', False)
        atom.SetBoolProp('head_neighbor', False)
        atom.SetBoolProp('tail_neighbor', False)
        if (atom.GetSymbol() == "H" and atom.GetIsotope() == 3) or atom.GetSymbol() == "*":
            atom.SetBoolProp('linker', True)
            if not flag:
                mol.SetIntProp('head_idx', atom.GetIdx())
                mol.SetIntProp('tail_idx', atom.GetIdx())
                flag = True
            else:
                if reverse:
                    mol.SetIntProp('head_idx', atom.GetIdx())
                else:
                    mol.SetIntProp('tail_idx', atom.GetIdx())

    if not flag: return False

    mol_head_idx = mol.GetIntProp('head_idx')
    mol.GetAtomWithIdx(mol_head_idx).SetBoolProp('head', True)

    mol_tail_idx = mol.GetIntProp('tail_idx')
    mol.GetAtomWithIdx(mol_tail_idx).SetBoolProp('tail', True)

    head_ne_idx = mol.GetAtomWithIdx(mol_head_idx).GetNeighbors()[0].GetIdx()
    mol.SetIntProp('head_ne_idx', head_ne_idx)
    mol.GetAtomWithIdx(head_ne_idx).SetBoolProp('head_neighbor', True)

    tail_ne_idx = mol.GetAtomWithIdx(mol_tail_idx).GetNeighbors()[0].GetIdx()
    mol.SetIntProp('tail_ne_idx', tail_ne_idx)
    mol.GetAtomWithIdx(tail_ne_idx).SetBoolProp('tail_neighbor', True)

    return True

def set_terminal_idx(mol):

    count = 0
    for atom in mol.GetAtoms():
        resinfo = atom.GetPDBResidueInfo()
        if resinfo is None: continue
        resname = resinfo.GetResidueName()

        if resname == 'TU0':
            for na in atom.GetNeighbors():
                if na.GetPDBResidueInfo() is None: continue
                elif na.GetPDBResidueInfo().GetResidueName() != 'TU0':
                    mol.SetIntProp('terminal_idx1', atom.GetIdx())
                    mol.SetIntProp('terminal_ne_idx1', na.GetIdx())
                    count += 1

        elif resname == 'TU1':
            for na in atom.GetNeighbors():
                if na.GetPDBResidueInfo() is None: continue
                elif na.GetPDBResidueInfo().GetResidueName() != 'TU1':
                    mol.SetIntProp('terminal_idx2', atom.GetIdx())
                    mol.SetIntProp('terminal_ne_idx2', na.GetIdx())
                    count += 1
    
    return count

def set_mainchain_flag(mol):

    for atom in mol.GetAtoms():
        atom.SetBoolProp('main_chain', False)
    
    linker_result = set_linker_flag(mol)
    terminal_result = set_terminal_idx(mol)
    if not linker_result and terminal_result == 0:
        return False

    if terminal_result > 0:
        if mol.GetIntProp('terminal_ne_idx1') == mol.GetIntProp('terminal_ne_idx2'):
            path = [mol.GetIntProp('terminal_ne_idx1')]
        else:
            path = Chem.GetShortestPath(mol, mol.GetIntProp('terminal_ne_idx1'), mol.GetIntProp('terminal_ne_idx2'))
    elif mol.GetIntProp('head_idx') == mol.GetIntProp('tail_idx'):
        path = Chem.GetShortestPath(mol, mol.GetIntProp('head_idx'), mol.GetIntProp('head_ne_idx'))
    else:
        path = Chem.GetShortestPath(mol, mol.GetIntProp('head_idx'), mol.GetIntProp('tail_idx'))

    for idx in path:
        atom = mol.GetAtomWithIdx(idx)
        atom.SetBoolProp('main_chain', True)
        for batom in atom.GetNeighbors():
            if batom.GetTotalDegree() == 1:  # Expect -H, =O, =S, -F, -Cl, -Br, -I
                batom.SetBoolProp('main_chain', True)
    
    rings = Chem.GetSymmSSSR(mol)
    m_rings = []
    for ring in rings:
        dup = list(set(path) & set(ring))
        if len(dup) > 0:
            m_rings.append(ring)
            for idx in ring:
                atom = mol.GetAtomWithIdx(idx)
                atom.SetBoolProp('main_chain', True)
                for batom in atom.GetNeighbors():
                    if batom.GetTotalDegree() == 1:  # Expect -H, =O, =S, -F, -Cl, -Br, -I
                        batom.SetBoolProp('main_chain', True)

    for m_ring in m_rings:
        for ring in rings:
            dup = list(set(m_ring) & set(ring))
            if len(dup) > 0:
                for idx in ring:
                    atom = mol.GetAtomWithIdx(idx)
                    atom.SetBoolProp('main_chain', True)
                    for batom in atom.GetNeighbors():
                        if batom.GetTotalDegree() == 1:  # Expect -H, =O, =S, -F, -Cl, -Br, -I
                            batom.SetBoolProp('main_chain', True)

    return True


#vander Waals volume
#V_vdW(Å3/molecule))=∑allatomcontributions-5.92NN_B-14.7R_A-3.8R_NR
def cal_VDW (smiles):
    mol= rdkit.MolFromSmiles(smiles)
    mol= rdkit.AddHs(mol)
    num_C=Count_number(mol,Index_atom('C'))
    num_H=Count_number(mol,Index_atom('H'))
    num_N=Count_number(mol,Index_atom('N'))
    num_O=Count_number(mol,Index_atom('O'))
    num_Cl=Count_number(mol,Index_atom('Cl'))
    num_Br=Count_number(mol,Index_atom('Br'))
    num_F=Count_number(mol,Index_atom('F'))
    num_I=Count_number(mol,Index_atom('I'))
    num_S=Count_number(mol,Index_atom('S'))
    num_P=Count_number(mol,Index_atom('P'))
    num_As=Count_number(mol,Index_atom('As'))
    num_B=Count_number(mol,Index_atom('B'))
    num_Aromatic=Count_Aromatic(mol)
    num_Noaromatic=Count_Noaromatic(mol)
    All_Bonds=num_C*20.58+num_H*7.24+num_N*15.60+num_O*14.71+num_Cl*22.45+num_Br*26.52+num_F*13.31+num_I*32.52+num_S*24.43+num_P*24.43+num_As*26.52+num_B*40.48
    N_B=num_C+num_H+num_N+num_O+num_Cl+num_Br+num_F+num_I+num_S+num_P+num_As+num_B+num_Aromatic+num_Noaromatic-1
    R_A=num_Aromatic
    R_NR=num_Noaromatic
    V_vdw=All_Bonds-5.92*N_B-14.7*R_A-3.8*R_NR    
    return np.round (V_vdw,2)

#Molecular weight
def cal_MW (smiles):
    mol= rdkit.MolFromSmiles(smiles)
    mol= rdkit.AddHs(mol)
    num_C=Count_number(mol,Index_atom('C'))
    num_H=Count_number(mol,Index_atom('H'))
    num_N=Count_number(mol,Index_atom('N'))
    num_O=Count_number(mol,Index_atom('O'))
    num_Cl=Count_number(mol,Index_atom('Cl'))
    num_Br=Count_number(mol,Index_atom('Br'))
    num_F=Count_number(mol,Index_atom('F'))
    num_I=Count_number(mol,Index_atom('I'))
    num_S=Count_number(mol,Index_atom('S'))
    num_P=Count_number(mol,Index_atom('P'))
    num_As=Count_number(mol,Index_atom('As'))
    num_B=Count_number(mol,Index_atom('B'))
    MW=num_C*12.01+num_H*1.01+num_N*14.01+num_O*16.00+num_Cl*35.45+num_Br*79.90+num_F*19.00+num_I*126.90+num_S*32.07+num_P*30.97+num_As*74.92+num_B*10.81
    return MW

#Calculation of monomer length
def Cal_length (smi,number):
    polymer=Line_polymer(Group(smi),number)
    polymer=polymer_Optimize(polymer)
    secong_index=Get_Connect(polymer)
    if secong_index!=0:
        index=secong_index
    else:
        index=Get_NoHydrogen_Connect(smi)        
    length=Cal_distance(polymer,1,index)
    return np.round (length,2)

# extract mainchain (Referred from https://github.com/RadonPy/RadonPy/blob/29a5b1c33da68a826151c9637258125c235094cd/radonpy/core/poly.py)
def extract_mainchain(smiles):
    main_smi = None
    if smiles.count('*') != 2:
        print('Illegal number of connecting points in SMILES. %s' % smiles)
        return main_smi

    smi = smiles.replace('[*]', '[3H]')
    smi = smi.replace('*', '[3H]')

    try:
        mol = Chem.MolFromSmiles(smi)
        mol = Chem.AddHs(mol)
    except:
        print('Cannot convert to Mol object from %s' % smiles)
        return main_smi

    set_mainchain_flag(mol)

    for atom in mol.GetAtoms():
        if atom.GetBoolProp('main_chain'):
            for na in atom.GetNeighbors():
                if not na.GetBoolProp('main_chain'):
                    bidx = mol.GetBondBetweenAtoms(atom.GetIdx(), na.GetIdx()).GetIdx()
                    mol = Chem.FragmentOnBonds(mol, [bidx], addDummies=False)

    RDLogger.DisableLog('rdApp.*')

    try:
        fsmi = [Chem.MolToSmiles(Chem.MolFromSmiles(x)) for x in Chem.MolToSmiles(mol).split('.')]
    except:
        print('Cannot convert to fragmented Mol')
        RDLogger.EnableLog('rdApp.*')
        return main_smi

    RDLogger.EnableLog('rdApp.*')

    for s in fsmi:
        if '[3H]' in s:
            try:
                main_smi = Chem.MolToSmiles(Chem.MolFromSmiles(s.replace('[3H]', '*')))
            except:
                print('Cannot convert to canonical SMILES from %s' % smi)

    return main_smi

#Molecular weight of Backbone/Molecular weight of Monomer
def cal_MWratio (smi):
    Backbone_smi=extract_mainchain(smi)
    MW_ratio=cal_MW (Backbone_smi)/cal_MW (smi)
    return np.round (MW_ratio,3)   


# In[12]:


if __name__ == '__main__':
    path="./dataset/"
    filename = sys.argv[1]
    df=pd.read_csv(path+str(filename))
    SMI=df["SMILES"]
    df['VDW'] = df['SMILES'].map(lambda x: cal_VDW (x)) 
    df['MW'] = df['SMILES'].map(lambda x: cal_MW (x))
    df['Monomer_length'] = df['SMILES'].map(lambda x: Cal_length (x,2))  
    df['MW_ratio'] = df['SMILES'].map(lambda x: cal_MWratio (x))
    df.to_csv(path+"Des_monomer.csv",index=None)




