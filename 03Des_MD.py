#!/usr/bin/env python
# coding: utf-8

from radonpy.core import poly, utils 
from radonpy.ff.gaff2_mod import GAFF2_mod
from radonpy.ff.descriptor import FF_descriptor
import numpy as np
import pandas as pd
import sys


def ff_descriptor(smiles):
    #smiles = '[*]=C1C(=O)Nc2ccc(-c3ccc4c(c3)C(=c3ccc(=c5ccc(=[*])s5)s3)C(=O)N4)cc21' 
    mol= poly.make_cyclicpolymer(smiles, n=2, return_mol=True)
    ff = GAFF2_mod() 
    result = ff.ff_assign(mol, charge='gasteiger') 
    desc=FF_descriptor(ff)
    all_desc=desc.get_param_list(mol)



    mass_max=np.max(all_desc[0])
    mass_min=np.min(all_desc[0])
    mass_ave=np.average(all_desc[0])

    charge_max=np.max(all_desc[1])
    charge_min=np.min(all_desc[1])
    charge_ave=np.average(all_desc[1])

    epsilon_max=np.max(all_desc[2])
    epsilon_min=np.min(all_desc[2])
    epsilon_ave=np.average(all_desc[2])

    sigma_max=np.max(all_desc[3])
    sigma_min=np.min(all_desc[3])
    sigma_ave=np.average(all_desc[3])

    k_bond_max=np.max(all_desc[4])
    k_bond_min=np.min(all_desc[4])
    k_bond_ave=np.average(all_desc[4])

    r0_max=np.max(all_desc[5])
    r0_min=np.min(all_desc[5])
    r0_ave=np.average(all_desc[5])

    k_ang_max=np.max(all_desc[6])
    k_ang_min=np.min(all_desc[6])
    k_ang_ave=np.average(all_desc[6])

    theta0_max=np.max(all_desc[7])
    theta0_min=np.min(all_desc[7])
    theta0_ave=np.average(all_desc[7])

    k_dih_max=np.max(all_desc[8])
    k_dih_min=np.min(all_desc[8])
    k_dih_ave=np.average(all_desc[8])
    
    return mass_max,mass_min,mass_ave,charge_max,charge_min,charge_ave,epsilon_max,epsilon_min,epsilon_ave,sigma_max,sigma_min,sigma_ave,\
        k_bond_max,k_bond_min,k_bond_ave,r0_max,r0_min, r0_ave, k_ang_max,k_ang_min, k_ang_ave, theta0_max,theta0_min,theta0_ave,k_dih_max,k_dih_min,k_dih_ave



if __name__ == '__main__':
    path="./data/"
    filename=sys.argv[1]
    dataframe = pd.read_csv(path+filename)
    SMILES=dataframe['SMILES']
    ID=dataframe['ID']
    Mass_max=[]
    Mass_min=[]
    Mass_ave=[]
    Charge_max=[]
    Charge_min=[]
    Charge_ave=[]
    Epsilon_max=[]
    Epsilon_min=[]
    Epsilon_ave=[]
    Sigma_max=[]
    Sigma_min=[]
    Sigma_ave=[]
    K_bond_max=[]
    K_bond_min=[]
    K_bond_ave=[]
    R0_max=[]
    R0_min=[]
    R0_ave=[]
    K_ang_max=[]
    K_ang_min=[]
    K_ang_ave=[]
    Theta0_max=[]
    Theta0_min=[]
    Theta0_ave=[]
    K_dih_max=[]
    K_dih_min=[]
    K_dih_ave=[]
    NAME=[]
    SMILE=[]
    success=[]
    for i in range(0,len(dataframe)):
        smi=SMILES[i]
        id_=ID[i]
        try:
            mass_max,mass_min,mass_ave,charge_max,charge_min,charge_ave,epsilon_max,epsilon_min,epsilon_ave,sigma_max,sigma_min,sigma_ave,\
                k_bond_max,k_bond_min,k_bond_ave,r0_max,r0_min, r0_ave, k_ang_max,k_ang_min, k_ang_ave, theta0_max,theta0_min,theta0_ave,\
                    k_dih_max,k_dih_min,k_dih_ave=ff_descriptor(smi)
            NAME.append (id_)
            SMILE.append(smi)    
            Mass_max.append(mass_max)
            Mass_min.append(mass_min)
            Mass_ave.append(mass_ave)
            Charge_max.append(charge_max)
            Charge_min.append(charge_min)
            Charge_ave.append(charge_ave)
            Epsilon_max.append(epsilon_max)
            Epsilon_min.append(epsilon_min)
            Epsilon_ave.append(epsilon_ave)
            Sigma_max.append(sigma_max)
            Sigma_min.append(sigma_min)
            Sigma_ave.append(sigma_ave)
            K_bond_max.append(k_bond_max)
            K_bond_min.append(k_bond_min)
            K_bond_ave.append(k_bond_ave)
            R0_max.append(r0_max)
            R0_min.append(r0_min)
            R0_ave.append(r0_ave)
            K_ang_max.append(k_ang_max)
            K_ang_min.append(k_ang_min)
            K_ang_ave.append(k_ang_ave)
            Theta0_max.append(theta0_max)
            Theta0_min.append(theta0_min)
            Theta0_ave.append(theta0_ave)
            K_dih_max.append(k_dih_max)
            K_dih_min.append(k_dih_min)
            K_dih_ave.append(k_dih_ave)
            success.append ("success")
        except:
            success.append ("fail")
            pass

    a = [x for x in Mass_max]
    b = [x for x in Mass_min]
    c = [x for x in Mass_ave]
    d = [x for x in Charge_max]
    e = [x for x in Charge_min]
    f = [x for x in Charge_ave]
    g = [x for x in Epsilon_max]
    h = [x for x in Epsilon_min]
    i = [x for x in Epsilon_ave]
    j = [x for x in Sigma_max]
    k = [x for x in Sigma_min]
    l = [x for x in Sigma_ave]
    m = [x for x in K_bond_max]
    n = [x for x in K_bond_min]
    o = [x for x in K_bond_ave]
    p = [x for x in R0_max]
    q = [x for x in R0_min]
    r = [x for x in R0_ave]
    s = [x for x in K_ang_max]
    t = [x for x in K_ang_min]
    u = [x for x in K_ang_ave]
    v = [x for x in Theta0_max]
    w = [x for x in Theta0_min]
    x = [x for x in Theta0_ave]
    y = [x for x in K_dih_max]
    z = [x for x in K_dih_min]
    aa = [x for x in K_dih_ave]
    name=[x for x in NAME]
    smiles=[x for x in SMILE]

    dataframe = pd.DataFrame({'ID':name,'Smiles':smiles,'Mass_max':a,'Mass_min':b,'Mass_ave':c,'Charge_max':d,'Charge_min':e,\
        'Charge_ave':f,'Epsilon_max':g,'Epsilon_min':h,'Epsilon_ave':i,'Sigma_max':j,'Sigma_min':k,'Sigma_ave':l,'K_bond_max':m,'K_bond_min':n,'K_bond_ave':o,\
            'R0_max':p,'R0_min':q,'R0_ave':r,'K_ang_max':s,'K_ang_min':t,'K_ang_ave':u,'Theta0_max':v,'Theta0_min':w,'Theta0_ave':x,'K_dih_max':y,\
                'K_dih_min':z,'K_dih_ave':aa})
    data_name="Des_MD.csv"
    #judge_name='judge_Des_MD.txt"
    dataframe.to_csv(path+data_name, index=False)
    #np.savetxt(judge_name,success,fmt='%s')