#!/usr/bin/env python
# coding: utf-8

import pandas as pd

path="./dataset/"
file0="Des_monomer.csv"
file1="Des_MD.csv"
file2="Des_Mordred.csv"

df=pd.read_csv(path+file0)
df1=pd.read_csv(path+file1)
df2=pd.read_csv(path+file2)
df11=df1.drop(df1.columns[0:2], axis=1)
df21=df2.drop(df2.columns[0:3], axis=1)
df_des=df.join(df11)
df_des=df_des.join(df21)
df_des.to_csv(path+"Des_optimized.csv",index=None)
