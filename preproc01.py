import numpy as np
import pandas as pd
from pathlib import Path
from datetime import datetime
import nmrglue as ng
import ray
from multiprocessing import cpu_count
import os
from pathlib import Path
import json
import pandas as pd

rootpath = Path(rootpath)

df_mi = pd.read_csv(rootpath/'cmuh/Minocycline.csv')
df_lvx = pd.read_csv(rootpath/'cmuh/Levofloxacin.csv')
df_caz = pd.read_csv(rootpath/'cmuh/Ceftazidime.csv')
df_sxt = pd.read_csv(rootpath/'cmuh/Trimethoprim-Sulfamethoxazole.csv')


## 2018-2019
directory = Path(path-to-file)
df_filepath = pd.read_csv(directory/'mbt_0318.csv')
# drop QC files
df_filepath.drop(df_filepath[df_filepath['accnum'].str.len() < 10].index, axis=0, inplace=True)
df_filepath.reset_index(drop=True, inplace=True)
# extract accnum
temp_test = df_filepath['accnum'].str.split('#', expand=True)
df_filepath = pd.concat([df_filepath, temp_test], 1)
df_filepath.columns = ['date', 'accum', 'Accession #', 'case #']
# drop more than one isolates in a single accnum
df_filepath.drop_duplicates(subset=['Accession #'], keep=False, inplace=True)
# # include dataset before 2020
df_filepath['date'] = pd.to_numeric(df_filepath['date'], errors='coerce', downcast="integer")
df_filepath.dropna(subset=['date'], inplace=True)
df_filepath.drop(df_filepath[df_filepath['date'] > 109000000].index, 0, inplace=True)
df_filepath.reset_index(drop=True, inplace=True)
df_filepath.to_csv(rootpath/'decoding_18_19.csv', index=False)
df_filepath


## 2020-2021
IndexDir = []
directory = Path(path-to-file)
for fidx in os.listdir(directory):
    if os.path.isfile(directory/fidx/'runInfo.json'):
        f = open(directory/fidx/'runInfo.json', 'r')
        data = json.load(f)
        for i in range(len(data["Analytes"])):
            IndexDir.append([data["ProjectUid"], data["Analytes"][i]["AnalyteUid"], 
                            data["ProjectName"], data["Analytes"][i]["AnalyteId"]])
        f.close()

decoding_2020 = pd.DataFrame(IndexDir, columns=['foldername', 'filename', 'date', 'accnum'])

# clean up unrelated data
decoding_2020.drop(decoding_2020[decoding_2020['accnum'].str.len() < 10].index, axis=0, inplace=True)
decoding_2020.reset_index(drop=True, inplace=True)
decoding_2020

IndexDir = []
directory = Path(path-to-file)
for fidx in os.listdir(directory):
    if os.path.isfile(directory/fidx/'runInfo.json'):
        f = open(directory/fidx/'runInfo.json', 'r')
        data = json.load(f)
        for i in range(len(data["Analytes"])):
            IndexDir.append([data["ProjectUid"], data["Analytes"][i]["AnalyteUid"], 
                            data["ProjectName"], data["Analytes"][i]["AnalyteId"]])
        f.close()

decoding_2021 = pd.DataFrame(IndexDir, columns=['foldername', 'filename', 'date', 'accnum'])

decoding_2021.drop(decoding_2021[decoding_2021['accnum'].str.len() < 10].index, axis=0, inplace=True)
decoding_2021.reset_index(drop=True, inplace=True)
decoding_2021


## meraging 2020-2021 dataset
decoding_df = pd.concat([decoding_2020, decoding_2021], 0)
decoding_df.drop_duplicates(inplace=True)
decoding_df.reset_index(drop=True, inplace=True)
decoding_df['Accession #'] = decoding_df['accnum'].str.split('#', expand=True)[0]
decoding_df.to_csv(rootpath/'decoding_20_21.csv', index=False)


def mergedata(df1, df2):
    temp_df1 = df2.merge(df1.astype('str'), how='inner', on=['Accession #'])
    temp_df1.drop_duplicates(subset=['Accession #'], keep=False, inplace=True)
    temp_df1.reset_index(drop=True, inplace=True)
    return temp_df1
  
  
## 2020 and 2021
min_df1 = mergedata(df_mi, decoding_df)
lvx_df1 = mergedata(df_lvx, decoding_df)
caz_df1 = mergedata(df_caz, decoding_df)
sxt_df1 = mergedata(df_sxt, decoding_df)


## 2019 and before
min_df2 = mergedata(df_mi, df_filepath)
lvx_df2 = mergedata(df_lvx, df_filepath)
caz_df2 = mergedata(df_caz, df_filepath)
sxt_df2 = mergedata(df_sxt, df_filepath)

for idx in [min_df2, lvx_df2, caz_df2, sxt_df2]:
    idx['foldername'] = idx['date'].astype('int')
    idx['filename'] = idx['accum']
    
    
## merge all 4 years
min_df0 = pd.concat([min_df1, min_df2], axis=0, join='inner')
lvx_df0 = pd.concat([lvx_df1, lvx_df2], axis=0, join='inner')
caz_df0 = pd.concat([caz_df1, caz_df2], axis=0, join='inner')
sxt_df0 = pd.concat([sxt_df1, sxt_df2], axis=0, join='inner')


def getfileaccnum(rootpath):
    # collect info from files
    path_all = []
    directory = rootpath
    for filepath in Path(directory).glob('**/*'):
        getpath = str(filepath.absolute())
        if getpath.endswith('1SLin/pdata/1'):
            path_all.append(filepath.absolute())

    path_all = pd.DataFrame(path_all, columns=['fullpath'])

    datapath = ((path_all['fullpath']).astype('str').str.replace(str(directory), ''))\
                .str.replace('/1SLin/pdata/1', '')\
                .str.split('/', expand=True)
    return pd.concat([path_all, datapath.iloc[:, 1:3]], 1)
  
## 2019 and before
check2019 = getfileaccnum(Path(file-to-path00))
check2019.columns = ['fullpath', 'foldername', 'filename']
check2019.drop(check2019[check2019['filename'].str.len() < 10].index, axis=0, inplace=True)
check2019.reset_index(drop=True, inplace=True)


## 2020 and 2021
check2021 = pd.DataFrame([])
for pathidx in [file-to-path01, file-to-path02]:
    check2021 = pd.concat([check2021, getfileaccnum(pathidx)], 0)

check2021.columns = ['fullpath', 'foldername', 'filename']


## merge all fullpath
check18_21 = pd.concat([check2019, check2021], 0)

## save fullpath file to csv
check18_21.to_csv(rootpath/'msdata_fullpath18to21_25jul22.csv', index=False)


## min
min_df =  min_df0.astype('str').merge(check18_21, how='inner', on=['foldername', 'filename'])
min_df.drop_duplicates(subset=['Accession #'], inplace=True)
min_df.reset_index(drop=True, inplace=True)

## lvx
lvx_df = lvx_df0.astype('str').merge(check18_21, how='inner', on=['foldername', 'filename'])
lvx_df.drop_duplicates(subset=['Accession #'], inplace=True)
lvx_df.reset_index(drop=True, inplace=True)

## caz
caz_df = caz_df0.astype('str').merge(check18_21, how='inner', on=['foldername', 'filename'])
caz_df.drop_duplicates(subset=['Accession #'], inplace=True)
caz_df.reset_index(drop=True, inplace=True)

## sxt
sxt_df = sxt_df0.astype('str').merge(check18_21, how='inner', on=['foldername', 'filename'])
sxt_df.drop_duplicates(subset=['Accession #'], inplace=True)
sxt_df.reset_index(drop=True, inplace=True)


min_df.to_csv(rootpath/'sm_min_clean25jul22.csv', index=False)
lvx_df.to_csv(rootpath/'sm_lvx_clean25jul22.csv', index=False)
caz_df.to_csv(rootpath/'sm_caz_clean25jul22.csv', index=False)
sxt_df.to_csv(rootpath/'sm_sxt_clean25jul22.csv', index=False)
