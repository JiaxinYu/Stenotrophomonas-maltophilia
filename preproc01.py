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
  
  
min_df1 = mergedata(df_mi, decoding_df)
lvx_df1 = mergedata(df_lvx, decoding_df)
caz_df1 = mergedata(df_caz, decoding_df)
sxt_df1 = mergedata(df_sxt, decoding_df)


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
  
  
check01 = pd.DataFrame([])
for pathidx in [file-to-path01, file-to-path02]:
    check01 = pd.concat([check01, getfileaccnum(pathidx)], 0)

check01.columns = ['fullpath', 'foldername', 'filename']
check01


## min
min_df2 = min_df1.merge(check01, how='inner', on=['foldername', 'filename'])
min_df2.drop_duplicates(subset=['Accession #'], inplace=True)
min_df2.reset_index(drop=True, inplace=True)

## lvx
lvx_df2 = lvx_df1.merge(check01, how='inner', on=['foldername', 'filename'])
lvx_df2.drop_duplicates(subset=['Accession #'], inplace=True)
lvx_df2.reset_index(drop=True, inplace=True)

## caz
caz_df2 = caz_df1.merge(check01, how='inner', on=['foldername', 'filename'])
caz_df2.drop_duplicates(subset=['Accession #'], inplace=True)
caz_df2.reset_index(drop=True, inplace=True)

## sxt
sxt_df2 = sxt_df1.merge(check01, how='inner', on=['foldername', 'filename'])
sxt_df2.drop_duplicates(subset=['Accession #'], inplace=True)
sxt_df2.reset_index(drop=True, inplace=True)


check01.to_csv(rootpath/'msdata_fullpath_18jul22.csv', index=False)
min_df2.to_csv(rootpath/'sm_min_clean18jul22.csv', index=False)
lvx_df2.to_csv(rootpath/'sm_lvx_clean18jul22.csv', index=False)
caz_df2.to_csv(rootpath/'sm_caz_clean18jul22.csv', index=False)
sxt_df2.to_csv(rootpath/'sm_sxt_clean18jul22.csv', index=False)
