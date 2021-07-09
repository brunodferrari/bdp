# -*- coding: utf-8 -*-
"""
Created on Wed Jun 30 02:37:23 2021

@author: Bruno Ferrari
"""

_path="C:/Users/Bruno Ferrari/Documents/Bruno/2019/2s/MC/artigos revisão/Artigos Mes/GD/"

import bgraph
import time
import pandas as pd
import numpy as np

from tabu_search import ts
from concurrent.futures import ThreadPoolExecutor

# inicio = datetime.now()
# inicio_time = inicio.strftime("%H:%M:%S")

def carrega_grafo(path):
    
    G = bgraph.BGraph()
    with open(path) as arq:
        dados = arq.readlines()
    
    G.v1(eval(dados[0]))
    G.v2(eval(dados[1]))
    G.edges(eval(dados[2]))
    
    return G

df_results = pd.concat([pd.read_excel(_path+"bdp/dbdp_instances/metafeat.xlsx", sheet_name="results"),
                        #pd.read_excel(_path+"bdp/dbdp_instances/metafeat2.xlsx", sheet_name="results"),
                        #pd.read_excel(_path+"bdp/dbdp_instances/metafeat3.xlsx", sheet_name="results")
                        ], axis=0).reset_index(drop=True)
df_results = df_results.set_index('Instance')
df_results['Crossing'] = np.nan
df_results['Time'] = np.nan

TS = carrega_grafo((_path+'bdp/dbdp_instances/GraphData/{name}.txt').format(name=df_results.index[1]))
ts(TS, verbose=0, max_it=5)
  
for i, inst in enumerate(df_results.index[162:]):
    TS = carrega_grafo((_path+'bdp/dbdp_instances/GraphData/{name}.txt').format(name=inst))
    inicio = time.time()
    ts(TS, verbose=0, max_it=5)
    fim = time.time()
    
    df_results.loc[inst, 'Crossing'] = TS.n_cross()
    df_results.loc[inst, 'Time'] = (fim - inicio)
    print(i)
    if not(i % 5):
        with pd.ExcelWriter(_path+"bdp/dbdp_instances/meta_tabu_3.xlsx") as writer:
            df_results.dropna().reset_index().to_excel(writer, sheet_name="results", index=False)