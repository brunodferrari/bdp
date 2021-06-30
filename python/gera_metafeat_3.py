# -*- coding: utf-8 -*-
"""
Created on Mon Jun 28 23:36:19 2021

@author: Bruno Ferrari
"""
_path="C:/Users/Bruno Ferrari/Documents/Bruno/2019/2s/MC/artigos revisão/Artigos Mes/GD/"

import pandas as pd
import numpy as np
import bgraph

from concurrent.futures import ThreadPoolExecutor

def gera_txt(nome, G):
    with open((_path+'bdp/dbdp_instances/GraphData/{inst}.txt').format(inst=nome), 'w') as arq:
        arq.write(str(G.v1()) + "\n")
        arq.write(str(G.v2()) + "\n")
        arq.write(str(G.edges()))
        
def gera_grafo(path):
    graph_data = pd.read_csv(path)
    
    graph_data_exp = graph_data.iloc[:,0].str.split(" ",expand=True)
    n_1 = int(graph_data_exp.iloc[0,0])
    n_2 = int(graph_data_exp.iloc[0,1])
    
    graph_edges = []
    for key, row in enumerate(graph_data.iloc[1:n_1+1, 0]):
        for adj in row.split(" ")[2:]:
            #if int(adj) not in graph_adj_nodes_incre: 
            graph_edges.append((key, int(adj)))
            
    order_1 = dict(zip(range(0,n_1),graph_data_exp.iloc[1:n_1+1,1].astype(int).values + 1))
    order_2 = dict(zip(range(n_1,n_1+n_2),graph_data_exp.iloc[n_1+1:n_1+1+n_2,1].astype(int).values + 1))
    
    New = bgraph.BGraph()
    New.v1(list(order_1.keys()))
    New.v2(list(order_2.keys()))
    New.edges(graph_edges)
    
    return New

# Descosidera vertices adicionais (ie com linha == 0)
def gera_grafo_2(path):

    graph_data = pd.read_csv(path)
    
    graph_data_exp = graph_data.iloc[:,0].str.split(" ",expand=True)
    n_1 = int(graph_data_exp.iloc[0,0])
    
    aux = (graph_data_exp[1:n_1+1].set_index(0).loc['1'])
    graph_edges = [(int(v), int(u)) for v in aux.iloc[:,0] for u in aux.set_index(1).loc[v] if u != None]
    
    New = bgraph.BGraph()
    New.v1(np.unique(np.array(np.matrix(graph_edges)[:,0])).tolist())
    New.v2(np.unique(np.array(np.matrix(graph_edges)[:,1])).tolist())
    New.edges(graph_edges)
    
    return New

def gera_features(output):
    for o in output:
        _aux_deg = list(o.degree().values())
        yield (o.n_v(),                         #vertices
               o.n_edge(),                      #arestas
               abs(o.n_v1() - o.n_v2()),        #diff nos
               np.mean(_aux_deg),               #mean deg
               np.std(_aux_deg),                #std deg
               np.median(_aux_deg),             #median deg                  
               o.density(),                     #graph density
               o.n_v1(), o.n_v2(), np.min(_aux_deg)) #aux feat


#df_results = pd.read_excel(_path+"bdp/dbdp_instances/stallman_reduced.xlsx", sheet_name="Sheet1")
df_results = pd.read_excel(_path+"bdp/dbdp_instances/instance.xlsx", sheet_name="instances").replace(regex=".txt", value="")
output_list = []
for i in range(130)[::10]:
    with ThreadPoolExecutor(6) as ex:
         output = list(ex.map(lambda x: gera_grafo_2((_path+"bdp/dbdp_instances/instances/{name}.txt").format(name=x)), list(df_results['Instance'][i-10:i])))
    with ThreadPoolExecutor(6) as ex:
        ex.map(lambda x: gera_txt(x[0], x[1]), zip(df_results['Instance'][i-10:i].replace(regex="inc", value=""), output))
    output_list.append(pd.DataFrame(gera_features(output)))
    if not(i % 20):
        with pd.ExcelWriter(_path+"bdp/dbdp_instances/metafeat3.xlsx") as writer:
            pd.concat([df_results.replace(regex="inc", value=""), pd.concat(output_list, axis=0).reset_index(drop=True)],axis=1).to_excel(writer, sheet_name="results", index=False)
            

