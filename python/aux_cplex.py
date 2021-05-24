# -*- coding: utf-8 -*-
"""
Created on Tue May  4 23:01:21 2021

@author: Bruno Ferrari
"""

import pandas as pd
import numpy as np
import networkx as nx

from python import bgraph
from python.grasp import grasp as gs

df_results = pd.read_excel("./dbdp_instances/instances.xlsx", sheet_name="instances").replace(regex=".txt", value="")

for i, inst in enumerate(df_results['Instance']):
    path="./dbdp_instances/instances/{name}.txt".format(name=inst)
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
    gs(New, 0.9, 1) 

    r2 = []
    while New.n_cross() >= 42:
        gs(New, 0.9, 1) 
        New.order_v1()
        New.order_v2()
        r2.append(New.n_cross())

'''
gerar dados para o cplex. grafo sera dado pela matrix de adjacencia.

1o teste: 
matrix de adjacencia gerada da seguinte maneira: coluna 1 do dados define a posicao do vertice no layer, o restante eh a lista d adj.

aux=graph_data_exp.iloc[1:,1:]
pd.unique(aux.values.ravel('K'))
t=np.delete(a, np.where(a == None)).astype(int)
'''

madj = np.zeros((31*2,2*31))
for a in graph_edges:
    #print(a[0],a[1])
    madj[a[0],a[1]]=1
    madj[a[1],a[0]]=1
    
text_file = open("sample.txt", "w")
n = text_file.write(str(madj.astype(int).tolist()))
text_file.close()
'''
2o teste
matrix de adjacencia gerada da seguinte maneira: coluna 1 do dados o vertice inicial da lista de adj.
'''


