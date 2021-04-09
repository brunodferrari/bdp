# -*- coding: utf-8 -*-
"""
Created on Thu Apr  8 23:57:08 2021

@author: bferrari

#diferencas tabelas e conjunto de testes
"""

import bgraph
import pandas as pd
import numpy as np

df_results = pd.read_excel("C:/Users/bferrari/Desktop/pessoal/bdp/dbdp_instances/dbgdp_results.xlsx", sheet_name="intersec")
df_results['Crossing_init'] = np.nan

for i, row in enumerate(df_results['Instance']):
    path="C:/Users/bferrari/Desktop/pessoal/bdp/dbdp_instances/stallman_reduced/{name}.txt".format(name=row)
    graph_data = pd.read_csv(path)
    
    graph_data_exp = graph_data.iloc[:,0].str.split(" ",expand=True)
    n_1 = int(graph_data_exp.iloc[0,0])
    n_2 = int(graph_data_exp.iloc[0,1])
    
    graph_edges = []
    for key, row in enumerate(graph_data.iloc[1:n_1+1, 0]):
        for adj in row.split(" ")[2:]:
            #if int(adj) not in graph_adj_nodes_incre: 
            graph_edges.append(( key, int(adj)))
            
    order_1 = dict(zip(range(0,n_1),graph_data_exp.iloc[1:n_1+1,1].astype(int).values + 1))
    order_2 = dict(zip(range(n_1,n_1+n_2),graph_data_exp.iloc[n_1+1:n_1+1+n_2,1].astype(int).values + 1))
    
    New = bgraph.BGraph()
    New.v1(list(order_1.keys()))
    New.v2(list(order_2.keys()))
    New.pi_1 = order_1
    New.pi_2 = order_2
    New.order_v1()
    New.order_v2()
    New.edges(graph_edges)
    
    df_results['Crossing_init'][i] = New.n_cross()
    
df_results['diff'] = df_results['Crossings'] - df_results['Crossing_init']

with pd.ExcelWriter("C:/Users/bferrari/Desktop/pessoal/bdp/dbdp_instances/diff.xlsx") as writer:
    df_results.to_excel(writer, sheet_name="intersec", index=False)