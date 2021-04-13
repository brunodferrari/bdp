# -*- coding: utf-8 -*-
"""
Created on Mon Apr 12 23:54:49 2021

@author: bferrari

automate test
"""

import bgraph
import time
import pandas as pd
import numpy as np

from grasp import grasp as gs


# inicio = datetime.now()
# inicio_time = inicio.strftime("%H:%M:%S")


df_results = pd.read_excel("C:/Users/bferrari/Desktop/pessoal/bdp/dbdp_instances/stallman_reduced.xlsx", sheet_name="Sheet1")[136:]
df_results['Crossing'] = np.nan
df_results['Time'] = np.nan


for i, inst in enumerate(df_results['Instance']):
    path="C:/Users/bferrari/Desktop/pessoal/bdp/dbdp_instances/stallman_reduced/{name}.txt".format(name=inst)
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
    
    print("\n" + inst + " : " + str(i))
    inicio = time.time()
    gs(New, alpha=1)
    fim = time.time()
    
    df_results['Crossing'][i] = New.n_cross()
    df_results['Time'][i] = (fim - inicio)
    
    
with pd.ExcelWriter("C:/Users/bferrari/Desktop/pessoal/bdp/dbdp_instances/grasph_2.xlsx") as writer:
    df_results.dropna()[:-1].to_excel(writer, sheet_name="results", index=False)