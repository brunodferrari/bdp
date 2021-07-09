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
from python.tabu_search import ts
from python.vns import VNS

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
    
    
    import time
    inicio = time.time()
    %%snakeviz
    GS = New.copy()
    %%snakeviz
    gs(GS, 0.9, 0) 
    fim = time.time()
    (fim - inicio)
    
    %%snakeviz
    #t = []
    #for _ in range(50):
        %%snakeviz
        TS = New.copy()
        inicio = time.time()
        %%snakeviz
        ts(TS, verbose=0, max_it=5)
        print(TS.n_cross())
        fim = time.time()
        t.append(fim - inicio)

    t = []
    for _ in range(20):
        TS = New.copy()
        inicio = time.time()
        VNS(TS, alpha=1.0, verbose=0, k_max=5)
        print(TS.n_cross())
        fim = time.time()
        t.append(fim - inicio)
    
    param_grid = list(product(np.arange(6,10,1)/10, [5,10,15,20]))
    res_dict = {}
    for i, param in enumerate(param_grid):
        c = []
        t = []
        for _ in range(10):
            inicio = time.time()
            gs(GS, param[0], max_it=param[1], verbose=0)
            fim = time.time()
            t.append(fim - inicio)    
            c.append(GS.n_cross())
        print(i)
        res_dict[i] = (c,t)
    
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


'''

'''
idx=np.argsort(list(cplex.pi_1.keys()))
x=np.tri(len(cplex.pi_1)).T 
np.fill_diagonal(x,0)
idx=np.argsort(list(cplex.pi_1.keys()))
x = x[idx,:][:,idx]

idx=np.argsort(list(cplex.pi_2.keys()))
y=np.tri(len(cplex.pi_2)).T 
np.fill_diagonal(y,0)
idx=np.argsort(list(cplex.pi_2.keys()))
y = y[idx,:][:,idx]

madj = np.zeros((cplex.n_v(),cplex.n_v()))
for a in cplex.edges():
    #print(a[0],a[1])
    madj[a[0],a[1]-min(cplex.v2()) + max(cplex.v1()) + 1]=1
    madj[a[1]-min(cplex.v2()) + max(cplex.v1()) + 1,a[0]]=1
    
text_file = open("sample1.txt", "w")
n = text_file.write(str(madj.astype(int).tolist()))
text_file.close()

text_file = open("x.txt", "w")
n = text_file.write(str(x.astype(int).tolist()))
text_file.close()

text_file = open("y.txt", "w")
n = text_file.write(str(y.astype(int).tolist()))
text_file.close()