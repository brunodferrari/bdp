# -*- coding: utf-8 -*-
"""


@author: Bruno Ferrari
"""
%cd C:\Users\Bruno Ferrari\Documents\Bruno\2019\2s\MC\artigos revisão\Artigos Mes\GD\bdp

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import networkx as nx
from python import bgraph

B = bgraph.BGraph()
B.v1([5,3,2,4])
B.v2(['a','b','c'])
B.edges([(5, "a"), (5, "b"), (2, "b"), (2, "c"), (3, "c"), (4, "a")])
bgraph.plotBGraph(B)
B.n_cross()

graph_data = pd.read_csv("./dbdp_instances/instances/incgraph_25_25_0.3_0.2_1.txt")

graph_adj_nodes_incre = graph_data.iloc[-6:,0].str.split(" ",expand=True).iloc[:, 1].astype(int).to_list()

graph_edges = []
for row in graph_data.iloc[1:26, 0]:
    for adj in row.split(" ")[2:]:
        if int(adj) not in graph_adj_nodes_incre: 
            graph_edges.append(( int(row.split(" ")[1]), int(adj)))

graph_edges2 = []
for row in graph_data.iloc[1:32, 0]:
    for adj in row.split(" ")[2:]:
        if int(adj) not in graph_adj_nodes_incre: 
            graph_edges2.append(( int(row.split(" ")[1]), int(adj)))

New = bgraph.BGraph()
New.v2(np.unique(np.array(np.matrix(graph_edges)[:,1])).tolist())
New.v1(np.unique(np.array(np.matrix(graph_edges)[:,0])).tolist())

stop = 0
i = 0 
while stop==0:
    New.edges(graph_edges[:i])
    bgraph.plotBGraph(New)
    plt.show()
    print('N crossing:', New.n_cross())
    i+=2
    stop=int(input())

New.edges(graph_edges)
bgraph.plotBGraph(New)
print('N crossing:', New.n_cross())

stop = 0
while stop==0:
    New.v1(bary_sort(np.vectorize(bary)(New, New.v1(), 1), New.v1()))
    bgraph.plotBGraph(New)
    plt.show()
    print('N crossing:', New.n_cross())
    
    New.v2(bary_sort(np.vectorize(bary)(New, New.v2(), 2), New.v2()))
    bgraph.plotBGraph(New)
    plt.show()
    print('N crossing:', New.n_cross())
    
    stop=int(input())


#New.v1(order_i)
#bgraph.plotBGraph(New)