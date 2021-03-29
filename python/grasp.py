# -*- coding: utf-8 -*-
"""
Created on Sat Mar 27 17:24:59 2021

@author: Bruno Ferrari
"""

'''
The GRASPH META-HEURISTC for bdp works in the two main phases:
    -> Construction Phase:
        Create 2 list of unassigned vertices, with U_1 + U_2 = U = V
        assing pi(v) = 0 for all v in V
    
    -> Improvement Phase


'''

#from python import bgraph

import networkx as nx
import random
from networkx.classes.reportviews import DegreeView


def greedy_selection(G, alpha=1.0, nodelist = None, subgraph = None):
    
    deg_dict = G.degree(nodelist, subgraph)
    deg_max = max(deg_dict.values())
    mask = np.array(list(deg_dict.values()), dtype=int) >= alpha * deg_max
    idx = np.random.choice( np.where(mask)[0] )
    v = list(deg_dict.keys())[idx]
    
    return v

U_1 = New.v1().copy()
U_2 = New.v2().copy()
V = U_1 + U_2

#U_list = New.v1() + New.v2()    
New.pi_1 = dict(zip(U_1, New.n_v1() * [0]))
New.pi_2 = dict(zip(U_2, New.n_v2() * [0]))

v = greedy_selection(New)
    
try:
    U_1.remove(v)
    New.pi_1[v] = 1
except ValueError:
    U_2.remove(v)
    New.pi_2[v] = 1

a = 1.0 #alpha values
while len(U_1 + U_2) > 0:
    v = greedy_selection(New, a, U_1+U_2, list(set(V) - set(U_1+U_2))  )    

    try:
        U_1.remove(v)
        bc_v = New.bc(v, k=1)
        
        fbc_v = np.floor(bc_v).astype(int)
        cbc_v = np.ceil(bc_v).astype(int)
        
        pos_assing = list(New.pi_2.values())
        
        if fbc_v not in pos_assing:
            New.pi_1[v] = fbc_v
            f_ncross = New.n_cross()
        elif fbc_v - 1 > 0:
            fbc_v = fbc_v
            New.pi_1[v] = fbc_v
            f_ncross = New.n_cross()
            
        if cbc_v not in pos_assing:
            New.pi_1[v] = cbc_v
            c_ncross = New.n_cross()
        elif cbc_v + 1 < New.n_v1():
            cbc_v = cbc_v + 1
            New.pi_1[v] = cbc_v
            c_ncross = New.n_cross()
            
        New.pi_1[v] = min( [(c_ncross, cbc),(f_ncross, fbc)]  )[1]
    except ValueError:
        U_2.remove(v)
        bc_v = New.bc(v, k=2)
        
        fbc_v = np.floor(bc_v).astype(int)
        cbc_v = np.ceil(bc_v).astype(int)
        
        pos_assing = list(New.pi_2.values())
        
        if fbc_v not in pos_assing:
            New.pi_2[v] = fbc_v
            f_ncross = New.n_cross()
        elif fbc_v - 1 > 0:
            fbc_v = fbc_v - 1
            New.pi_2[v] = fbc_v
            f_ncross = New.n_cross()
            
        if cbc_v not in pos_assing:
            New.pi_2[v] = cbc_v
            c_ncross = New.n_cross()
        elif cbc_v + 1 < New.n_v2():
            cbc_v = cbc_v + 1
            New.pi_2[v] = cbc_v
            c_ncross = New.n_cross()
            
        New.pi_2[v] = min( [(c_ncross, cbc_v),(f_ncross, fbc_v)]  )[1]


    New.degree(U_linha, [2])

G = nx.bipartite.gnmk_random_graph(3,5,10, seed=111)

%%timeit -n 10
mask = np.array(list(New.degree().values()), dtype=int) >= 10
v = list(New.degree().keys())[ np.random.choice(np.where(mask)[0])  ]

## aux
New2=nx.Graph()
New2.add_nodes_from(New.v1()+New.v2())
New2.add_edges_from( graph_edges )



aux = DegreeView(G)

nx.Graph()


np.array(list(New.degree().items()), dtype=int)