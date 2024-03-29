# -*- coding: utf-8 -*-
"""
Created on Sat Mar 27 17:24:59 2021

@author: Bruno Ferrari

The GRASPH META-HEURISTC for bdp works in the two main phases:
    -> Construction Phase:
        Create 2 list of unassigned vertices, with U_1 + U_2 = U = V
        assing pi(v) = 0 for all v in V
        select random vertice those vertices with maximum degree
        U - v (remove vertice)
        
        Select a vertice in U_linha =  v in U | deg(v) > alpha deg_max, where degree is calc. with respect V - U
        
        
    -> Improvement Phase


"""
#%cd C:\Users\bferrari\Desktop\pessoal\bdp\python

import random
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import networkx as nx

#from python 
#import bgraph

from datetime import datetime
inicio = datetime.now()
inicio_time = inicio.strftime("%H:%M:%S")
print("Current Time =", inicio_time)


def greedy_selection(G, alpha=1.0, nodelist=None, subgraph=None):
    
    deg_dict = G.degree(nodelist, subgraph)
    deg_max = max(deg_dict.values())
    mask = np.array(list(deg_dict.values()), dtype=int) >= alpha * deg_max
    idx = np.random.choice( np.where(mask)[0] )
    v = list(deg_dict.keys())[idx]
    
    return v

def construction_phase(G, U_1, U_2, V, alpha=1.0):
   
    v = greedy_selection(G)
    try:
        U_1.remove(v)
        G.pi_1[v] = 1
    except ValueError:
        U_2.remove(v)
        G.pi_2[v] = 1
    
    a = alpha #alpha values
    while len(U_1 + U_2) > 0:
        v = greedy_selection(G, a, U_1+U_2, list(set(V) - set(U_1+U_2))  )    
    
        try:
            U_1.remove(v)
            bc_v = G.bc(v, k=1)
            
            fbc_v = np.floor(bc_v).astype(int)
            cbc_v = np.ceil(bc_v).astype(int)
            
            pos_assing = list(G.pi_1.values())
            
            f_ncross = np.inf
            c_ncross = np.inf
            
            if fbc_v not in pos_assing:
                G.pi_1[v] = fbc_v
                f_ncross = G.n_cross()
            elif fbc_v - 1 > 0:
                fbc_v = fbc_v - 1
                G.pi_1[v] = fbc_v
                f_ncross = G.n_cross()
                
            if cbc_v not in pos_assing:
                G.pi_1[v] = cbc_v
                c_ncross = G.n_cross()
            elif cbc_v + 1 < G.n_v1():
                cbc_v = cbc_v + 1
                G.pi_1[v] = cbc_v
                c_ncross = G.n_cross()
            #verificar se é necessaria essa parte
            pos = min( [(c_ncross, cbc_v),(f_ncross, fbc_v)]  )[1]
            if pos in pos_assing: pos = max(pos_assing) + 1
            G.pi_1[v] = pos
            
        except ValueError:
            U_2.remove(v)
            bc_v = G.bc(v, k=2)
            
            fbc_v = np.floor(bc_v).astype(int)
            cbc_v = np.ceil(bc_v).astype(int)
            
            pos_assing = list(G.pi_2.values())
            
            f_ncross = np.inf
            c_ncross = np.inf
            if fbc_v not in pos_assing:
                G.pi_2[v] = fbc_v
                f_ncross = G.n_cross()
            elif fbc_v - 1 > 0:
                fbc_v = fbc_v - 1
                G.pi_2[v] = fbc_v
                f_ncross = G.n_cross()
                
            if cbc_v not in pos_assing:
                G.pi_2[v] = cbc_v
                c_ncross = G.n_cross()
            elif cbc_v + 1 < G.n_v2():
                cbc_v = cbc_v + 1
                G.pi_2[v] = cbc_v
                c_ncross = G.n_cross()
            #verificar se é necessaria essa parte
            pos = min( [(c_ncross, cbc_v),(f_ncross, fbc_v)]  )[1]
            if pos in pos_assing: pos = max(pos_assing) + 1
            G.pi_2[v] = pos

def improvement_phase(G, U_1, U_2, V, verbose=0):
    Pr_list = list(G.degree().values()) / np.array(list(G.degree().values())).sum()
    Pr_dict = dict(zip(G.degree().keys(),  Pr_list) ) #vetor de probabilidades
    
    while len(U_1 + U_2) > 0:
        
        v = random.choices(list(Pr_dict.keys()),weights=list(Pr_dict.values()))[0]
        Pr_dict.pop(v)
        cross_i = G.n_cross()
        
        #if verbose: 
            #print("move v: ", v, " n_cross: ", cross_i)
        
        try:
            U_1.remove(v)
            bc_v = G.bc(v, k=1)
               
            fbc_v = np.floor(bc_v).astype(int)
            cbc_v = np.ceil(bc_v).astype(int)
            pi_aux = G.pi_1.copy()
        
            pi_f = G.move_v1(v, fbc_v)
            G.pi_1 = pi_f
            f_ncross =  G.n_cross()
            
            pi_c = G.move_v1(v, cbc_v)
            G.pi_1 = pi_c
            c_ncross = G.n_cross()
            
            f_1_ncross = np.inf
            pi_f_1 = pi_aux
            if fbc_v - 1 > 0:
                pi_f_1 = G.move_v1(v, fbc_v-1)
                G.pi_1 = pi_f_1
                f_1_ncross = G.n_cross()
            
            c_1_ncross = np.inf
            pi_c_1 = pi_aux
            if cbc_v + 1 < G.n_v1():
                pi_c_1 = G.move_v1(v, cbc_v + 1)
                G.pi_1 = pi_c_1
                c_1_ncross = G.n_cross()
            
            G.pi_1 = min( [(c_ncross, pi_c),(f_ncross, pi_f), \
                             (f_1_ncross, pi_f_1), (c_1_ncross, pi_c_1), (cross_i, pi_aux)], key=lambda x: x[0])[1]
        except ValueError:
            U_2.remove(v)
            bc_v = G.bc(v, k=2)
               
            fbc_v = np.floor(bc_v).astype(int)
            cbc_v = np.ceil(bc_v).astype(int)
            pi_aux = G.pi_2.copy()
        
            pi_f = G.move_v2(v, fbc_v)
            G.pi_2 = pi_f
            f_ncross =  G.n_cross()
            
            pi_c = G.move_v2(v, cbc_v)
            G.pi_2 = pi_c
            c_ncross = G.n_cross()
            
            f_1_ncross = np.inf
            pi_f_1 = pi_aux
            if fbc_v - 1 > 0:
                pi_f_1 = G.move_v2(v, fbc_v-1)
                G.pi_2 = pi_f_1
                f_1_ncross = G.n_cross()
            
            c_1_ncross = np.inf
            pi_c_1 = pi_aux
            if cbc_v + 1 < G.n_v2():
                pi_c_1 = G.move_v2(v, cbc_v + 1)
                G.pi_2 = pi_c_1
                c_1_ncross = G.n_cross()
            
            G.pi_2 = min( [(c_ncross, pi_c), (f_ncross, pi_f), \
                             (f_1_ncross, pi_f_1), (c_1_ncross, pi_c_1), (cross_i, pi_aux)], key=lambda x: x[0] )[1]

def grasp(G, alpha=1.0, max_it=5, verbose=0):
    
    V = G.v1() + G.v2()
    min_c = np.inf
    it = 1
    while it < max_it:
        
        U_1 = G.v1().copy()
        U_2 = G.v2().copy()
        construction_phase(G, U_1, U_2, V, alpha)
        G.order_v1()
        G.order_v2()
        U_1 = G.v1().copy()
        U_2 = G.v2().copy()
        improvement_phase(G, U_1, U_2, V, verbose)
        C = G.n_cross()
        if verbose: 
            print(it, ":", str(min_c), "::", C, ":::", G.n_cross())
        
        if C < min_c:
            G_min = G.copy()
            min_c = C
        else:
            it+=1
    
    G.pi_1 = G_min.pi_1
    G.pi_2 = G_min.pi_2
###################################    
# def grasp(G, alpha=1.0, verbose=0):
    
#     U_1 = G.v1().copy()
#     U_2 = G.v2().copy()
#     V = U_1 + U_2
    
#     construction_phase(G, U_1, U_2, V, alpha)
#     #G.order_v1()
#     #G.order_v2()

#     min_cross_ant = G.n_cross()
#     min_cross = min_cross_ant - 1
#     while min_cross < min_cross_ant:
    
#         U_1 = G.v1().copy()
#         U_2 = G.v2().copy()   
        
#         improvement_phase(G, U_1, U_2, V, verbose)
#       #   G.order_v1()
#       #  G.order_v2()
        
#         min_cross_ant = min_cross
#         min_cross = G.n_cross()
        
#     #G.order_v1()
#     #G.order_v2()
######################################      

if __name__=='__main__':
    #graph_data = pd.read_csv("C:/Users/bferrari/Desktop/pessoal/bdp/dbdp_instances/instances/incgraph_25_25_0.3_0.2_1.txt")
    graph_data = pd.read_csv("C:/Users/bferrari/Desktop/pessoal/bdp/dbdp_instances/stallman_reduced/G_00_02_scr_0002_30.txt")
    graph_adj_nodes_incre = graph_data.iloc[-6:,0].str.split(" ",expand=True).iloc[:, 1].astype(int).to_list()
    
    graph_data_exp = graph_data.iloc[:,0].str.split(" ",expand=True)
    n_1 = int(graph_data_exp.iloc[0,0])
    n_2 = int(graph_data_exp.iloc[0,1])
    
    graph_edges = []
    for key, row in enumerate(graph_data.iloc[1:n_1+1, 0]):
        for adj in row.split(" ")[2:]:
            #if int(adj) not in graph_adj_nodes_incre: 
            graph_edges.append(( key, int(adj)))
    
    graph_edges2 = []
    for row in graph_data.iloc[n_1+1:32, 0]:
        for adj in row.split(" ")[2:]:
            if int(adj) not in graph_adj_nodes_incre: 
                graph_edges2.append(( int(row.split(" ")[1]), int(adj)))
    
    order_1 = dict(zip(range(0,n_1),graph_data_exp.iloc[1:n_1+1,1].astype(int).values + 1))
    order_2 = dict(zip(range(n_1,n_1+n_2),graph_data_exp.iloc[n_1+1:n_1+1+n_2,1].astype(int).values + 1))
    
    
    New = bgraph.BGraph()
    #New.v2(np.unique(np.array(np.matrix(graph_edges)[:,1])).tolist())
    #New.v1(np.unique(np.array(np.matrix(graph_edges)[:,0])).tolist())
    
    New.v1(list(order_1.keys()))
    New.v2(list(order_2.keys()))
    New.pi_1 = order_1
    New.pi_2 = order_2
    New.order_v1()
    New.order_v2()
    New.edges(graph_edges)
    
    bgraph.plotBGraph(New, size=20/2, height=150)
    plt.title("N crossing: "+ str(bgraph.crossing(New)) )
    plt.show()
    
    U_1 = New.v1().copy()
    U_2 = New.v2().copy()
    V = U_1 + U_2
    
    #U_list = New.v1() + New.v2()    
    New.pi_1 = dict(zip(U_1, New.n_v1() * [0]))
    New.pi_2 = dict(zip(U_2, New.n_v2() * [0]))
    
    construction_phase(New, U_1, U_2, V, 1)
    New.order_v1()
    New.order_v2()
    bgraph.plotBGraph(New)
    plt.title("N crossing: "+ str(bgraph.crossing(New)) )
    plt.show()
    
    
    
    
    min_cross_ant = New.n_cross()
    min_cross = min_cross_ant - 1
    while min_cross < min_cross_ant:
        
        U_1 = New.v1().copy()
        U_2 = New.v2().copy()   
        
        improvement_phase(New, U_1, U_2, V)
        New.order_v1()
        New.order_v2()
        
        min_cross_ant = min_cross
        min_cross = bgraph.crossing(New)
        
        
        bgraph.plotBGraph(New)
        plt.title("N crossing: "+ str(min_cross) )
        plt.show()
    
# New.degree(U_linha, [2])

# G = nx.bipartite.gnmk_random_graph(3,5,10, seed=111)

# %%timeit -n 10
# mask = np.array(list(New.degree().values()), dtype=int) >= 10
# v = list(New.degree().keys())[ np.random.choice(np.where(mask)[0])  ]

# ## aux
# New2=nx.Graph()
# New2.add_nodes_from(New.v1()+New.v2())
# New2.add_edges_from( graph_edges )



# aux = DegreeView(G)

# nx.Graph()


# np.array(list(New.degree().items()), dtype=int)