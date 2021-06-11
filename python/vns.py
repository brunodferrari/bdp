# -*- coding: utf-8 -*-
"""
Created on Fri Jun 11 00:05:43 2021

@author: Bruno Ferrari
"""
import random
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import networkx as nx

from heapq import nlargest as maxq
#import heapq as max_q

from python import bgraph

from datetime import datetime
inicio = datetime.now()
inicio_time = inicio.strftime("%H:%M:%S")
print("Current Time =", inicio_time)

global n_cross
n_cross = np.inf

def greedy_selection(G, alpha=1.0, nodelist=None, subgraph=None):
    
    deg_dict = G.degree(nodelist, subgraph)
    deg_max = maxq(1, deg_dict.values())[0]
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

def __auxplot(G_aux, v, sign=" +"):
    #print(v)
    G_aux.plot(order=1)
    plt.title(str(v) + sign)
    plt.show()

def _make_neighborhood(G, neighlist, k=1, verbose=0):
    
    if verbose:
        G_aux = G.copy()
        G_pi_1_org = G.pi_1
        G_pi_2_org = G.pi_2
    
    trsh_crossing = np.inf
    idx_min = np.nan
    for i, v in enumerate(neighlist):
        if v[0] in G.v1(): #Layer 1
            if (G.pi_1[v[0]] + k <= G.n_v1()) and np.sign(v[1])==1: 
            #if (m_plus is not None):
                chg_pos = G.find_pos(v[0], G.pi_1[v[0]]+k)
                K_vj = G.K(v[0], chg_pos) 
                K_jv = G.K(chg_pos, v[0])
                v[1] = k
                v[2] = v[2] + (K_jv - K_vj)
                if (v[2] < trsh_crossing):
                    idx_min = i
                    trsh_crossing = v[2]
                if verbose: 
                    G_aux.pi_1 = G.move_v1(v[0], G.pi_1[v[0]]+k)
                    __auxplot(G_aux, v[0])
            
            if (G.pi_1[v[0]] - k > 0) and np.sign(v[1])==-1:
            #if (m_minus is not None):               
                chg_pos = G.find_pos(v[0], G.pi_1[v[0]]-k)
                K_vj = G.K(v[0], chg_pos) 
                K_jv = G.K(chg_pos, v[0])
                v[1] = -k
                v[2] = v[2] + (K_vj - K_jv)
                if (v[2]  < trsh_crossing):
                    idx_min = i
                    trsh_crossing = v[2] 
                if verbose: 
                    G_aux.pi_1 = G.move_v1(v[0], G.pi_1[v[0]]-k)
                    __auxplot(G_aux, v[0], " -")
            
        else:           #Layer 2
            if (G.pi_2[v[0]] + k <= G.n_v2()) and np.sign(v[1])==1: 
            #if (m_plus is not None):
                chg_pos = G.find_pos(v[0], G.pi_2[v[0]]+k)
                K_vj = G.K(v[0], chg_pos) 
                K_jv = G.K(chg_pos, v[0])
                v[1] = k
                v[2] = v[2] + (K_jv - K_vj)
                if (v[2] < trsh_crossing):
                    idx_min = i
                    trsh_crossing = v[2]
                if verbose: 
                    G_aux.pi_2 = G.move_v1(v[0], G.pi_2[v]+k)
                    __auxplot(G_aux, v)
            
            if (G.pi_2[v[0]] - k > 0) and np.sign(v[1])==-1:
            #if (m_minus is not None):               
                chg_pos = G.find_pos(v[0], G.pi_2[v[0]]-k)
                K_vj = G.K(v[0], chg_pos) 
                K_jv = G.K(chg_pos, v[0])
                v[1] = -k
                v[2] = v[2] + (K_vj - K_jv)
                if (v[2]  < trsh_crossing):
                    idx_min = i
                    trsh_crossing = v[2] 
                if verbose: 
                    G_aux.pi_2 = G.move_v1(v[0], G.pi_2[v]-k)
                    __auxplot(G_aux, v[0], " -")
            
        if verbose:
            G_aux.pi_1 = G_pi_1_org
            G_aux.pi_2 = G_pi_2_org
    
    neighlist.append(neighlist.pop(idx_min))
    
def improvement_phase(G, V, k, verbose=0):
    
    global n_cross
    
    min_cross = n_cross #G.n_cross()

    G_pi_1_min = G.pi_1
    G_pi_2_min = G.pi_2
    
    i=0
    neigh_list = V
    neigh_set = list(_make_neighborhood(G, neigh_list))
    while neigh_list:
        
        neigh_tup = neigh_set[-1]
        tabu_dict[neigh_tup[0]] = i
        cross = neigh_tup[2]
        n_cross = cross
        
        v = neigh_tup[0]
        if v in G.v1():
            G.pi_1 = G.move_v1(v, G.pi_1[v] + neigh_tup[1])
        else:
            G.pi_2 = G.move_v2(v, G.pi_2[v] + neigh_tup[1])
        
        if cross < min_cross:
            min_cross = cross
            G_pi_1_min = G.pi_1
            G_pi_2_min = G.pi_2
        
        i = i + 1        
        neigh_list = [v for v in V if (i - tabu_dict[v]) > ternure ]
        neigh_set = list(_make_neighborhood(G, neigh_list))
        
        
    G.pi_1 = G_pi_1_min
    G.pi_2 = G_pi_2_min
    n_cross = min_cross
    return min_cross

def VNS(G, alpha=1.0, k_max=5, verbose=0):
    
    U_1 = G.v1().copy()
    U_2 = G.v2().copy()
    V = U_1 + U_2
    
    construction_phase(G, U_1, U_2, V, alpha)
    G.order_v1()
    G.order_v2()
    global n_cross
    n_cross = G.n_cross()
    min_c = n_cross
    
    it = 0
    while it < max_it:
        C = improvement_phase(G, V, len(V))
        if verbose: 
            print(it, ":", str(min_c), "::", C)#, ":::", G.n_cross())
        if C < min_c:
            it = 0
            min_c = C
        else: 
            it = it + 1