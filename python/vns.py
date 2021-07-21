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

#from python import bgraph

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
'''
def _make_neighborhood(G, neighlist, k=1, verbose=0):
    
    if verbose:
        G_aux = G.copy()
        G_pi_1_org = G.pi_1
        G_pi_2_org = G.pi_2
    
    trsh_crossing = np.inf
    idx_min = np.nan
    n=len(neighlist)
    v=neighlist.copy()
    for i in range(n):
        if v[i][0] in G.v1(): #Layer 1
            if (G.pi_1[v[i][0]] + k <= G.n_v1()) and np.sign(v[i][1])==1: 
            #if (m_plus is not None):
                chg_pos = G.find_pos(v[i][0], G.pi_1[v[i][0]]+k)
                K_vj = G.K(v[i][0], chg_pos) 
                K_jv = G.K(chg_pos, v[i][0])
                v[i][1] = k
                v[i][2] = v[i][2] + (K_jv - K_vj)
                if (v[i][2] < trsh_crossing):
                    idx_min = i
                    trsh_crossing = v[i][2]
                if verbose: 
                    G_aux.pi_1 = G.move_v1(v[i][0], G.pi_1[v[i][0]]+k)
                    __auxplot(G_aux, v[i][0])
            elif (G.pi_1[v[i][0]] + k > G.n_v1()) and np.sign(v[i][1])==1: 
                print(":1", neighlist.pop(i - (len(v)-len(neighlist)) ))
                
            if (G.pi_1[v[i][0]] - k > 0) and np.sign(v[i][1])==-1:
            #if (m_minus is not None):               
                chg_pos = G.find_pos(v[i][0], G.pi_1[v[i][0]]-k)
                K_vj = G.K(v[i][0], chg_pos) 
                K_jv = G.K(chg_pos, v[i][0])
                v[i][1] = -k
                v[i][2] = v[i][2] + (K_vj - K_jv)
                if (v[i][2]  < trsh_crossing):
                    idx_min = i
                    trsh_crossing = v[i][2] 
                if verbose: 
                    G_aux.pi_1 = G.move_v1(v[i][0], G.pi_1[v[i][0]]-k)
                    __auxplot(G_aux, v[i][0], " -")
            elif (G.pi_1[v[i][0]] - k <= 0) and np.sign(v[i][1])==-1: 
                print(":2", neighlist.pop(i - (len(v)-len(neighlist)) ))
        else:           #Layer 2
            if (G.pi_2[v[i][0]] + k <= G.n_v2()) and np.sign(v[i][1])==1: 
            #if (m_plus is not None):
                chg_pos = G.find_pos(v[i][0], G.pi_2[v[i][0]]+k)
                K_vj = G.K(v[i][0], chg_pos) 
                K_jv = G.K(chg_pos, v[i][0])
                v[i][1] = k
                v[i][2] = v[i][2] + (K_jv - K_vj)
                if (v[i][2] < trsh_crossing):
                    idx_min = i
                    trsh_crossing = v[i][2]
                if verbose: 
                    G_aux.pi_2 = G.move_v1(v[i][0], G.pi_2[v]+k)
                    __auxplot(G_aux, v)
            
            elif (G.pi_2[v[i][0]] + k > G.n_v2()) and np.sign(v[i][1])==1: 
                print(":3", neighlist.pop(i - (len(v)-len(neighlist)) ))
            
            if (G.pi_2[v[i][0]] - k > 0) and np.sign(v[i][1])==-1:
            #if (m_minus is not None):               
                chg_pos = G.find_pos(v[i][0], G.pi_2[v[i][0]]-k)
                K_vj = G.K(v[i][0], chg_pos) 
                K_jv = G.K(chg_pos, v[i][0])
                v[i][1] = -k
                v[i][2] = v[i][2] + (K_vj - K_jv)
                if (v[i][2]  < trsh_crossing):
                    idx_min = i
                    trsh_crossing = v[i][2] 
                if verbose: 
                    G_aux.pi_2 = G.move_v1(v[i][0], G.pi_2[v]-k)
                    __auxplot(G_aux, v[i][0], " -")
            
            elif (G.pi_2[v[i][0]] - k <= 0) and np.sign(v[i][1])==-1: 
                print(":4", neighlist.pop(i - (len(v)-len(neighlist)) ))
            
        if verbose:
            G_aux.pi_1 = G_pi_1_org
            G_aux.pi_2 = G_pi_2_org
    return idx_min
    #neighlist.append(neighlist.pop(idx_min - ((len(v)-len(neighlist))) ))
''' 



def _eval_neighborhood(G, neighlist, k=1, verbose=0):

    if verbose:
        G_aux = G.copy()
        G_pi_1_org = G.pi_1
        G_pi_2_org = G.pi_2
    
    trsh_crossing = np.inf
    idx_min = np.nan
    
    for i, v in enumerate(neighlist):
        if v[0] in G.v1(): #Layer 1
            pi = G.pi_1
        else:              #Layer 2 
            pi = G.pi_2
         
        if (pi[v[0]] + k <= len(pi)) and np.sign(v[1])==1: 
        #if (m_plus is not None):
            chg_pos = G.find_pos(v[0], pi[v[0]]+k)
            K_vj = G.K(v[0], chg_pos) 
            K_jv = G.K(chg_pos, v[0])
            if ((K_jv - K_vj) <= 0):
                v[1] = k
                v[2] = v[2] + (K_jv - K_vj)
                if (K_jv - K_vj) < trsh_crossing: 
                    idx_min = i
                    trsh_crossing = v[2]
        
        if (pi[v[0]] - k > 0) and np.sign(v[1])==-1:
            chg_pos = G.find_pos(v[0], pi[v[0]]-k)
            K_vj = G.K(v[0], chg_pos) 
            K_jv = G.K(chg_pos, v[0])
            if ( ((K_vj - K_jv)) <= 0):
                v[1] = -k
                v[2] = v[2] + (K_vj - K_jv)
                if (K_vj - K_jv) < trsh_crossing: 
                    idx_min = i
                    trsh_crossing = v[2] 
    
    r = neighlist.copy()        
    try: 
        r.append(r.pop(idx_min))
        return [x for x in r if abs(x[1]) == k]
    except: 
        return [x for x in r if abs(x[1]) == k]
             
 
def improvement_phase(G, V, k, k_max, verbose=0):
    
    global n_cross
    
    min_cross = n_cross #G.n_cross()

    G_pi_1_min = G.pi_1
    G_pi_2_min = G.pi_2
    
    if k==1:
        neigh_list = [[v,i,0] for v in G._nodes for i in [np.inf,-np.inf]]
    
    neigh_list = _eval_neighborhood(G, neigh_list, k=1)
    while neigh_list and k<k_max:
        
        neigh_tup = neigh_list[-1]
        
        if n_cross > n_cross + neigh_tup[2]:
            v = neigh_tup[0]
            if v in G.v1():
                G.pi_1 = G.move_v1(v, G.pi_1[v] + neigh_tup[1])
            else:
                G.pi_2 = G.move_v2(v, G.pi_2[v] + neigh_tup[1])
            G_pi_1_min = G.pi_1
            G_pi_2_min = G.pi_2
            min_cross = n_cross + neigh_tup[2]
            n_cross = min_cross
            break
        
        k+=1
        neigh_list = _eval_neighborhood(G, neigh_list, k)
    if not(neigh_list) and k<k_max: 
        k = k_max
    return min_cross, k


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
    
    it = 1
    while it < k_max:
        C, k_aux = improvement_phase(G, V, 1, k_max, verbose)
        if verbose: 
            print(k_aux, ":", str(min_c), "::", C, ":::", G.n_cross())
        if C < min_c:
            it = 1
            min_c = C
        else: 
            it = k_aux