# -*- coding: utf-8 -*-
"""
Created on Thu Mar 25 22:57:24 2021

@author: Bruno Ferrari
"""
import numpy as np
import matplotlib.pyplot as plt
import networkx as nx

def pi(setlist, i):
        try:
            return np.int(np.where(np.array(setlist) == i )[0])
        except FutureWarning:
            print(setlist, i)
        except TypeError:
            return -1    
        
def crossing(G):
    aux = G.edges().copy()
    c = 0
    while len(aux) > 0:
        e1 = aux.pop(0)
        i = e1[0]
        k = e1[1]
        for e2 in aux:
            j = e2[0]
            l = e2[1]
            
            if pi(G.v1(), i) < pi(G.v1(), j) and pi(G.v2(), k) > pi(G.v2(), l):
                c+=1
            elif pi(G.v1(), i) > pi(G.v1(), j) and pi(G.v2(), k) < pi(G.v2(), l):
                c+=1  
    return c


def bdp_lyt(G): ## Formata lyt adequado para o plot do grafo
    
    import numpy as np
    #G, center = _process_params(G, center=center, dim=2)
    #if len(G) == 0:
     #   return {}
     
    #center = np.zeros(2)

    top = G.v1()[::-1]
    bottom = G.v2()[::-1]

    height = 100
    width = (4/3) * height
    offset = (width/2, height/2)

    nodes = top + bottom

    left_xs = np.repeat(0, len(top))
    right_xs = np.repeat(width, len(bottom))
    left_ys = np.linspace(0, height, len(top))
    right_ys = np.linspace(0, height, len(bottom))
    
    top_pos = np.column_stack([left_xs, left_ys]) - offset
    bottom_pos = np.column_stack([right_xs, right_ys]) - offset
    
    pos = np.concatenate([top_pos, bottom_pos])
    #pos = rescale_layout(pos, scale=scale) + center
    pos = dict(zip(nodes, pos))
    return pos


#plot utilizando o lyt adequado
def plotBGraph(G):
    
    B = nx.Graph()
    B.add_nodes_from(G.v1(), bipartite=1)
    B.add_nodes_from(G.v2(), bipartite=2)
    B.add_edges_from(G.edges())
    
    pos = bdp_lyt(G)    
    nx.draw(B, pos)
    nx.draw_networkx_labels(B, pos)
    #plt.savefig("test.pdf")    
    #plt.show()

def bary_sort(barylist, nodelist):
    
    aux = [(pos, v) for (pos, v) in zip(barylist, nodelist)]
    aux.sort(key=lambda tup: tup[0])
    
    return list(np.int0(np.array(aux)[:,1]))

#encontra baricentro do vertice
def bary(G, v, v_layer = None):
    
    if v_layer == None:
        return
    elif v_layer == 1: 
        pi_k = G.v2()
        K = [x for x in pi_k if (v, x) in G.edges()] #encontra os viznho do vertice v na 2a camada
        return G.perm_v2(K).mean()    
    elif v_layer == 2:
        pi_k = G.v1()
        K = [x for x in pi_k if (x, v) in G.edges()] #encontra os viznho do vertice v na 1a camada
        return G.perm_v1(K).mean()

class BGraph:
    
    
    def __init__(self):
        self.set_v1 = []
        self.set_v2 = []
        self.set_edges = []
            
    def edges(self, edgelist = None):
        if edgelist != None:
            self.set_edges = []
            self.set_edges = edgelist
        return self.set_edges
    
    def v1 (self, setlist = None):
        if setlist != None:
            self.set_v1 = []
            self.set_v1 = setlist
        return self.set_v1
        
    def v2 (self, setlist = None):
        if setlist != None:
            self.set_v2 = []
            self.set_v2 = setlist
        return self.set_v2
    
    def perm_v1(self, pos = None):
        if pos != None:
            return np.vectorize(lambda i: pi(self.v1(), i))(pos) + 1
        else:
            return np.vectorize(lambda i: pi(self.v1(), i))(self.v1()) + 1
    
    def perm_v2(self, pos = None):
        if pos != None:
            return np.vectorize(lambda i: pi(self.v2(), i))(pos) + 1
        else:
            return np.vectorize(lambda i: pi(self.v2(), i))(self.v2()) + 1
    
    def add_v1(self, i, pos):
        if pos != -1: 
            self.set_v1 = self.v1()[:pos] + [i] + self.v1()[pos:]
        else:
            self.set_v1 = self.v1()[:] + [i]
        
    def add_v2(self, i, pos):
        if pos != -1: 
            self.set_v2 = self.v2()[:pos] + [i] + self.v2()[pos:]
        else:
            self.set_v2 = self.v2()[:] + [i]
    
    def n_v1(self):
        #self.n_v1 = len(self.v1)
        return len(self.set_v1)
        
    def n_v2(self):
        #self.n_v2 = len(self.v2)
        return len(self.set_v2)
            
    def n_edge(self, n):
        self.n_edge = n
        
    def density(self):
        return len(self.set_edges) / (len(self.set_v1)*len(self.set_v2))
        
    def n_cross(self):
        return crossing(self)