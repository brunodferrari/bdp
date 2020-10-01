# -*- coding: utf-8 -*-
"""
Created on Sun Sep 27 19:59:50 2020

@author: Bruno Ferrari
"""

import numpy as np
import matplotlib as plt
import networkx as nx


G = nx.bipartite.gnmk_random_graph(3,5,10, seed=111)

top = nx.bipartite.sets(G)[0]
pos = nx.bipartite_layout(G, top)

nx.draw(G, pos)
nx.draw_networkx_labels(G, pos)

def pi(setlist, i):
    try:
        return np.int(np.where(np.array(setlist) == i )[0])
    except TypeError:
        return -1    

def plotBGraph(G):
    
    B = nx.Graph()
    B.add_nodes_from(G.v1(), bipartite=1)
    B.add_nodes_from(G.v2(), bipartite=2)
    B.add_edges_from(G.edges())
    
    pos = bdp_lyt(G)    
    nx.draw(B, pos)
    nx.draw_networkx_labels(B, pos)

def bdp_lyt(G):
    
    import numpy as np
    #G, center = _process_params(G, center=center, dim=2)
    #if len(G) == 0:
     #   return {}
     
    #center = np.zeros(2)

    top = G.v1()[::-1]
    bottom = G.v2()[::-1]

    height = 1
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
        
    def adj_v1(self, G):
        self.adj_v1 = G
        
B = BGraph()
B.v1([5,3,2,4])
B.v2(['a','b','c'])
B.edges([(5, "a"), (5, "b"), (2, "b"), (2, "c"), (3, "c"), (4, "a")])
plotBGraph(B)



C = nx.Graph()

# Add nodes with the node attribute "bipartite"

C.add_nodes_from(B.v1(), bipartite=1)
C.add_nodes_from(B.v2(), bipartite=2)

# Add edges only between nodes of opposite node sets

C.add_edges_from([(5, "a"), (5, "b"), (2, "b"), (2, "c"), (3, "c"), (4, "a")])
