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



B = nx.Graph()

# Add nodes with the node attribute "bipartite"

B.add_nodes_from([5, 3, 2, 4], bipartite=1)
B.add_nodes_from(["a", "b", "c"], bipartite=2)

# Add edges only between nodes of opposite node sets

B.add_edges_from([(5, "a"), (5, "b"), (2, "b"), (2, "c"), (3, "c"), (4, "a")])


def pi(setlist, i):
    try:
        return np.int(np.where(np.array(setlist) == i )[0])
    except TypeError:
        return -1    

def plotBGraph(G):
    
    top = list(nx.bipartite.sets(G)[0])[::-1]
    #top = nx.bipartite.sets(G)[0]
    pos = nx.bipartite_layout(G, top)
    
    nx.draw(G, pos)
    nx.draw_networkx_labels(G, pos)


class Graph:
    
    #def __init__(self):
     #   self.set_v1 = []

    def set_v1(self):
        self.set_v1 = []
    
    def set_v2(self):
        self.set_v2 = []
    
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

C = Graph()
C.v1([10,3,4,5,61, 20])