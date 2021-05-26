# -*- coding: utf-8 -*-
"""
Created on Wed May 26 01:45:16 2021

@author: Bruno Ferrari
"""
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
            
            if (G.pi_1[i] < G.pi_1[j]) and (G.pi_2[k] > G.pi_2[l]) and (G.pi_1[i] * G.pi_1[j] * G.pi_2[k] * G.pi_2[l]):
                c = c + 1
            elif (G.pi_1[i] > G.pi_1[j]) and (G.pi_2[k] < G.pi_2[l]) and (G.pi_1[i] * G.pi_1[j] * G.pi_2[k] * G.pi_2[l]):
                c = c + 1
    return c



from collections.abc import Sequence
from collections.abc import Mapping

class BGraph:

    """ aux """
    def __init__(self):
        self.set_v1 = []
        self.set_v2 = []
        self.set_edges = []
        self.pi_1 = {}
        self.pi_2 = {}
        self._adj = {}
        self._nodes = {}
        
    @property
    def adj(self):
        return AdjacencyView(self._adj)
    
    def __getitem__(self, n):
        return self.adj[n]
        
        
    def edges(self, edgelist = None):
        if edgelist != None:
            self.set_edges = []
            self.set_edges = edgelist
            
            for e in edgelist:
                u, v = e
                
                if u not in self._nodes:
                     self._adj[u] = {}
                     self._nodes[u] = {}
                if v not in self._nodes:
                     self._adj[v] = {}
                     self._nodes[v] = {}
                
                datadict = self._adj[u].get(v, {})
                self._adj[u][v] = datadict
                self._adj[v][u] = datadict
        else:
            return self.set_edges
        
class AtlasView(Mapping):
    
    __slots__ = ("_atlas")
    
    def __init__(self, d):
        self._atlas = d
    
    def __len__(self):
        return len(self._atlas)    
    
    def __iter__(self):
        return iter(self._atlas)
    
    def __getitem__(self, key):
        return self._atlas[key]


class AdjacencyView(AtlasView):

    __slots__ = ()  # Still uses AtlasView slots names _atlas

    def __getitem__(self, name):
        return AtlasView(self._atlas[name])

    def copy(self):
        return {n: self[n].copy() for n in self._atlas}


class MyClass(Sequence):
    def __init__(self, L):
        self.L = L
        super().__init__()
    def __getitem__(self, i):
        return self.L[i]
    def __len__(self):
        return len(self.L)

# Let's test it:
myobject = MyClass([1,2,3])
try:
    for idx,_ in enumerate(myobject):
        print(myobject[idx])
except Exception:
    print("Gah! No good!")
    raise
