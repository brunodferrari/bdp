# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt
import networkx as nx
import pandas as pd
import copy

from numba import njit
from numba.typed import Dict, List
from numba.core import types

from concurrent.futures import ThreadPoolExecutor

np.seterr(over='ignore')
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
            
            if (G.pi_1[i] < G.pi_1[j]) and (G.pi_2[k] > G.pi_2[l]) and (G.pi_1[i] * G.pi_1[j] * G.pi_2[k] * G.pi_2[l]):
                c = c + 1
            elif (G.pi_1[i] > G.pi_1[j]) and (G.pi_2[k] < G.pi_2[l]) and (G.pi_1[i] * G.pi_1[j] * G.pi_2[k] * G.pi_2[l]):
                c = c + 1
    return c

def _cross_n(G, e1, e2):
    i = e1[0]
    k = e1[1]
    j = e2[0]
    l = e2[1]
    if (G.pi_1[i] < G.pi_1[j]) and (G.pi_2[k] > G.pi_2[l]) and (G.pi_1[i] * G.pi_1[j] * G.pi_2[k] * G.pi_2[l]):
        return 1
    elif (G.pi_1[i] > G.pi_1[j]) and (G.pi_2[k] < G.pi_2[l]) and (G.pi_1[i] * G.pi_1[j] * G.pi_2[k] * G.pi_2[l]):
        return 1
    
    return 0

def _cross_w(G, edgeslist):
    e1 = edgeslist[0]
    #with ThreadPoolExecutor(6) as ex:
    output = list(map(lambda x: _cross_n(G, e1, x), edgeslist[1:]))
    
    #output=0
    #for e2 in edgeslist[1:]:
    #    output = output + _cross_n(G, e1, e2)
    
    return sum(output)
    
def crossing2(G):
    edgeslist = G.edges()
    c = 0
    with ThreadPoolExecutor(6) as ex:
        output = list(ex.map(lambda x: _cross_w(G, x), [edgeslist[i:] for i in range(len(edgeslist))]))
        
    return c + int(sum(output))

@njit
def _numba_cross(pi_1, pi_2, edgeslist):
    c = 0
    for s, e1 in enumerate(edgeslist):
        i = e1[0]
        k = e1[1]
        for e2 in edgeslist[s+1:]:
            j = e2[0]
            l = e2[1]
        
            if (pi_1[i] < pi_1[j]) and (pi_2[k] > pi_2[l]) and (pi_1[i] * pi_1[j] * pi_2[k] * pi_2[l]):
                c = c + 1
            elif (pi_1[i] > pi_1[j]) and (pi_2[k] < pi_2[l]) and (pi_1[i] * pi_1[j] * pi_2[k] * pi_2[l]):
                c = c + 1
        
    return c

def crossing3(G):
    edgeslist = G.edges()
    c = 0
    
    pi_1 = Dict.empty(
        key_type=types.int64,
        value_type=types.int64,
                    )
    pi_2 = Dict.empty(
        key_type=types.int64,
        value_type=types.int64,
                    )
    pi_1.update(G.pi_1)
    pi_2.update(G.pi_2)
    
    output = _numba_cross(pi_1, pi_2, List(edgeslist))
    
    return c + output  
    
    
def bdp_lyt(G, size = 4/3, height = 100): ## Formata lyt adequado para o plot do grafo
    
    import numpy as np
    #G, center = _process_params(G, center=center, dim=2)
    #if len(G) == 0:
     #   return {}
     
    #center = np.zeros(2)

    top = G.v1()[::-1]
    bottom = G.v2()[::-1]

    height = 100
    width = size * height
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
def plotBGraph(G, size = 4/3, height=100):
    
    B = nx.Graph()
    B.add_nodes_from(G.v1(), bipartite=1)
    B.add_nodes_from(G.v2(), bipartite=2)
    B.add_edges_from(G.edges())
    
    pos = bdp_lyt(G, size, height)    
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
    b = 0
    if v_layer == None:
        return
    elif v_layer == 1: 
        pi_k = G.v2()
        K = [x for x in pi_k if (((v, x) in G.edges()) or ((x, v) in G.edges())) and G.pi_2[x] > 0] #encontra os viznho do vertice v na 2a camada
        if len(K) > 0:
            b = G.perm_v2(K).mean()    
    elif v_layer == 2:
        pi_k = G.v1()
        #K = [x for x in pi_k if (x, v) in G.edges()] #encontra os viznho do vertice v na 1a camada
        K = [x for x in pi_k if (((v, x) in G.edges()) or ((x, v) in G.edges())) and G.pi_1[x] > 0]
        if len(K) > 0:
            b = G.perm_v1(K).mean()

    return b

@njit
def _deg(nodelist, subgraph, edges):
    deg = Dict.empty(
        key_type=types.int64,
        value_type=types.int64,
                    )
    for v in nodelist:
            K = [x for x in subgraph if ((v, x) in edges) or ((x, v) in edges)]
            deg[v] = len(K)
    
    return deg

class BGraph:

    """ aux """
    def __init__(self):
        self._set_v1 = []
        self._set_v2 = []
        self._set_edges = []
        self.pi_1 = {}
        self.pi_2 = {}
        self._adj = {}
        self._nodes = {}
        
    @property
    def adj(self):
        return self._adj
    
    def __getitem__(self, n):
        return self.adj[n]
    
    def edges(self, edgelist = None):
        if edgelist != None:
            self._set_edges = []
            self._set_edges = edgelist
            
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
            return self._set_edges
    
    def v1 (self, setlist = None):
        if setlist != None:
            self._set_v1 = []
            self._set_v1 = setlist
            self.pi_1 = dict(zip(self._set_v1,self.perm_v1()))
            
            for u in setlist:
                if u not in self._nodes: 
                    self._nodes[u] = {}
                    self._adj[u] = {}
        else:
            return self._set_v1
        
    def v2 (self, setlist = None):
        if setlist != None:
            self._set_v2 = []
            self._set_v2 = setlist
            self.pi_2 = dict(zip(self._set_v2, self.perm_v2()))
            
            for u in setlist:
                if u not in self._nodes: 
                    self._nodes[u] = {}
                    self._adj[u] = {}
        else:
            return self._set_v2
    
    def perm_v1(self, pos = None):
        if pos != None:
            return np.vectorize(lambda i: self.pi_1[i])(pos)
        else:
            return np.vectorize(lambda i: pi(self._set_v1, i))(self._set_v1) + 1

    def perm_v2(self, pos = None):
        if pos != None:
            return np.vectorize(lambda i: self.pi_2[i])(pos)
        else:
            return np.vectorize(lambda i: pi(self._set_v2, i))(self._set_v2) + 1
    
    def order_v1(self):
        #aux = [(pos, v) for (pos, v) in zip(self.pi_1, self. set_v1)]
        aux = list( self.pi_1.items() )
        aux.sort(key=lambda tup: tup[1])
        self.v1(list(np.int0(np.array(aux)[:,0])))
        return
        
    def order_v2(self):
        #aux = [(pos, v) for (pos, v) in zip(self.pi_2,self. set_v2)]
        aux = list( self.pi_2.items() )
        aux.sort(key=lambda tup: tup[1])
        self.v2(list(np.int0(np.array(aux)[:,0])))
        return
    
    def add_v1(self, i, pos):
        if pos != -1: 
            self._set_v1 = self.v1()[:pos] + [i] + self.v1()[pos:]
        else:
            self._set_v1 = self.v1()[:] + [i]
        
    def add_v2(self, i, pos):
        if pos != -1: 
            self._set_v2 = self.v2()[:pos] + [i] + self.v2()[pos:]
        else:
            self._set_v2 = self.v2()[:] + [i]
    
    def n_v1(self):
        #self.n_v1 = len(self.v1)
        return len(self._set_v1)
        
    def n_v2(self):
        #self.n_v2 = len(self.v2)
        return len(self._set_v2)
            
    def n_edge(self):
        return len(self._set_edges)
        
    def n_v(self):
        return len(self._nodes)
    
    def density(self):
        return len(self._set_edges) / (len(self._set_v1)*len(self._set_v2))
        
    def n_cross(self):
        return crossing3(self)
    
    def bc(self, v, k):
        return bary(self, v, k)
    
    #@jit(nopython=True)
    def degree(self, nodelist = None, subgraph = None):
        #deg = {}
        
        if nodelist is None:
            nodelist = (self.v1() + self.v2())
        
        if subgraph is None:
            subgraph = (self.v1() + self.v2())
        
        # for v in self.v1():
        #     K = [x for x in self.v2() if (v, x) in self.edges()]
        #     deg[v] = len(K)
        # for v in self.v2():
        #     K = [x for x in self.v1() if (x, v) in self.edges()]
        #     deg[v] = len(K)
        
        d = _deg(List(nodelist), List(subgraph), List(self.edges()))
            
    
        return d
    
    def move_v1(self, v, pos, inplace=False):
        
        aux = self.pi_1.copy()
        pos_v = aux.pop(v)
        aux = np.array(list(aux.items()))
        aux[ aux[:,1] >= pos_v, 1] = aux[ aux[:,1] >= pos_v, 1] - 1
        aux[ aux[:,1] >= pos, 1] = aux[ aux[:,1] >= pos, 1] + 1
        aux = dict(aux)
        aux[v] = pos   
                
        if not inplace:
            return aux
        else:
            self.pi_1 = aux
    
    def move_v2(self, v, pos, inplace=False):
        
        aux = self.pi_2.copy()
        pos_v = aux.pop(v)
        aux = np.array(list(aux.items()))
        aux[ aux[:,1] >= pos_v, 1] = aux[ aux[:,1] >= pos_v, 1] - 1
        aux[ aux[:,1] >= pos, 1] = aux[ aux[:,1] >= pos, 1] + 1
        aux = dict(aux)
        aux[v] = pos   
                
        if not inplace:
            return aux
        else:
            self.pi_2 = aux
    
    def plot(self, size=4/3, height=100, order=0):
        if order:
            self.order_v1()
            self.order_v2()
        plotBGraph(self, size=size, height=height)
        
    
    def K(self, u, v):
        i = u
        j = v
    
        #G.move_v1(i, j, True)
        #G.move_v1(j, i, True)
        
        c = 0
         
        if u in self._set_v1:
            pi = self.pi_2
        elif u in self._set_v2:
            pi = self.pi_1
        #nodes_between = [v for v in G.pi_1 if G.pi_1[v] >= G.pi_1[i] and G.pi_1[v] <= G.pi_1[j]]
            
        #while nodes_between:
        #    i = nodes_between.pop()
        #    for j in nodes_between:
        #print(pi)            
        for k in self._adj[i]:
            for l in self._adj[j]:
                if (pi[k] > pi[l]):
                    c = c + 1
                            
        return c
    
    def find_pos(self, u, pos):
        
        if u in self._set_v1:
            pi = self.pi_1
        elif u in self._set_v2:
            pi = self.pi_2
        
        return [u for u in pi if pi[u] == pos][0]
    
    def copy(self):
        return copy.deepcopy(self)