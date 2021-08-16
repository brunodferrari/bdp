# -*- coding: utf-8 -*-
"""
Created on Thu Aug  5 13:33:45 2021

@author: bferrari
"""

import xlwings as xw
import pandas as pd
import numpy as np

import scipy.cluster.hierarchy as clst
import matplotlib.pyplot as plt


def getHRP(df_historico_retornos, link=None, sortIx=None, by_name=1):
        # A distance matrix based on correlation, where 0<=d[i,j]<=1      
        # This is a proper distance metric   
        dist = np.sqrt((1-df_historico_retornos.corr())/2)  
        if link is None:
            link = clst.linkage(dist)
        
        if sortIx is None:
            sortIx=getQuasiDiag(link)
            if by_name: 
                sortIx=dist.index[sortIx].tolist() # recover labels     
        
        try:
            hrp=pd.Series(getRecBipart(df_historico_retornos.cov(),sortIx).values, index=dist.index[sortIx])
        except:
            hrp=getRecBipart(df_historico_retornos.cov(),sortIx)
            
        return hrp



def getRecBipart(cov, sortIx):
    # Compute HRP alloc 
    print(sortIx)
    w = pd.Series(1,index=sortIx)
    cItems = [sortIx] # initialize all items in one cluster 
    while len(cItems) > 0:
        cItems=[i[j:k] for i in cItems for j,k in ((0,int(len(i)/2)), (int(len(i)/2),len(i))) if len(i)>1]    # bi-section         
        for i in range(0,len(cItems),2): # parse in pairs
            cItems0 = cItems[i]       # cluster 1             
            cItems1 = cItems[i+1]     # cluster 2 
            cVar0 = getClusterVar(cov,cItems0)
            cVar1 = getClusterVar(cov,cItems1)
            alpha = 1 - cVar0/(cVar0+cVar1)
            w[cItems0]*= alpha      # weight 1 
            w[cItems1]*= 1-alpha    # weight 2 

    return w

def getClusterVar(cov,cItems): 
    # Compute variance per cluster 
    try:
        cov_ = cov.loc[cItems,cItems] # matrix slice 
    except KeyError:
        cov_ = cov.iloc[cItems,cItems]
        
    w_ = getIVP(cov_).reshape(-1,1)
    cVar = ((w_.T@cov_) @ w_).values[0,0] 
    return cVar

def getIVP(cov,**kwargs): 
    # Compute the inverse-variance portfolio 
    ivp= 1./np.diag(cov)
    ivp/=ivp.sum() # equivalent to np.trace on 1/cov ( this is not the inverse matrix !!)
    return ivp

def plotDendrogram(link, labels=None, leaf_rotation=45):
    plt.clf()
    clst.dendrogram(link, labels=labels, leaf_rotation=leaf_rotation)
    plt.title('Dendrograma')
    fig = plt.gcf()
    plt.close()
    return fig
    
def getQuasiDiag(link): 
    # Sort clustered items by distance 
    link=link.astype(int)     
    sortIx=pd.Series([link[-1,0],link[-1,1]])     
    numItems=link[-1,3] # number of original items 
    while sortIx.max()>=numItems:         
        sortIx.index=range(0,sortIx.shape[0]*2,2) # make space 
        df0=sortIx[sortIx>=numItems] # find clusters 
        i=df0.index;j=df0.values-numItems         
        sortIx[i]=link[j,0] # item 1         
        df0=pd.Series(link[j,1],index=i+1)         
        sortIx=sortIx.append(df0) # item 2 
        sortIx=sortIx.sort_index() # re-sort         
        sortIx.index=range(sortIx.shape[0]) # re-index 

    return sortIx.tolist()

def plotCorrMatrix(corr, labels=None, cmap='RdBu'):
    
    if labels is None: labels=[]
    plt.pcolormesh(corr, cmap=cmap)
    plt.colorbar()
    plt.yticks(np.arange(.5,corr.shape[0]+.5),labels)
    plt.xticks(np.arange(.5,corr.shape[0]+.5),labels, rotation=90)
    plt.title('Matriz de Correlação')
    fig = plt.gcf()
    plt.close()
    return fig

def getClusters(df_historico_retornos):
    dist = np.sqrt(1-df_historico_retornos.corr()/2)
    link = clst.linkage(dist)
    
    return link

import random
def generateData(nObs,size0,size1,sigma1): 
    # Time series of correlated variables    
    #1) generating some uncorrelated data 
    np.random.seed(seed=12345); random.seed(12345)
    x = np.random.normal(0,1,size=(nObs,size0)) # each row is a variable 
    #2) creating correlation between the variables 
    cols = [random.randint(0,size0-1) for i in range(size1)]    
    y=x[:,cols]+np.random.normal(0,sigma1,size=(nObs,len(cols)))    
    x=np.append(x,y,axis=1)     
    x=pd.DataFrame(x,columns=range(1,x.shape[1]+1))     
    return x,cols

#def getHRP(cov,corr): 
#    # Construct a hierarchical portfolio 
#    corr,cov=pd.DataFrame(corr),pd.DataFrame(cov)
#    dist=np.sqrt(1-corr/2)
#    link=clst.linkage(dist,'single')
#    sortIx=getQuasiDiag(link)
#    sortIx=corr.index[sortIx].tolist() # recover labels 
#   hrp=getRecBipart(cov,sortIx) 
#    return hrp.sort_index() 

def hrpMC(numIters=1e3,nObs=520,size0=5,size1=5,mu0=0,sigma0=1e-2, sigma1F=.25,sLength=260,rebal=22): 
    # Monte Carlo experiment on HRP 
    methods = [getIVP,getHRP]
    stats = {i.__name__:pd.Series() for i in methods}
    numIter = 0
    pointers = range(sLength,nObs,rebal)
    while numIter<numIters:
        print(numIter)
        #1) Prepare data for one experiment 
        x,cols = generateData(nObs,size0,size1,sigma1F)
        r = {i.__name__:pd.Series() for i in methods} 
        #2) Compute portfolios in-sample 
        for pointer in pointers:
            x_ = x[pointer-sLength:pointer]
            cov_, corr_ = np.cov(x_,rowvar=0), np.corrcoef(x_,rowvar=0)
            x_= x[pointer:pointer+rebal]
            for func in methods:                 
                w_ = func(cov=cov_, corr=corr_) # callback 
                r_= pd.Series(np.dot(x_,w_))
                r[func.__name__]=r[func.__name__].append(r_)
            #4) Evaluate and store results 
            for func in methods:             
                r_=r[func.__name__].reset_index(drop=True)
                p_=(1+r_).cumprod()
                stats[func.__name__].loc[numIter]=p_.iloc[-1]-1 # terminal return 
            numIter+=1
        #5) Report results 
        stats=pd.DataFrame.from_dict(stats,orient='columns')
        stats.to_csv('stats.csv')
        df0, df1= stats.std(),stats.var()
        print(pd.concat([df0,df1,df1/df1['getHRP']-1],axis=1))
        return 