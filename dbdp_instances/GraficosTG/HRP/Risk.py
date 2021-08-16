# -*- coding: utf-8 -*-
"""
Created on Wed Apr 28 13:16:51 2021

@author: bferrari
"""

import numpy as np
import pandas as pd
import scipy.stats as st
import scipy.cluster.hierarchy as clst

from itertools import combinations

class Risk:
    
    def vol_correl(self, df_weight, df_vols, df_correl):
        "Calcula Vol utilizando peso e vols e a matrix de correlação dos ativos"
        return self.risk_contribution_corr(df_weight, df_vols, df_correl).sum(axis=1)
    
    
    def param_var(self, vol, days=1, conf_level=0.95):
        "Calcula o VaR parametrico utilizando a vol, o periodo e o nivel de confianca"
        
        return st.norm.ppf(conf_level) * np.sqrt(days/252) * vol
    
    def hist_var(self, df_weight, df_historico_retornos, conf_level=0.95):
        "Calcula o VaR parametrico utilizando os retornos historicos e o nivel de confianca"
        r_hist = df_historico_retornos @ df_weight
        return -r_hist.quantile(1-conf_level)
        
    def risk_contribution_corr(self, df_weight, df_vols, df_correl):
        
        "Calcula o RC utilizando peso, vols e a matriz de correlação dos ativos"
        
        cov = (np.diag(df_vols) @ df_correl) @ np.diag(df_vols)                       
        sigma_w = df_weight.values @ cov
        vol = np.sqrt(sigma_w @ df_weight.values.T) 
        rc = (sigma_w / vol) * df_weight.values
        rc = pd.DataFrame(rc.values.reshape(1,-1), index = [df_weight.name], columns = df_weight.index)
    
        return rc
    
    def vol_cov(self, df_weight, df_historico_retornos, freq=252):
        "Calcula Vol utilizando peso e historico de retornos dos ativos via Covariancia"
        
        cov = df_historico_retornos.loc[:df_weight.index[0], df_weight.columns][-freq:].cov() * 252
        sigma_w = df_weight @ cov
        vol = np.sqrt(sigma_w @ df_weight.T)
        rc = (sigma_w / vol.values) * df_weight
        
        return rc
    
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
    
    def getIVP(cov): 
        # Compute the inverse-variance portfolio 
        ivp= 1./np.diag(cov)
        ivp/=ivp.sum() # equivalent to np.trace on 1/cov ( this is not the inverse matrix !!)
        return ivp

    def getIVP(cov): 
        # Compute the inverse-variance portfolio 
        
        
        
        
        ivp= 
        ivp/=ivp.sum() # equivalent to np.trace on 1/cov ( this is not the inverse matrix !!)
        return ivp


    def getClusterVar(cov,cItems): 
        # Compute variance per cluster 
        try:
            cov_ = cov.loc[cItems,cItems] # matrix slice 
        except KeyError:
            cov_ = cov.iloc[cItems,cItems]
            
        w_ = getIVP(cov_).reshape(-1,1)
        cVar = ((w_.T@cov_) @ w_).values[0,0] 
        return cVar
    
    def getRecBipart(cov, sortIx):
        # Compute HRP alloc 
        w=pd.Series(1,index=sortIx)
        cItems=[sortIx] # initialize all items in one cluster 
        while len(cItems)>0:
            cItems=[i[j:k] for i in cItems for j,k in ((0,int(len(i)/2)), (int(len(i)/2),len(i))) if len(i)>1]    # bi-section         
            for i in range(0,len(cItems),2): # parse in pairs
                cItems0=cItems[i]       # cluster 1             
                cItems1=cItems[i+1]     # cluster 2 
                cVar0=getClusterVar(cov,cItems0)
                cVar1=getClusterVar(cov,cItems1)
                alpha= 1 - cVar0/(cVar0+cVar1)
                w[cItems0]*= alpha      # weight 1 
                w[cItems1]*= 1-alpha    # weight 2 

        return w
        
    def getHRP(df_historico_retornos, by_name=1):
        # A distance matrix based on correlation, where 0<=d[i,j]<=1      
        # This is a proper distance metric   
        dist = np.sqrt((1-df_historico_retornos.corr())/2)  
        link = clst.linkage(dist)
        
        sortIx=getQuasiDiag(link)[:7] 
        if by_name: 
            sortIx=dist.index[sortIx].tolist() # recover labels     
        
        try:
            hrp=pd.Series(getRecBipart(df_historico_retornos.cov(),sortIx).values, index=dist.index[sortIx])
        except:
            hrp=getRecBipart(df_historico_retornos.cov(),sortIx)
            
        return hrp

    def covHRP(cov, df_weight):
        X_Y = list(combinations(df_weight.columns,2))
        
        HRP_cor = pd.DataFrame(index=df_weight.columns, columns=df_weight.columns).replace(np.nan,1)
        
        HRP_cov = HRP_cor*np.diag(cov.loc[df_weight.columns,df_weight.columns])
        for xy in X_Y:
            cVar = getClusterVar(cov, xy)
            cxy = (cVar - ((df_weight[xy[0]][0]**2)*cov.loc[xy[0], xy[0]] + 
                           (df_weight[xy[1]][0]**2)*cov.loc[xy[1], xy[1]] )) / 2*df_weight[xy[0]][0]*df_weight[xy[1]][0]
            
            cxy = cxy/(df_weight[xy[0]][0]*df_weight[xy[1]][0])
            HRP_cov.loc[xy[0], xy[1]] = cxy
            HRP_cov.loc[xy[1], xy[0]] = cxy
            
            corr_xy = cxy / np.sqrt(cov.loc[xy[0], xy[0]]*cov.loc[xy[1], xy[1]])
            HRP_cor.loc[xy[0], xy[1]] = corr_xy
            HRP_cor.loc[xy[1], xy[0]] = corr_xy
    
    df_weight = pd.DataFrame()
    aux_w = reg_dict_sadj['21.470.947/0001-36']['beta'][-1:]
    for col in aux_w:
        if aux_w[col][0] < 0:
            df_weight[col + '_inv'] = -aux_w[col].values
        else:
            df_weight[col] = aux_w[col].values
        
    