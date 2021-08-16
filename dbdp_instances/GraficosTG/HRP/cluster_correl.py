# -*- coding: utf-8 -*-
"""
Created on Wed Jul 28 16:49:51 2021

@author: bferrari
"""

from datetime import datetime
import sys
sys.path.insert(1, 'P:/quant/projetos/Projeto 45/')

import pandas as pd
import numpy as np

from Positioning import sel_fundos, ajusta_cotas
from Positioning import (
    get_tipo_retorno, 
    get_fundos_data, 
    get_historico_cotas, 
    get_risk_data
)

from sklearn.decomposition import PCA

import matplotlib.pyplot as plt

import scipy.cluster.hierarchy as sch
import networkx as nx
import plotly.express as px
import plotly.io as pio
pio.renderers.default='browser'

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


def main():
    inicio = datetime.now()
    inicio_time = inicio.strftime("%H:%M:%S")
    print("Current Time =", inicio_time)
    
    
    dt_i = '2005-01-04'
    dt_i = '2007-01-03'
    dt_f = None
    
    ts_i = pd.date_range(end=pd.Timestamp.today(), periods=3, freq="M")[0]
    ts_f = pd.Timestamp.today()
    if dt_f == None:
        dt_f = ts_f.strftime('%Y-%m-%d')
        
    #dt_i = ts_i.strftime('%Y-%m-%d')
    
    risk_fact = [
        'IBOV',
        'SPX',
        'DI 1y', 
        'USDBRL',
                 
        'EURUSD',
        'AUDUSD',
        'USDCNH',
        'USDCLP',
        'USDMXN',
        'USDZAR',
   
        'MXN_1y',
        'USGT10',
        'IHFA',
        'CDI Acum'
        ]
    
    fx_global = [
        'EURUSD',
        'AUDUSD',
        'USDCNH',
        'USDCLP',
        'USDMXN',
        'USDZAR'
        ]
    
    bskt_w = [-0.214, -0.172, 0.214, 0.157, 0.138, 0.105] 
    
    df_tipo_ret = get_tipo_retorno(risk_fact, table="ProdutoGestao")
    df_tipo_ret = df_tipo_ret.set_index('produto', drop=True)
    df_tipo_ret = df_tipo_ret.rename({'DI_1Y': 'DI 1y'}, axis=0)
    df_tipo_ret = df_tipo_ret.rename({'tipo_ret': 'tipo_retorno'}, axis=1)
    
    ret_diff = df_tipo_ret.index[df_tipo_ret.tipo_retorno == 'diff'].to_list()
    ret_perc = [prod for prod in risk_fact if prod not in ret_diff]
    
    df_produto_risco = pd.read_excel('P:/quant/projetos/Projeto 45/PositioningFundos.xlsb', sheet_name='ProdutoRisco', engine='pyxlsb')
    df_produto_risco = df_produto_risco.set_index('produto_fundo', drop=True)
    df_produto_risco = df_produto_risco.merge(df_tipo_ret, how='left', left_on='risk_fact', right_index=True )
    
    df_fundos = get_fundos_data()
    fundos_sel = sel_fundos(df_fundos)
    dict_cnpj = df_fundos.set_index('cnpj').loc[fundos_sel, 'fundo'].to_dict()
    
    # %% Tratando Dados
    dados_cotas = get_historico_cotas(fundos_sel, dt_i, dt_f=dt_f)
    dados_risk = get_risk_data(risk_fact, dt_i=dt_i, dt_f=dt_f)
    
    dados_risk['dt_referencia'] = pd.to_datetime(dados_risk['dt_referencia'], format="%Y-%m-%d")
    dados_cotas['dt_referencia'] = pd.to_datetime(dados_cotas['dt_referencia'], format="%Y-%m-%d")
    
    dados_risk_pivot=dados_risk.pivot(index='dt_referencia', columns='produto', values='valor').ffill()
    dados_cotas_pivot=dados_cotas.pivot(index='dt_referencia', columns='cnpj', values='valor_cota').replace(0, np.nan).ffill() 
    
    dados_cotas_pivot=ajusta_cotas(dados_cotas_pivot, n=5)
    
    # %%Calculando Retornos
    dados_cotas_ret = dados_cotas_pivot.ffill().pct_change()#.dropna()
    df_diff = (dados_risk_pivot.loc[dados_cotas_pivot.index, ret_diff]/100).diff()#.dropna()
    df_perc = dados_risk_pivot.loc[dados_cotas_pivot.index, ret_perc].pct_change()#.dropna()
    
    df_rf = df_perc['CDI Acum'].copy()
    df_ihfa = df_perc[['IHFA']].copy()
    
    w_fx = pd.DataFrame({'fx': fx_global,'weight': bskt_w})
    
    df_perc['FX-Global'] = (df_perc[fx_global] @ w_fx.set_index('fx')) / w_fx.set_index('fx').sum()
    
    
    df_perc.drop(['IHFA', 'CDI Acum'] + fx_global, axis=1, inplace=True)
    
    dados_risk_ret = pd.concat([df_perc, df_diff], axis=1)
    dados_cotas_ret['IHFA'] = df_ihfa
    #dados_cotas_ret =  dados_cotas_ret - df_rf.to_numpy().reshape(-1,1)
    
    dados_risk_ret.drop(index=dados_risk_ret.index[0], axis=0, inplace=True)
    dados_cotas_ret.drop(index=dados_cotas_ret.index[0], axis=0, inplace=True)
    #dados_risk_ret.rename({'SPX': 'S&P', 'MXN_1y': 'MXN 1y'}, axis=1, inplace=True)
    
    dist = np.sqrt((1-pd.concat([dados_risk_ret],axis=1)[-126:].corr())*0.5)
    plt.matshow(corr)
    plt.colorbar()
    plt.yticks(np.arange(.5,corr.shape[0]+.5),labels)
    plt.xticks(np.arange(.5,corr.shape[0]+.5),labels)
    
    link = sch.linkage(dist)
    dist_sort = dist.iloc[getQuasiDiag(link), getQuasiDiag(link)]
    plt.matshow(dist_sort)
    plt.colorbar()

    #np.diag(np.fliplr(df_historico_retornos[['IBOV','SPX']].cov())).mean()/ np.diag(df_historico_retornos[['IBOV','SPX']].cov()).cumprod()[-1]**0.5

    pca = PCA(n_components=3)
    df_pca = pd.DataFrame(pca.fit_transform(dist_sort))
    
    sch.dendrogram(link)
    
    size = []
    text = []
    thsh_fundos = df_fundos.set_index('cnpj').loc[dados_cotas_ret.T.drop('IHFA').index,'vl_patrim_liq'].mean()
    for f in dist_sort:
        try:
            size.append(df_fundos.set_index('cnpj').loc[f,'vl_patrim_liq']/1e6)
            if df_fundos.set_index('cnpj').loc[f,'vl_patrim_liq'] > thsh_fundos:
                text.append(" ".join(
                    dict_cnpj[f].split(" ")[:2]))
            else:
                text.append('')
        except:
            size.append(1e10/1e6)
            text.append(f)
    #size.insert(dist_sort.index.get_loc('IHFA'),1e10)
    dict_apelido = df_fundos[['cnpj','fundo']].set_index('cnpj').to_dict()['apelido']
    
    dict_apelido2 = df_fundos[['cnpj','fundo']].set_index('cnpj')['fundo'].str.split(" ").apply(
        lambda x: " ".join(x[:2])).to_dict()
    px.imshow(dist.rename(dict_apelido2).rename(dict_apelido2, axis=1))
    px.imshow(dist_sort.rename(dict_apelido2).rename(dict_apelido2, axis=1))
    px.scatter(df_pca,
           x=0,y=1,#z=2,
           hover_name=dist_sort.rename(dict_cnpj).index, 
           color=sch.fcluster(link,6,criterion='maxclust').astype(str), 
                              size=size,
                              text=text)    
     px.scatter_3d(df_pca,
           x=0,y=1,z=2,
           hover_name=dist_sort.rename(dict_cnpj).index, 
           color=sch.fcluster(link,6,criterion='maxclust').astype(str), 
                              size=size,
                              text=text)    
