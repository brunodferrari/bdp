# -*- coding: utf-8 -*-
"""
Created on Thu Aug  5 16:47:08 2021

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
    get_risk_data,
    gera_regressoes
)

from sklearn.decomposition import PCA

import matplotlib.pyplot as plt

import scipy.cluster.hierarchy as sch
import networkx as nx
import plotly.express as px
import plotly.io as pio
pio.renderers.default='browser'

from HRP import *

#Gera pos
if __name__ == "__main__":
    
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
    dados_cotas_ret =  dados_cotas_ret - df_rf.to_numpy().reshape(-1,1)
    
    dados_risk_ret.drop(index=dados_risk_ret.index[0], axis=0, inplace=True)
    dados_cotas_ret.drop(index=dados_cotas_ret.index[0], axis=0, inplace=True)
    
    reg_dict_sadj, fundos_erros2 = gera_regressoes(dados_cotas_ret.drop('IHFA',axis=1), dados_risk_ret, 21)
    reg_ihfa, _ = gera_regressoes(dados_cotas_ret[['IHFA']], dados_risk_ret, 21)
    reg_ihfa_curta, _ = gera_regressoes(dados_cotas_ret[['IHFA']], dados_risk_ret, 10)

results = []
for i in range(252*10):
    df_historico_retornos=pd.concat([dados_risk_ret],axis=1)[-(126+i):][:126]
    hrp=getHRP(df_historico_retornos)
    results.append(pd.DataFrame({df_historico_retornos.index[-1]: hrp}))

df_results = pd.concat(results,axis=1).T.sort_index()

for col in df_results:
    (df_results[col]*100).plot(title=col)
    plt.show()
    
std = df_results.rolling(252*5, min_periods=252).std()
df_results / std

concat_ix = pd.concat([reg_ihfa['IHFA']['beta'], df_results],axis=1).dropna().index
concat_ix = pd.concat([reg_dict_sadj['21.470.947/0001-36']['beta'], df_results],axis=1).dropna().index
for col in df_results:
    df_results.loc[concat_ix, col][-252:].plot(label='HRP',title=col)
    reg_ihfa['IHFA']['beta'].loc[concat_ix, col][-252:].plot(label='IHFA',title=col, secondary_y=True)
    #reg_dict_sadj['21.470.947/0001-36']['beta'].loc[concat_ix, col][-252:].plot(label='Beta Vertex',title=col, secondary_y=True)
    plt.legend()
    plt.show()
    
    
results2 = []
for i in range(252*1):
    df_historico_retornos=pd.concat([dados_risk_ret, -dados_risk_ret.add_suffix('_inv')],axis=1)[-(126+i):][:126]    
    
    if not(i % 21):
        hrp=getHRP(df_historico_retornos)
        idx = hrp.index
    else:
        hrp=getHRP(df_historico_retornos[idx])
        s
    results2.append(pd.DataFrame({df_historico_retornos.index[-1]: hrp}))
    
df2 = pd.concat(results2,axis=1)

aux = df2[df2.index.str.contains("_inv")].T.sort_index()
aux2 = df2[~df2.index.str.contains("_inv")].T.sort_index()]
    
aux.columns = aux.columns.str.replace("_inv","")
aux = -aux
for col in aux:
    (aux[col]).plot(title=col)
    aux2[col.replace("_inv","")].plot(title=col)
    plt.show()

df_results2=pd.DataFrame()
for col in aux:
    df_results2[col] = pd.concat([aux[col].dropna(), aux2[col].dropna()],axis=0)

for col in df_results2:
    df_results2.sort_index()[col].plot(title=col)
    plt.show()
    
    
#####################################
results = []
for i in range(252*10):
    df_historico_retornos=pd.concat([dados_risk_ret,dados_cotas_ret],axis=1)[-(126+i):][:126]
    hrp=getHRP(df_historico_retornos)
    results.append(pd.DataFrame({df_historico_retornos.index[-1]: hrp}))
'''
for i in range(21,252,21):
    plt.scatter(pca.fit_transform(dados_cotas_ret.corr()[-i:])[0], pca.fit_transform(dados_cotas_ret.corr()[-i:])[1])
    plt.show()
