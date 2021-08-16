# -*- coding: utf-8 -*-
"""
Created on Sat Aug  7 02:01:25 2021

@author: bferrari
"""
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from sklearn.svm import SVC
from itertools import combinations

from xgboost import XGBClassifier
from sklearn.pipeline import Pipeline
from sklearn.model_selection import cross_val_score, cross_val_predict, GridSearchCV, ParameterGrid
from sklearn.svm import SVC
from sklearn.ensemble import RandomForestClassifier
from sklearn.tree import DecisionTreeClassifier
from sklearn.neighbors import KNeighborsClassifier
from sklearn.metrics import confusion_matrix, recall_score, f1_score, ConfusionMatrixDisplay

import pylab as pl


dados_org = pd.read_excel('final_results.xlsx')
dados_scatter = pd.concat([dados_org.iloc[:,1:11]],axis=1)


#Graficos TG #nos
dados_org[dados_org['Instance'].str.contains('scr')][['V', 'V1', 'V2', 'V_diff']].hist(figsize=(8,6))
pl.suptitle("Histograma # Nós - Conjunto 1")

dados_org[1000:1120][['V', 'V1', 'V2', 'V_diff']].hist(figsize=(8,6))
pl.suptitle("Histograma # Nós - Conjunto 2")

dados_org[1120:1240][['V', 'V1', 'V2', 'V_diff']].hist(figsize=(8,6))
pl.suptitle("Histograma # Nós - Conjunto 3")

dados_org[['V', 'V1', 'V2', 'V_diff']].hist(figsize=(8,6))
pl.suptitle("Histograma # Nós - Geral")

################
dados_org[dados_org['Instance'].str.contains('scr')][['V', 'V1', 'V2', 'V_diff']].boxplot(figsize=(8,6))
pl.suptitle("Boxplot # Nós - Conjunto 1")

dados_org[1000:1120][['V', 'V1', 'V2', 'V_diff']].boxplot(figsize=(8,6))
pl.suptitle("Boxplot # Nós - Conjunto 2")

dados_org[1120:1240][['V', 'V1', 'V2', 'V_diff']].boxplot(figsize=(8,6))
pl.suptitle("Boxplot # Nós - Conjunto 3")

dados_org[['V', 'V1', 'V2', 'V_diff']].boxplot(figsize=(8,6))
pl.suptitle("Boxplot # Nós - Geral")

#########################################################



#Densidade
pd.concat([dados_org[dados_org['Instance'].str.contains('scr')][['density']].rename({'density':"Conjunto 1"}, axis=1),
            dados_org[dados_org['Instance'].str.contains('inc')][['density']].rename({'density':"Conjunto 2"}, axis=1),
            dados_org[(~dados_org['Instance'].str.contains('scr')) & 
                      (~dados_org['Instance'].str.contains('inc'))][['density']].rename({'density':"Conjunto 3"}, axis=1),
            dados_org[['density']].rename({'density':"Geral"}, axis=1)],
            axis=1).boxplot(figsize=(8,6))

pl.suptitle("Boxplot Densidades")


pd.concat([dados_org[dados_org['Instance'].str.contains('scr')][['density']].rename({'density':"Conjunto 1"}, axis=1),
            dados_org[dados_org['Instance'].str.contains('inc')][['density']].rename({'density':"Conjunto 2"}, axis=1),
            dados_org[(~dados_org['Instance'].str.contains('scr')) & 
                      (~dados_org['Instance'].str.contains('inc'))][['density']].rename({'density':"Conjunto 3"}, axis=1),
            dados_org[['density']].rename({'density':"Geral"}, axis=1)],
            axis=1).hist(figsize=(8,6))
pl.suptitle("Histograma Densidades")


################################ edges
pd.concat([dados_org[dados_org['Instance'].str.contains('scr')][['E']].rename({'E':"Conjunto 1"}, axis=1),
            dados_org[dados_org['Instance'].str.contains('inc')][['E']].rename({'E':"Conjunto 2"}, axis=1),
            dados_org[(~dados_org['Instance'].str.contains('scr')) & 
                      (~dados_org['Instance'].str.contains('inc'))][['E']].rename({'E':"Conjunto 3"}, axis=1),
            dados_org[['E']].rename({'E':"Geral"}, axis=1)],
            axis=1).boxplot(figsize=(8,6))

pl.suptitle("Boxplot # Arestas")


pd.concat([dados_org[dados_org['Instance'].str.contains('scr')][['E']].rename({'E':"Conjunto 1"}, axis=1),
            dados_org[dados_org['Instance'].str.contains('inc')][['E']].rename({'E':"Conjunto 2"}, axis=1),
            dados_org[(~dados_org['Instance'].str.contains('scr')) & 
                      (~dados_org['Instance'].str.contains('inc'))][['E']].rename({'E':"Conjunto 3"}, axis=1),
            dados_org[['E']].rename({'E':"Geral"}, axis=1)],
            axis=1).hist(figsize=(8,6))
pl.suptitle("Histograma # Arestas")




# graus
dados_org[dados_org['Instance'].str.contains('scr')][dados_org.columns[dados_org.columns.str.contains('deg')]].boxplot(figsize=(8,6))
pl.title('Boxplot Grau- Conjunto 1')

dados_org[1000:1120][dados_org.columns[dados_org.columns.str.contains('deg')]].boxplot(figsize=(8,6))
pl.title('Boxplot Grau- Conjunto 2')

dados_org[1120:1240][dados_org.columns[dados_org.columns.str.contains('deg')]].boxplot(figsize=(8,6))
pl.title('Boxplot Grau - Conjunto 3')

dados_org[:][dados_org.columns[dados_org.columns.str.contains('deg')]].boxplot(figsize=(8,6))
pl.title('Boxplot Grau - Geral')


fig, ax = plt.subplots(2, 2, figsize=(8,6))
dados_org[dados_org['Instance'].str.contains('scr')][['deg_mean']].hist(figsize=(8,6), ax=ax[0,0])
dados_org[dados_org['Instance'].str.contains('scr')][['deg_median']].hist(figsize=(8,6), ax=ax[0,1])
dados_org[dados_org['Instance'].str.contains('scr')][['deg_std']].hist(figsize=(8,6), ax=ax[1,0])
dados_org[dados_org['Instance'].str.contains('scr')][['deg_mean', 'deg_std', 'deg_median']].boxplot(figsize=(8,6),ax=ax[1,1])
plt.suptitle('Grau - Conjunto 1')


fig, ax = plt.subplots(2, 2, figsize=(8,6))
dados_org[1000:1120][['deg_mean']].hist(figsize=(8,6), ax=ax[0,0])
dados_org[1000:1120][['deg_median']].hist(figsize=(8,6), ax=ax[0,1])
dados_org[1000:1120][['deg_std']].hist(figsize=(8,6), ax=ax[1,0])
dados_org[1000:1120][['deg_mean', 'deg_std', 'deg_median']].boxplot(figsize=(8,6),ax=ax[1,1])
plt.suptitle('Grau - Conjunto 2')

fig, ax = plt.subplots(2, 2, figsize=(8,6))
dados_org[1120:1240][['deg_mean']].hist(figsize=(8,6), ax=ax[0,0])
dados_org[1120:1240][['deg_median']].hist(figsize=(8,6), ax=ax[0,1])
dados_org[1120:1240][['deg_std']].hist(figsize=(8,6), ax=ax[1,0])
dados_org[1120:1240][['deg_mean', 'deg_std', 'deg_median']].boxplot(figsize=(8,6),ax=ax[1,1])
plt.suptitle('Grau - Conjunto 3')



fig, ax = plt.subplots(2, 2, figsize=(8,6))
dados_org[:][['deg_mean']].hist(figsize=(8,6), ax=ax[0,0])
dados_org[:][['deg_median']].hist(figsize=(8,6), ax=ax[0,1])
dados_org[:][['deg_std']].hist(figsize=(8,6), ax=ax[1,0])
dados_org[:][['deg_mean', 'deg_std', 'deg_median']].boxplot(figsize=(8,6),ax=ax[1,1])
plt.suptitle('Grau - Geral')

##################
dados_org.loc[:,dados_org.columns.str.contains('Crossing')].rename({'Crossing_vns': 'VND',
                                   'Crossing_ts': 'TABU', 
                                   'Crossing_gs_vns': 'GRASP'},axis=1).dropna().drop('Crossing_gs', axis=1).plot(
                                       logy=True, title='# Cruzamento de Arestas',figsize= (8,6))
pl.xlabel('Instâncias')   
pl.ylabel('Cruzamentos')

dados_org.loc[:,dados_org.columns.str.contains('Time')].rename({'Time_vns': 'VND',
                                   'Time_ts': 'TABU', 
                                   'Time_gs_vns': 'GRASP'},axis=1).dropna().drop('Time_gs', axis=1).plot(
                                       logy=True, title='Tempo de Execução (s)',figsize= (8,6))
pl.xlabel('Instâncias')   
pl.ylabel('S')
################################ PAIRWISE PLOT ################################

X_Y = list(combinations(X.columns, 2))

for cols in X_Y[:]:
    ax=pd.concat([X[list(cols)], best_mh_labeled.replace(1,2)], axis=1).rename({0:'MH'},axis=1)[:].plot(kind='scatter', x=cols[0], y=cols[1], c='MH', cmap="viridis", sharex=False,
                                                                          title='', label=['TABU', 'GRASP'])
    plt.legend(handles = ax.legend_.legendHandles[0].legend_elements()[0], 
               labels=['TABU', 'GRASP'])


fig, ax = plt.subplots(7, 3,figsize=(15,25),constrained_layout=False)
ij = [(i,j) for i in range(3) for j in range(7)]
for cols, pos in zip(X_Y[:], ij):
    aux = pd.concat([X[list(cols)], best_mh_labeled.replace(1,2)], axis=1).rename({0:'MH'},axis=1)
    
    #fig2 = plt.figure()
    ax[pos[1],pos[0]].scatter(aux.iloc[:,0].loc[aux['MH']==2], aux.iloc[:,1].loc[aux['MH']==2], c = 'black', marker='o', label='VND/TABU')
    ax[pos[1],pos[0]].scatter(aux.iloc[:,0].loc[aux['MH']==3], aux.iloc[:,1].loc[aux['MH']==3], c = 'orange', marker='.', label='GRASP')
    
    ax[pos[1],pos[0]].set_xlabel(aux.columns[0])
    ax[pos[1],pos[0]].set_ylabel(aux.columns[1])

    ax[pos[1],pos[0]].legend()
    
    #ax[pos[1],pos[0]]=plt.gcf()
    #plt.show()                                                          
