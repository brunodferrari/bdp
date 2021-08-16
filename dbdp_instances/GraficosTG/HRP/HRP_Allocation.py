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

from HRP import *

def main():
    wb = xw.Book("HRP_Allocation.xlsm")
    df_historico_retornos = wb.sheets['input'].range('A1').expand().options(pd.DataFrame).value
    
    cov, corr = df_historico_retornos.cov(), df_historico_retornos.corr()
    
    link = getClusters(df_historico_retornos)
    sortIx = getQuasiDiag(link)
    
    labels_true =  df_historico_retornos.columns.tolist() 
    labels_sorted = df_historico_retornos.columns[sortIx].tolist() 
    
    
    fig_dendo = plotDendrogram(link, labels=labels_true, leaf_rotation=90)
    fig_corr = plotCorrMatrix(corr.iloc[sortIx,sortIx], labels=labels_sorted)
    
    sht = wb.sheets['output_dashboard']
    sht.pictures.api.Delete()
    sht.pictures.add(fig_corr, update=False, left=sht.range('A1').left, top=sht.range('A1').top)
    sht.pictures.add(fig_dendo, update=False, left=sht.range('G1').left, top=sht.range('G1').top)
    
    df_hrp = getHRP(df_historico_retornos, link=link, sortIx=sortIx)
    
    sht.range('M1').expand('table').clear_contents()
    sht.range('M1').expand('table').options(index=False).value = df_hrp.reset_index().rename({'index':'Ativo', 0:'Peso'},axis=1)

if __name__ == "__main__":
    main()
