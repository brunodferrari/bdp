# -*- coding: utf-8 -*-
"""
Created on Sat Aug  7 13:38:07 2021

@author: bferrari
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from sklearn.svm import SVC
from itertools import combinations

from xgboost import XGBClassifier
from sklearn.pipeline import Pipeline
from sklearn.model_selection import cross_val_score, cross_val_predict, GridSearchCV, ParameterGrid, train_test_split
from sklearn.ensemble import RandomForestClassifier
from sklearn.tree import DecisionTreeClassifier
from sklearn.neighbors import KNeighborsClassifier
from sklearn.metrics import confusion_matrix, recall_score, f1_score, ConfusionMatrixDisplay, plot_confusion_matrix
from sklearn.decomposition import PCA
from sklearn.linear_model import LogisticRegressionCV
from sklearn.preprocessing import StandardScaler

from sklearn.multiclass import OneVsRestClassifier
from sklearn.metrics import roc_curve, auc
from sklearn.preprocessing import label_binarize


dados_org = pd.read_excel('final_results.xlsx')
dados_scatter = pd.concat([dados_org.iloc[:,1:11]],axis=1)

X = dados_scatter.drop(['V1', 'V2', 'deg_min'],axis=1)

crossing_mh = ['Crossing_vns', 'Crossing_ts', 'Crossing_gs_vns']
best_mh = dados_org[crossing_mh].idxmin(axis=1).fillna('Crossing_gs_vns')
dummy_mh = pd.get_dummies(best_mh)

best_mh_labeled = best_mh.replace({'Crossing_vns': 1,
                                   'Crossing_ts': 1, 
                                   'Crossing_gs_vns': 3})

###################### SVM ##########################
poly_kernel_svm_clf = Pipeline([ 
        ("scaler", StandardScaler()), 
        ("model", SVC(probability=True,random_state=42) 
    ]) 

tuned_parameters = [{'model__kernel': ['rbf'], 'model__gamma': [1e4, 1e3, 1e2, 1e1, 1, 1e-3, 1e-4],
                     'model__C': [1e3, 1e2, 1e1, 1, 1e-3, 1e-4]},
                    {'model__kernel': ['poly'], 'model__degree': [1,2], 'model__C': [1e3, 1e2, 1e1, 1, 1e-3, 1e-4]}]

svc_grid = GridSearchCV(poly_kernel_svm_clf, tuned_parameters, n_jobs=-1, verbose=1)

best_svm = svc_grid.fit(X, best_mh_labeled)

acc_svm = cross_val_score(best_svm.best_estimator_, X, best_mh_labeled, n_jobs=5, verbose=1, scoring='accuracy')
print('Accuracy')
print(round(acc_svm.mean(), 2))
print(round(acc_svm.std(), 2))
roc_svm = cross_val_score(best_svm.best_estimator_, X, best_mh_labeled, n_jobs=5, verbose=1, scoring='roc_auc_ovr')
print('ROC')
print(round(roc_svm.mean(), 2))
print(round(roc_svm.std(), 2))
svm_hat = cross_val_predict(best_svm.best_estimator_, X, best_mh_labeled, n_jobs=-1, verbose=1)

###################### AD ##########################
Xtr, Xval, Ytr, Yval = train_test_split(X, best_mh_labeled)
from sklearn.tree import plot_tree

tuned_parameters = [{'model__max_features': [None, 'sqrt'], 'model__max_depth': [2,3,5,8,13,21,34],
                     'model__min_samples_leaf':[3,5,8,13,21,34]}]

ad = Pipeline([ 
        ("model", DecisionTreeClassifier(max_depth=2, random_state=42))
    ]) 

ad.fit(Xtr,Ytr)
print(ad.score(Xval,Yval))
import pylab as pl
pl.figure(figsize=(20,16))
plot_tree(ad.named_steps['model'], feature_names=X.columns)


ad_grid = GridSearchCV(ad, tuned_parameters, n_jobs=5, verbose=1)
best_ad = ad_grid.fit(X, best_mh_labeled)
acc_ad = cross_val_score(best_ad.best_estimator_, X, best_mh_labeled, n_jobs=5, verbose=1, scoring='accuracy')
print('Accuracy')
print(round(acc_ad.mean(), 2))
print(round(acc_ad.std(), 2))
print('ROC')
roc_ad = cross_val_score(best_ad.best_estimator_, X, best_mh_labeled, n_jobs=5, verbose=1, scoring='roc_auc_ovr')
print(round(roc_ad.mean(), 2))
print(round(roc_ad.std(), 2))
ad_hat = cross_val_predict(best_ad.best_estimator_, X, best_mh_labeled, n_jobs=-1, verbose=1)


###################### RF ##########################
Xtr, Xval, Ytr, Yval = train_test_split(X, best_mh_labeled)

tuned_parameters = [{'model__max_features': [None, 'sqrt'], 
                     'model__min_samples_leaf': [1,2,3,5,8,13,21,34],
                     'model__n_estimators': [10, 100, 1000]}]

rf = Pipeline([ 
        ("model", RandomForestClassifier(random_state=42,n_jobs=-1))
    ]) 

rf.fit(Xtr,Ytr)
print(rf.score(Xval,Yval))



rf_grid = GridSearchCV(rf, tuned_parameters, n_jobs=5, verbose=1)
best_rf = rf_grid.fit(X, best_mh_labeled)
acc_rf = cross_val_score(best_rf.best_estimator_, X, best_mh_labeled, n_jobs=5, verbose=1, scoring='accuracy')
print('Accuracy')
print(round(acc_rf.mean(), 2))
print(round(acc_rf.std(), 2))
print('ROC')
roc_rf  = cross_val_score(best_rf.best_estimator_, X, best_mh_labeled, n_jobs=5, verbose=1, scoring='roc_auc_ovr')
print(round(roc_rf.mean(), 2))
print(round(roc_rf.std(), 2))
rf_hat = cross_val_predict(best_rf.best_estimator_, X, best_mh_labeled, n_jobs=-1, verbose=1)


###################### KNN ##########################
Xtr, Xval, Ytr, Yval = train_test_split(X, best_mh_labeled)

tuned_parameters = [{'scaler': ['passthrough', StandardScaler()],
                    'model__n_neighbors': [3,5,7,13,15,21], 
                     'model__p': [1,2,3,5]}]
                    

knn = Pipeline([ 
    ("scaler", StandardScaler()), 
    ("model", KNeighborsClassifier(n_neighbors=15, p=5, n_jobs=-1))
    ]) 

knn.fit(Xtr,Ytr)
print(knn.score(Xval,Yval))

knn_grid = GridSearchCV(knn, tuned_parameters, n_jobs=5, verbose=1)
best_knn = knn_grid.fit(X, best_mh_labeled)

acc_knn = cross_val_score(best_knn.best_estimator_, X, best_mh_labeled, n_jobs=5, verbose=1, scoring='accuracy')
print('Accuracy')
print(round(acc_knn.mean(), 2))
print(round(acc_knn.std(), 2))
print('ROC')
roc_knn  = cross_val_score(best_knn.best_estimator_, X, best_mh_labeled, n_jobs=5, verbose=1, scoring='roc_auc_ovr')
print(round(roc_knn.mean(), 2))
print(round(roc_knn.std(), 2))
knn_hat = cross_val_predict(best_knn.best_estimator_, X, best_mh_labeled, n_jobs=-1, verbose=1)

###################### CM ##########################
#fig, ax = plt.subplots(2, 2, figsize=(12,12))
#pos = [(i,j) for i in range(2) for j in range(2)]
for var in globals():
    #print(var)
    if str(var).find("_hat") > 0:
        #pos_ = pos.pop(0)
        ConfusionMatrixDisplay(confusion_matrix(best_mh_labeled, globals()[var]), ['VND/TABU', 'GRASP']).plot()
        pl.title(var.split("_")[0].upper())
        pl.xlabel('Valor Previsto')
        pl.ylabel('Valor Atribuído')                                                                                                        
        