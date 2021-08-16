# -*- coding: utf-8 -*-
"""
Created on Sun Jul 18 22:39:53 2021

@author: Bruno Ferrari
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




dados_org = pd.read_excel('final_results.xlsx')
dados_scatter = pd.concat([dados_org.iloc[:,1:11]],axis=1)

X = dados_scatter.drop(['V1', 'V2', 'deg_min'],axis=1)

crossing_mh = ['Crossing_vns', 'Crossing_ts', 'Crossing_gs', 'Crossing_gs_vns']

best_mh = dados_org[crossing_mh].idxmin(axis=1)#.fillna('Crossing_gs_vns')
dummy_mh = pd.get_dummies(best_mh)
best_mh_labeled = best_mh.replace({'Crossing_vns': 1,
                                   'Crossing_ts': 2, 
                                   'Crossing_gs': 3,
                                   'Crossing_gs_vns': 4})#.fillna(1).astype(int)

z_cross = ((dados_org[crossing_mh] - 
            dados_org[crossing_mh].mean(axis=1).values.reshape(-1,1)) / 
            dados_org[crossing_mh].std(axis=1,ddof=0).values.reshape(-1,1))

time_mh = ['Time_vns', 'Time_ts','Time_gs_vns']
z_time = ((dados_org[time_mh] - 
           dados_org[time_mh].mean(axis=1).values.reshape(-1,1)) / 
           dados_org[time_mh].std(axis=1,ddof=0).values.reshape(-1,1))

w_time=1.5
wbest_mh = pd.DataFrame((z_cross.values + w_time*z_time.values)/2,
                        columns= ['Crossing_vns', 'Crossing_ts', 'Crossing_gs_vns']).idxmin(axis=1).fillna('Crossing_gs_vns')

wbest_mh_labeled = wbest_mh.replace({'Crossing_ts': 3, 'Crossing_vns': 2, 'Crossing_gs_vns': 1})
dummy_wmh = pd.get_dummies(wbest_mh)

for col in dummy_mh:
    dummy_mh[col].plot(style='*', figsize=(8,6), title=col)
    plt.show()
dummy_mh.plot(style='*', figsize=(8,6), title=col)    

for col in dummy_wmh:
    dummy_wmh[col].plot(style='*', figsize=(8,6), title=col)
    plt.show()
dummy_wmh.plot(style='*', figsize=(8,6), title=col)

for col in dados_scatter:
    pd.plotting.scatter_matrix(pd.concat([dados_scatter[col], best_mh_labeled], axis=1))
    plt.show()

# sensibilidade do w_time pro mh
for mh in crossing_mh[2:3]:
    dict_w={}
    for w_time in np.logspace(0, np.log10(2)):
        dict_w[w_time] = pd.DataFrame((z_cross.values + w_time*z_time.values)/2, 
                     columns=crossing_mh).idxmin(axis=1).fillna('Crossing_gs_vns').value_counts()[mh]
    pd.Series(dict_w).plot(title=mh)
    plt.show()

# nuvem plot mh
X_Y = list(combinations(dados_scatter.columns, 2))

for cols in X_Y:
    ax=pd.concat([dados_scatter[list(cols)], best_mh_labeled], axis=1).plot(kind='scatter', x=cols[0], y=cols[1], c=0, cmap="viridis", sharex=False,
                                                                          title='Time Constrained', label=crossing_mh)
    plt.legend(handles = ax.legend_.legendHandles[0].legend_elements()[0], 
               labels=crossing_mh)
  #  plt.show()    

from sklearn.svm import LinearSVC
from sklearn.decomposition import PCA
from sklearn.metrics import plot_confusion_matrix

svc_mh = LinearSVC(max_iter=1e8)
svc_mh = SVC(gamma=5)

svc_mh.fit(X, best_mh_labeled)
y_hat_mh = svc_mh.predict(X)
print(best_mh_labeled.value_counts(normalize=1))
(best_mh_labeled==y_hat_mh).mean()
plot_confusion_matrix(svc_mh, X, best_mh_labeled)

w_time=1.5
wbest_mh = pd.DataFrame((z_cross.values + w_time*z_time.values)/2,
                        columns=crossing_mh).idxmin(axis=1).fillna('Crossing_gs_vns')
wbest_mh_labeled = wbest_mh.replace({'Crossing_ts': 3, 'Crossing_vns': 2, 'Crossing_gs_vns': 1}).astype(int)
print(wbest_mh.value_counts())
print(wbest_mh_labeled.value_counts(normalize=1))
print(wbest_mh_labeled.replace({3:2}).value_counts(normalize=1))

svc_wmh = SVC(kernel="rbf", gamma=5)
svc_wmh.fit(X, 
            wbest_mh_labeled.replace({3:2})
)
y_hat_wmh = svc_wmh.predict(X)
(wbest_mh_labeled.replace({3:2})==y_hat_wmh).mean()



###########################
#https://medium.com/swlh/the-hyperparameter-cheat-sheet-770f1fed32ff
results = svc_grid.fit(X, best_mh_labeled)
print(results.best_score_)
mod = cross_val_score(svc_grid, X, best_mh_labeled, n_jobs=-1, verbose=1)
print(mod, mod.mean())
y_hat = cross_val_predict(svc_grid, X, best_mh_labeled, n_jobs=-1, verbose=1)
ConfusionMatrixDisplay(confusion_matrix(best_mh_labeled, y_hat)).plot()


wresults = svc_grid.fit(X, wbest_mh_labeled)
print(wresults.best_score_)
wmod = cross_val_score(svc_grid, X, wbest_mh_labeled, n_jobs=-1, verbose=1)
print(wmod, wmod.mean())
print(cross_val_predict(svc_grid, X, wbest_mh_labeled, n_jobs=-1, verbose=1))


wwresults = svc_grid.fit(X, wbest_mh_labeled.replace({3:2}))
print(wwresults.best_score_)
wwmod = cross_val_score(svc_grid, X, wbest_mh_labeled.replace({3:2}), n_jobs=-1, verbose=1)
print(wwmod, wmod.mean())
print(cross_val_predict(svc_grid, X, wbest_mh_labeled.replace({3:2}), n_jobs=-1, verbose=1))

###########################


tuned_parameters = [{'model__kernel': ['rbf'], 'model__gamma': [1e4, 1e3, 1e2, 1e1, 1, 1e-3, 1e-4],
                     'model__C': [1e4, 1e3, 1e2, 1e1, 1, 1e-3, 1e-4]},
                    {'model__kernel': ['poly'], 'model__degree': [1,2, 3], 'model__C': [1e4, 1e3, 1e2, 1e1, 1, 1e-3, 1e-4]}]

from sklearn.preprocessing import StandardScaler
poly_kernel_svm_clf = Pipeline([ 
        ("scaler", StandardScaler()), 
        ("model", SVC()) 
    ]) 

paramgrid=[
    {'model':[XGBClassifier()],
     'model__max_depth':[2,3,5,7,10,15],
     'model__n_estimators':[50,100,150,200],
     'model__booster': ['gbtree', 'gblinear','dart']}
    ]

svc_grid = GridSearchCV(poly_kernel_svm_clf, tuned_parameters, n_jobs=-1, verbose=1)




poly_kernel_svm_clf.fit(X, wbest_mh_labeled)

poly_kernel_svm_clf.predict(X)

pipeline = Pipeline([('model', XGBClassifier())])

paramgrid=[
    {'model':[XGBClassifier()],
     'model__max_depth':[2,3,5,7,10,15],
     'model__n_estimators':[50,100,150,200],
     'model__booster': ['gbtree', 'gblinear','dart']}
    ]

#grid = GridSearchCV(pipeline, paramgrid, n_jobs=-1, scoring='f1_macro')
grid = GridSearchCV(pipeline, paramgrid, n_jobs=-1)
results = grid.fit(X, best_mh_labeled)


desc={}
for mh in crossing_mh:
    svc = SVC(gamma=5)
    svc.fit(X, dummy_mh[mh])
    y_hat = svc.predict(X)
    desc[mh] = svc.decision_function(X)
    print(dummy_mh[mh].value_counts(normalize=1))
    print( (dummy_mh[mh]==y_hat).mean())

######################### PCA plots

pca = PCA()
scaling_pca = Pipeline([("scaler", StandardScaler()), 
                        ("pca", PCA()),
                        ])

px.scatter_matrix(pd.DataFrame(
                    np.hstack(
                        [best_mh.values.reshape(-1,1),
                        scaling_pca.fit_transform(X)]
                        )
                    ), color=0 )

##################### SVC with scaler and PCA

pipe = Pipeline([("scaler", StandardScaler()), 
                 ("pca", PCA()),
                 ("model", SVC())
                 ])

pipe.fit(X, best_mh_labeled)
y_hat_mh = pipe.predict(X)
print(best_mh_labeled.value_counts(normalize=1))
print((best_mh_labeled==y_hat_mh).mean())
plot_confusion_matrix(pipe, X, best_mh_labeled)



tuned_parameters = [{
                     'model__kernel': ['rbf'], 
                     'model__gamma': [1e5, 1e4, 1e3, 1e2, 1e1, 1, 1e-3, 1e-4],
                     'model__C': [1e4, 1e3, 1e2, 1e1, 1, 1e-3, 1e-4],
                     },
                    {'model__kernel': ['poly'], 
                     'model__degree': [1,2, 3], 
                     'model__C': [1e4, 1e3, 1e2, 1e1, 1, 1e-3, 1e-4]
                     },
                     {'model':[XGBClassifier()],
                      'model__max_depth':[2,3,5,7,10,15],
                      'model__n_estimators':[50,100,150,200],
                      'model__booster': ['gbtree', 'gblinear','dart']}               
                    {'mode'}
                    ]
mod_grid = GridSearchCV(pipe, tuned_parameters, n_jobs=-1, verbose=1)

results = mod_grid.fit(X, best_mh_labeled)
print(results.best_score_)
mod = cross_val_score(mod_grid, X, best_mh_labeled, n_jobs=-1, verbose=1)
print(mod, mod.mean())
y_hat = cross_val_predict(mod_grid, X, best_mh_labeled, n_jobs=-1, verbose=1)
print(cross_val_predict(mod_grid, X, best_mh_labeled, n_jobs=-1, verbose=1))


ConfusionMatrixDisplay(confusion_matrix(best_mh_labeled, y_hat)).plot()
plot_confusion_matrix(pipe, X, best_mh_labeled)

import plotly.express as px
import plotly.io as pio

pio.renderers.default='browser'

px.scatter_matrix(pd.DataFrame(np.hstack([best_mh_labeled.astype(str).values.reshape(-1,1),scaling_pca.fit_transform(X)])),
                  color=0)

from sklearn.decomposition import KernelPCA 
rbf_pca = KernelPCA(n_components = 2, kernel="rbf", gamma=1e-1) 
X_reduced = rbf_pca.fit_transform(X)

px.scatter_matrix(pd.DataFrame(np.hstack([best_mh_labeled.astype(str).values.reshape(-1,1),X_reduced])),
                  color=0)
