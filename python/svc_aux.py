# -*- coding: utf-8 -*-
"""
Created on Sun Jul 18 22:40:52 2021

@author: Bruno Ferrari
"""
import pandas as pd
import numpy as np

from sklearn.datasets import make_classification
from sklearn.svm import SVC
X, y = make_classification(n_samples=1000,
                           n_features=10, 
                           n_informative=5, 
                           n_redundant=5,
                           n_classes=3, 
                           random_state=1)

X = pd.DataFrame(X)
y = pd.Series(y)

y_dumm = pd.get_dummies(y)

svc = SVC(verbose=1)
svc.fit(X,y)
y_hat = svc.predict(X)

#acuracia OvR sklearn
(y==y_hat).mean()

svc_0 = SVC()
svc_0.fit(X,y_dumm[0])
desc_0 = pd.DataFrame(svc_0.decision_function(X))

svc_1 = SVC()
svc_1.fit(X,y_dumm[1])
desc_1 = pd.DataFrame(svc_1.decision_function(X))

svc_2 = SVC()
svc_2.fit(X,y_dumm[2])
desc_2 = pd.DataFrame(svc_2.decision_function(X))

##### acc feita com modelos individuais
(np.argmax(np.hstack([desc_0, desc_1, desc_2]),axis=1)==y_hat).mean()
