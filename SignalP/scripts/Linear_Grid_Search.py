# -*- coding: utf-8 -*-
"""
Created on Tue Feb 21 21:00:20 2017

@author: Revant
"""
# General Libraries

import time
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
from collections import OrderedDict as oD
#import warnings
#warnings.filterwarnings("ignore", category=DeprecationWarning)
    
# Sklearn modules

from sklearn.preprocessing import normalize
from sklearn.ensemble import BaggingClassifier, RandomForestClassifier
from sklearn.multiclass import OneVsRestClassifier
from sklearn.svm import LinearSVC
from sklearn.svm import SVC
from sklearn.metrics import confusion_matrix
from sklearn.metrics import precision_recall_fscore_support as score
from sklearn.grid_search import GridSearchCV
from sklearn.metrics import f1_score

#Custom modules

import dense_data_parser as ddp

# Set parameters

window_size = [2,5,8,10,13]

filepath = '''/home/u2196/Desktop/KB8024/KB8024/data/globular_signal_tm_3state.txt'''
output = '''/home/u2196/Desktop/KB8024/KB8024/SignalP/output/Linear_grid_search/'''

# Starting script

start = time.time()

final_list = oD()

for windows in window_size:
 
    data = ddp.pre_vec_parser(filepath, windows)

    clf = LinearSVC(class_weight = 'balanced', cache_size = 2000)

    X, Y = ddp.skl_parser(data, windows)

    parameters = {"C": [1,2,4,8] }

    model_tunning = GridSearchCV(clf, param_grid=parameters, score_func=f1_score)

    model_tunning.fit(X,Y)

    s = model_tunning.best_score_
    p = model_tunning.best_params_
    
    final_list[str(windows)] = [s,p]       

end = time.time()

best_table = pd.DataFrame.from_dict(final_list, orient='columns')
best_table.to_csv(output+'Linear_grid_search.txt')

print("Script took %f time"%(end-start))
'''
print("Scoring over cross validation sets...")

start = time.time()

score = cross_val_score(LinearSVC(), X, Y, cv=3)

end = time.time()

print ("Time taken for cross validation scoring : %f"%(end-start))
print("Individual scores: ",score)
print("Overall Accuracy: %0.2f (+/- %0.2f)" % (score.mean(), score.std() * 2))
print('\n')

print("One vs Rest Classification using SVC")
start = time.time()
clf = OneVsRestClassifier(SVC(kernel='linear', probability=True, class_weight='auto'))
clf.fit(X, Y)
end = time.time()
print ("Single SVC", end - start, clf.score(X,Y))
proba = clf.predict_proba(X)
print('\n')

print("One vs Rest SVC with Bagging")
n_estimators = 10
start = time.time()
clf = OneVsRestClassifier(BaggingClassifier(SVC(kernel='linear', probability=True, class_weight='auto'), max_samples=1.0 / n_estimators, n_estimators=n_estimators))
clf.fit(X, Y)
end = time.time()
print ("Bagging SVC", end - start, clf.score(X,Y))
proba = clf.predict_proba(X)
print('\n')

print("Random forest")
start = time.time()
clf = RandomForestClassifier(min_samples_leaf=20)
clf.fit(X, Y)
end = time.time()
print ("Random Forest", end - start, clf.score(X,Y))
proba = clf.predict_proba(X)
'''


