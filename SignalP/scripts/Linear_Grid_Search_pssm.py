# -*- coding: utf-8 -*-
"""
Created on Tue Feb 21 21:00:20 2017

@author: Revant
"""
# General Libraries

import time
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt
from collections import OrderedDict as oD
#import warnings
#warnings.filterwarnings("ignore", category=DeprecationWarning)
    
# Sklearn modules

from sklearn.svm import LinearSVC
from sklearn.grid_search import GridSearchCV
from sklearn.metrics import f1_score
from sklearn.metrics import make_scorer

#Custom modules

import pssm_data_parser as pdp

# Set parameters

window_size = [25]

filepath = '''/home/u2196/Desktop/KB8024/KB8024/data/globular_signal_tm_3state.txt'''
output = '''/home/u2196/Desktop/KB8024/KB8024/SignalP/output/Linear_grid_search/'''
pssm_type = 'sub'
pssm_loc = '''/home/u2196/Desktop/KB8024/KB8024/SignalP/input/pssms/'''
# Starting script

start = time.time()

f1_scorer = make_scorer(f1_score, labels=[-1,1], average='macro')

final_list = oD()

for windows in window_size:

    clf = LinearSVC(class_weight = 'balanced')
 
    data = pdp.pre_vec_parser(filepath, windows)

    p_data = pdp.pssm_parser(data, windows, pssm_loc)    

    X, Y = pdp.skl_pssm_parser(p_data, windows, pssm_type=pssm_type)
    
    print("Data Parsing Done...")

    parameters = {"C": [1,2,4,8,10,12,14,16] }

    model_tunning = GridSearchCV(clf, param_grid=parameters, scoring=f1_scorer)

    model_tunning.fit(X,Y)

    s = model_tunning.best_score_
    p = model_tunning.best_params_
    
    final_list[str(windows)] = [s,p['C']]       

end = time.time()

best_table = pd.DataFrame.from_dict(final_list, orient='index')
best_table.columns= ['F1Score', 'C']
best_table.to_csv(output+'Linear_grid_search_'+pssm_type+'_4.csv')

fig = plt.figure() # Create matplotlib figure

ax = fig.add_subplot(111) # Create matplotlib axes
ax2 = ax.twinx() # Create another axes that shares the same x-axis as ax.

best_table.plot(x=best_table.index, kind='bar', secondary_y=best_table.columns[1], mark_right=False)

ax.set_ylabel("F1 Score for Signal and Transmembrane Domains")
ax2.set_ylabel(" C parameter")

plt.xlabel("+/- frames around target residue")
plt.figure(figsize=(20,10)).savefig(output+'graph_'+pssm_type+'_4.png')
plt.close()
print("Script took %f seconds"%(end-start))
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


