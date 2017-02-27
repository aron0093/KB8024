# -*- coding: utf-8 -*-
"""
Created on Tue Feb 21 21:00:20 2017

@author: Revant
"""

import SS_parser
import time
import numpy as np
from sklearn.ensemble import BaggingClassifier, RandomForestClassifier
from sklearn.multiclass import OneVsRestClassifier
from sklearn.svm import LinearSVC
from sklearn.svm import SVC
from sklearn.model_selection import cross_val_score

filepath = '''/KB8024/KB8024/data/globular_signal_tm_3state.txt'''
window_size = 3

print("Parsing data...")

X, Y = SS_parser.skl_parser(filepath, window_size)

print("Scoring over cross validation sets...")

start = time.time()

score = cross_val_score(LinearSVC(), X, Y, cv=3)

end = time.time()

print ("Time taken for cross validation scoring : %f"%(end-start))
print("Individual scores: ",score)
print("Overall Accuracy: %0.2f (+/- %0.2f)" % (score.mean(), score.std() * 2))

#start = time.time()
#clf = OneVsRestClassifier(SVC(kernel='linear', probability=True, class_weight='auto'))
#clf.fit(X, y)
#end = time.time()
#print "Single SVC", end - start, clf.score(X,y)
#proba = clf.predict_proba(X)

#n_estimators = 10
#start = time.time()
#clf = OneVsRestClassifier(BaggingClassifier(SVC(kernel='linear', probability=True, class_weight='auto'), max_samples=1.0 / n_estimators, #n_estimators=n_estimators))
#clf.fit(X, y)
#end = time.time()
#print "Bagging SVC", end - start, clf.score(X,y)
#proba = clf.predict_proba(X)

#start = time.time()
#clf = RandomForestClassifier(min_samples_leaf=20)
#clf.fit(X, y)
#end = time.time()
#print "Random Forest", end - start, clf.score(X,y)
#proba = clf.predict_proba(X)



