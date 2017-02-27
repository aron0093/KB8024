# -*- coding: utf-8 -*-
"""
Created on Tue Feb 21 21:00:20 2017

@author: Revant
"""

import SS_parser
import time

filepath = '''/KB8024/KB8024/data/globular_signal_tm_3state.txt'''
window_size = 3

X, Y = SS_parser.skl_parser(filepath, window_size)

print("Data parsing successful")

from sklearn.svm import SVC
from sklearn.model_selection import cross_val_score

print("Starting cross validation")

start = time.time()

score = cross_val_score(SVC(kernel ='linear'), X, Y, cv=2)

end = time.time()

print ("Time taken for cross validation %f"%(end-start))
print('\n')
print(score)
print('\n')
print("Accuracy: %0.2f (+/- %0.2f)" % (score.mean(), score.std() * 2))


