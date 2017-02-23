# -*- coding: utf-8 -*-
"""
Created on Tue Feb 21 21:00:20 2017

@author: Revant
"""

import data_parser

filepath = '''/KB8024/KB8024/data/globular_signal_tm_3state.txt'''
outpath = '''/KB8024/KB8024/SignalP/input/Window_3'''
window_size = 2

X, Y = data_parser.sklearn_parser(filepath, window_size)

from sklearn.svm import SVC
from sklearn.model_selection import cross_val_score

score = cross_val_score(SVC(), X, Y, cv=10)

print(score)
print("Accuracy: %0.2f (+/- %0.2f)" % (score.mean(), score.std() * 2))


