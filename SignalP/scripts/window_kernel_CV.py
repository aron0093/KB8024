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

#Custom modules

import dense_data_parser as ddp

# Set parameters

window_size = [2,5,8,10,13]

kernel = str(input("Input kernel type as linear, rbf, poly or sigmoid: "))

filepath = '''/home/u2196/Desktop/KB8024/KB8024/data/globular_signal_tm_3state_30_slice.txt'''
output = '''/home/u2196/Desktop/KB8024/KB8024/SignalP/output/window_kernel/'''
cv_sets = 5
structure_dic = {'S':-1, 'M':1, 'G':0}

# Starting script

start = time.time()

final_list = []

for windows in window_size:
 
    data = ddp.pre_vec_parser(filepath, windows)

    clf = OneVsRestClassifier(SVC( kernel = kernel, class_weight = 'balanced'), n_jobs=-1)

    scores = oD()
    scores['labels'] = np.array([-1,0,1])
    p = np.zeros(3)
    r = np.zeros(3)
    f = np.zeros(3)
    s = np.zeros(3)
    ft = np.zeros(3)

    for train_data, test_data in ddp.cv_data_gen(data, cv_sets, randomise=False):

        X_train, Y_train = ddp.skl_parser(train_data, windows)
        X_test, Y_test = ddp.skl_parser(test_data, windows)
        
        fit_start = time.time()
        clf.fit(X_train,Y_train)
        fit_end = time.time()
        
        fit_time = fit_end - fit_start
        
        predicted = clf.predict(X_test)
        
        precision, recall, fscore, support = score(Y_test, predicted)
        
        p = p + precision
        r = r + recall
        f = f + fscore
        s = s + support
        ft = ft + np.array([fit_time])
        
    scores['window']= windows
    scores['precision'] = normalize(p/cv_sets).reshape(3,)
    scores['recall'] = normalize(r/cv_sets).reshape(3,)
    scores['fscore'] = normalize(f/cv_sets).reshape(3,)
    scores['support'] = normalize(s/cv_sets).reshape(3,)
    scores['fit_time'] = normalize(ft/cv_sets).reshape(3,)        
          
    print("Done window %s..."%(str(windows)))
     
    scores_table = pd.DataFrame.from_dict(scores, orient='columns')
    
    final_list.append(scores_table)

final_table = pd.concat(final_list, ignore_index = True)
final_table.to_csv(output+str(kernel)+'.csv')

for clas in set(final_table['labels']):
    temp = final_table.loc[final_table['labels']==clas]
    temp.plot(x='window', y = ['precision', 'recall', 'fscore', 'support', 'fit_time'], kind='bar', colormap='Pastel1')
    plt.xlabel("+/- frames around target residue")
    plt.ylabel("Normalised Scores")
    plt.title(str(cv_sets)+' CV score for windows for class '+structure_dic[clas]+' using kernel '+kernel)
    plt.savefig(output+str(kernel)+structure_dic[clas]+'.png')         

end = time.time()

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


