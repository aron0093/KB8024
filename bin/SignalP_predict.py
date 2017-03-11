# -*- coding: utf-8 -*-
"""
Demonstration of a SVM based predictor of signal and membrane domains created for the course 'PROJECT IN MOLECULAR LIFE SCIENCE (KB8024) 2017'. 

Run this script to get predictions the model. To train a new model see the 'model_builder.py' in 'SignalP/scripts/'.

@author: Revant Gupta

"""
###### Enter the following parameters ######

repo_loc = '''/home/u2196/Desktop/KB8024/KB8024/''' # Enter local repo location. See current version at https://github.com/aron0093/KB8024 for updates.
test_data_gen = bool(input("You actually want to use this to predict stuff??? Enter True for madness, enter False to check my work!: ")) # Enter if test data is to be parsed.
use_pssm = True
###### Review following parameters ######

filepath = repo_loc+'''data/TOPDB_50.txt''' # Location of test raw_data
inpath = repo_loc+'SignalP/output/model/' # Location of models
outpath = repo_loc+'SignalP/output/test_stats/' # Location of results
database = '/local_uniref/uniref/uniref50/uniref50.db' 
###### Importing the required libraries and modules ######

# General Libraries

import os
import time
import pandas as pd
import numpy as np
from collections import OrderedDict as oD
from matplotlib import pyplot as plt
import itertools

#Sklearn Modules

from sklearn.externals import joblib
from sklearn.metrics import precision_recall_fscore_support as score
from sklearn.metrics import classification_report
from sklearn.metrics import confusion_matrix
from sklearn.preprocessing import normalize
#Custom Modules

import pssm_data_parser as pdp
import dense_data_parser_basicGrade as ddp

# Custom functions

def plot_confusion_matrix(cm, classes, normalize=False, title='Confusion matrix', cmap=plt.cm.Blues):

    from matplotlib import pyplot as plt
    
    plt.imshow(cm, interpolation='nearest', cmap=cmap)
    plt.title(title)
    plt.colorbar()
    tick_marks = np.arange(len(classes))
    plt.xticks(tick_marks, classes, rotation=45)
    plt.yticks(tick_marks, classes)

    if normalize:
        cm = cm.astype('float') / cm.sum(axis=1)[:, np.newaxis]
        print("Normalized confusion matrix")
    else:
        print('Confusion matrix, without normalization')

    print(cm)

    thresh = cm.max() / 2.
    for i, j in itertools.product(range(cm.shape[0]), range(cm.shape[1])):
        plt.text(j, i, round(cm[i, j]*100, 2),
                 horizontalalignment="center",
                 color="green" if cm[i, j] > thresh else "black")

    plt.tight_layout()
    plt.ylabel('True label')
    plt.xlabel('Predicted label')
    
    return

###### Model selection ######
print("\n")
print(str(os.listdir(inpath)))
print("\n")

model_name = str(input("Enter Model Name from the list above: ")) # Enter the desired model file
print("\n")
window_size = int(input("Enter Window Size associated with model: ")) # Fixed and unchangeable depending on model

###### Data preprocessing and storage. ######

# Import data as pandas dataframe
if test_data_gen == True:

    
    print("Parsing data...")
    print("\n")

    data = pdp.pre_vec_parser(filepath, window_size)

    # Vectorise data

    print("Vectorising data...")
    print("\n")

    if use_pssm == False:
    
        X_test, Y_test =  ddp.skl_parser(data, window_size)
    else:
        
        pssm_loc = repo_loc+'SignalP/input/test_pssms/'
        data_pssm = pdp.pssm_parser(data, window_size, pssm_loc)
        
        X_test, Y_test = pdp.skl_pssm_parser(data, window_size, pssm_type='freq')
        
    # Test model
  
    scores = oD()
    scores['labels'] = np.array([-1,0,1])

    model = joblib.load(inpath+model_name)

    print("Testing this thing...")
    print("\n")

    Y_pred = model.predict(X_test)

    precision, recall, fscore, support = score(Y_test, Y_pred)
    
    cm = confusion_matrix(Y_test, Y_pred)
    
    classes = ['Signal', 'Globular', 'TransMembrane']
    
    scores['precision'] = precision
    scores['recall'] = recall
    scores['fscore'] = fscore
    scores['norm_support'] = normalize(support).reshape(3,)

    scores_table = pd.DataFrame.from_dict(scores, orient='columns')
    scores_table.to_csv(outpath+model_name+'.csv')
    scores_table.plot(x='labels', y = ['precision', 'recall', 'fscore', 'norm_support'], kind='bar', colormap='Pastel1')
    plt.figure(1)
    plt.xlabel("Classes -1:Signal 0:Globular, 1:TransMembrane")
    plt.ylabel("Scores")
    plt.title("Are you not entertained...")

    print("Test results...")
    print("\n")
    print(classification_report(Y_test, Y_pred, target_names=['Signal', 'Globular', 'Transmembrane']))
    print("\n")
    plt.figure(1).savefig(outpath+model_name+'data_plot.png')
    plt.figure(1).show()
    plt.figure(2)
    plot_confusion_matrix(cm, classes, normalize=True)
    plt.figure(2).savefig(outpath+model_name+'confusion.png')
    plt.figure(2).show()
else:

    print("\n")
    
    X_test_loc = input("Enter location of target sequence in fasta format")
    
    with open(X_test_loc, 'r') as f:
        out = open(repo_loc+'SignalP/output/temp/temp.txt', 'w')
        for lines in f:
            out.write(lines)
            a = len(lines)
        out.write('G'*a)
        out.close()
     
     
    data = pdp.pre_vec_parser(repo_loc+'SignalP/output/temp/temp.txt', window_size)
    
    have_pssm = bool(input("Do you have psiblast pssm in csv format?"))
    
    if have_pssm==True:
        X_pssm_loc = input("Enter location of psiblast pssm in csv format")
        
        X_test, Y_test = pdp.skl_pssm_parser(data, window_size, 'freq')
    else:
    
        X_test, Y_test = ddp.skl_parser(data, window_size)
            
    Y_pred = model.predict(X_test)
    
    structure_dic = {-1:'S', 1:'M', 0:'G'}
    
    prediction_list = [structure_dic[item] for item in Y_pred]
    
    predcition = str(prediction_list)
    
    print(prediction)
      







