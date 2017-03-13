# -*- coding: utf-8 -*-
"""
Demonstration of a SVM based predictor of signal and membrane domains created for the course 'PROJECT IN MOLECULAR LIFE SCIENCE (KB8024) 2017'. 

Run this script to get predictions the model. To train a new model see the 'model_builder.py' in 'SignalP/scripts/'.

@author: Revant Gupta

"""
###### Enter the following parameters ######

repo_loc = '''/home/u2196/Desktop/KB8024/KB8024/''' # Enter local repo location. See current version at https://github.com/aron0093/KB8024 for updates.
test_data_gen = input("You actually want to use this to predict stuff??? Enter True for subpar predictions, enter False to check my work on 50 proteins from TOPDB!: ") # Enter if test data is to be parsed.
use_pssm = True # Use pssm for test data or simple encoding, pssm prefered.

###### Review following parameters ######

filepath = repo_loc+'''data/TOPDB_50.txt''' # Location of test raw_data
inpath = repo_loc+'SignalP/output/model/' # Location of models
outpath = repo_loc+'SignalP/output/test_stats/' # Location of results
database = '/local_uniref/uniref/uniref50/uniref50.db' # Psiblast database

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
import dense_data_parser as ddp

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

# Testing predictor on TOPDB_50 data

# Import data as pandas dataframe
if test_data_gen == 'False' or test_data_gen == 'false':

    
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
    plt.figure(1, figsize=(20,10)).savefig(outpath+model_name+'data_plot.png')
    plt.figure(1, figsize=(20,10)).show()
    plt.figure(2)
    plot_confusion_matrix(cm, classes, normalize=True)
    plt.figure(2, figsize=(20,10)).savefig(outpath+model_name+'confusion.png')
    plt.figure(2, figsize=(20,10)).show()

elif test_data_gen == 'True' or test_data_gen == 'true':

# Generate predictions for novel proteins

    print("\n")
    
    import re
    import subprocess
    from subprocess import Popen, PIPE
    
    # Import sequence to be predicted
    
    X_test_loc = input("Enter location of target sequences in fasta format")    
    with open(X_test_loc, 'r') as f:
        len_f = len(f.readlines())
    
    collect = ''
    len_collector = []
    prot_collector = []
    lol = 0
    a = 0
    
    # Generate temp files to store sequence, generate base annotation and psiblast prep
    
    with open(X_test_loc, 'r') as f:
        
        out = open(repo_loc+'SignalP/output/temp/temp.txt', 'w')
        
        for i,lines in enumerate(f):
            if lines.startswith('>'):
                
                if lol == 1:
                    a = a + len(collect)
                    len_collector.append(a)
                    prot_collector.append(lines)
                    out.write(collect+'\n')
                    out_.write(collect+'\n')
                    out_.close()
                    out.write(('G'*len(collect))+'\n')
                    collect = ''
                out.write(lines)
                out_ = open(repo_loc+'SignalP/output/temp/prots/'+re.sub('[^A-Za-z0-9]+', '', lines[1:10])+'.txt', 'w')           
                out_.write(lines)
                lol = 1
            
            else:
                collect = collect + lines.strip()
                if i == (len_f-1):
                    a = a + len(collect)
                    len_collector.append(a)
                    prot_collector.append(lines)
                    out.write(collect+'\n')
                    out_.write(collect+'\n')
                    out_.close()
                    out.write(('G'*len(collect))+'\n')
                    collect = ''             
        out.close()
    print("File input complete")
     
    data = pdp.pre_vec_parser(repo_loc+'SignalP/output/temp/temp.txt', window_size)
        
    pssm_loc = repo_loc+'SignalP/output/temp/pssms/'
    
    print("Running psiblast")
    
    # Run psiblast as using subprocess 
        
    script = open(repo_loc+'SignalP/output/temp/prots/'+'script.sh','w')
        
    script.write('#!/bin/bash'+'\n'+'for files in '+repo_loc+'SignalP/output/temp/prots/'+'*.txt'+'\n'+'do'+'\n'+'psiblast -db '+database+' -query $files -num_threads 3 -num_iterations 4 -outfmt 10 -out_ascii_pssm '+pssm_loc+'$files.csv'+'\n'+'done')
        
    script.close()
    
    query_text = 'bash '+repo_loc+'SignalP/output/temp/prots/'+'script.sh'
    
    
    now = subprocess.Popen(query_text, stdout=PIPE, stderr=PIPE, shell=True)
    stdout, stderr = now.communicate()
  
    print("Vectorising data")
    
    pssm_data = pdp.pssm_parser(data, window_size, pssm_loc)
    
    X_test, Y_test = pdp.skl_pssm_parser(pssm_data, window_size, 'sub')    
            
    model = joblib.load(inpath+model_name)
    
    print(" Generating predictions")
    
    # Output predictions by protein
    
    Y_pred = model.predict(X_test)
    
    structure_dic = {-1:'S', 1:'M', 0:'G'}
    
    prediction_list = [structure_dic[item] for item in Y_pred]
    
    print(" Here are your results!")
    
    for i in range(0,len(prot_collector)):
        print(prot_collector[i])
        print(prediction_list[len_collector[i-1]:len_collector[i]])
      
else:

    print("Sorry, thats not an option! Try again!")






