# -*- coding: utf-8 -*-
"""
Demonstration of a SVM based predictor of signal and membrane domains created for the course 'PROJECT IN MOLECULAR LIFE SCIENCE (KB8024) 2017'. 

Run this script to get predictions the model. To train a new model see the 'model_builder.py' in 'SignalP/scripts/'.

@author: Revant Gupta

"""
###### Enter the following parameters ######

repo_loc = '''/home/u2196/Desktop/KB8024/KB8024/''' # Enter local repo location. See current version at https://github.com/aron0093/KB8024 for updates.
test_data_gen = bool(input("You actually want to use this to predict stuff??? Enter True for madness, enter False to check my work!: ")) # Enter if test data is to be parsed.
use_pssm = False
###### Review following parameters ######

filepath = repo_loc+'''data/globular_signal_tm_3state_100_slice.txt''' # Location of test raw_data
inpath = repo_loc+'SignalP/output/model/' # Location of models
database = '/local_uniref/uniref/uniref50/uniref50.db' 
###### Importing the required libraries and modules ######

# General Libraries

import os
import time
import pandas as pd
import numpy as np
from collections import OrderedDict as oD
from matplotlib import pyplot as plt

#Sklearn Modules

from sklearn.externals import joblib
from sklearn.metrics import precision_recall_fscore_support as score
from sklearn.preprocessing import normalize
#Custom Modules

import pssm_data_parser as pdp
import dense_data_parser as ddp

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
        data_pssm = pdp.pssm_parser(data, pssm_loc)
        X_test, Y_test = pdp.skl_pssm_parser(data_pssm, window_size)
        
    # Test model
  
    scores = oD()
    scores['labels'] = np.array([-1,0,1])

    model = joblib.load(inpath+model_name)

    print("Testing this thing...")
    print("\n")

    Y_pred = model.predict(X_test)

    precision, recall, fscore, support = score(Y_test, Y_pred)

    scores['precision'] = precision
    scores['recall'] = recall
    scores['fscore'] = fscore
    scores['norm_support'] = normalize(support).reshape(3,)

    scores_table = pd.DataFrame.from_dict(scores, orient='columns')
    scores_table.plot(x='labels', y = ['precision', 'recall', 'fscore', 'norm_support'], kind='bar', colormap='Pastel1')

    plt.xlabel("Classes -1:Signal 0:Globular, 1:Transmembrane")
    plt.ylabel("Scores")
    plt.title("Are you not entertained...")

    print("Test results...")
    print("\n")
    print(scores_table)
    print("\n")
    plt.show()

else:

    print("\n")
    
    X_test_loc = input("Enter location of target sequence in fasta format")
    
    target_seq = pd.read_csv(X_test_loc)
    
    from Bio.Blast.Applications import NcbipsiblastCommandline
    
    out_pssm = repo_loc+'/SignalP/output/temp.csv'
    psi_cline = NcbipsiblastCommandline('psiblast', db = database, query = inp_fasta, num_threads = num_thr, num_iterations = num_iter , outfmt = 10, out_ascii_pssm = out_pssm)
    psi_cline()
    
    #Code to read PSSM and encode protein
    
    #X_test =
    
    Y_pred = model.predict(X_test)
    
    print(Y_pred)
      







