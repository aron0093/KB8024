# -*- coding: utf-8 -*-
"""
Demonstration of a SVM based predictor of signal and membrane domains created for the course 'PROJECT IN MOLECULAR LIFE SCIENCE (KB8024) 2017'. 

Run this script to train the model. To obtain predictions run the 'SignalP_predict.py' file in '/bin'.

CAUTION: Do not re run script once the model is built.

@author: Revant Gupta

"""
###### Enter the following parameters ######

repo_loc = '''/home/u2196/Desktop/KB8024/KB8024/''' # Enter local repo location. See current version at https://github.com/aron0093/KB8024 for updates.
window_size = 10 # Create a frame with +/- window_size residues around the target residue. Frame size will be (1+2*window_size).
use_pssm = True # If True then check wether PSSMs are avilable in SignalP/input/pssms else generate using pssm_script.py
n_estimators = 10 # Number of bagging models

###### Review following parameters ######

filepath = repo_loc+'''data/globular_signal_tm_3state.txt''' # Location of raw_data
outpath = repo_loc+'SignalP/output/model/' # Location of output model
pssm_loc = repo_loc+'SignalP/input/pssms/'

###### Importing the required libraries and modules ######

# General Libraries

import os
import time

#Sklearn Modules

from sklearn.svm import LinearSVC
from sklearn.ensemble import RandomForestClassifier
from sklearn.externals import joblib

#Custom Modules

import pssm_data_parser as pdp
import dense_data_parser as ddp

###### Data preprocessing and storage. ######

# Import data as pandas dataframe

print("Parsing data...")
print("\n")

data = pdp.pre_vec_parser(filepath, window_size)

# Vectorise data

print("Vectorising data...")
print("\n")

if use_pssm == False:

    X, Y =  ddp.skl_parser(data, window_size)
else:
    data_pssm = pdp.pssm_parser(data, window_size, pssm_loc)
    X, Y = pdp.skl_pssm_parser(data_pssm, window_size, pssm_type='freq')
  
# Train model

model = RandomForestClassifier(class_weight='balanced', n_estimators=n_estimators)

print("Training this magnificient model...")
print("\n")

fit_start = time.time()

model.fit(X,Y)

fit_end = time.time()

print("Training completed in %f seconds"%(fit_end - fit_start))
print("\n")

model_name = 'SignalP_RandF_wind_'+str(window_size)

joblib.dump(model, outpath+model_name)

print(" Model saved at %s as %s"%(outpath, model_name))
print("\n")
print("It's gonna be great...immortal words by D.J.T.")





