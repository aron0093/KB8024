# -*- coding: utf-8 -*-
"""
Demonstration of a SVM based predictor of signal and membrane domains created for the course 'PROJECT IN MOLECULAR LIFE SCIENCE (KB8024) 2017'. 

Run this script to train the a Random Forest for comparison. To obtain predictions run the 'SignalP_predict.py' file in '/bin'.

CAUTION: Do not re run script once the model is built.

@author: Revant Gupta

"""
###### Enter the following parameters ######

repo_loc = '''/home/u2196/Desktop/KB8024/KB8024/''' # Enter local repo location. See current version at https://github.com/aron0093/KB8024 for updates.
window_size = 25 # Create a frame with +/- window_size residues around the target residue. Frame size will be (1+2*window_size).
use_pssm = True # If True then check wether PSSMs are avilable in SignalP/input/pssms else generate using pssm_script.py
n_estimators = 128 # Number of Trees

###### Review following parameters ######

filepath = repo_loc+'''data/globular_signal_tm_3state_30_slice.txt''' # Location of raw_data
outpath = repo_loc+'SignalP/output/model/' # Location of output model
pssm_loc = repo_loc+'SignalP/input/pssms/'

###### Importing the required libraries and modules ######

# General Libraries

import time
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt
from collections import OrderedDict as oD

#Sklearn Modules

from sklearn.preprocessing import normalize
from sklearn.ensemble import RandomForestClassifier
from sklearn.externals import joblib
from sklearn.model_selection import cross_val_predict
from sklearn.metrics import precision_recall_fscore_support as score
from sklearn.metrics import confusion_matrix

#Custom Modules

import pssm_data_parser as pdp
import dense_data_parser as ddp

# Function to format confusion_matrix

def plot_confusion_matrix(cm, classes, normalize=False, title='Confusion matrix', cmap=plt.cm.Blues):

    from matplotlib import pyplot as plt
    import itertools
    
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

model = RandomForestClassifier(n_estimators = n_estimators, class_weight='balanced')

print("Training this magnificient model...")
print("\n")

fit_start = time.time()

scores = oD()
scores['labels'] = np.array([-1,0,1])


Y_pred = cross_val_predict(model, X, Y)

precision, recall, fscore, support = score(Y, Y_pred)

cm = confusion_matrix(Y, Y_pred)
    
classes = ['Signal', 'Globular', 'TransMembrane']

model.fit(X,Y)

scores['precision'] = precision
scores['recall'] = recall
scores['fscore'] = fscore
scores['norm_support'] = normalize(support).reshape(3,)
scores_table = pd.DataFrame.from_dict(scores, orient='columns')
scores_table.to_csv(outpath+'model_stats/'+'SignalP_Rand_forest_wind_final_bagging'+str(window_size)+'data.csv')    
scores_table.plot(x='labels', y = ['precision', 'recall', 'fscore', 'norm_support'], kind='bar', colormap='Pastel1')
plt.figure(1)
plt.xlabel("Classes -1:Signal 0:Globular, 1:TransMembrane")
plt.ylabel("Scores")
plt.title("Are you not entertained...")

plt.figure(1).savefig(outpath+'model_stats/'+'SignalP_Rand_forest_wind_final_bagging'+str(window_size)+'data_plot.png')
  
plt.figure(2)
plot_confusion_matrix(cm, classes, normalize=True)
plt.figure(2).savefig(outpath+'model_stats/'+'SignalP_Rand_forest_wind_final_bagging'+str(window_size)+'confusion_plot.png')


fit_end = time.time()

print("Training completed in %f seconds"%(fit_end - fit_start))
print("\n")

model_name = 'SignalP_Rand_forest_wind_final_bagging'+str(window_size)

joblib.dump(model, outpath+model_name)

print(" Model saved at %s as %s"%(outpath, model_name))
print("\n")
print("It's gonna be great...immortal words by D.J.T.")





