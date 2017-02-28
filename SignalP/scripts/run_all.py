# -*- coding: utf-8 -*-
"""
SIGMEMB

A SVM based predictor of signal and membrane domains created for the course 'PROJECT IN MOLECULAR LIFE SCIENCE (KB8024) 2017'. 

@author: Revant Gupta

"""
###### Enter the following parameters ######

repo_loc = '''/home/u2196/Desktop/KB8024/''' # Enter local repo location. See current version at https://github.com/aron0093/KB8024 for updates.
window_size = 1






###### Review following parameters ######

filepath = repo_loc+'''/KB8024/data/globular_signal_tm_3state.txt''' # 
outpath = repo_loc+'''/KB8024/KB8024/SignalP/input/Window_'''+str(window_size+'/'



###### Importing the required libraries and modules ######

import os

###### Data preprocessing and storage. ######

if pre_process == True:

    # Creating folder to store parsed data as sequence-wise files.
    
    try: 
        os.makedirs(path)
    except OSError:
        if not os.path.isdir(path):
            raise







