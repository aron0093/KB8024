# -*- coding: utf-8 -*-
"""
Created on Tue Feb 27 13:11:47 2017

@author: Revant Gupta
"""

#Input raw data and window_size and store pre vectorised data as pandas Dataframe

def pre_vec_parser(filepath, window_size):
    
    # Import data as a pandas dataframe and pivot it to wide form

    import pandas as pd

    raw_data = pd.read_csv(filepath, header=None)

    headers = ['Title', 'Sequence', 'Structure']
    new_index = []
    for i in range(0,(int(len(raw_data[0])/3))):
        new_index.extend([i,i,i])

    raw_data[1] = pd.Series(headers*(int(len(raw_data[0])/3)))
    raw_data[2] = pd.Series(new_index)

    data = pd.DataFrame()
    data = raw_data.pivot(index = 2, columns= 1, values = 0)

    pre_suf = ['B']*window_size

    data['Sequence_windowed'] = ''

    # Use ''.join for efficiency instead of +. Also use nested loops instead of storing variables.

    for i in range(0,len(data['Sequence'])):
        data['Sequence_windowed'][i] = ''.join([''.join(pre_suf), data['Sequence'][i], ''.join(pre_suf)])

    return data

# Function to parse pssms

def pssm_parser(data, window_size, pssm_loc):
    
    import pandas as pd
    import numpy as np
    import math
    import os

    # Sigmoid function to normalise substitution matrix
    def sigmoid(x):
        return (1 / (1 + math.exp(-x)))

    data['PSSM_sub']=''
    data['PSSM_freq']=''
    
    # Selecting PSSMs

    for fil in os.listdir(pssm_loc):
        for index, prot in data['Title'].iteritems():
            if prot[1:] == fil.partition('.')[0]:
                
                with open(pssm_loc+fil, 'r') as f:
                     raw_file = pd.read_csv(f)
                     raw_pssm = raw_file[raw_file.columns[0]]
                     
                     # Expected vector shape for numpy arrays
                     
                     sub = np.zeros(((len(raw_pssm)-6+(2*window_size)),20), dtype=float)
                     freq = np.zeros(((len(raw_pssm)-6+(2*window_size)),20), dtype=float)
                     
                     # Storing substituion and frequency matrices in separate cells for each sequence
                     
                     for i in range(1,(len(raw_pssm)-5)):
                        sub[i-1+(2*window_size)] = np.array([sigmoid(float(x)) for x in raw_pssm[i].split()[2:22]], dtype=float)
                        freq[i-1+(2*window_size)] = (np.array((raw_pssm[i].split()[22:42]), dtype=float))/(100.0)
                        
                     data['PSSM_sub'][index]= sub 
                     data['PSSM_freq'][index] = freq  
    return data


def cv_data_gen(data, K, randomise=False):

    # Protein level division of data for cross validation

    import pandas as pd
    import numpy as np
    from sklearn.utils import shuffle

    if randomise:
        data = shuffle(data)
        
        
    splits = np.array_split(data, K)
                         
    for k in range(K):
        train_sets = [item for i,item in enumerate(splits) if i!=k ]
        train_data = pd.concat(train_sets, ignore_index=True)
        test_data = pd.DataFrame(splits[k]).reset_index(drop =True)    
        
        yield train_data, test_data



# Generates input arrays for sklearn using pssm.
    
def skl_pssm_parser(data, window_size, pssm_type='freq'):      
        
    import numpy as np
    import pandas as pd
    import re
    
    # PSSM type
    
    if pssm_type == 'sub':
        pssm = 'PSSM_sub'
    elif pssm_type == 'freq':
        pssm = 'PSSM_freq'
    else:
        raise ValueError("Invalid pssm type. Expected one of: %s" % (['sub', 'freq']))

    # Create formatted input file

    structure_dic = {'S':-1, 'M':1, 'G':0}

    frame = (2*window_size)+1
    # Creating separate files for each peptide to avoid working with entire data     

    for i in range(0,len(data['Sequence'])):

        ######## Create vector dictionary from PSSM

        # Using numpy arrays instead of lists to save memory.
        X_ = np.zeros([(len(data['Sequence_windowed'][i])-2*window_size),20
        *frame])
        Y_ = np.zeros(len(data['Structure'][i])) 
            
        for j in range(window_size, (len(data['Sequence_windowed'][i])-window_size)):

            temp_vector = list()
            
            for a in [data[pssm][i][k] for k in range(j-window_size, j+window_size+1)]: 
                temp_vector.extend(a)

            X_[j-window_size] = np.array(temp_vector) 
            
        for m in range(0, len(data['Structure'][i])):
            Y_[m]= float(structure_dic[data['Structure'][i][m]])

        assert len(X_) == len(Y_)

        if i==0:
            X = X_
            Y = Y_
        else:
            X = np.concatenate((X,X_), axis=0)
            Y = np.concatenate((Y,Y_), axis=0)
     
    return X, Y


