# -*- coding: utf-8 -*-
"""
Created on Tue Feb 21 13:11:47 2017

@author: Revant Gupta
"""
def svmlight_len3_parser(filepath, outpath):

    # Import data as a pandas dataframe and pivot it to wide form
    
    import pandas as pd
    import re
    
    raw_data = pd.read_csv(filepath, header=None)
    
    headers = ['Title', 'Sequence', 'Structure']
    new_index = []
    for i in range(0,(int(len(raw_data[0])/3))):
        new_index.extend([i,i,i])
    
    raw_data[1] = pd.Series(headers*(int(len(raw_data[0])/3)))
    raw_data[2] = pd.Series(new_index)
    
    data = pd.DataFrame()
    data = raw_data.pivot(index = 2, columns= 1, values = 0)
    
    # Create formatted input file
    aa_dic = {'A':1, 'C':2, 'D':3, 'E':4, 'F':5, 'G':6, 'H':7, 'I':8, 'K':9, 'L':10, 'M':11, 'N':12, 
              'P':13, 'Q':14, 'R':15, 'S':16, 'T':17, 'U':18, 'V':19, 'W':20, 'Y':21 }
    structure_dic = {'S':-1, 'M':1, 'G':0}
    
    for i in range(0,len(data['Sequence'])):
        out = open(outpath+'/'+re.sub("[^a-zA-Z0-9]+", '.', data['Title'][i][1:])+'.txt', 'w')
        out.close
        out = open(outpath+'/'+re.sub("[^a-zA-Z0-9]+", '.', data['Title'][i][1:])+'.txt', 'a')
        out.write('#'+str(aa_dic)[1:-1])
        for j in range(0,(len(data['Sequence'][i])-2)):
            qid = str(structure_dic[data['Structure'][i][j+1]])
            aa_1 = str(aa_dic[data['Sequence'][i][j]])
            aa_2 = str(aa_dic[data['Sequence'][i][j+1]])
            aa_3 = str(aa_dic[data['Sequence'][i][j+2]])
            out.write(qid+' '+aa_1+':1 '+aa_2+':1 '+aa_3+':1')
        out.close
    return

def sklearn_parser(filepath, outpath, window_size):
    
# Import data as a pandas dataframe and pivot it to wide form
    
    import pandas as pd
    import numpy as np
    import re
    
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
    
    for i in range(0,len(data['Sequence'])):
        data['Sequence'][i] = ''.join([''.join(pre_suf), data['Sequence'][i], ''.join(pre_suf)])
    
# Create formatted input file
    
# Create dictionaries
    aa_dic = {	'A':[1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
		'C':[0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0], 
		'D':[0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
		'E':[0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
		'F':[0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
		'G':[0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0], 
		'H':[0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0], 
		'I':[0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0],
		'K':[0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0],
		'L':[0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0],
		'M':[0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0],
		'N':[0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0],
		'P':[0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0],
		'Q':[0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0],
		'R':[0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0],
		'S':[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0],
		'T':[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0],
		'V':[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0],
		'W':[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0],
		'Y':[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0],
		'B':[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0] }
    
    structure_dic = {'S':-1, 'M':1, 'G':0}
    
# Creating separate files for each peptide to avoid working with entire data     
    for i in range(0,len(data['Sequence'])):
        out = open(outpath+'/'+re.sub("[^a-zA-Z0-9]+", '.', data['Title'][i][1:])+'_'+str(window_size)+'.txt', 'w')
        out.close
        out = open(outpath+'/'+re.sub("[^a-zA-Z0-9]+", '.', data['Title'][i][1:])+'_'+str(window_size)+'.txt', 'a')
    
        X=np.empty([len(data['Sequence'][i]),20*((2*window_size)+1)])
        Y=np.empty([len(data['Sequence'][i])])
    
        for j in (0, len(data['Sequence'][i])):
            X[j][:20*((2*window_size)+1)] = temp.extend([aa_dic[data['Sequence'][i][k]] for k in range(-window_size, window_size)
    return X
