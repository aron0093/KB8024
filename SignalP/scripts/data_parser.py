# -*- coding: utf-8 -*-
"""
Created on Tue Feb 21 13:11:47 2017

@author: Revant Gupta
"""

#Stores data as svmligh sparse format with window length 3 at designated location

def svmL_inp_gen_len3(filepath, outpath):

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

#Stores data for each sequence in separate files at a designated location

def skl_inp_gen(filepath, outpath, window_size, single_file=True):
    
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

    # Use ''.join for efficiency instead of +. Also use nested loops instead of storing variables.

    for i in range(0,len(data['Sequence'])):
        data['Sequence'][i] = ''.join([''.join(pre_suf), data['Sequence'][i], ''.join(pre_suf)])

    # Create formatted input file

    # Create dictionaries
    aa_dic = {  'A':[1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
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

    frame = (2*window_size)+1
    # Creating separate files for each peptide to avoid working with entire data     

    for i in range(0,len(data['Sequence'])):

        # Using numpy arrays instead of lists to save memory.
        X = np.zeros([(len(data['Sequence'][i])-2*window_size),21*frame])
        Y = np.zeros(len(data['Structure'][i])) 
            
        for j in range(window_size, (len(data['Sequence'][i])-window_size)):

            temp_vector = list()
            
            for a in [aa_dic[data['Sequence'][i][k]] for k in range(j-window_size, j+window_size+1)]: 
                temp_vector.extend(a)

            X[j-window_size] = np.array(temp_vector) 
            
        for m in range(0, len(data['Structure'][i])):
            Y[m]= float(structure_dic[data['Structure'][i][m]])

        if single_file==False:
            
            out1 = open(outpath+'/'+re.sub("[^a-zA-Z0-9]+", '.', data['Title'][i][1:])+'_'+str(window_size)+'_Vectors'+'.gz', 'w')
            out2 = open(outpath+'/'+re.sub("[^a-zA-Z0-9]+", '.', data['Title'][i][1:])+'_'+str(window_size)+'_Labels'+'.gz', 'w')
            out1.close
            out2.close
            np.savetxt(outpath+'/'+re.sub("[^a-zA-Z0-9]+", '.', data['Title'][i][1:])+'_'+str(window_size)+'_Vectors'+'.gz', X)
            np.savetxt(outpath+'/'+re.sub("[^a-zA-Z0-9]+", '.', data['Title'][i][1:])+'_'+str(window_size)+'_Labels'+'.gz', Y)

        elif single_file==True:
            if i==0:
                out1 = open(outpath+'/'+str(window_size)+'_Vectors'+'.gz', 'w')
                out2 = open(outpath+'/'+str(window_size)+'_Labels'+'.gz', 'w')
                out1.close
                out2.close
                np.savetxt(open(outpath+'/'+str(window_size)+'_Vectors'+'.gz', 'a'), X)
                np.savetxt(open(outpath+'/'+str(window_size)+'_Labels'+'.gz', 'a'), Y)
                
            else:            
                np.savetxt(open(outpath+'/'+str(window_size)+'_Vectors'+'.gz', 'a'), X)
                np.savetxt(open(outpath+'/'+str(window_size)+'_Labels'+'.gz', 'a'), Y)
            
        else:
            print("Specify output type as Single or Sequence-wise file(s)")
    return
           

# Generates input arrays for sklearn

def skl_parser(filepath, window_size):
    
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

    # Use ''.join for efficiency instead of +. Also use nested loops instead of storing variables.

    for i in range(0,len(data['Sequence'])):
        data['Sequence'][i] = ''.join([''.join(pre_suf), data['Sequence'][i], ''.join(pre_suf)])

    # Create formatted input file

    # Create dictionaries
    aa_dic = {  'A':[1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
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

    frame = (2*window_size)+1
    # Creating separate files for each peptide to avoid working with entire data     

    for i in range(0,len(data['Sequence'])):

        # Using numpy arrays instead of lists to save memory.
        X_ = np.zeros([(len(data['Sequence'][i])-2*window_size),21*frame])
        Y_ = np.zeros(len(data['Structure'][i])) 
            
        for j in range(window_size, (len(data['Sequence'][i])-window_size)):

            temp_vector = list()
            
            for a in [aa_dic[data['Sequence'][i][k]] for k in range(j-window_size, j+window_size+1)]: 
                temp_vector.extend(a)

            X_[j-window_size] = np.array(temp_vector) 
            
        for m in range(0, len(data['Structure'][i])):
            Y_[m]= float(structure_dic[data['Structure'][i][m]])

        if i==0:
            X = X_
            Y = Y_
        else:
            X = np.concatenate((X,X_), axis=0)
            Y = np.concatenate((Y,Y_), axis=0)
     
    return X, Y
