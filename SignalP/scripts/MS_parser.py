# -*- coding: utf-8 -*-
"""
Created on Tue Feb 21 13:11:47 2017

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



def pssm_gen(filepath, db_address, inp_address, out_address):

       
    # Default dictionary
    
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

    aa_list = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y', 'B']

    # Code to run PSIBLAST and generate PSSM

    from Bio.Blast.Applications import NcbipsiblastCommandline
    #import subprocess
    
    raw_data = open(filepath, 'r')
    fasta_in = open(inp_address+'fasta_form.fasta', 'w')
    fasta_out = open(out_address+'pssm_out.csv', 'w')
    
    for i, lines in enumerate(raw_data.readlines()):
        if (i+1) % 3 !=0:
            fasta_in.write(lines)
        
    psi_cline = NcbipsiblastCommandline('psiblast', db = db_address, query = inp_address+'fasta_form.fasta', num_iterations = 4 , outfmt = 10, out_pssm = out_address+'pssm_out.csv', save_pssm_after_last_round=True)

    #p = subprocess.Popen(str(psi_cline),stdout=subprocess.PIPE,stderr=subprocess.PIPE,shell=True)
    #blastParser(p.stdout)
    
    str(psi_cline)
    
    psi_cline()
    
    print("PSSM stored at %spssm_out"%(out_address))
    
    return 

    # Replace default dic with pssm
    
    for vect in pssm:
        for aa in aa_list[:-1]:
            aa_dic[aa] = list(vect)
    return aa_dic    

# Generates input arrays for sklearn

def skl_pssm_parser(filepath, window_size):

    import numpy as np
    import re

    # Import data as a pandas dataframe and pivot it to wide form
    
    data = pre_vec_parser(filepath, window_size)

    # Create formatted input file

    structure_dic = {'S':-1, 'M':1, 'G':0}

    frame = (2*window_size)+1
    # Creating separate files for each peptide to avoid working with entire data     

    for i in range(0,len(data['Sequence'])):

        # Create vector dictionary from PSSM

        # Using numpy arrays instead of lists to save memory.
        X_ = np.zeros([(len(data['Sequence_windowed'][i])-2*window_size),21*frame])
        Y_ = np.zeros(len(data['Structure'][i])) 
            
        for j in range(window_size, (len(data['Sequence_windowed'][i])-window_size)):

            temp_vector = list()
            
            for a in [aa_dic[data['Sequence_windowed'][i][k]] for k in range(j-window_size, j+window_size+1)]: 
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

# Generates input arrays for sklearn

def skl_pssm_inp_gen(filepath, outpath, window_size, single_file=True):

    import numpy as np
    import re

    # Import data as a pandas dataframe and pivot it to wide form
    
    data = pre_vec_parser(filepath, window_size)

    # Create formatted input file

    structure_dic = {'S':-1, 'M':1, 'G':0}

    frame = (2*window_size)+1
    # Creating separate files for each peptide to avoid working with entire data     

    for i in range(0,len(data['Sequence'])):

        # Create vector dictionary from PSSM

        # Using numpy arrays instead of lists to save memory.
        X = np.zeros([(len(data['Sequence_windowed'][i])-2*window_size),21*frame])
        Y = np.zeros(len(data['Structure'][i])) 
            
        for j in range(window_size, (len(data['Sequence_windowed'][i])-window_size)):

            temp_vector = list()
         
            for a in [aa_dic[data['Sequence_windowed'][i][k]] for k in range(j-window_size, j+window_size+1)]: 
                temp_vector.extend(a)

            X[j-window_size] = np.array(temp_vector) 
            
        for m in range(0, len(data['Structure'][i])):
            Y[m]= float(structure_dic[data['Structure'][i][m]])

        assert len(X) == len(Y)

        if single_file==False:
            
            out1 = open(outpath+re.sub("[^a-zA-Z0-9]+", '.', data['Title'][i][1:])+'_'+str(window_size)+'_Vectors'+'.gz', 'w')
            out2 = open(outpath+re.sub("[^a-zA-Z0-9]+", '.', data['Title'][i][1:])+'_'+str(window_size)+'_Labels'+'.gz', 'w')
            out1.close
            out2.close
            np.savetxt(outpath+re.sub("[^a-zA-Z0-9]+", '.', data['Title'][i][1:])+'_'+str(window_size)+'_Vectors'+'.gz', X)
            np.savetxt(outpath+re.sub("[^a-zA-Z0-9]+", '.', data['Title'][i][1:])+'_'+str(window_size)+'_Labels'+'.gz', Y)

        elif single_file==True:
            if i==0:
                np.savetxt(open(outpath+str(window_size)+'_Vectors'+'.gz', 'w'), X)
                np.savetxt(open(outpath+str(window_size)+'_Labels'+'.gz', 'w'), Y)
                
            else:            
                np.savetxt(open(outpath+str(window_size)+'_Vectors'+'.gz', 'a'), X)
                np.savetxt(open(outpath+str(window_size)+'_Labels'+'.gz', 'a'), Y)
            
        else:
            print("Specify output type as Single or Sequence-wise file(s)")
    return
           



