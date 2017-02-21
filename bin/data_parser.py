# -*- coding: utf-8 -*-
"""
Created on Tue Feb 21 13:11:47 2017

@author: Revant Gupta
"""
def signalP_parser(filepath, outpath, window_size):

    # Import data as a pandas dataframe and pivot it to wide form
    
    import pandas as pd
    import re
    
    raw_data = pd.read_csv(filepath, header=None)
    
    headers = ['Title', 'Sequence', 'Structure']
    new_index = []
    for i in range(0,(len(raw_data[0])/3)):
        new_index.extend([i,i,i])
    
    raw_data[1] = pd.Series(headers*(len(raw_data[0])/3))
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