# -*- coding: utf-8 -*-
"""
Created on Tue Feb 21 21:00:20 2017

@author: Revant
"""

import pssm_gen
import time

# pssm_ssh

outpath = '''/home/u2196/Desktop/KB8024/KB8024/SignalP/input/pssms/'''
database = '/local_uniref/uniref/uniref50/uniref50.db'
username = 'u2196'
password = 'Aim7ete9g'
num_list = ['01', '02', '03', '04', '05', '06', '07', '08', '10', '11', '13', '14', '17', '18', '19', '20', '21', '24', '26', '27', '29', '30', '32', '33', '34', '35', '36', '37']

# data_divide for pssm

filepath = '''/home/u2196/Desktop/KB8024/KB8024/data/globular_signal_tm_3state_30_slice.txt'''
pssm_split_loc = '/home/u2196/Desktop/KB8024/KB8024/data/pssm_split/'
divisions = len(num_list)

# Dividing data

pssm_gen.data_divide(filepath, pssm_split_loc, divisions)

# Generating lists for psiblast on ssh

server_list = []
file_list =[]

for i in range(divisions):
    
    server_list.extend(['stud'+num_list[i]+'.dbb.su.se'])
    file_list.extend([pssm_split_loc+'raw_data_'+str(i+1)+'.txt'])

# Purge processes

pssm_gen.process_purge(server_list, username, password)

# Start psiblast on ssh

start = time.time()

pssm_gen.pssm_gen_ssh(database, server_list, file_list, outpath, username, password)

end = time.time()

print("Time taken was %f"%(end-start))


           

