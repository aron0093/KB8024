# -*- coding: utf-8 -*-
"""
Created on Tue Feb 21 21:00:20 2017

@author: Revant
"""

import pssm_gen
import time

# pssm_ssh

outpath = '''/home/u2196/Desktop/KB8024/KB8024/SignalP/input/pssms/all/'''
database = '/local_uniref/uniref/uniref50/uniref50.db'
username = 'u2196'
password = '/home/u2196/Desktop/KB8024/justafile.txt'
num_list = ['01', '02', '03', '04', '05', '06', '07', '08', '10', '11', '13', '17', '18', '19', '20', '21', '24', '26', '29', '32', '34', '35', '36', '37', '38'] # 14, 27, 30, 33

# data_divide for pssm

filepath = '''/home/u2196/Desktop/KB8024/KB8024/data/globular_signal_tm_3state.txt'''
pssm_split_loc = '/home/u2196/Desktop/KB8024/KB8024/data/pssm_split/'
divisions = len(num_list)

# Dividing data

pssm_gen.data_divide(filepath, pssm_split_loc, divisions)

# Generating lists for psiblast on ssh

server_list = []
file_loc =[]

for i in range(divisions):
    
    server_list.extend(['stud'+num_list[i]+'.dbb.su.se'])
    file_loc.extend([pssm_split_loc+'raw_data_'+str(i+1)+'/'])

assert divisions == len(server_list) == len(file_loc)
# Purge processes

pssm_gen.process_purge(server_list, username, password)

# Start psiblast on ssh

pssm_gen.pssm_gen_ssh(database, server_list, file_loc, outpath, username, password)


print("Psiblast is running on your servers, Good luck!")


           

