# -*- coding: utf-8 -*-
"""
Created on Tue Feb 21 21:00:20 2017

@author: Revant
"""

import pandas as pd
import numpy as np
import dense_data_parser as dpp
import pssm_data_parser as pdp
import pssm_func as pf
from Bio.Blast.Applications import NcbipsiblastCommandline
import subprocess
import time
import paramiko
import os

# data_parsing + general
filepath = '''/home/u2196/Desktop/KB8024/KB8024/data/globular_signal_tm_3state.txt'''
#outpath = '''/home/u2196/Desktop/KB8024/KB8024/SignalP/input/Window_1/'''
window_size = 1
single_file = True

#cv_set_gen
K = 5

# pssm_gen
db_address = '/local_uniref/uniref/uniref50/uniref50.db'
inp_address = '/home/u2196/Desktop/KB8024/KB8024/SignalP/input/'
out_address = '/home/u2196/Desktop/KB8024/KB8024/SignalP/output/'
num_iter = 4
num_thr = 8

# data_divide for pssm
outpath = '/home/u2196/Desktop/KB8024/KB8024/data/pssm_split/'
#divisions = 28

# Pssm location
pssm_loc = '/home/u2196/Desktop/KB8024/KB8024/SignalP/input/test/'

start = time.time()


data = dpp.pre_vec_parser(filepath, window_size)

pssm_data = pdp.pssm_parser(data, window_size, pssm_loc)

#s = set()

#for seq in data['Sequence']:
#    for char in seq:
#        s.add(char)
           
#for a, b in dense_data_parser.cv_data_gen(data, 5, randomise=False):
#    print(len(a), len(b))

X, Y = pdp.skl_pssm_parser(pssm_data, window_size, pssm_type='sub')

#####for a,b,c,d in cv_set_gen.cv_set_gen(X, Y, K, randomise=False):
#####    print(len(a), len(b), len(c), len(d))

end = time.time()

print("Time taken was %f"%(end-start))


           

