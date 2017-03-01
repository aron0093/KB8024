# -*- coding: utf-8 -*-
"""
Created on Tue Feb 21 21:00:20 2017

@author: Revant
"""
import dense_data_parser
#import cv_set_gen
import pssm_gen
from Bio.Blast.Applications import NcbipsiblastCommandline
import subprocess
import time
import paramiko
import os

# data_parsing + general
filepath = '''/home/u2196/Desktop/KB8024/KB8024/data/globular_signal_tm_3state.txt'''
outpath = '''/home/u2196/Desktop/KB8024/KB8024/SignalP/input/Window_1/'''
window_size = 1
single_file = True

#cv_set_gen
K = 2

# pssm_gen
db_address = '/local_uniref/uniref/uniref50/uniref50.db'
inp_address = '/home/u2196/Desktop/KB8024/KB8024/SignalP/input/'
out_address = '/home/u2196/Desktop/KB8024/KB8024/SignalP/output/'
num_iter = 4
num_thr = 8

# data_divide for pssm
pssm_split_loc = '/home/u2196/Desktop/KB8024/KB8024/data/pssm_split/'
divisions = 20

start = time.time()



#pssm_gen.data_divide(filepath, pssm_split_loc, divisions)

#p = subprocess.Popen(str(psi_cline),stdout=subprocess.PIPE,stderr=subprocess.PIPE, shell=True)
#def bb(db_address, inp_address, out_address, num_iter, num_thr):
#    psi_cline = NcbipsiblastCommandline('psiblast', db = db_address, query = inp_address+'fasta_form.fasta', num_threads = num_thr, num_iterations = num_iter , outfmt = 5, out_ascii_pssm = out_address+'pssm_out.xml')
#    psi_cline()
#    return
    
#bb(db_address, inp_address, out_address, num_iter, num_thr)
#X, Y = SS_parser.skl_parser(filepath, window_size)

#for a,b,c,d in cv_set_gen.cv_set_gen(X, Y, K, randomise=False):
#    print(len(a), len(b), len(c), len(d))

end = time.time()

print("Time taken was %f"%(end-start))


           

