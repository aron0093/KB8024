# -*- coding: utf-8 -*-
"""
Created on Tue Feb 21 21:00:20 2017

@author: Revant
"""
import SS_parser
import MS_parser
import cv_set_gen

filepath = '''C:\Users\Goodbaccha\Desktop\Masters\Semesters\Spring 17\KB8024\KB8024\data\globular_signal_tm_3state_slice.txt'''
outpath = '''C:\Users\Goodbaccha\Desktop\Masters\Semesters\Spring 17\KB8024\KB8024\SignalP\input\Window_1/'''
window_size = 1
single_file = True
K = 2
db_address = ''
inp_address = ''
out_address = ''

MS_parser.pssm_gen(filepath, db_address, inp_address, out_address)

#X, Y = SS_parser.skl_parser(filepath, window_size)

#for a,b,c,d in cv_set_gen.cv_set_gen(X, Y, K, randomise=False):
#    print(len(a), len(b), len(c), len(d))




           

