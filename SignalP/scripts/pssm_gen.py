# -*- coding: utf-8 -*-
"""
Created on Tue Feb 21 13:11:47 2017

@author: Revant Gupta
"""

# Function to generate PSSM
 
def pssm_gen(database, input_fasta, out_pssm, num_iter, num_thr):

    from Bio.Blast.Applications import NcbipsiblastCommandline
    
    psi_cline = NcbipsiblastCommandline('psiblast', db = database, query = inp_fasta, num_threads = num_thr, num_iterations = num_iter , outfmt = 5, out_ascii_pssm = out_pssm)
    psi_cline()
    return out
 
def pssm_add(filepath, db_address, inp_address, out_address, num_iter, num_thr):
       
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

    raw_data = open(filepath, 'r')
    fasta_in = open(inp_address+'fasta_form.fasta', 'w')
       
    for i, lines in enumerate(raw_data.readlines()):
        if (i+1) % 3 !=0:
            fasta_in.write(lines)
            
    from Bio.Blast import NCBIXML
    
        
    pssm_gen(db_address, inp_address, out_address, num_iter, num_thr)
       
    print("PSSM will be stored at %spssm_out"%(out_address)) 
            
    return




