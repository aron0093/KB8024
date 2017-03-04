# -*- coding: utf-8 -*-
"""
Created on Tue Feb 21 13:11:47 2017

@author: Revant Gupta
"""

# Function to divide raw data into n parts. Should I generalise? Not required for the project!

def data_divide(filepath, outpath, divisions):
    
    
    import pandas as pd
    import os

    raw_data = pd.read_csv(filepath, header=None)
    
    drop_list = []
    
    for l in range(len(raw_data)):
        if (l+1) % 3 == 0:
            drop_list.extend([l])
    
        
            
    raw_data.drop(raw_data.index[drop_list], inplace=True)
    raw_data.reset_index(drop = True, inplace = True)
    
    prot_count = len(raw_data)/2
    
    if  prot_count % divisions == 0:
    
        i = 0
        
        for k in range(divisions):
            temp = raw_data[(int(i*(len(raw_data)/divisions))): int(((i+1)*(len(raw_data)/divisions)))]
            try: 
                os.makedirs(outpath + 'raw_data_'+str(k+1)+'/')
            except OSError:
                if not os.path.isdir(outpath + 'raw_data_'+str(k+1)+'/'):
                    raise
            for s in range(0,len(temp),2):
                temp.iloc[s:s+2].to_csv(outpath+'raw_data_'+str(k+1)+'/'+str(temp.iloc[s].to_string(index=False, header=False))+'.txt', index=False, header=None)  
            i =+ 1
    
    else:
        
        i = 0
        
        residue = 2*(prot_count % divisions)

        for k in range(divisions):

            if (k+1) == divisions:
                temp = pd.concat([raw_data[(int(i*((len(raw_data)-residue)/divisions))): int(((i+1)*((len(raw_data)-residue)/divisions)))], raw_data[int(-residue):]], ignore_index =True)
            else:
                temp = raw_data[(int(i*((len(raw_data)-residue)/divisions))): int(((i+1)*((len(raw_data)-residue)/divisions)))]
            try: 
                os.makedirs(outpath + 'raw_data_'+str(k+1)+'/')
            except OSError:
                if not os.path.isdir(outpath + 'raw_data_'+str(k+1)+'/'):
                    raise
            for s in range(0,len(temp),2):
                temp.iloc[s:s+2].to_csv(outpath+'raw_data_'+str(k+1)+'/'+str(temp.iloc[s].to_string(index=False, header=False))+'.txt', index=False, header=None)  
            i += 1
            
    return
    

# Function to generate PSSM
 
def pssm_gen_ncbi(database, input_fasta, out_pssm, num_iter, num_thr):

    from Bio.Blast.Applications import NcbipsiblastCommandline
    
    psi_cline = NcbipsiblastCommandline('psiblast', db = database, query = inp_fasta, num_threads = num_thr, num_iterations = num_iter , outfmt = 5, out_ascii_pssm = out_pssm)
    psi_cline()
    
    return

# Function to execute psiblast over multiple servers
    
def pssm_gen_ssh(database, server_list, file_loc, outpath, username, password_loc):
    
    import paramiko
    import os
    import time
    import re
    
    password = (open(password_loc, 'r').read()).replace('\n','')
        
    for i in range(len(file_loc)):
    
        server = server_list[i]
        
        script = open(file_loc[i]+'script.sh','w')
        
        script.write('#!/bin/bash'+'\n'+'for files in '+file_loc[i]+'*.txt'+'\n'+'do'+'\n'+'psiblast -db '+database+' -query $files -num_threads 8 -num_iterations 4 -outfmt 10 -out_ascii_pssm '+outpath+'${files#*>}_'+str(i+1)+'.csv'+'\n'+'done')
        
        script.close()
  
        #query_text_a= 'chmod 755'+file_loc[i]+'script.sh'
        query_text= 'bash '+file_loc[i][12:]+'script.sh'
        #print(query_text)
       
        #query_text = 'psiblast -db '+database+' -query '+file_list[i]+' -num_threads 4 -num_iterations 4 -outfmt 10 -out_ascii_pssm '+outpath+'pssm_'+str(i+1)+'.csv'
    
        # Running the queries over ssh
        
        ssh = paramiko.SSHClient()
        ssh.set_missing_host_key_policy(paramiko.AutoAddPolicy())
        
        #In case the server's key is unknown, we will be adding it automatically to the list of known hosts.
        
        ssh.load_host_keys(os.path.expanduser(os.path.join("~", ".ssh", "known_hosts")))
        
        #Loads the user's local known host file.
        
        print("Connecting to server: ", server)
        
        ssh.connect(server, username=username, password=password)
          
        #ssh_stdin, ssh_stdout, ssh_stderr = ssh.exec_command(query_text_a)
        ssh_stdin, ssh_stdout, ssh_stderr = ssh.exec_command(query_text)
        ssh_stdin.close()
        
        time.sleep(1)
        #session = ssh.invoke_shell()
        #session.send("\n")

        #session.send(query_text+"\n")
        
        #print ("output", ssh_stdout.read()) #Reading output of the executed command
        #error = ssh_stderr.read()
        
        #Reading the error stream of the executed command
        #print ("err", error, len(error))

    return

# Purge processes on multiple servers

def process_purge(server_list, username, password_loc):
    
    import paramiko
    import os
    import time
    
    password = (open(password_loc, 'r').read()).replace('\n','')
    
     
    for i in range(len(server_list)):
    
        server = server_list[i]
        e_mc_2 = 'pkill -u '+username 
    
        # Running the queries over ssh
        
        ssh = paramiko.SSHClient()
        ssh.set_missing_host_key_policy(paramiko.AutoAddPolicy())
        
        #In case the server's key is unknown, we will be adding it automatically to the list of known hosts.
        
        ssh.load_host_keys(os.path.expanduser(os.path.join("~", ".ssh", "known_hosts")))
        
        #Loads the user's local known host file.
        
        print("Purging processes on: ", server)
        
        ssh.connect(server, username=username, password=password)
        ssh_stdin, ssh_stdout, ssh_stderr = ssh.exec_command(e_mc_2)
        ssh_stdin.close()
        
        time.sleep(1)
        
    return
