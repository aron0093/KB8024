#!/bin/bash
for files in /home/u2196/Desktop/KB8024/KB8024/data/pssm_split/raw_data_15/*.txt
do
psiblast -db /local_uniref/uniref/uniref50/uniref50.db -query $files -num_threads 3 -num_iterations 4 -outfmt 10 -out_ascii_pssm /home/u2196/Desktop/KB8024/KB8024/SignalP/input/pssms/${files#*>}_15.csv
done