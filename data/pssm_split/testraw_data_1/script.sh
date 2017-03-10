#!/bin/bash
for files in /home/u2196/Desktop/KB8024/KB8024/data/pssm_split/testraw_data_1/*.txt
do
psiblast -db /local_uniref/uniref/uniref50/uniref50.db -query $files -num_threads 3 -num_iterations 4 -outfmt 10 -out_ascii_pssm /home/u2196/Desktop/KB8024/KB8024/SignalP/input/test_pssms/${files#*>}.csv
done
