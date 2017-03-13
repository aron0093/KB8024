#!/bin/bash
for files in /home/u2196/Desktop/KB8024/KB8024/SignalP/output/temp/prots/*.txt
do
psiblast -db /local_uniref/uniref/uniref50/uniref50.db -query $files -num_threads 3 -num_iterations 4 -outfmt 10 -out_ascii_pssm /home/u2196/Desktop/KB8024/KB8024/SignalP/output/temp/pssms/$files.csv
done