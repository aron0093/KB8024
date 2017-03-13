#!/bin/bash
for files in /home/u2196/Desktop/KB8024/KB8024/data/pssm_split/dimitri/*.fasta
do
psiblast -db /local_uniref/uniref/uniref90/uniref90.db -query $files -num_threads 8 -num_iterations 3 -evalue 0.01 -out_ascii_pssm /home/u2196/Desktop/KB8024/KB8024/SignalP/input/dimitri/$files.pssm
done
