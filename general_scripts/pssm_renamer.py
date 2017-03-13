# Script to rename all pssms

inpath = '/home/u2196/Desktop/KB8024/KB8024/SignalP/input/pssms_old/'
outpath ='/home/u2196/Desktop/KB8024/KB8024/SignalP/input/pssms/'


import os
import re
import shutil

for fil in os.listdir(inpath):
    out_name = re.sub('[^A-Za-z0-9]+','', fil.partition('.')[0])+''.join(fil.partition('.')[1:])
    if not os.path.exists(outpath):
        os.makedirs(outpath)
    shutil.copy(inpath+fil, outpath+out_name)


