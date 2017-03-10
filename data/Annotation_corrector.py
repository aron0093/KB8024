 
import pandas as pd

filepath = 'TOPDB_50.txt'

raw_data = pd.read_csv(filepath, header=None)
    
drop_list = []
    
for l in range(len(raw_data)):
    if (l+1) % 3 == 0:
        raw_data[0][l]=raw_data[0][l].replace('I','G').replace('T', 'M')
        
            

raw_data.to_csv(filepath+'_corrected.txt', header=None, index=None )
