 
import pandas as pd

filepath = 'globular_signal_tm_3state.txt'

raw_data = pd.read_csv(filepath, header=None)
    
drop_list = []
    
for l in range(len(raw_data)):
    if (l+1) % 3 == 0:
        drop_list.extend([l])
    
        
            
raw_data.drop(raw_data.index[drop_list], inplace=True)
raw_data.reset_index(drop = True, inplace = True)
    
raw_data.to_csv(filepath+'.fasta', header=None, index=None )
