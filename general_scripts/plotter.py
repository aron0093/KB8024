# Read data into pandas dataframe and generate stastistics

import pandas as pd


file_loc = ''

data = pd.read_csv(file_loc, headers=True, index=True)

data.plot(x=, y= [], kind='bar', colormap = 'Pastel1')
plt.figure(1)
plt.xlabel("")
plt.ylabel("")
plt.title("")
plt.figure(1, figsize=(20,10)).savefig("")
