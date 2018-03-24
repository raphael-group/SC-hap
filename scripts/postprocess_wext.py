## The wext merged datafile
import sys

data_file = sys.argv[1]
output_file = sys.argv[2]

import pandas as pd
from sklearn.metrics import precision_recall_curve
from random import random
import math
from scipy.stats import chi2
import numpy as np

import logging
#filename = 
logging.basicConfig(filename="log/"+ sys.argv[0].split('/')[0] + ".log", level=logging.INFO)

#logging.getLogger().setLevel('DEBUG')

logging.debug("Data file: "+data_file)
logging.debug("Output file: "+output_file)

data = pd.read_csv(data_file, delimiter='\t')

logging.debug(data.columns)
data = data[(data['Index1'] != "Index1")] # get rid of extra headers from cat
data['Distance'] = data['Distance'].apply(int)
data['P-Value'] = data['P-Value'].apply(float)
data = data[data['Distance'] > 0]
data = data.drop_duplicates()

logging.debug("Imported data. DF size:"+str(data.shape))

data_sorted = data.fillna(1).replace(float('-inf'), 1).sort_values(by = 'P-Value')
data_sorted = data_sorted[data_sorted['P-Value'] > 0] 
data_sorted['Allele1'] = data_sorted['Gene1'].apply(lambda x: x[-1])
data_sorted['Allele2'] = data_sorted['Gene2'].apply(lambda x: x[-1])
#data_sorted['True'] = data_sorted.apply(lambda x: x['Gene1'][-1] == x['Gene2'][-1], axis=1)
data_sorted['True'] = data_sorted['Allele1'] == data_sorted['Allele2']
data_sorted['Pair'] = data_sorted.apply(lambda x: x['Gene1'][:-1] + ':' +  x['Gene2'][:-1], axis = 1)
data_sorted['Index1'] = data_sorted['Index1'].apply(int)
data_sorted['Index2'] = data_sorted['Index2'].apply(int)
data_sorted['log p'] = data_sorted['P-Value'].apply(math.log)

logging.debug("Sorted data. DF size:"+str(data_sorted.shape))

def p_value(x):
    k = 4
    v = chi2.logsf(-2*x, k)
    return v    

def is_correct(x):
    if x != 0:
        return x < 0
    else:
        return random() < 0.5
 
diff = data_sorted
diff = diff[['Pair', 'Distance',  'log p', 'True']].groupby(['Pair', 'Distance', 'True']).sum()
diff.reset_index(inplace=True)
diff['log p joint'] = diff['log p'].apply(p_value)
diff['log p sum'] = diff.apply(lambda x: x['log p joint'] if x['True'] else -x['log p joint'], axis=1)

logging.debug("Created diff. Diff size:" + str(diff.shape))

diff.to_csv(output_file)

