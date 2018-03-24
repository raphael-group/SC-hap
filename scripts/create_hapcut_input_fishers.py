## The wext merged datafile
import sys
input_file = sys.argv[1]
data_file = sys.argv[2]
output_file = sys.argv[3]
cutoff = float(sys.argv[4])
#cutoff = 5

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
logging.debug("cuttoff: "+str(cutoff))

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
data_sorted['True'] = data_sorted.apply(lambda x: x['Gene1'][-1] == x['Gene2'][-1], axis=1)
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

div = diff[['Pair', 'Distance', 'log p sum']].groupby(by = ['Pair','Distance']).sum()
div = div[div['log p sum'] != 0]
div['log p abs'] = div['log p sum'].apply(abs)
div['correct'] = div['log p sum'].apply(lambda x: x < 0)
div['correct2'] = div['log p sum'].apply(is_correct)
div.reset_index(inplace=True)

logging.debug("Created div. Div size:" + str(div.shape))
logging.debug("Score span: " + str(div['log p sum'].max()) +"\t"+ str(div['log p sum'].min()))

import vcf

pos_index_map = {}


data = pd.read_csv(input_file)    

pos_index_map = {}
def make_pos_index_map(x):
    pos_index_map[x['loc']] = x.name + 1
data.apply(make_pos_index_map, axis=1)

def write_fragment(cutoff, output_file):
    with open(output_file,'w') as out:
        dclip = div[div['log p abs'] > cutoff]
        #dclip = div
        for i in range(len(dclip)):
            d = dclip.iloc[i]
            v1 = min(map(int, d['Pair'].split(':'))) 
            v2 = max(map(int, d['Pair'].split(':')))

            i1 = pos_index_map[v1]
            i2 = pos_index_map[v2]
            a1 = "0"
            a2 = "0"

            if d['log p sum'] < 0:
                if random() > 0.5: allele = "00"
                else: allele = "11"
            else:
                if random() > 0.5: allele = "10"
                else: allele = "01"
                
            if abs(i1 - i2) == 1:
                q = '++'
                line = "1 {0} {1} {2} {3}\n".format(d['Pair'], i1, allele, q)
            else:
                q = '++'
                line = "2 {0} {1} {2} {3} {4} {5}\n".format(d['Pair'], i1, allele[0], i2, allele[1], q)
            out.write(line)
            
write_fragment(cutoff, output_file)   
