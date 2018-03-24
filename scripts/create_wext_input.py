# Event file. Tab-separated "adjacency list" format. Each line lists a sample in the first column, with each other column containing the name of an alteration in that sample (e.g. gene names). Warning: alteration names that appear in both the MAF and event files can overwrite one another.

# Create input in adjacency list format with each line being a sample and containing name of samples and list of positions
import sys
filename = sys.argv[1]
output = sys.argv[2]
import pandas as pd

data = pd.read_csv(filename)


samples = {}
for column in [d for d in data.columns if d.startswith('a_') or d.startswith('b_')]:
    sample = column.split('_')[-1]
    allele = column.split('_')[0]
    if sample not in samples:
        samples[sample]=[]

    vs = [str(v)+allele for v in list(data[data[column] > 0]['loc'])]
    samples[sample] += vs

with open(output, 'w') as out:
    for i,s in enumerate(samples):
        values = samples[s] 
        out.write(str(s) +'\t'+ '\t'.join(values)+'\n')

