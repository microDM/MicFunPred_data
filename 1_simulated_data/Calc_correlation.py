import pandas as pd
import os 
from collections import defaultdict

folders = ['4_predicted_97']
datasets = set([x.split('_')[0] for x in os.listdir(folders[0])])

#create metagenome dictionary
met_dict = {}
for i in datasets:
    try:
        df = pd.read_csv('3_metagenome/' + i + '_met.txt.gz',sep='\t',index_col=0,compression='gzip')
        met_dict[i] = df
    except:
        pass

out = open('correlations.txt','w')
out.write('Sample\tDataset\tTool\tPercent_identity\tCorrelation\n')

for i in folders:
    pident = i.split('_')[-1]
    for j in os.listdir(i):
        d = j.split('_')[0]
        tool = j.split('_')[1]
        df = pd.read_csv(os.path.join(i,j),sep='\t',index_col=0,compression='gzip')
        cordf = df.corrwith(met_dict[d],method='spearman').dropna()
        cordf = dict(zip(cordf.index,cordf.to_list()))
        for sam,corr in cordf.items():
            out.write(sam + '\t' + d + '\t' + tool + '\t' + str(pident) + '\t' + str(corr) + '\n')
            out.flush()