"""
################### Linux ##################
What it does?
1. Extract specified variable region from 16S using vxtactor.pl
2. Multiplex data by inserting n mutations at random sites in sequence
# Make otu table
3. Make dirichlet distribution of 100 samples for n number of otus
4. Insert 40-60% random zeros in distribution 
# Extract gene contents
5. Extract KO, pfam and other gene contents
6. Normalize it with 16S rRNA gene copy numbers
# Make metagenome
7. Make metagenome

Author: Dattatray Mongad
Date: 11/07/2020

"""

from Bio import SeqIO
import os
import pandas as pd
import numpy as np
import subprocess
import argparse
import random

parser = argparse.ArgumentParser('Prepare data')
parser.add_argument('-i',help='16S sequence directory',dest='indir')

args = parser.parse_args()
indir = args.indir

random.seed(100)

def run_vxtractor(infile,outfile,region):
    cmd = 'perl /mnt/d/tools/vxtractor/vxtractor.pl -h /mnt/d/tools/vxtractor/HMMs/SSU/bacteria/ -r ' + region + ' -o ' + outfile + ' ' + infile
    print(cmd)
    sub = subprocess.Popen(cmd,stdin=subprocess.PIPE,stdout=subprocess.PIPE,shell=True)
    out,err = sub.communicate()
    #print(out,err)

def readVarFasta(fastaFile):
    d = dict()
    orgList = []
    for record in SeqIO.parse(fastaFile,'fasta'):
        org = record.description.split(' ')[1].split('_')[0]
        if(org not in orgList):
            d[org] = str(record.seq)
            orgList.append(org)
        else:
            d[org+'_'+str(orgList.count(org))] = str(record.seq)
            orgList.append(org)
    return d

def mutateSequencesMultipleTimes(d):
    d_to_return = d.copy()
    for k,v in d.items():
        # how many times a sequence has to mutated
        for i in range(1,random.choice(range(1,20))):
            seq = list(v)
            # select 3% random sites and mutate them
            for j in random.sample(k=round(0.03*len(seq)),population=range(0,len(seq))):
                bases = ['A','T','G','C','N']
                bases.remove(seq[j])
                seq[j] = random.choice(bases)
            d_to_return[str(k)+'_'+str(i)] = ''.join(seq) # save in dictionary
    return d_to_return

def justMutate(d):
    d_to_return = d.copy()
    for k,v in d.items():
        # how many times a sequence has to mutated
        seq = list(v)
        # select 3% random sites and mutate them
        for j in random.sample(k=round(0.03*len(seq)),population=range(0,len(seq))):
            bases = ['A','T','G','C','N']
            bases.remove(seq[j])
            seq[j] = random.choice(bases)
        d_to_return[k] = ''.join(seq) # save in dictionary
    return d_to_return

def writeFasta(d,outfile):
    out = open(outfile,'w')
    for k,v in d.items():
        out.write('>'+str(k)+'\n'+str(v)+'\n')
    out.close()

def simulateDirichlet(n):
    dist = list(np.random.lognormal(0.01,1,n))
    # 10% singletons
    n1 = round(0.1*len(dist))
    for i in range(0,n1):
        dist[i] = 1
    # 60% zeros
    n2 = n1 + round(0.6*len(dist))
    for i in range(n1+1,n2):
        dist[i] = 0
    for i in range(len(dist)):
        dist[i] = dist[i]*random.randint(1,1000)
    return dist


################### MAIN #####################
for i in os.listdir(indir):
    data = i.split('_')[0]
    # extract variable region
    if(not os.path.isdir('v3-v5')):
        os.mkdir('v3-v5')
    varfile = os.path.join('v3-v5',data+'_v3-v5.fasta')
    infile = os.path.join(indir,i)
    run_vxtractor(infile,varfile,'.V3-V5.')
    # mutate sequence
    d = readVarFasta(varfile)
    #d = mutateSequencesMultipleTimes(d)
    #d = justMutate(d)
    #varfile_mutated = varfile.split('.')[0] + '_mutated.fasta'
    #writeFasta(d,varfile_mutated)
    varfile_reformated = varfile.split('.')[0] + '_reformated.fasta'
    writeFasta(d,varfile_reformated)
    # dirichlet abundance table
    df_abund = pd.DataFrame(index=d.keys(),columns=['sample_'+str(x) for x in range(0,101)])
    for sam in range(101): 
        sampleName = 'sample_' + str(sam)
        dist = simulateDirichlet(len(d))
        random.shuffle(dist)
        df_abund[sampleName] = dist
    if(not os.path.isdir('OTU_table')):
        os.mkdir('OTU_table')
    df_abund.to_csv(os.path.join('OTU_table',data+'v3-v5.tsv'),sep='\t')
    # read 16S rRNA cp
    otu_ids = df_abund.index
    df_cp = pd.read_csv(os.path.join('1_new','files',data+'_16S_cp.txt'),sep='\t',index_col=0)
    df_cp = df_cp.reindex([int(x.split('_')[0]) for x in otu_ids])
    df_cp.index = otu_ids
    df_ko = pd.read_csv(os.path.join('1_new','files',data+'_ko.txt'),sep='\t',index_col=0)
    df_ko = df_ko.reindex([int(x.split('_')[0]) for x in otu_ids])
    df_ko.index = otu_ids
    # normalize with respect to 16S rRNA cp
    df_ko = df_ko.div(df_cp['0'].to_dict(),axis=0)
    if(not os.path.isdir('ko')):
        os.mkdir('ko')
    df_ko.to_csv(os.path.join('ko',data+'_ko.txt.gz',compression='gzip'),sep='\t')
    # multiply with abundance
    df_met = df_ko.transpose().dot(df_abund)
    if(not os.path.isdir('metagenome')):
        os.mkdir('metagenome')
    df_met.to_csv(os.path.join('metagenome',data+'_met.txt.gz',compression='gzip'),sep='\t')

