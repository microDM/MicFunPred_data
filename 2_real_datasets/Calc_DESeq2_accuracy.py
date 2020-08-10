import pandas as pd
import glob
import subprocess
import os

def run_deseq2(infile,grpfile,outfile):
    cmd = 'Rscript run_deseq2.R ' + infile + ' ' + grpfile + ' ' + outfile
    deseq2_run = subprocess.Popen(cmd,stderr=subprocess.PIPE,stdout=subprocess.PIPE,shell=True)
    output, error = deseq2_run.communicate()
    print(output,error)


#run deseq2
folders = ['2_humann2_metagenome','3_predicted_97','4_predicted_98','5_predicted_99','6_predicted_100']
"""
for i in folders:
    for j in os.listdir(i):
        data = j.split('_')[0]
        infile = os.path.join(i,j)
        grpfile = os.path.join('group',data+'.txt')
        outfile = os.path.join(i,j.replace('.tsv.gz','_deseq2.tsv'))
        print(i,j)
        run_deseq2(infile,grpfile,outfile)
"""

out = open('accuracy_deseq2.csv','w')
out.write(','.join(['Tool','Dataset','Perrcent_identity','TP','TN','FP','FN','TPR','FPR','Precision','F1 score','Balanced accuracy'])+'\n')
predDfList = []

folders.remove('2_humann2_metagenome')

for i in folders:
    for j in os.listdir(i):
        if(j.endswith('deseq2.tsv')):
            pident = i.split('_')[-1]
            tool = j.split('_')[1]
            d = j.split('_')[0]
            # prepare predictions
            pred_deseq2 = pd.read_csv(os.path.join(i,j),sep='\t',index_col=0)
            pred_deseq2 = pred_deseq2.loc[pred_deseq2['padj']<=0.05].dropna()
            met_deseq2 = pd.read_csv(os.path.join('2_humann2_metagenome',d+'_met_deseq2.tsv'),sep='\t',index_col=0)
            met_deseq2 = met_deseq2.loc[met_deseq2['padj']<=0.05].dropna()
            # possible kos
            possible_kos_tool = list(pd.read_csv(os.path.join('possible_ko',tool+'.txt'),sep='\t',index_col=0,header=None).index)
            possible_kos_met = list(pd.read_csv(os.path.join('possible_ko','humann2.txt'),sep='\t',index_col=0,header=None).index)
            # calculation
            tp = len(set(pred_deseq2.index).intersection(set(met_deseq2.index)))
            fp = len(set(pred_deseq2.index).difference(set(met_deseq2.index)))
            not_predicted_tool = set(possible_kos_tool).difference(set(pred_deseq2.index))
            not_predicted_met = set(possible_kos_met).difference(set(met_deseq2.index))
            tn = len(set(not_predicted_tool).intersection(not_predicted_met))
            fn = len(set(not_predicted_tool).intersection(set(met_deseq2.index)))
            # calculation of measures
            tpr = (tp/(tp+fn))*100
            fpr = (fp/(fp+tn))*100
            try:
                precision = (tp/(tp+fp))*100
            except ZeroDivisionError:
                precision = 0
            try:
                f1_score = 2*((tpr*precision)/(tpr+precision))
            except ZeroDivisionError:
                f1_score = 0
            try:
                balanced_acc = (0.5*((tp/(tp+fp)) + (tn/(tn+fn))))*100
            except ZeroDivisionError:
                balanced_acc = 0
            # write
            write_list = [str(x) for x in [tool,d,pident,tp,tn,fp,fn,tpr,fpr,precision,f1_score,balanced_acc]]
            out.write(','.join(write_list)+'\n')