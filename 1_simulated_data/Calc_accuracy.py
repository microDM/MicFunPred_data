import pandas as pd
import os
import csv

def calcAccuracy(refDf,testDf,tool,d,pident,writer,pos_tool,pos_met):
    for i in set(testDf.columns).intersection(refDf.columns):
        ref_pred = refDf[i].loc[(refDf[i]>0)].index
        test_pred = testDf[i].loc[(testDf[i]>0)].index
        tp = len(set(ref_pred).intersection(test_pred))
        fp = len(set(test_pred).difference(ref_pred))
        ref_not_pred = list(refDf[i].loc[(refDf[i]==0)].index)
        ref_not_pred.extend(list(set(pos_met).difference(ref_pred)))
        test_not_pred = list(testDf[i].loc[(testDf[i]==0)].index)
        test_not_pred.extend(list(set(pos_tool).difference(test_pred)))
        tn = len(set(test_not_pred).intersection(ref_not_pred))
        fn = len(set(ref_pred).difference(set(test_not_pred)))
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
        temp = [tool,d,pident,i,tp,fp,tn,fn,tpr,fpr,precision,f1_score,balanced_acc]
        writer.writerow(temp)

################################ MAIN ######################################

# open file for performance
csvfile = open('accuracy.csv','w',newline='\n')
csvwriter = csv.writer(csvfile,delimiter=',')
csvwriter.writerow(['Tool','Data','Percent_identity','Sample','TP','FP','TN','FN','TPR','FPR','Precision','F1_score','BA'])

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

for i in folders:
    pident = i.split('_')[-1]
    for j in os.listdir(i):
        d = j.split('_')[0]
        tool = j.split('_')[1]
        df_pred = pd.read_csv(os.path.join(i,j),sep='\t',index_col=0,compression='gzip')
        df_met = met_dict[d]
        possible_ko_tool = list(pd.read_csv(os.path.join('../real_datasets/possible_ko',tool+'.txt'),sep='\t',index_col=0,header=None).index)
        possible_ko_met = list(pd.read_csv(os.path.join('../real_datasets/possible_ko','humann2.txt'),sep='\t',index_col=0,header=None).index)
        perform = calcAccuracy(df_met,df_pred,tool,d,pident,csvwriter,possible_ko_tool,possible_ko_met)
        csvfile.flush()
