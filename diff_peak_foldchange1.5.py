# -*- coding: utf-8 -*-
"""
Created on Wed Mar 21 16:57:46 2018

@author: xxli
"""

from __future__ import division
import numpy as np
from scipy import stats
import math


cell1 = 'NT6' ; cell2 = 'NT5'
tad1 = 'C:\\Users\\xxli\\Desktop\\diff_peaks_new\\' + cell1 + '-MboI-40K-tworeps.txt'
tad2 = 'C:\\Users\\xxli\\Desktop\\diff_peaks_new\\' + cell2 + '-MboI-40K-tworeps.txt'
peak1 = 'C:\\Users\\xxli\\Desktop\\diff_peaks_new\\' + cell1 +'_promoter_euchr.narrowPeak'
peak2 = 'C:\\Users\\xxli\\Desktop\\diff_peaks_new\\' + cell2 +'_promoter_euchr.narrowPeak'
w1 = 'C:\\Users\\xxli\\Desktop\\diff_peaks_new\\results\\' + cell1 +'_' + cell1 + cell2 +'_diffpeaks_promoter_euchr.narrowPeak'
#w2 = 'C:\\Users\\xxli\\Desktop\\diff_peaks_new\\results\\' + cell2 +'_' + cell1 + cell2 +'_diffpeaks_promoter_euchr.narrowPeak'

chrs = ['1' ,'2' ,'3' ,'4' ,'5' ,'6' ,'7' ,'8' ,'9' ,'10' ,'11' ,'12' 
        ,'13' ,'14' ,'15' ,'16' ,'17' ,'18' ,'19' ,'X' , 'Y']
itype = np.dtype({'names':['chr','start','end','name' ,'score' , 'strand' , 'value','p-value','q-value','peak'],
                      'formats':['<S20' , np.int , np.int , '<S20', np.float , '<S2' , np.float, np.float, np.float , np.int]})
z_score = np.dtype({'names':['chr','start','end','name' ,'score' , 'strand' , 'value','p-value','q-value','peak','z_score'],
                      'formats':['<S20' , np.int , np.int , '<S20', np.float , '<S2' , np.float, np.float, np.float , np.int , np.float]})
t = np.dtype({'names':['start','end'],
                  'formats':[np.int , np.int]})
def Load_TAD(tad):
    TAD ={}
    f = open(tad,'r')
    for i in f:
        i = i.split()
        if i[-1] == '0':
            if i[0] not in TAD.keys():
                TAD[i[0]] = []
                TAD[i[0]].append((i[1],i[2]))
            else:
                TAD[i[0]].append((i[1],i[2]))
        else:
            continue
    for c,v in TAD.items():
        v = np.array(v,dtype = t)
        TAD[c] = v
    f.close()
    return TAD

def Load_Peaks(peak_file):
    '''
    '''
    Peak_In = {}
    f = open(peak_file,'r')
    for line in f:
        line = line.strip().split()
        chro = line[0].lstrip('chr')
        if chro not in chrs:
            continue
        else:
            if chro not in Peak_In.keys():
                Peak_In[chro] = []
                Peak_In[chro].append((line[0],line[1],line[2],line[3],line[4],line[5],line[6],line[7],line[8],line[9]))
            else:
                Peak_In[chro].append((line[0],line[1],line[2],line[3],line[4],line[5],line[6],line[7],line[8],line[9]))
    for c,v in Peak_In.items():
        v = np.array(v,dtype = itype)
        Peak_In[c] = v
    f.close()
    
    return Peak_In   

def diff_peak_inTAD(tad1,tad2,peak1,peak2,euchr,w):
    w1 = open(w,'w')
#    w2 = open(w2,'w')
    TAD1 = Load_TAD(tad1)
    TAD2 = Load_TAD(tad2)
    peak1 = Load_Peaks(peak1)
    peak2 = Load_Peaks(peak2)
#--------------------------------------merge_TAD--------------------------------------------------------
    if euchr == 'euchr':
        chro = chrs[:19]
    else:
        chro = ['X']
#    a =[]
    for g in chro:
        tmp_1 = TAD1[g]
        tmp_2 = TAD2[g]
        for p in tmp_1:
            mask = (tmp_2['start'] < p['end']) & (tmp_2['end'] > p['start'])
            overlap = tmp_2[mask]
            if overlap.size != 0:
                mask_1 = (peak1[g]['start']<p['end']) & (peak1[g]['end']>p['start'])
                overlap_1 = peak1[g][mask_1]
#                a.append(overlap_1.size)
                if overlap_1.size == 0:
                    continue
                

                score_1 = stats.zscore(overlap_1['score'])
                overlap1 = []
                for i in range(len(score_1)):
                    j = list(overlap_1[i])
                    if score_1[i] > 0:
                        j.append(score_1[i])
                        overlap1.append(tuple(j))
                    else:
                        continue   
                if overlap1 ==[]:
                    continue
                overlap_1 = np.array(overlap1 , dtype = z_score)
#                a_arg = overlap_1.argsort(order='z_score')
#                overlap_1 = overlap_1[a_arg][-5:]
                overlap_2 = []
                for i in overlap:
                    mask_2 = (peak2[g]['start']<i['end']) & (peak2[g]['end']>p['start'])
                    overlap_2_1 = peak2[g][mask_2]
                    if overlap_2_1.size == 0:
                        continue
                    score_2 = stats.zscore(overlap_2_1['score'])
                    overlap2 =[]
                    for j in range(len(score_2)):
                        k = list(overlap_2_1[j])
                        if score_2[j]>0 :
                            k.append(score_2[j])
                            overlap2.append(tuple(k))
                        else:
                            continue
                    if overlap2 ==[]:
                        continue
                    overlap_2_1 = np.array(overlap2 , dtype = z_score)
                    a_arg = overlap_2_1.argsort(order='z_score')
                    overlap_2_1 = overlap_2_1[a_arg][-5:]
                if overlap_2 ==[]:
                    overlap_2 = overlap_2_1
                else:
                    overlap_2 = np.concatenate((overlap_2,overlap_2_1),axis = 0)
                        
                for j in overlap_1:
                    mask_3 = (overlap_2['start']<j['end']) & (overlap_2['end']>j['start'])
                    overlap_3 = overlap_2[mask_3]
                    fold = []
                    if overlap_3.size != 0:
#                        a.append(overlap_3)
                        for k in overlap_3:
                            fold_change = math.log(j['z_score']/k['z_score'],2)
                            fold.append(fold_change)
                        fold_change = max(fold)
                    else:
                        fold_change =100
                        
                    if fold_change > 1.5 or fold_change<-1.5:
                        j = np.array(list(j),dtype = np.str)
                        w1.writelines('\t'.join(j[:10]) + '\n')
                    else:
                        continue
            else:
                continue
    w1.close()
diff_peak_inTAD(tad1,tad2,peak1,peak2,'euchr',w)           
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
#def diff_peak_inTAD(tad1,tad2,peak1,peak2,euchr,w):
#    w1 = open(w,'w')
#    TAD1 = Load_TAD(tad1)
#    TAD2 = Load_TAD(tad2)
#    TAD ={}
#    peak1 = Load_Peaks(peak1)
#    peak2 = Load_Peaks(peak2)
#    #-----------------merge_TAD---------------------
#    if euchr == 'euchr':
#        chro = chrs[:19]
#    else:
#        chro = ['X']
#    for g in chro:
#        tmp_1 = TAD1[g]
#        tmp_2 = TAD2[g]
#        TAD[g] = []
#        for p in tmp_1:
#            mask = (tmp_2['start'] < p['end']) & (tmp_2['end'] > p['start'])
#            overlap = tmp_2[mask]
#            if overlap.size != 0:
##                min_number = min((min(p),min([min(m) for m in overlap])))
##                max_number = max((max(p),max([max(m) for m in overlap])))
#                mask_1 = (peak1[g]['start']<p['end']) & (peak1[g]['end']>p['start'])
#                overlap_1 = peak1[g][mask_1]
#                score_1 = stats.zscore(overlap_1['score'])
#                overlap1 = []
#                for i in range(len(score_1)):
#                    j = list(overlap_1[i])
#                    j.append(score_1[i])
#                    overlap1.append(tuple(j))
#                overlap_1 = np.array(overlap1 , dtype = z_score)
#                overlap2 =[]
#                for i in overlap:
#                    mask_2 = (peak2[g]['start']<i['end']) & (peak2[g]['end']>p['start'])
#                    overlap_2 = peak2[g][mask_2]
#                    score_2 = stats.zscore(overlap_2['score'])
#                    for j in range(len(score_2)):
#                        k = list(overlap_2[j])
#                        k.append(score_2[j])
#                        overlap2.append(tuple(k))
#                overlap_2 =np.array(overlap2 , dtype = z_score)
#                for j in overlap_1:
#                    mask_3 = (overlap_2['start']<j['end']) & (overlap_2['end']>j['start'])
#                    overlap_3 = overlap_2[mask_3]
#                    fold = []
#                    if overlap_3.size != 0:
#                        for k in overlap_3:
#                            if ((j['z_score']) > 0 and (k['z_score']>0)) or ((j['z_score']) > 0 and (k['z_score']<0)):
#                                fold_change = abs(j['z_score'])/abs(k['z_score'])
#                                fold.append(fold_change)
#                                fold_change = max(fold)
#                            elif ((j['z_score']) < 0 and (k['z_score'] < 0)):
#                                fold_change = abs(k['z_score'])/abs(j['z_score'])
#                                fold.append(fold_change)
#                                fold_change = max(fold)
#                            else:
#                                fold_change = 0
#                    else:
#                        fold_change = 1000
#                    if fold_change > 1.5:
#                        j = np.array(list(j),dtype = np.str)
#                        w1.writelines('\t'.join(j[:10]) + '\n')
#                    else:
#                        continue
#            else:
#                mask_1 = (peak1[g]['start']<p['end']) & (peak1[g]['end']>p['start'])
#                overlap_1 = peak1[g][mask_1]
#                score_1 = stats.zscore(overlap_1['score'])
#                overlap1 = []
#                for i in range(len(score_1)):
#                    j = list(overlap_1[i])
#                    j.append(score_1[i])
#                    overlap1.append(tuple(j))
#                overlap_1 = np.array(overlap1 , dtype = z_score)
#                mask_2 = (peak2[g]['start']<p['end']) & (peak2[g]['end']>p['start'])
#                overlap_2 = peak2[g][mask_2]
#                score_2 = stats.zscore(overlap_2['score'])
#                overlap2 = []
#                for i in range(len(score_2)):
#                    j = list(overlap_2[i])
#                    j.append(score_2[i])
#                    overlap2.append(tuple(j))
#                overlap_2 = np.array(overlap2 , dtype = z_score)
#                for j in overlap_1:
#                    mask_3 = (overlap_2['start']<j['end']) & (overlap_2['end']>j['start'])
#                    overlap_3 = overlap_2[mask_3]
#                    fold = []
#                    if overlap_3.size != 0:
#                        for k in overlap_3:
#                            if ((j['z_score']) > 0 and (k['z_score']>0)) or ((j['z_score']) > 0 and (k['z_score']<0)):
#                                fold_change = abs(j['z_score'])/abs(k['z_score'])
#                                fold.append(fold_change)
#                                fold_change = max(fold)
#                            elif ((j['z_score']) < 0 and (k['z_score'] < 0)):
#                                fold_change = abs(k['z_score'])/abs(j['z_score'])
#                                fold.append(fold_change)
#                                fold_change = max(fold)
#                            else:
#                                fold_change = 0
#                    else:
#                        fold_change = 1000
#                    if fold_change > 1.5:
#                        j = np.array(list(j),dtype = np.str)
#                        w1.writelines('\t'.join(j[:10]) + '\n')
#                    else:
#                        continue
#                    
#n = diff_peak_inTAD(tad1,tad2,peak1,peak2,'euchr',w)
#                TAD[g].append((min_number,max_number))
#            else:
#                TAD[g].append(p)
#        for i in TAD[g]:
#            
#        for p in tmp_2:
#            mask = (tmp_1['start'] < p['end']) & (tmp_1['end'] > p['start'])
#            overlap = tmp_1[mask]
#            if overlap.size != 0:
#                continue
#            else:
#                tad.append(p)
#        tad.sort()
#        tmp = ()
#        for i in tad:
#            if tmp == ():
#                tmp = i
#            else:
#                if i[0] <= tmp[1]:
#                    tmp =(min(tmp[0],i[0]),max(tmp[1],i[1]))
#                else:
#                    TAD[g].append(tmp)
#                    tmp = i
#                if i ==tad[-1]:
#                    TAD[g].append(tmp)
#                else:
#                    continue
#                
#                    
                
            