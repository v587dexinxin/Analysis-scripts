# -*- coding: utf-8 -*-
"""
Created on Mon Mar 12 19:50:20 2018

@author: hluo
"""

from __future__ import division
import numpy as np

chrs = ['1' ,'2' ,'3' ,'4' ,'5' ,'6' ,'7' ,'8' ,'9' ,'10' ,'11' ,'12' 
        ,'13' ,'14' ,'15' ,'16' ,'17' ,'18' ,'19' ,'X' , 'Y']

itype = np.dtype({'names':['chr','start','end','name' ,'score' , 'strand' , 'value','p-value','q-value','peak'],
                      'formats':['<S8' , np.int , np.int , '<S64', np.float , '<S5' , np.float, np.float, np.float , np.int]})
def LoadPeaks(peak_file):
    '''
    '''
    itype = np.dtype({'names':['chr','start','end','name' ,'score' , 'strand' , 'value','p-value','q-value','peak'],
                      'formats':['<S8' , np.int , np.int , '<S64', np.float , '<S5' , np.float, np.float, np.float , np.int]})
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


def Find_Common_Peak(peak1,peak2):
    '''
    '''
    Common_Peak_1 = {}
    Common_Peak_2 = {}
    Peak1 = LoadPeaks(peak1)
    Peak2 = LoadPeaks(peak2)
    chro = Peak1.keys()
    n = 0
    with open('C:\\Users\\xxli\\Desktop\\IDR\\' + cell + '_Common_peaks.txt','w') as out:
        for c in chro:
            tmp_1 = Peak1[c]
            tmp_2 = Peak2[c]
            Common_Peak_1[c] = []
            Common_Peak_2[c] = []
            for p in tmp_1:
                mask = (tmp_2['start'] < p['end']) & (tmp_2['end'] > p['start'])
                overlap = tmp_2[mask]
                if overlap.size == 0:
                    continue
                else:
                    if overlap.size > 1:
                        n += 1
                        print overlap.size
                    out.writelines(p['chr'] + '\t' + str(p['start']) + '\t' + str(p['end']) + '\t' + str(overlap[0]['start'])+ '\t' + str(overlap[0]['end']) + '\t' + str(p['score']) + '\t' + str(overlap[0]['score']) + '\t' + str(p['p-value']) + '\t' + str(overlap[0]['p-value']) + '\n')
                    i = overlap[0]
                    Common_Peak_1[c].append((p['chr'],p['start'],p['end'],p['name'] ,p['score'] , p['strand'] , p['value'],p['p-value'],p['q-value'],p['peak'],i['p-value']))
                    Common_Peak_2[c].append((i['chr'],i['start'],i['end'],i['name'] ,i['score'] , i['strand'] , i['value'],i['p-value'],i['q-value'],i['peak'],p['p-value']))
#    with open('C:\\Users\\xxli\\Desktop\\IDR\\' + cell + '_R1_Common_peaks.txt','w') as out:
#        
#        for c in Common_Peak_1.keys():
#            for j in Common_Peak_1[c]:
#                j = np.array(list(j),dtype = str)
#                out.writelines('\t'.join(j)+'\n')
#    out.close()
#    with open('C:\\Users\\xxli\\Desktop\\IDR\\' + cell + '_R2_Common_peaks.txt','w') as out:
#        for c in Common_Peak_2.keys():
#            for j in Common_Peak_2[c]:
#                j = np.array(list(j),dtype = str)
#                out.writelines('\t'.join(j)+'\n')
#    out.close()
    
    return Common_Peak_1 ,Common_Peak_2,n

def Find_Common3_Peak(peak1,peak2,peak3):
    '''
    '''
    Common_Peak_1 = {}
    Common_Peak_2 = {}
    Common_Peak_3 = {}
    Peak1 = LoadPeaks(peak1)
    Peak2 = LoadPeaks(peak2)
    Peak3 = LoadPeaks(peak3)
    n = 0
    chro = Peak1.keys()
#    with open('C:\\Users\\xxli\\Desktop\\chrX\\' + cell + '_Common_peaks.txt','w') as out:
    for c in chro:
        tmp_1 = Peak1[c]
        tmp_2 = Peak2[c]
        tmp_3 = Peak3[c]
        Common_Peak_1[c] = []
        Common_Peak_2[c] = []
        Common_Peak_3[c] = []
        for p in tmp_1:
            mask = (tmp_2['start'] < p['end']) & (tmp_2['end'] > p['start'])
            overlap = tmp_2[mask]
            if overlap.size == 0:
                continue
            else:
#                    out.writelines(p['chr'] + '\t' + str(p['start']) + '\t' + str(p['end']) + '\t' + str(overlap[0]['start'])+ '\t' + str(overlap[0]['end']) + '\t' + str(p['score']) + '\t' + str(overlap[0]['score']) + '\t' + str(p['p-value']) + '\t' + str(overlap[0]['p-value']) + '\n')
                for i in overlap:
                    mask1 = (tmp_3['start'] < i['end']) & (tmp_3['end'] > i['start'])
                    overlap1 = tmp_3[mask1]
                    if overlap1.size == 0:
                        continue
                    else:
                        for j in overlap1:
                            Common_Peak_1[c].append((p['chr'],p['start'],p['end'],p['name'] ,p['score'] , p['strand'] , p['value'],p['p-value'],p['q-value'],p['peak'],i['p-value']))
                            Common_Peak_2[c].append((i['chr'],i['start'],i['end'],i['name'] ,i['score'] , i['strand'] , i['value'],i['p-value'],i['q-value'],i['peak'],p['p-value']))
                            Common_Peak_3[c].append((j['chr'],j['start'],j['end'],j['name'] ,j['score'] , j['strand'] , j['value'],j['p-value'],j['q-value'],j['peak'],i['p-value']))
                        n += 1
    return n 
def Find_diff_Peak(f1,f2,f3,euchr):
    '''
    '''
    Peak1 = LoadPeaks(f1)
    Peak2 = LoadPeaks(f2)
    diff_1 = 0
    diff_2 = 0
    common_1 = 0
    common_2 = 0
    common_all = 0
    Peak3 = LoadPeaks(f3)
    common = {}
#    n = 0
#    cell1num = 0
#    common_all = 0
#    itype = np.dtype({'names':['chr','start','end','name' ,'score' , 'strand' , 'value','p-value','q-value','peak'],
#                      'formats':['<S8' , np.int , np.int , '<S20', np.float , '<S2' , np.float, np.float, np.float , np.int]})
#    w1 = open('C:\\Users\\xxli\\Desktop\\chrX\\allreplicates\\' + cell1 + '_' + cell1 + cell2 + '_' + euchr + '_diff_peaks.txt','w')
#    w2 = open('C:\\Users\\xxli\\Desktop\\chrX\\allreplicates\\' + cell2 + '_' + cell1 + cell2 + '_' + euchr + '_diff_peaks.txt','w')
#    w_common = open('C:\\Users\\xxli\\Desktop\\diff_peak\\' + cell1 + cell2 + cell3 + '_' + euchr + '_common_peaks.txt','w')
    if euchr == 'euchr':
        chrom = chrs[:19]
    else:
        chrom = ['X']
        
    for c in chrom:
        x = []
        tmp_1 = Peak1[c]
        tmp_2 = Peak2[c]
        tmp_3 = Peak3[c]
        for p in tmp_1:
#            cell1num += 1
            mask = (tmp_2['start'] < p['end']) & (tmp_2['end'] > p['start'])
            overlap = tmp_2[mask]
            if overlap.size == 0:
                diff_1 += 1
#                w1.writelines('\t'.join([str(m) for m in p]) + '\n')
            else:
                common_1 += 1
                x.append (p)
        for p in tmp_2:
            mask = (tmp_1['start'] < p['end']) & (tmp_1['end'] > p['start'])
            overlap = tmp_1[mask]
            if overlap.size == 0:
                diff_2 += 1
#                w2.writelines('\t'.join([str(m) for m in p]) + '\n')
            else:
                common_2 += 1
#        w1.close()
#        w2.close()
#    print cell1num , n 
        
        common[c] = np.array(x,dtype = itype)
        for p in tmp_3:
            mask = (common[c]['start'] < p['end']) & (common[c]['end'] > p['start'])
            overlap = common[c][mask]
            if overlap.size == 0:
                continue
            else:
                common_all += 1
#                c2 = np.array(list(p),dtype = str)
#                w_common.writelines('\t'.join(c2)+'\n')
    return diff_1,common_1,diff_2,common_2,common_all
        

def Recalling_IDR_peaks(Common_1,IDR):
    '''
    '''
    IDR_type = np.dtype({'names':['P1','P2'],
                         'formats':[np.float,np.float]})
    IDR_P = np.loadtxt(IDR,dtype = IDR_type,usecols=(0,1))
    f1 = open(Common_1,'r')
#    f2 = open(Common_2,'r')
    with open('C:\\Users\\xxli\\Desktop\\IDR\\selected_peaks\\idr_select_peaks_chr\\' + cell + '_IDR_Selected.narrowPeak','w') as out:
        for line in f1:
            l = line.strip().split()
            mask = (IDR_P['P1'] == float(l[7]))&(IDR_P['P2'] == float(l[-1]))
            tmp = IDR_P[mask]
            if tmp.size != 0:
                out.writelines('chr' + l[0] + '\t' + '\t'.join(l[1:10]) + '\n')
            else:
                continue
#    with open('C:\\Users\\xxli\\Desktop\\chrX\\' + cell + '_R2_IDR_Selected.narrowPeak','w') as out:
#        for line in f2:
#            l = line.strip().split()
#            mask = (IDR_P['P1'] == float(l[-1]))&(IDR_P['P2'] == float(l[7]))
#            tmp = IDR_P[mask]
#            if tmp.size != 0:
#                out.writelines('\t'.join(l[:10]) + '\n')
#            else:
#                continue
#
#cell1 = 'NT6' ; cell2 = 'NT5' ; cell3 = 'NT6' ; euchr = 'chrX'
#f1 = 'C:\\Users\\xxli\\Desktop\\IDR\\selected_peaks\\idr_select_peaks_chr\\' + cell1 + '_IDR_Selected.narrowPeak'
#f2 = 'C:\\Users\\xxli\\Desktop\\IDR\\selected_peaks\\idr_select_peaks_chr\\' + cell2 + '_IDR_Selected.narrowPeak'
#f3 = 'C:\\Users\\xxli\\Desktop\\IDR\\selected_peaks\\idr_select_peaks_chr\\' + cell3 + '_IDR_Selected.narrowPeak'

def Find_diff_Peaks(f1,f2):
    '''
    '''
    Peak1 = LoadPeaks(f1)
    Peak2 = LoadPeaks(f2)
    diff_1 = 0
    diff_2 = 0
    for c in chrs:
        tmp_1 = Peak1[c]
        tmp_2 = Peak2[c]
        for p in tmp_1:
            mask = (tmp_2['start'] < p['end']) & (tmp_2['end'] > p['start'])
            overlap = tmp_2[mask]
            if overlap.size == 0:
                diff_1 += 1
                w1.writelines('\t'.join([str(m) for m in p]) + '\n')
            else:
                continue
        for p in tmp_2:
            mask = (tmp_1['start'] < p['end']) & (tmp_1['end'] > p['start'])
            overlap = tmp_1[mask]
            if overlap.size == 0:
                diff_2 += 1
                w2.writelines('\t'.join([str(m) for m in p]) + '\n')
            else:
                continue
    w1.close()
    w2.close()
#    print cell1num , n 
    return diff_1,diff_2



cell = 'CCS'  
f1 = 'C:\\Users\\xxli\\Desktop\\IDR\\' + 'uniq_' + cell + '_R1_peaks.narrowPeak'
f2 = 'C:\\Users\\xxli\\Desktop\\IDR\\' + 'uniq_' + cell + '_R3_peaks.narrowPeak'
(common_1,common_2,n) = Find_Common_Peak(f1,f2)
##
f3 = 'C:\\Users\\xxli\\Desktop\\IDR\\' + cell + '_R1_Common_peaks.txt'
#f4 = 'C:\\Users\\xxli\\Desktop\\chrX\\' + cell + '_R2_Common_peaks.txt'
f5 = 'C:\\Users\\xxli\\Desktop\\IDR\\' + cell + '_R1_R3_selected.txt'
Recalling_IDR_peaks(f3,f5)

cell1 = 'fESC' ; cell2 = 'NT6' 
f1 = 'C:\\Users\\xxli\\Desktop\\diff_peaks_new\\all\\' + cell1 + '_IDR_Selected.narrowPeak'
f2 = 'C:\\Users\\xxli\\Desktop\\diff_peaks_new\\all\\' + cell2 + '_IDR_Selected.narrowPeak'
w1 = open('C:\\Users\\xxli\\Desktop\\diff_peaks_new\\results\\' + cell1 + '_' + cell1 + cell2 + '.narrowPeak' , 'w')
w2 = open('C:\\Users\\xxli\\Desktop\\diff_peaks_new\\results\\' + cell2 + '_' + cell1 + cell2 + '.narrowPeak' , 'w')
Find_diff_Peaks(f1,f2)

#n = 1
#for i in f:
#    i = i.split()
#    w1.writelines(i[0] + '\t' + i[1] + '\t' + i[2] + '\t' + cell + '_peak_' + str(n) + '\t' + i[5] + '\t' + '.' + '\t' + '.' + '\t' + i[7] + '\t.\t.\n')
#    w2.writelines(i[0] + '\t' + i[3] + '\t' + i[4] + '\t' + cell + '_peak_' + str(n) + '\t' + i[6] + '\t' + '.' + '\t' + '.' + '\t' + i[8] + '\t.\t.\n')
#    n += 1
#f.close()
#w1.close()
#w2.close()
                
            
            
#cell1 = 'ESC' ; cell2 = 'NT5';cell3 = 'NT6'
#f1 = 'C:\\Users\\xxli\\Desktop\\chrX\\allreplicates\\' + cell1 + '_peaks.narrowPeak'
#f2 = 'C:\\Users\\xxli\\Desktop\\chrX\\allreplicates\\' + cell2 + '_peaks.narrowPeak'
#f3 = 'C:\\Users\\xxli\\Desktop\\chrX\\allreplicates\\' + cell3 + '_peaks.narrowPeak'
#diff_1,common_1,diff_2,common_2 ,common_all= Find_diff_Peak(f1,f2,f3,'chrX')

#cell1 = 'ESC' ; cell2 = 'NT5'
#f1 = 'C:\\Users\\xxli\\Desktop\\chrX\\diffpeaks\\' + cell1 + '_R1_IDR_Selected.narrowPeak'
#f2 = 'C:\\Users\\xxli\\Desktop\\chrX\\diffpeaks\\' + cell2 + '_R1_IDR_Selected.narrowPeak'
#(common_1,common_2,n) = Find_Common_Peak(f1,f2)

#cell1 = 'ESC' ; cell2 = 'NT5' ; cell3 = 'NT6'
#f1 = 'C:\\Users\\xxli\\Desktop\\chrX\\diffpeaks\\' + cell1 + '_R1_IDR_Selected.narrowPeak'
#f2 = 'C:\\Users\\xxli\\Desktop\\chrX\\diffpeaks\\' + cell2 + '_R1_IDR_Selected.narrowPeak'
#f3 = 'C:\\Users\\xxli\\Desktop\\chrX\\diffpeaks\\' + cell3 + '_R1_IDR_Selected.narrowPeak'
#n = Find_Common3_Peak(f1,f2,f3)

#cell1 = 'ESC' ; cell2 = 'NT5'
#f1 = 'C:\\Users\\xxli\\Desktop\\chrX\\diffpeaks\\' + cell1 + '_R1_IDR_Selected.narrowPeak'
#f2 = 'C:\\Users\\xxli\\Desktop\\chrX\\diffpeaks\\' + cell2 + '_R1_IDR_Selected.narrowPeak'
#diff_1,common_1,diff_2,common_2 = Find_diff_Peak(f1,f2,'chrX')
