# -*- coding: utf-8 -*-
"""
Created on Tue Sep 04 19:22:20 2018

@author: xxli
"""

import numpy as np


cells = ['CCS','NT5','NT6','fESC']
datatype = ({'names':['chr' ,	'start' , 'end' , 'M_C'	, 'P_C' , 'QuantileRatio' , 'log2(FC)' ,	'stat' ,	'P-value' , 'adj_P_value'],
             'formats':['S8' , np.int , np.int , np.int , np.int , 'S64' , 'S64' , 'S64' , 'S64' ,'S64']})
datatype_1 = ({'names':['chr' ,	'start' , 'end' , 'score'],
             'formats':['S8' , np.int , np.int , np.float]})
gtf_type = ({'names':['chr' , 'feature' , 'start' , 'end' , 'strand' , 'gene_name'],
             'formats':['<S8' , 'S64' , np.int , np.int , 'S4' , 'S64']})
gtf_type_1 = ({'names':['gene_name' , 'chr' , 'start' , 'end'],
             'formats':['<S64' , 'S8' , np.int , np.int]})


gtf = np.loadtxt('E:\\data\\genome\\gencode.vM15.chr_patch_hapl_scaff.annotation.gtf' , usecols = ((0, 2 , 3, 4 , 6 , 13)) , skiprows = 5 , dtype = gtf_type)
gtf_1 = []
for i in gtf:
    if i['feature'] == 'gene':
        if i['strand'] == '+':
            gtf_1.append((i['gene_name'].split("\"")[1] , i['chr'] , i['start'] , i['end']))
        elif i['strand'] == '-':
            gtf_1.append((i['gene_name'].split("\"")[1] , i['chr'] , i['end'] , i['start']))
        else:
            continue
    else:
        continue
    
gtf = np.array(gtf_1 , dtype = gtf_type_1)

for c in cells:
    f1 = np.loadtxt('C:\\Users\\xxli\\Desktop\\allelic_data\\Raw\\ATAC\\Sort_New_' + c + '_coverage.bed', skiprows = 1, dtype = datatype)
    f2 = np.loadtxt('C:\\Users\\xxli\\Desktop\\IDR\\selected_peaks\\idr_select_peaks_chr\\' + c + '_IDR_Selected.narrowPeak' , usecols = (0,1,2,4) ,dtype = datatype_1)    
    o = open('C:\\Users\\xxli\\Desktop\\allelic_data\\Raw\\ATAC\\Select_Sort_New_' + c + '_enhancer.bed','w')
    o.writelines('\t'.join(['chr' ,	'start' , 'end' , 'M_C'	, 'P_C' , 'QuantileRatio' , 'log2(FC)' ,	'stat' ,	'P-value' , 'adj_P_value' , 'peak_score']) + '\n')

    for i in f1:
        chro = 'chr' + i['chr']
        mask_gtf = (i['start'] <= (gtf[gtf['chr'] == chro]['start'] + 500)) & (i['end'] >= (gtf[gtf['chr'] == chro]['start'] - 500))
        overlap1 = gtf[gtf['chr'] == chro][mask_gtf]
        if overlap1.size != 0:
            continue
        mask = (i['start'] == f2[f2['chr'] == chro]['start']) & (i['end'] == f2[f2['chr'] == chro]['end'])
        overlap =  f2[f2['chr'] == chro][mask]
        if overlap.size != 0:
            o.writelines('\t'.join([str(x) for x in i]) + '\t' + str(overlap[0]['score']) + '\n')
        else:
            continue
    
    o.close()   
        
        