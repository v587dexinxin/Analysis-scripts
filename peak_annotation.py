# -*- coding: utf-8 -*-
"""
Created on Tue Mar 20 19:34:35 2018

@author: xxli
"""

import numpy as np
from itertools import islice

chrs = ['1' ,'2' ,'3' ,'4' ,'5' ,'6' ,'7' ,'8' ,'9' ,'10' ,'11' ,'12' 
        ,'13' ,'14' ,'15' ,'16' ,'17' ,'18' ,'19','X','Y']
peaktype = np.dtype({'names':['chr','start','end','name' ,'score' , 'strand' , 'value','p-value','q-value','peak'],
                      'formats':['<S8' , np.int , np.int , '<S20', np.float , '<S2' , np.float, np.float, np.float , np.int]})

cell ='CCS' 
#cell2 = '' 

f1 = np.loadtxt('C:\\Users\\xxli\\Desktop\\IDR\\selected_peaks\\idr_select_peaks_chr\\' + cell + '_IDR_Selected.narrowPeak' , dtype = peaktype)
f2 = 'E:\\data\\genome\\gencode.vM15.chr_patch_hapl_scaff.annotation.gtf'

w1 = 'C:\\Users\\xxli\\Desktop\\IDR\\selected_peaks\\idr_select_peaks_chr\\' + cell + '_IDR_Selected_promoter.narrowPeak'
w2 = 'C:\\Users\\xxli\\Desktop\\IDR\\selected_peaks\\idr_select_peaks_chr\\' + cell + '_IDR_Selected_enhancer.narrowPeak'


def Load_annotation(annotation_file):
    '''
    '''
    itype = np.dtype({'names':['chr' , 'start' , 'end' , 'strand' , 'gene_name'],
                      'formats':['<S8' , np.int , np.int , '<S2' , '<S30']})
    anno = {}
    f = open(annotation_file,'r')
    for line in islice(f,5,None):
        i = line.strip().split()
        j = line.split("gene_name")[1].split(";")[0]
        chro = i[0].lstrip('chr')
        if chro not in chrs:
            continue
        if i[2] != 'gene':
            continue
        if chro not in anno.keys():
            anno[chro] = []
            if i[6] == '-':
                anno[chro].append((i[0],i[4],i[3],i[6],j))
            else:
                anno[chro].append((i[0],i[3],i[4],i[6],j))
        else:
            if i[6] == '-':
                anno[chro].append((i[0],i[4],i[3],i[6],j))
            else:
                anno[chro].append((i[0],i[3],i[4],i[6],j))
    for c,v in anno.items():
        v = np.array(v,dtype = itype)
        anno[c] = v
    f.close()
    
    return anno

def annotation_peak(promoter,enhancer):
    '''
    '''
    n = 0
    anno = Load_annotation(f2)
    w1 = open(promoter,'w')
    w2 = open(enhancer,'w')
#    if euchr == 'euchr':
#        chro = chrs[:19]
#    else:
#        chro = ['X']
    for g in chrs:
        tmp_1 = f1[f1['chr'] == 'chr' + g]
        tmp_2 = anno[g]
        for p in tmp_1:
            n += 1
            mask = (tmp_2['start'] - 500 <= p['end']) & (tmp_2['start'] + 500 > p['start'])
            overlap = tmp_2[mask]
            p = np.array(list(p),dtype = np.str)
            if overlap.size != 0:
                w1.writelines('\t'.join(p) + overlap[0]['gene_name'] + '\n')
            else:
                w2.writelines('\t'.join(p) + '\n')
    w1.close()
    w2.close()
    return n 


n = annotation_peak(w1,w2)


                   

       
                   
                   
                   
                   
