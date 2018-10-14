# -*- coding: utf-8 -*-
"""
Created on Wed Aug 08 20:34:54 2018

@author: xxli
"""

from __future__ import division
import numpy as np
import csv
from itertools import islice
from sklearn.cluster import KMeans
from scipy import stats
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from scipy.cluster.hierarchy import dendrogram, linkage
from scipy.cluster import hierarchy
from scipy import cluster 
import seaborn as sns
import copy
import scipy
import scipy.cluster.hierarchy as sch
from scipy.cluster.vq import vq,kmeans,whiten


o = open('D:\\xxli\\data_analysis\\BDF1\\RNA_seq\\gene_expression\\all_genes.txt','w')

cells = ['CCS_1','CCS_2','CCS_3','fESC_1','fESC_2','fESC_3','NT5_1','NT5_2','NT5_3','NT5_4','NT6_1','NT6_2','NT6_3']
datatype = ({'names':['gene_id','gene_name','ref','strand','start','end','CCS_1','CCS_2','CCS_3','fESC_1','fESC_2','fESC_3','NT5_1','NT5_2','NT5_3','NT5_4','NT6_1','NT6_2','NT6_3'],
             'formats':['S64','S64','S32','S2','S32','S32',np.float,np.float,np.float,np.float,np.float,np.float,np.float,np.float,np.float,np.float,np.float,np.float,np.float]})
datatype_1 = ({'names':['gene_id','gene_name','ref','strand','start','end','FPKM'],
             'formats':['S64','S64','S32','S2','S32','S32',np.float]})

al_geneID = set([])

#----------------------------------------FPKM_union---------------------------------------------------------------------------------
for c in cells:
    f = np.loadtxt('D:\\xxli\\data_analysis\\BDF1\\RNA_seq\\gene_expression\\' + c + '.txt' ,usecols = (0),dtype = 'S64',skiprows = 1)
    f = set(f)
    al_geneID = al_geneID.union(f)
    
al_gene = {}

for i in al_geneID:
    al_gene[i] = np.array((0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),dtype = datatype)
    
for c in cells:
    f = np.loadtxt('D:\\xxli\\data_analysis\\BDF1\\RNA_seq\\gene_expression\\' + c + '.txt' ,usecols = (0,1,2,3,4,5,7),dtype = datatype_1,skiprows = 1)
    if c == 'CCS_1':
        for i in f:
            al_gene[i['gene_id']]['gene_id'] = i['gene_id']
            al_gene[i['gene_id']]['gene_name'] = i['gene_name']
            al_gene[i['gene_id']]['ref'] = i['ref']
            al_gene[i['gene_id']]['strand'] = i['strand']
            al_gene[i['gene_id']]['start'] = i['start']
            al_gene[i['gene_id']]['end'] = i['end']
            al_gene[i['gene_id']][c] += i['FPKM']
    else:
        for i in f:
            al_gene[i['gene_id']][c] += i['FPKM']
        
o.writelines('\t'.join(['gene_id','gene_name','ref','strand','start','end','CCS_1','CCS_2','CCS_3','fESC_1','fESC_2','fESC_3','NT5_1',
                        'NT5_2','NT5_3','NT5_4','NT6_1','NT6_2','NT6_3']) + '\n')
for k,v in al_gene.items():
    o.writelines('\t'.join([str(v['gene_id']),str(v['gene_name']),str(v['ref']),str(v['strand']),str(v['start']),str(v['end']),
                            str(v['CCS_1']),str(v['CCS_2']),str(v['CCS_3']),str(v['fESC_1']),str(v['fESC_2']),str(v['fESC_3']),
                            str(v['NT5_1']),str(v['NT5_2']),str(v['NT5_3']),str(v['NT5_4']),str(v['NT6_1']),str(v['NT6_2']),
                            str(v['NT6_3'])]) + '\n')
    
    
o.close()    





