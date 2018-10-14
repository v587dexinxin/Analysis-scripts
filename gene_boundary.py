# -*- coding: utf-8 -*-
"""
Created on Thu Jul 26 15:18:56 2018

@author: xxli
"""

import csv
from itertools import islice
import numpy as np
import matplotlib.pyplot as plt
import os
import codecs

cell = 'NT6'

datatype = ({'names':['chr','gene_name','start','end'],
             'formats':['S2','S32',np.int,np.int]})
datatype_1 = ({'names':['chr','pos','CCS','fESC','NT5','NT6'],
             'formats':['S2',np.int,np.int,np.int,np.int,np.int]})
gene_M = np.loadtxt('C:\\Users\\xxli\\Desktop\\allele-specific\\RNA_seq\\RNA_newSNP\\Select_Sort_New_' + cell + '_coverage_M.txt',dtype = datatype,usecols=(0,1,2,3),skiprows = 1)
gene_P = np.loadtxt('C:\\Users\\xxli\\Desktop\\allele-specific\\RNA_seq\\RNA_newSNP\\Select_Sort_New_' + cell + '_coverage_P.txt',dtype = datatype,usecols=(0,1,2,3),skiprows = 1)
gene_B = np.loadtxt('C:\\Users\\xxli\\Desktop\\allele-specific\\RNA_seq\\RNA_newSNP\\Select_Sort_New_' + cell + '_coverage_B.txt',dtype = datatype,usecols=(0,1,2,3),skiprows = 1)

boundary_data = csv.reader(open('C:\\Users\\xxli\\Desktop\\allele-specific\\RNA_seq\\boundary_bias2gene_bias\\states.5_specific difference.csv'))

#---------------------------------------------------boundary---------------------------------------------------------

boundary = {}
for line in islice(boundary_data, 1, None):
    chro = line[0]
    pos = int(line[1])
    if line[2] == line[3]:
        CCS = 0
    elif line[2] != line[3] :
        CCS = 1
    if line[4] == line[5]:
        fESC = 0
    elif line[4] != line[5]:
        fESC = 1
    if line[6] == line[7]:
        NT5 = 0
    elif line[6] != line[7]:
        NT5 = 1
    if line[8] == line[9]:
        NT6 = 0
    elif line[8] != line[9]:
        NT6 = 1
    if chro not in boundary.keys():
        boundary[chro] = []
        boundary[chro].append((chro,pos,CCS,fESC,NT5,NT6))
    else:
        boundary[chro].append((chro,pos,CCS,fESC,NT5,NT6))
        
for k,v in boundary.items():
    v = np.array(v,dtype = datatype_1)
    boundary[k] = v
    
    
#--------------------------------------gene---------------------------------------------------------------------        

a = {'CCS':2,'fESC':3,'NT5':4,'NT6':5}

n = 0
m = 0
m1 = 0
m2 = 0

m_1 = []
m_2 = []
#for i in gene_B:
#    if i['chr'] == '9' or i['chr'] == 'Y':
#        continue
#    start = i['start']
#    mask = ((start-200000) <= boundary[i['chr']]['pos']) & ((start+200000) >= boundary[i['chr']]['pos'])  
#    overlap = boundary[i['chr']][mask]
#    if overlap.size == 0 :
#        m += 1
#        continue
#    elif overlap.size != 0:
#        distance = {}
#        for j in overlap:
#            d = abs(start-j['pos'])
#            distance[d] = j
#        d = min(distance.keys())
#        b = distance[d][a[cell]]
#    m_1.append(0)
#    m_2.append(b)
#    n += 1
#        

for i in gene_B:
    if i['chr'] == '9' or i['chr'] == 'Y':
        continue
    start = i['start']
    distance = {}
    for j in boundary[i['chr']]:
        d = abs(start-j['pos'])
        distance[d] = j
    d = min(distance.keys())
    b = distance[d][a[cell]]
    if b == 0:
        m1 += 1
    elif b == 1:
        m2 += 1
    m_1.append(1)
    m_2.append(b)
    n += 1

for i in gene_M:
    if i['chr'] == '9' or i['chr'] == 'Y':
        continue
    start = i['start']
    distance = {}
    for j in boundary[i['chr']]:
        d = abs(start-j['pos'])
        distance[d] = j
    d = min(distance.keys())
    b = distance[d][a[cell]]
    if b == 0:
        m1 += 1
    elif b != 0:
        m2 += 1
    m_1.append(1)
    m_2.append(b)
    n += 1

for i in gene_P:
    if i['chr'] == '9' or i['chr'] == 'Y':
        continue
    start = i['start']
    distance = {}
    for j in boundary[i['chr']]:
        d = abs(start-j['pos'])
        distance[d] = j
    d = min(distance.keys())
    b = distance[d][a[cell]]
    if b == 0:
        m1 += 1
    elif b != 0:
        m2 += 1
    m_1.append(1)
    m_2.append(b)
    n += 1    
    
#for i in gene_P:
#    if i['chr'] == '9' or i['chr'] == 'Y':
#        continue
#    start = i['start']
#    mask = ((start-100000000) <= boundary[i['chr']]['pos']) & ((start+100000000) >= boundary[i['chr']]['pos'])  
#    overlap = boundary[i['chr']][mask]
#    if overlap.size == 0:
#        m += 1
#        continue
#    elif overlap.size != 0:
#        distance = {}
#        for j in overlap:
#            d = abs(start-j['pos'])
#            distance[d] = j
#        d = min(distance.keys())
#        b = distance[d][a[cell]]
#    m_1.append(2)
#    m_2.append(b)
#    n += 1    
    

        
x = np.array(m_1)
x = x.reshape(x.shape[0],1)
y = np.array(m_2)
y = y.reshape(y.shape[0],1)

matrix = np.hstack((x,y))
plt.figure(figsize=(4, 10))   
plt.imshow(matrix, cmap = 'Blues', vmax = 2, vmin = 0,aspect = 'auto')
plt.title(cell,fontsize = 25)
plt.xticks([0,1],['Gene','Boundary'],fontsize = 20)



    