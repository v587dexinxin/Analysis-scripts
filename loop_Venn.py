# -*- coding: utf-8 -*-
"""
Created on Sat Apr 28 13:15:46 2018

@author: xxli
"""

from __future__ import division
import numpy as np
import sys, cPickle
import os
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

loopFolder = 'C:\\Users\\xxli\\Desktop\\loop'
chrom = ['1', '2', '3', '4', '5', '6', '7', '8', '9','10','11', '12', '13', '14', '15', '16', '17', '18', '19', 'X']

loopType = np.dtype({'names':['chr','S','E'],
                     'formats':['S2', np.int, np.int]})
uniontype = np.dtype({'names':['chr','S','E','ESC','NT5','NT6','CCS'],
                          'formats':['S8' , np.int , np.int, np.float,np.float,np.float,np.float]})

#---------------------------------union_loop--------------------------------
cell = ['ESC','NT5','NT6','CCS']
Len = 10000
union = {}
loops = {}
standard = 'ESC'
if standard == 'ESC':
    order = [1,2,3]
    o = 0
elif standard == 'NT5':
    order = [0,2,3]
    o = 1
elif standard == 'NT6':
    ororder = [0,1,3]
    o = 2
elif standard == 'CCS':
    order = [0,1,2]
    o = 3
s = os.path.join(loopFolder, standard + '_cluster_filter_nums_0.2-10K.txt')
loops[o] = np.loadtxt(s, usecols = (0,1,2) , dtype = loopType, skiprows = 1 )

for i in loops[0]:
    if i[0] not in union.keys():
        union[i[0]] = []
        union[i[0]].append(i)
    else:
        union[i[0]].append(i)
for k,v in union.items():
    v = np.array(v,loopType)
    union[k] = v
    
n1 = 0
n2 = 0
n3 = 0

for i in order:
    loopFil = cell[i] + '_cluster_filter_nums_0.2-10K.txt'
    loopSource = os.path.join(loopFolder, loopFil)
    loops[i] = np.loadtxt(loopSource, usecols = (0,1,2) , dtype = loopType, skiprows = 1 )           
    for j in loops[i]:
        mask = (j['S'] >=  union[j['chr']]['S'] - Len) & (j['S'] <=  union[j['chr']]['S'] + Len) & \
               (j['E'] >=  union[j['chr']]['E'] - Len) & (j['E'] <=  union[j['chr']]['E'] + Len)
        tmp = union[j['chr']][mask]
        if tmp.size == 0:
            t = list(union[j['chr']])
            t.append(j)
            union[j['chr']] = np.array(t , dtype = loopType)
            n3 += 1
        else:
          
            if tmp.size != 1:
                n2 += 1
#                print tmp
                
            elif tmp.size == 1:  
                n1 += 1
            
Loop = {}
for i in union:
    Loop[i] = []
    for j in range(len(union[i])):
        tmp = union[i][j]
        if tmp['E'] - tmp['S'] >= 150000:
            Loop[i].append(tmp)
union = Loop.copy()
for k,v in union.items():
    v = np.array(v,loopType)
    union[k] = v

Loop = {}
for i in range(4):
    Loop[i] = []
    for j in chrom:
        for k in range(len(loops[i][loops[i]['chr'] == j])):
            tmp = loops[i][loops[i]['chr'] == j][k]
            if tmp['E'] - tmp['S'] >= 150000:
                Loop[i].append(tmp)
    Loop[i] = np.array(Loop[i],loopType)
            
loops = Loop       
#-------------------------------loop_venn_figure---------------------------------------


ESC = loops[0]
NT5 = loops[1]
NT6 = loops[2]
CCS = loops[3]

single = {1:len(ESC) , 2:len(NT5) , 3:len(NT6) , 4:len(CCS)}
n1 = 0 ; n2 = 0 ; n3 = 0 ; n4 = 0
n12 = 0 ; n13 = 0 ; n14 = 0 ; n23 = 0 ; n24 = 0 ; n34 = 0
n123 = 0 ; n124 = 0 ; n134 = 0 ; n234 = 0
n1234 = 0
for i in union:
    for j in union[i]:
        mask_ESC = (j['S'] >=  ESC[j['chr']]['S'] - Len) & (j['S'] <=  ESC[j['chr']]['S'] + Len) & \
                   (j['E'] >=  ESC[j['chr']]['E'] - Len) & (j['E'] <=  ESC[j['chr']]['E'] + Len)
        mask_NT5 = (j['S'] >=  NT5[j['chr']]['S'] - Len) & (j['S'] <=  NT5[j['chr']]['S'] + Len) & \
                   (j['E'] >=  NT5[j['chr']]['E'] - Len) & (j['E'] <=  NT5[j['chr']]['E'] + Len)
        mask_NT6 = (j['S'] >=  NT6[j['chr']]['S'] - Len) & (j['S'] <=  NT6[j['chr']]['S'] + Len) & \
                   (j['E'] >=  NT6[j['chr']]['E'] - Len) & (j['E'] <=  NT6[j['chr']]['E'] + Len)
        mask_CCS = (j['S'] >=  CCS[j['chr']]['S'] - Len) & (j['S'] <=  CCS[j['chr']]['S'] + Len) & \
                   (j['E'] >=  CCS[j['chr']]['E'] - Len) & (j['E'] <=  CCS[j['chr']]['E'] + Len)
        ESC_c = ESC[j['chr']][mask_ESC] ; NT5_c = NT5[j['chr']][mask_NT5] 
        NT6_c = NT6[j['chr']][mask_NT6] ; CCS_c = CCS[j['chr']][mask_CCS]
        if ESC_c.size == 0 and NT5_c.size == 0 and NT6_c.size == 0 and CCS_c.size == 0:
            print 1
        elif ESC_c.size != 0 and NT5_c.size == 0 and NT6_c.size == 0 and CCS_c.size == 0:
            n1 += 1
        elif ESC_c.size == 0 and NT5_c.size != 0 and NT6_c.size == 0 and CCS_c.size == 0:
            n2 += 1
        elif ESC_c.size == 0 and NT5_c.size == 0 and NT6_c.size != 0 and CCS_c.size == 0:
            n3 += 1
        elif ESC_c.size == 0 and NT5_c.size == 0 and NT6_c.size == 0 and CCS_c.size != 0:
            n4 += 1
        elif ESC_c.size != 0 and NT5_c.size != 0 and NT6_c.size == 0 and CCS_c.size == 0:
            n12 += 1
        elif ESC_c.size != 0 and NT5_c.size == 0 and NT6_c.size != 0 and CCS_c.size == 0:
            n13 += 1
        elif ESC_c.size != 0 and NT5_c.size == 0 and NT6_c.size == 0 and CCS_c.size != 0:
            n14 += 1
        elif ESC_c.size == 0 and NT5_c.size != 0 and NT6_c.size != 0 and CCS_c.size == 0:
            n23 += 1
        elif ESC_c.size == 0 and NT5_c.size != 0 and NT6_c.size == 0 and CCS_c.size != 0:
            n24 += 1
        elif ESC_c.size == 0 and NT5_c.size == 0 and NT6_c.size != 0 and CCS_c.size != 0:
            n34 += 1
        elif ESC_c.size != 0 and NT5_c.size != 0 and NT6_c.size != 0 and CCS_c.size == 0:
            n123 += 1
            n12 += 1
            n23 += 1
            n13 += 1
        elif ESC_c.size != 0 and NT5_c.size != 0 and NT6_c.size == 0 and CCS_c.size != 0:
            n124 += 1
            n12 += 1
            n24 += 1
            n14 += 1
        elif ESC_c.size != 0 and NT5_c.size == 0 and NT6_c.size != 0 and CCS_c.size != 0:
            n134 += 1
            n13 += 1
            n14 += 1
            n34 += 1
        elif ESC_c.size == 0 and NT5_c.size != 0 and NT6_c.size != 0 and CCS_c.size != 0:
            n234 += 1
            n23 += 1
            n24 += 1
            n34 += 1
        elif ESC_c.size != 0 and NT5_c.size != 0 and NT6_c.size != 0 and CCS_c.size != 0:
            n1234 += 1
            n12 += 1 ; n13 += 1 ; n14 += 1 ; n23 += 1 ; n24 += 1 ; n34 += 1
            n123 += 1 ; n124 += 1 ; n134 += 1 ; n234 += 1
            
    


#n = 0
#for i in union:
#    for j in union[i]:
#        mask_ESC = (j['S'] >=  Loop_r[1][Loop_r[1]['chr'] == j['chr']]['S'] - 1) & (j['S'] <=  Loop_r[1][Loop_r[1]['chr'] == j['chr']]['S'] + 1) & \
#                   (j['E'] >=  Loop_r[1][Loop_r[1]['chr'] == j['chr']]['E'] - 1) & (j['E'] <=  Loop_r[1][Loop_r[1]['chr'] == j['chr']]['E'] + 1)
#        ESC_c = Loop_r[1][Loop_r[1]['chr'] == j['chr']][mask_ESC] 
#        if ESC_c.size != 0 :
#            n += 1
#
#n = 0            
#for i in ESC:
#    mask = (i['S'] >=  union[i['chr']]['S'] - 1) & (i['S'] <=  union[i['chr']]['S'] + 1) & \
#           (i['E'] >=  union[i['chr']]['E'] - 1) & (i['E'] <=  union[i['chr']]['E'] + 1)
#    ESC_c = union[i['chr']][mask]
#    if ESC_c.size != 0 :
#        n += 1
        
        
print 'area1=' + str(single[1])
print 'area2=' + str(single[2])
print 'area3=' + str(single[3])
print 'area4=' + str(single[4])
print 'n12=' + str(n12) 
print 'n13=' + str(n13)
print 'n14=' + str(n14)
print 'n23=' + str(n23)
print 'n24=' + str(n24)
print 'n34=' + str(n34)
print 'n123=' + str(n123)
print 'n124=' + str(n124)
print 'n134=' + str(n134)
print 'n234=' + str(n234)
print 'n1234=' + str(n1234)


