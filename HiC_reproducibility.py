# -*- coding: utf-8 -*-
"""
Created on Tue Mar 06 11:02:12 2018

@author: xxli
"""

from __future__ import division
import numpy as np
from scipy import sparse
import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt
from scipy.sparse import lil_matrix
from scipy.sparse.linalg import eigsh

chrs = ['1' ,'2' ,'3' ,'4' ,'5' ,'6' ,'7' ,'8' ,'9' ,'10' ,'11' ,'12' 
        ,'13' ,'14' ,'15' ,'16' ,'17' ,'18' ,'19' ,'X']

data_fold = 'E:\\data\\nuclear_transfer\\Raw_matrix\\matrix_1M\\'
cell = ['fESC','NT5','NT6','CCS']
gap_file = np.loadtxt('E:\\data\\nuclear_transfer\\gap.txt', usecols=(1,2,3,6), dtype=[('chr', '<S6'), ('S', '<i4'), ('E', '<i4'), ('lenth', '<i4')])
res = 1000000



def get_matrix(data): 
    '''
    from sparse to matrix
    '''
    matrix = {}
    lenth = {}
    for g in chrs:
        tuple_matrix = data[g]
        sizes = tuple_matrix['bin2'].max() + 1
        P1 = sparse.coo_matrix((tuple_matrix['IF'],(tuple_matrix['bin1'],tuple_matrix['bin2'])),shape=(sizes,sizes)).toarray()
        P2 = sparse.coo_matrix((tuple_matrix['IF'],(tuple_matrix['bin2'],tuple_matrix['bin1'])),shape=(sizes,sizes)).toarray()
        mat = np.zeros((sizes, sizes))
        mat[np.nonzero(P1)] = P1[np.nonzero(P1)]
        mat[np.nonzero(P2)] = P2[np.nonzero(P2)]
        matrix[g] = mat
        lenth[g] = sizes
    return matrix , lenth

def Dis_distribution(matrix,lenth,gaps):
    dict_disttribution = {}
    for i in range(1,61):
        dict_disttribution[i] = []
    gaps['S'] = gaps['S']/res
    gaps['E'] = gaps['E']/res
    gaps['lenth'] = gaps['lenth']/res
    gaps = gaps[gaps['lenth'] > 0 ]
    for g in chrs:
#        matrix_lenth = np.zeros((lenth[g], lenth[g]))
        c_gap = gaps[gaps['chr'] == 'chr'+g]
        for x in c_gap:
            gap_site = range(x[1], x[2])
            matrix[g][gap_site,:] = 0
            matrix[g][:,gap_site] = 0
        for i in range(lenth[g]):
            for j in range(lenth[g]):
                d = j - i
                if (d >=1 and d <= 60):
                    dict_disttribution[d].append(matrix[g][i,j])
    return dict_disttribution
    

corr = {}
for i in cell:
    f1 = np.load(data_fold + i + '-MboI-R1-filtered-1M-sparse.npz')
    f2 = np.load(data_fold + i + '-MboI-R2-filtered-1M-sparse.npz')
    (matrix_1,lenth_1) = get_matrix(f1)
    (matrix_2,lenth_2) = get_matrix(f2)
    x = []
    dict_disttribution_1 = Dis_distribution(matrix_1,lenth_1,gap_file)
    dict_disttribution_2 = Dis_distribution(matrix_2,lenth_2,gap_file)
    for j in range(1,61):
        minlenth_1 = min(len(dict_disttribution_1[j]),len(dict_disttribution_2[j]))
        x.append(np.corrcoef(dict_disttribution_1[j][:minlenth_1],dict_disttribution_2[j][:minlenth_1])[0,1])
    corr[i] = x



fig = plt.figure(figsize = (12, 6))
left, width, bottom, height = 0.1, 0.80, 0.1, 0.60
size_axes = [left, bottom, width, height]
ax = fig.add_axes(size_axes)
ax.plot(range(1,61), corr['fESC'], label = 'fESC')
#ax.plot(range(1,61), corr['ESC'], label = 'ESC')
ax.plot(range(1,61), corr['NT5'], label = 'NT5')
ax.plot(range(1,61), corr['NT6'], label = 'NT6')
ax.plot(range(1,61), corr['CCS'], label = 'CCS')

#ax.set_xticks(chrlen)
#ax.set_xticklabels(chromosome)
ax.set_xlabel('Mb')
#ax.set_yticks([0.5, 0.75, 0.9])
#ax.set_yticklabels(['0.5', '0.75', '0.9'])
ax.set_ylabel('Correlation')
ax.legend()
#plt.ylim(0.5,1)
plt.title('Correlation between Biological replicates_1M')
