# -*- coding: utf-8 -*-
"""
Created on Wed Aug 29 22:03:05 2018

@author: xxli
"""

from __future__ import division
import math
import numpy as np
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
from itertools import islice  
from sklearn.cluster import KMeans

cells = ['CCS_coverage_M','CCS_coverage_P','fESC_coverage_M','fESC_coverage_P','NT5_coverage_M','NT5_coverage_P','NT6_coverage_P','NT6_coverage_P']
cell = {'CCS':0,'fESC':1,'NT5':2,'NT6':3}
datatype = ({'names':['gene_name','M_C','P_C'],
             'formats':['<S64' , np.int , np.int]})

data = {}
union = set([])

for c in cells:
    f = np.loadtxt('C:\\Users\\xxli\\Desktop\\allelic_data\\Select\\RNA\\Select_Sort_New_' + c + '.txt' , usecols = ((1,4,5)) , skiprows = 1 , dtype = datatype)
    data[c] = f
    union = union.union(f['gene_name'])
union = list(union)
    

matrix = np.zeros((len(union),4))

for i in data:
    for j in data[i]:
        x = union.index(j['gene_name'])
        c = i.split("_")[0]
        y = cell[c]
        if j['M_C'] == 0 and j['P_C'] == 0:
            continue
        elif j['M_C'] == 0 and j['P_C'] != 0:
            matrix[x][y] = 15
        elif j['M_C'] != 0 and j['P_C'] == 0:
            matrix[x][y] = 11
        else:
            matrix[x][y] = math.log((j['M_C']/j['P_C']),2)
            
            
#matrix_1 = matrix        
matrix_1 = stats.zscore(matrix, axis=1, ddof=1)
matrix_1[np.isnan(matrix_1)] = 0

def Zeros_index(matrix):
    mask = (matrix[:,0] == 0) & (matrix[:,1] == 0) & (matrix[:,2] == 0) & (matrix[:,3] == 0)
    index = []
    n = 0
    for i in mask:
        if i == True:
            index.append(n)
        else:
            pass
        n += 1
    return index

index = Zeros_index(matrix_1)
matrix_1 = np.delete(matrix_1 , index , axis = 0)
            
#--------------------------------------K-means--------------------------------------------

def K_means_cluster(matrix_0,matrix,n):
    kmeans = KMeans(n_clusters=n, random_state=0).fit(matrix)
    a = kmeans.labels_
    a = a.reshape(len(a),1)

    matrix_1 = np.array(np.hstack((matrix_0,a)))
    matrix_1 = np.array(sorted(matrix_1, key=lambda x:x[-1]))
    return matrix_1

matrix_b = K_means_cluster(matrix_1,matrix_1,20)
    

#-------------------------------------plot-------------------------------------------------
left, bottom, width, height = 0.1, 0.1, 0.60, 0.80
size_axes = [left, bottom, width, height]
fig = plt.figure(figsize = (12, 12))
ax = fig.add_axes(size_axes)
#matrix_1 = np.hstack((matrix_b[:,:-1] , gene_loop_matrix_1))
matrix_1 = matrix_b[:,:-1]
vmax = np.percentile(matrix_1,95)
vmin = np.percentile(matrix_1,5)
im = ax.imshow(matrix_1,vmax=vmax,vmin=vmin,cmap='RdBu',aspect = 'auto')
plt.colorbar(im)
x = ['CCS','NT5','NT6','fESC']
ax.set_xticks(np.arange(len(x)))
ax.set_xticklabels(x,fontsize = 15)
plt.title('log2(M_C/P_C),gene_numbers:' + str(len(matrix_1)),fontsize = 20)