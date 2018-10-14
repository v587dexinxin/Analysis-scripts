# -*- coding: utf-8 -*-
"""
Created on Fri Sep 14 21:49:05 2018

@author: Administrator
"""

from __future__ import division
import math
import numpy as np
import csv
from itertools import islice
from sklearn.cluster import KMeans
from scipy import stats
from matplotlib.backends.backend_pdf import PdfPages
from scipy.cluster.hierarchy import dendrogram, linkage
from scipy.cluster import hierarchy
from scipy import cluster 
import seaborn as sns
import copy
import scipy
import scipy.cluster.hierarchy as sch
from itertools import islice  
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap

#------------------------------------------------------------------------------------------------------------------------------
## Matplotlib Settings
matplotlib.rcParams['xtick.direction'] = 'out'
matplotlib.rcParams['ytick.direction'] = 'out'


# Our Own Color Map
my_cmap = plt.get_cmap('bwr')
my_cmap.set_bad('#D3D3D3')




chrs = ['1' , '2' , '3' , '4' , '5' , '6' , '7' , '8' , '9' , '10' , '11' , '12' , '13' , '14' , '15' , '16' , '17' , '18' , '19' , 'X']

cells = ['CCS_coverage_M','CCS_coverage_P','NT5_coverage_M','NT5_coverage_P','NT6_coverage_M','NT6_coverage_P','fESC_coverage_M','fESC_coverage_P']
cell = {'CCS':0 , 'NT5':1 , 'NT6':2 , 'fESC':3 }


datatype = ({'names':['gene_name' , 'chr' , 'start' , 'end' , 'M_C' , 'P_C'],
             'formats':['<S64' , 'S8' , np.int , np.int , np.int , np.int]})
gtf_type = ({'names':['chr' , 'feature' , 'start' , 'end' , 'strand' , 'gene_name'],
             'formats':['<S8' , 'S64' , np.int , np.int , 'S4' , 'S64']})
gtf_type_1 = ({'names':['gene_name' , 'chr' , 'start' , 'end'],
             'formats':['<S64' , 'S8' , np.int , np.int]})
peak_type = ({'names':['chr' , 'start' , 'end' , 'M_C' , 'P_C'],
             'formats':['<S8' , np.int , np.int , np.int , np.int]})
peak_type_1 = ({'names':['chr' , 'start' , 'end' , 'M_C' , 'P_C' , 'score'],
             'formats':['<S8' , np.int , np.int , np.int , np.int , np.float]})
peak_type_2 = ({'names':['chr' , 'start' , 'end' , 'score'],
             'formats':['<S8' , np.int , np.int , np.float]})
loop_type = ({'names':['chr' , 'start' , 'end'],
             'formats':['S8' , np.int , np.int]})



def Sort(a):
    '''
    a: list of needed to be sorted
    '''
    a.sort(key = lambda x:(x[0],x[1]))
    return a

def union_interval(a,peak_type):
    '''
    a: np.array of peaks ,dtype = peak_type
    '''
    peak_new = []
    a = np.array(a , dtype = peak_type)
    for i in a:
        mask = (i['start'] <= a[a['chr'] == i['chr']]['end']) & (i['end'] >= a[a['chr'] == i['chr']]['start'])
        overlap = a[a['chr'] == i['chr']][mask]
        if overlap.size != 0:
            x = overlap['start'].min()
            y = overlap['end'].max()
            peak_new.append((i[0],x,y))
            a = list(a)
            for j in overlap:
                a.remove(j)
            a = np.array(a, dtype = peak_type_1)
        else:
            continue

    peak_new.sort() 
    return peak_new
        
#----------------------------------------------allelic-gene-------------------------------------------------------------------
data = {}
union_name = set([])

for c in cells:
    f = np.loadtxt('D:\\workspace\\allelic_data\\Select\\RNA\\Select_Sort_New_' + c + '.txt' , usecols = ((1 , 0 , 2 , 3 , 4 , 5)) ,
                   skiprows = 1 , dtype = datatype)
    data[c] = f
    union_name = union_name.union(f['gene_name'])
union_name = list(union_name)

gtf = np.loadtxt('G:\\data\\genome\\gencode.vM15.chr_patch_hapl_scaff.annotation.gtf' , usecols = ((0, 2 , 3, 4 , 6 , 13)) , skiprows = 5 , dtype = gtf_type)
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

union_gene = []
for i in union_name:
    mask = gtf['gene_name'] == i
    gene = gtf[mask]
    try:
        union_gene.append(gene[0])
    except:
        continue
    
union_gene = np.array(union_gene , dtype = gtf_type_1)    


#----------------------------------------------allelic_gene_matrix----------------------------------------------------------------------


matrix_1 = np.zeros((len(union_gene),4))
    
n = 0
for c in cell:
    f = np.loadtxt('D:\\workspace\\allelic_data\\Raw\\RNA\\Sort_New_' + c + '_coverage.txt' , usecols = (1 , 0 , 2 , 3 , 4, 5) , 
                   skiprows = 1 , dtype = datatype)
    for i in range(len(union_gene)):
        gene_name = union_gene[i]['gene_name']
        overlap = f[f['gene_name'] == gene_name]
        if overlap.size != 0:
            if overlap[0]['M_C'] == 0 and overlap[0]['P_C'] == 0:
                continue
            elif overlap[0]['M_C'] == 0 and overlap[0]['P_C'] != 0:
                matrix_1[i][cell[c]] = -11
            elif overlap[0]['M_C'] != 0 and overlap[0]['P_C'] == 0:
                matrix_1[i][cell[c]] = 15
            else:
                matrix_1[i][cell[c]] = math.log((overlap[0]['M_C']/overlap[0]['P_C']),2)
        else:
            continue
        
            
        



##--------------------------------------------Z-score 标准化------------------------------------------------------------------------            
#matrix_1 = stats.zscore(matrix_1, axis=1, ddof=1)
#matrix_1[np.isnan(matrix_1)] = 0 

##--------------------------------------------matrix_index-------------------------------------------------------------------------
matrix_index = np.zeros((len(union_gene),1))
for i in range(len(union_gene)):
    matrix_index[i][0] = i

matrix_1 = np.hstack((matrix_1,matrix_index))

##-------------------------------------------K-means cluster------------------------------------------------------------------------

def K_means_cluster(matrix_0,matrix,n):
    kmeans = KMeans(n_clusters=n, random_state=0).fit(matrix)
    a = kmeans.labels_
    a = a.reshape(len(a),1)

    matrix_1 = np.array(np.hstack((matrix_0,a)))
    matrix_1 = np.array(sorted(matrix_1, key=lambda x:x[-1]))
    return matrix_1

matrix_1 = K_means_cluster(matrix_1,matrix_1[:,:-1],6)


##-------------------------------------------Sort_uniongene------------------------------------------------------------------------------

index = matrix_1[:,4]
sort_uniongene = []
for i in index:
    sort_uniongene.append(union_gene[int(i)])
    



#-----------------------------------------------gene_expression-enhancer_peak---------------------------------------------------------------    
loops = {}
peaks = {}
allelic_M = {}
allelic_P = {}

for c in cell:
    f1 = np.loadtxt('F:\\xxli\\data_analysis\\BDF1\HiC\\Loop_newSNP\\Clustered_' + c + '_hiccups_loops.txt' , skiprows = 1,
                    usecols = (0 , 1 , 2) , dtype = loop_type)
    f2 = np.loadtxt('D:\\workspace\\IDR\\selected_peaks\\idr_select_peaks_chr\\' + c + '_IDR_Selected_enhancer.narrowPeak' ,
                    usecols = (0 , 1 , 2 , 4 ) , dtype = peak_type_2)
    f3 = np.loadtxt('D:\\workspace\\allelic_data\\Select\\ATAC\\Select_Sort_New_' + c + '_coverage_M.bed' , skiprows = 1 ,
                    usecols = (0 , 1 , 2 , 3 , 4) , dtype = peak_type)
    f4 = np.loadtxt('D:\\workspace\\allelic_data\\Select\\ATAC\\Select_Sort_New_' + c + '_coverage_P.bed' , skiprows = 1 ,
                    usecols = (0 , 1 , 2 , 3 , 4) , dtype = peak_type)
    
    allelic_M[cell[c]] = f3
    allelic_P[cell[c]] = f4
    loops[cell[c]] = f1
    peaks[cell[c]] = f2 
        
    
peak_all = {0:[] , 1:[]}

for c in cell:    
    for i in range(len(sort_uniongene)):
        classify = matrix_1[i][-1]
        chro = sort_uniongene[i]['chr']     
        s = min(sort_uniongene[i]['start'],sort_uniongene[i]['end'])
        e = max(sort_uniongene[i]['start'],sort_uniongene[i]['end'])
        loop = loops[cell[c]]
        peak = peaks[cell[c]]
        M = allelic_M[cell[c]]
        P = allelic_P[cell[c]]
        mask_loop = (s <= loop[loop['chr'] == chro.lstrip("chr")]['end']) & (e >= loop[loop['chr'] == chro.lstrip("chr")]['start'])
        overlap = loop[loop['chr'] == chro.lstrip("chr")][mask_loop]
        enhancer_loop = {}
        if overlap.size != 0:
            for j in overlap:
                if min(abs(int(j['start']) - int(sort_uniongene[i][2])),abs(int(j['end']) - int(sort_uniongene[i][3]))) > 100000:
                    continue
                dis_score = {}
                mask_enhancer = (j['start'] <=  peak[peak['chr'] == chro]['end']) & \
                                (j['end'] >=  peak[peak['chr'] == chro]['start'])
                overlap1 = peak[peak['chr'] == chro][mask_enhancer]
                if overlap1.size != 0:
                    for k in overlap1:
                        dis = min(abs(j['start'] - k['start']) , abs(j['end'] - k['start']))
                        if dis > 100000:
                            continue
                        dis_score[dis] = k['score']
                    try:
                        enhancer_loop[dis_score[min(dis_score.keys())]] = j
                    except:
                        enhancer_loop[0] = j
                else:
                    enhancer_loop[0] = j
                    
            
            
            if (classify == 0) or (classify == 2) or (classify == 3) or (classify == 5):
                if enhancer_loop != {}:
                    loop = enhancer_loop[max(enhancer_loop.keys())]
                    mask_1 = (loop['start'] <=  M[M['chr'] == chro.lstrip("chr")]['end']) & \
                             (loop['end'] >= M[M['chr'] == chro.lstrip("chr")]['start'])
                else:
                    mask_1 = (s <= M[M['chr'] == chro.lstrip("chr")]['end']) & \
                             (e >= M[M['chr'] == chro.lstrip("chr")]['start'])
                overlap_1 = M[M['chr'] == chro.lstrip("chr")][mask_1]
                if overlap_1.size != 0:  
                    for x in overlap_1:
                        peak_all[0].append(x)
            elif (classify == 1) or (classify == 4):
                if enhancer_loop != {}:
                    loop = enhancer_loop[max(enhancer_loop.keys())]
                    mask_1 = (loop['start'] <=  P[P['chr'] == chro.lstrip("chr")]['end']) & \
                             (loop['end'] >= P[P['chr'] == chro.lstrip("chr")]['start'])
                else:
                    mask_1 = (s <= P[P['chr'] == chro.lstrip("chr")]['end']) & \
                             (e >= P[P['chr'] == chro.lstrip("chr")]['start'])
                overlap_1 = P[P['chr'] == chro.lstrip("chr")][mask_1]
                if overlap_1.size != 0:
                    for x in overlap_1:
                        peak_all[1].append(x)


#--------------------------------------loop_enhancer_matrix-----------------------------------------------------
loop_enhancer_matrix = np.zeros((len(peak_new)))
for i in range(len(peak_new)):
    chro = peak_new[i][0]
    s = peak_new[i][1]
    e = peak_new[i][2]
    for c in cell:
        peak = []
        for j in peaks_all[cell[c]]:
            for k in peaks_all[cell[c]][j]:
                peak.append(k)
        peak = np.array(peak , dtype= peak_type_1)
        mask = (s <= peak[peak['chr'] == chro]['end']) & (e >= peak[peak['chr'] == chro]['start'])
        overlap = peak[peak['chr'] == chro][mask]
        if overlap.size != 0:
            if overlap[0]['M_C'] == 0 and overlap[0]['P_C'] == 0:
                ss = 0
            elif overlap[0]['M_C'] == 0 and overlap[0]['P_C'] != 0:
                ss = -11
            elif overlap[0]['M_C'] != 0 and overlap[0]['P_C'] == 0:
                ss = 15
            else:
                ss = math.log((overlap[0]['M_C']/overlap[0]['P_C']),2)
            loop_enhancer_matrix[i][cell[c]] = ss
            
    

f = np.loadtxt('D:\\workspace\\allele-specific\\ATAC_seq\\gene_expression_peaks\\union.bed' , dtype = loop_type)    
new = []
for c in f:
    for i in a:
        mask = (i['start'] <= a[a['chr'] == i['chr']]['end']) & (i['end'] >= a[a['chr'] == i['chr']]['start'])
        overlap = a[a['chr'] == i['chr']][mask]
        if overlap.size != 0:
            x = overlap['start'].min()
            y = overlap['end'].max()
            new.append((i[0],x,y))
            a = list(a)
            for j in overlap:
                a.remove(j)
            a = np.array(a, dtype = peak_type)
        else:
            continue

    
    
    
o = open('D:\\workspace\\allele-specific\\ATAC_seq\\gene_expression_peaks\\all_allelic_peaks_classify2.bed' , 'w')
for i in peak_new[1]:
    o.writelines('chr' + i[0] + '\t' + str(i[1]) + '\t' + str(i[2]) + '\n')
o.close()            
 


#---------------------------------------------------matrix------------------------------------------------------------------


def score(overlap):
    if overlap.size != 0:
        if overlap[0]['M_C'] == 0 and  overlap[0]['P_C']== 0:
            score = 0
        elif overlap[0]['M_C'] == 0 and  overlap[0]['P_C']!= 0:
            score = -11
        elif overlap[0]['M_C'] != 0 and  overlap[0]['P_C']== 0:
            score = 15
        else:
            score = math.log((overlap[0]['M_C']/overlap[0]['P_C']),2)
    else:
        score = 0
    return score
        
        
        
        
matrix_2 = np.zeros((len(peaks) , 4))        
for i in range(len(peaks)):
    chro = peaks[i][0]
    s = peaks[i][1]
    e = peaks[i][2]
    mask1 = (s <= allelic_peaks[0][allelic_peaks[0]['chr'] == chro]['end']) & (e >= allelic_peaks[0][allelic_peaks[0]['chr'] == chro]['start'])
    mask2 = (s <= allelic_peaks[1][allelic_peaks[1]['chr'] == chro]['end']) & (e >= allelic_peaks[1][allelic_peaks[1]['chr'] == chro]['start'])
    mask3 = (s <= allelic_peaks[2][allelic_peaks[2]['chr'] == chro]['end']) & (e >= allelic_peaks[2][allelic_peaks[2]['chr'] == chro]['start'])
    mask4 = (s <= allelic_peaks[3][allelic_peaks[3]['chr'] == chro]['end']) & (e >= allelic_peaks[3][allelic_peaks[3]['chr'] == chro]['start'])
    overlap1 = allelic_peaks[0][allelic_peaks[0]['chr'] == chro][mask1]
    overlap2 = allelic_peaks[1][allelic_peaks[1]['chr'] == chro][mask2]
    overlap3 = allelic_peaks[2][allelic_peaks[2]['chr'] == chro][mask3]
    overlap4 = allelic_peaks[3][allelic_peaks[3]['chr'] == chro][mask4]
    if overlap1.size > 1 or overlap2.size > 1 or overlap3.size > 1 or overlap4.size >1:
        print overlap1.size , overlap2.size , overlap3.size , overlap4.size 
    
    score1 = score(overlap1)
    score2 = score(overlap2)
    score3 = score(overlap3)
    score4 = score(overlap4)
    matrix_2[i][0] = score1
    matrix_2[i][1] = score2
    matrix_2[i][2] = score3
    matrix_2[i][3] = score4

           
#----------------------------------plot-----------------------------------------------
#n = range(16)
#pp = PdfPages('C:\\Users\\xxli\\Desktop\\gene_expression\\gene_expression_individual_classify_16.pdf')



                
left, bottom, width, height = 0.1, 0.1, 0.60, 0.80
size_axes = [left, bottom, width, height]
fig = plt.figure(figsize = (12, 12))
ax = fig.add_axes(size_axes)
#matrix_1 = np.hstack((matrix_b[:,:-1] , gene_loop_matrix_1))
matrix = matrix_1[:,:-2]
vmax = np.percentile(matrix,95)
vmin = np.percentile(matrix,5)
im = ax.imshow(matrix,vmax=3,vmin=-3,cmap=my_cmap,aspect = 'auto')
plt.colorbar(im)
x = ['CCS','NT5','NT6','fESC']
ax.set_xticks(np.arange(len(x)))
ax.set_xticklabels(x,fontsize = 15)
plt.title('log2(M_C/P_C),gene_numbers:' + str(len(matrix)),fontsize = 20)
#pp.savefig(fig)
#plt.close()

