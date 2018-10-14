# -*- coding: utf-8 -*-
"""
Created on Sat Aug 11 10:15:32 2018

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
from itertools import islice  

#from scipy.cluster.vq import vq,kmeans,whiten

chrs = ['1' , '2' , '3' , '4' , '5' , '6' , '7' , '8' , '9' , '10' , '11' , '12' 
        , '13' , '14' , '15' , '16' , '17' , '18' , '19' , 'X']

cells = ['CCS_1','CCS_2','CCS_3','fESC_1','fESC_2','fESC_3','NT5_1','NT5_2','NT5_3','NT5_4','NT6_1','NT6_2','NT6_3']
datatype = ({'names':['gene_id','gene_name','chr','strand','start','end','CCS_1','CCS_2','CCS_3','fESC_1','fESC_2','fESC_3','NT5_1','NT5_2','NT5_3','NT5_4','NT6_1','NT6_2','NT6_3'],
             'formats':['S64','S64','S32','S2','S32','S32',np.float,np.float,np.float,np.float,np.float,np.float,np.float,np.float,np.float,np.float,np.float,np.float,np.float]})
#datatype_1 = ({'names':['gene_name','chr','start','CCS','fESC','NT5','NT6'],
#             'formats':['S64','S2','S32',np.float,np.float,np.float,np.float]})
peaktype = ({'names':['chr','start','end','score'],
             'formats':['<S8' , np.int , np.int , np.float]})
looptype = ({'names':['chr' , 'start' , 'end' , 'CCS'  , 'NT5' , 'NT6' , 'fESC'],
             'formats':['S64' , np.int , np.int , np.float , np.float , np.float , np.float]})
looptype_1 = ({'names':['chr' , 'start' , 'end'],
               'formats':['S64' , np.int , np.int]})
datatype_2 = ({'names':['CCS','NT5','NT6','fESC','CCS_peak','NT5_peak','NT6_peak','fESC_peak','classify'],
             'formats':[np.float,np.float,np.float,np.float,np.float,np.float,np.float,np.float,np.float]})

#----------------------------------------------raw_matrix----------------------------------------------------------------------

diff_fil = ['Filtered_new_fESC_vs_CCS.csv','Filtered_new_CCS_vs_NT5.csv','Filtered_new_CCS_vs_NT6.csv',
            'Filtered_new_fESC_vs_NT5.csv','Filtered_new_fESC_vs_NT6.csv','Filtered_new_NT5_vs_NT6.csv']

gene_id = []
for f in diff_fil:
    data = csv.reader(open('D:\\xxli\\data_analysis\\BDF1\\RNA_seq\\difference_expression\\fc_1.5\\' + f,'r'))
    for i in islice(data,1,None):
        gene_id.append(i[0])
        
gene_id = set(gene_id)

        
        
f = np.loadtxt('E:\\data\\genome\\gencode.vM15.chr_patch_hapl_scaff.annotation.gtf' ,dtype = datatype ,skiprows = 1)

diff_gene = []
for i in gene_id:
    diff_gene.append(f[f['gene_id'] == i])
    
matrix_a = np.zeros((len(diff_gene),3),dtype = 'S32')
matrix_b = np.zeros((len(diff_gene),4))
for i in range(len(diff_gene)):
    matrix_a[i][0] = diff_gene[i]['gene_name'][0]
    matrix_a[i][1] = diff_gene[i]['chr'][0]
    if diff_gene[i]['strand'] == '+':
        matrix_a[i][2] = diff_gene[i]['start'][0]
    elif diff_gene[i]['strand'] == '-':
        matrix_a[i][2] = diff_gene[i]['end'][0]
    else:
        print 1
        continue
for i in range(len(diff_gene)):
    matrix_b[i][0] = (float(diff_gene[i]['CCS_1']) + float(diff_gene[i]['CCS_2']) + float(diff_gene[i]['CCS_3']))/3
    matrix_b[i][1] = (float(diff_gene[i]['NT5_1']) + float(diff_gene[i]['NT5_2']) + float(diff_gene[i]['NT5_3']) + float(diff_gene[i]['NT5_4']))/4
    matrix_b[i][2] = (float(diff_gene[i]['NT6_1']) + float(diff_gene[i]['NT6_2']) + float(diff_gene[i]['NT6_3']))/3
    matrix_b[i][3] = (float(diff_gene[i]['fESC_1']) + float(diff_gene[i]['fESC_2']) + float(diff_gene[i]['fESC_3']))/3
    
#matrix = np.hstack((matrix_a,matrix_b))
#matrix = np.array(matrix,dtype = datatype_1)
    
#--------------------------------Z-score标准化--------------------------------------------------------------------
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

index = Zeros_index(matrix_b)
matrix_a = np.delete(matrix_a , index , axis = 0)
matrix_b = np.delete(matrix_b , index , axis = 0)

matrix_b = stats.zscore(matrix_b, axis=1, ddof=1)    
        

gene_fpkm = {}
index = []
for i in range(len(matrix_a)):
    if matrix_a[i][0] in gene_fpkm.keys():
        index.append(i)
        continue
    gene_fpkm[matrix_a[i][0]] = list(matrix_b[i])

matrix_a = np.delete(matrix_a , index , axis = 0)
matrix_b = np.delete(matrix_b , index , axis = 0)    
    
    


    
#------------------------------- -1，1标准化------------------------------------------------------------------
matrix_new = np.zeros((len(diff_gene),4))
for i in range(len(matrix)):
    if matrix[i].max() == matrix[i].min():
        continue
    for j in range(4):
        norm = 2 * (matrix[i][j]  - matrix[i].min())/(matrix[i].max()-matrix[i].min()) - 1
        matrix_new[i][j] = norm
matrix = matrix_new        

#------------------------------------K-means cluster-------------------------------------------------------      

def K_means_cluster(matrix_0,matrix,n):
    kmeans = KMeans(n_clusters=n, random_state=0).fit(matrix)
    a = kmeans.labels_
    a = a.reshape(len(a),1)

    matrix_1 = np.array(np.hstack((matrix_0,a)))
    matrix_1 = np.array(sorted(matrix_1, key=lambda x:x[-1]))
    return matrix_1

matrix_b = K_means_cluster(matrix_b,matrix_b,6)
#new_matrix_b = np.zeros((len(matrix_b),5))
#index =[0.0,1.0,3.0,5.0,2.0,4.0]
#n = 0
#for i in index:
#    for j in matrix_b:
#        if j[-1] == i:   
#            new_matrix_b[n] = j
#            n += 1
#        else:
#            pass
#matrix_b = new_matrix_b


#-----------------------------------hierarchical clustering---------------------------------------------------------------
sns.clustermap(matrix,method ='ward',cmap= 'bwr',metric='euclidean',xticklabels=['CCS','NT5','NT6','fESC'])



#---------------------------------fPKM_gene------------------------------------------

gene = list(gene_fpkm.keys())[list(gene_fpkm.values()).index(list(matrix_b[:,:-1][0]))]
mask = (matrix_a[:,0] == gene)
sort_gene = matrix_a[mask]
for i in islice(matrix_b[:,:-1],1,None):
    i = list(i)
    gene = list(gene_fpkm.keys())[list(gene_fpkm.values()).index(i)]
    mask = (matrix_a[:,0] == gene)
    sort_gene = np.vstack((sort_gene,matrix_a[mask]))
    
    
#----------------------------gene_expression-promoter_peak--------------------------------------------
    
cell = {'CCS':0 , 'NT5':1 , 'NT6':2 , 'fESC':3}
matrix_peaks = np.zeros((len(gene_fpkm),4))
for c in cell:
    f = np.loadtxt('C:\\Users\\xxli\\Desktop\\IDR\\' + c + '_peaks.narrowPeak' , usecols = (0,1,2,6) , dtype = peaktype)
    for i in range(len(sort_gene)):
        s = int(sort_gene[:,-1][i]) - 500
        e = int(sort_gene[:,-1][i]) + 500
        mask = (f['start'] < e) & (f['end'] > s)
        overlap = f[mask]
        if overlap.size != 0:
#            if overlap.size > 1:
#                print overlap.size
            score = sum(overlap['score'])/len(overlap)
            matrix_peaks[i][cell[c]] = score
        else:
            matrix_peaks[i][cell[c]] =  -1
    
matrix_peaks = stats.zscore(matrix_peaks, axis=1, ddof=1)
matrix_peaks[np.isnan(matrix_peaks)] = -1

#----------------------------gene_expression-enhancer_peak--------------------------------------------    

f_CCS = np.loadtxt('D:\\xxli\\data_analysis\\BDF1\\HiC\\Loop_newSNP\\Clustered_CCS_hiccups_loops.txt',skiprows = 1,usecols = ((0,1,2)),dtype = looptype_1)
f_NT5 = np.loadtxt('D:\\xxli\\data_analysis\\BDF1\\HiC\\Loop_newSNP\\Clustered_NT5_hiccups_loops.txt',skiprows = 1,usecols = ((0,1,2)),dtype = looptype_1)
f_NT6 = np.loadtxt('D:\\xxli\\data_analysis\\BDF1\\HiC\\Loop_newSNP\\Clustered_NT6_hiccups_loops.txt',skiprows = 1,usecols = ((0,1,2)),dtype = looptype_1)
f_fESC = np.loadtxt('D:\\xxli\\data_analysis\\BDF1\\HiC\\Loop_newSNP\\Clustered_fESC_hiccups_loops.txt',skiprows = 1,usecols = ((0,1,2)),dtype = looptype_1)

loops = {0:f_CCS , 1:f_NT5 , 2:f_NT6 , 3:f_fESC}
union = {}
union['Y'] = np.zeros((1,3),dtype = looptype_1)
Len = 20000
for i in f_CCS:
    if i['chr'] not in union.keys():
        union[i['chr']] = []
        union[i['chr']].append(i)
    else:
        union[i['chr']].append(i)
for k , v in union.items():
    v = np.array(v,dtype = looptype_1)
    union[k] = v
        
for i in [f_NT5,f_NT6,f_fESC]:
    for j in i:
        chro = j['chr']
        tmp = union[chro]
        mask = (j['start'] >= tmp['start'] - Len) & (j['start'] <= tmp['start'] + Len) & \
                       (j['end'] >= tmp['end'] - Len) & (j['end'] <= tmp['end'] + Len)
        overlap = tmp[mask]
        if overlap.size == 0:
            tmp = list(tmp)
            tmp.append(j)
            union[chro] = np.array(tmp,dtype = looptype_1)
        else:
            continue
            


peak_CCS = np.loadtxt('C:\\Users\\xxli\\Desktop\\IDR\\selected_peaks\\idr_select_peaks_chr\\CCS_IDR_Selected_enhancer.narrowPeak' , usecols = ((0,1,2,4)) , dtype = peaktype)
peak_NT5 = np.loadtxt('C:\\Users\\xxli\\Desktop\\IDR\\selected_peaks\\idr_select_peaks_chr\\NT5_IDR_Selected_enhancer.narrowPeak' , usecols = ((0,1,2,4)) , dtype = peaktype)
peak_NT6 = np.loadtxt('C:\\Users\\xxli\\Desktop\\IDR\\selected_peaks\\idr_select_peaks_chr\\NT6_IDR_Selected_enhancer.narrowPeak' , usecols = ((0,1,2,4)) , dtype = peaktype)
peak_fESC = np.loadtxt('C:\\Users\\xxli\\Desktop\\IDR\\selected_peaks\\idr_select_peaks_chr\\fESC_IDR_Selected_enhancer.narrowPeak' , usecols = ((0,1,2,4)) , dtype = peaktype)

peaks = {0:peak_CCS , 1:peak_NT5 , 2:peak_NT6 , 3:peak_fESC}


#(---------------------------------------loop_strength_matrix---------------------------------------------------)
#loop_matrix = np.zeros((len(f),),dtype = looptype)
#n = 0
#for i in f_CCS:
#    loop_matrix[n]['chr'] = i[0]
#    loop_matrix[n]['start'] = i[1]
#    loop_matrix[n]['end'] = i[2]
#    loop_matrix[n]['CCS'] = f_CCS[i]
#    loop_matrix[n]['NT5'] = f_NT5[i]
#    loop_matrix[n]['NT6'] = f_NT6[i]
#    loop_matrix[n]['fESC'] = f_fESC[i]
#    n += 1


gene_loop_matrix = np.zeros((len(sort_gene),4))    

for i in range(len(sort_gene)):
    index = list(matrix_b[:,:-1][i]).index(matrix_b[:,:-1][i].max())
    chro = sort_gene[i][1].lstrip('chr')
    gene_start = int(sort_gene[i][2])
#    mask = (gene_start >= loops[index][loops[index]['chr'] == chro]['start']) & (gene_start <= loops[index][loops[index]['chr'] == chro]['end'])
    mask = (gene_start >= union[chro]['start']) & (gene_start <= union[chro]['end'])
#    overlap = loops[index][loops[index]['chr'] == chro][mask]
    overlap = union[chro][mask]
    if overlap.size == 0:
        continue
    else:
        n = overlap.size 
        enhancer_loop = {}
        for j in overlap:
            if min(abs(int(j['start']) - int(sort_gene[i][-1])),abs(int(j['end']) - int(sort_gene[i][-1]))) > 100000:
                gene_loop_matrix[i] = np.array([-1,-1,-1,-1])
                continue
            mask1 = (j['start'] <=  peaks[index][peaks[index]['chr'] == 'chr' + chro]['end']) & \
                    (j['end'] >=  peaks[index][peaks[index]['chr'] == 'chr' + chro]['start'])
            overlap1 = peaks[index][peaks[index]['chr'] == 'chr' + chro][mask1]
            if overlap1.size == 0:
                enhancer_loop[0] = j
            else:
                dis_score = {}
                for k in overlap1:
                    dis = min(abs(j['start'] - k['start']) , abs(j['end'] - k['start']))
                    if dis > 100000:
                        continue
                    dis_score[dis] = k['score']
                        
                    
                if dis_score != {}:
                    enhancer_loop[dis_score[min(dis_score.keys())]] = j
#                        score = max(dis_score.values())
                else:
                    enhancer_loop[0] = j
        if enhancer_loop != {}:  
            loop = enhancer_loop[max(enhancer_loop.keys())]
            for c in peaks:
                mask2 = (loop['start'] <=  peaks[c][peaks[c]['chr'] == 'chr' + chro]['end']) & \
                        (loop['end'] >=  peaks[c][peaks[c]['chr'] == 'chr' + chro]['start'])
                overlap2 = peaks[c][peaks[c]['chr'] == 'chr' + chro][mask2]
                if overlap2.size != 0:
                    dis_score = {}
                    for m in overlap2:
                        dis = min(abs(loop['start'] - m['start']) , abs(loop['end'] - m['start']))
                        if dis > 100000:
                            continue
                        dis_score[dis] = m['score']
#                        dis_score[dis] = 0
#                    enhancer_score = 1
                    try:    
                        enhancer_score = dis_score[min(dis_score.keys())]
                    except:
                        enhancer_score = 0
                else:
#                    enhancer_score = 0
                    enhancer_score = 0
                gene_loop_matrix[i][c] = enhancer_score
                
        else:
            gene_loop_matrix[i] = np.array([-1,-1,-1,-1])
        
   
gene_loop_matrix = stats.zscore(gene_loop_matrix, axis=1, ddof=1)
gene_loop_matrix[np.isnan(gene_loop_matrix)] = -1


    
    
#
#matrix_all = np.hstack((matrix_b[:,:-1],matrix_peaks))    

#----------------------------gene_expression-loop_strength--------------------------------------------

f_CCS = np.load('C:\\Users\\xxli\\Desktop\\loop\\loop_strength_20K_hluo\\CCS_loop_ave_strength.npz')
f_NT5 = np.load('C:\\Users\\xxli\\Desktop\\loop\\loop_strength_20K_hluo\\NT5_loop_ave_strength.npz')
f_NT6 = np.load('C:\\Users\\xxli\\Desktop\\loop\\loop_strength_20K_hluo\\NT6_loop_ave_strength.npz')
f_fESC = np.load('C:\\Users\\xxli\\Desktop\\loop\\loop_strength_20K_hluo\\fESC_loop_ave_strength.npz')

#for i in [f_CCS,f_NT5,f_NT6,f_fESC]:
#    for j in i:
#        try:
#            i[j] = float(i[j])
#        except:
#            print i[j]
#            i[j] = 0

peak_CCS = np.loadtxt('C:\\Users\\xxli\\Desktop\\IDR\\selected_peaks\\idr_select_peaks_chr\\CCS_IDR_Selected_enhancer.narrowPeak' , usecols = ((0,1,2,4)) , dtype = peaktype)
peak_NT5 = np.loadtxt('C:\\Users\\xxli\\Desktop\\IDR\\selected_peaks\\idr_select_peaks_chr\\NT5_IDR_Selected_enhancer.narrowPeak' , usecols = ((0,1,2,4)) , dtype = peaktype)
peak_NT6 = np.loadtxt('C:\\Users\\xxli\\Desktop\\IDR\\selected_peaks\\idr_select_peaks_chr\\NT6_IDR_Selected_enhancer.narrowPeak' , usecols = ((0,1,2,4)) , dtype = peaktype)
peak_fESC = np.loadtxt('C:\\Users\\xxli\\Desktop\\IDR\\selected_peaks\\idr_select_peaks_chr\\fESC_IDR_Selected_enhancer.narrowPeak' , usecols = ((0,1,2,4)) , dtype = peaktype)

peaks = {0:peak_CCS , 1:peak_NT5 , 2:peak_NT6 , 3:peak_fESC}



loop_matrix = np.zeros((len(f_CCS),),dtype = looptype)
n = 0
for i in f_CCS:
    loop_matrix[n]['chr'] = i[0]
    loop_matrix[n]['start'] = i[1]
    loop_matrix[n]['end'] = i[2]
    loop_matrix[n]['CCS'] = f_CCS[i]
    loop_matrix[n]['NT5'] = f_NT5[i]
    loop_matrix[n]['NT6'] = f_NT6[i]
    loop_matrix[n]['fESC'] = f_fESC[i]
    n += 1
    
gene_loop_matrix = np.zeros((len(sort_gene),4))    

for i in range(len(sort_gene)):
    index = list(matrix_b[:,:-1][i]).index(matrix_b[:,:-1][i].max())
    chro = sort_gene[i][1].lstrip('chr')
    gene_start = int(sort_gene[i][2])
    mask = (gene_start >= loop_matrix[loop_matrix['chr'] == chro]['start']) & (gene_start <= loop_matrix[loop_matrix['chr'] == chro]['end'])
    overlap = loop_matrix[loop_matrix['chr'] == chro][mask]
    if overlap.size == 0:
        continue
    else:
        n = overlap.size 
        peak_loop = {}
        for j in overlap:
            if min(abs(int(j['start']) - int(sort_gene[i][-1])),abs(int(j['end']) - int(sort_gene[i][-1]))) > 500000:
                continue
            mask1 = (j['start'] <=  peaks[index][peaks[index]['chr'] == 'chr' + chro]['end']) & \
                    (j['end'] >=  peaks[index][peaks[index]['chr'] == 'chr' + chro]['start'])
            overlap1 = peaks[index][peaks[index]['chr'] == 'chr' + chro][mask1]
            if overlap1.size == 0:
                peak_loop[0] = j
            else:
                dis_score = {}
                for k in overlap1:
                    dis = min(abs(j['start'] - k['start']) , abs(j['end'] - k['start']))
                    dis_score[dis] = k['score']
                dis = min(dis_score.keys())
                score = dis_score[dis]
                peak_loop[score] = j
        if peak_loop != {}:  
            loop = peak_loop[max(peak_loop.keys())]
            gene_loop_matrix[i][0] = loop['CCS']
            gene_loop_matrix[i][1] = loop['NT5']
            gene_loop_matrix[i][2] = loop['NT6']
            gene_loop_matrix[i][3] = loop['fESC']
        else:
            gene_loop_matrix[i][0] = 0
            gene_loop_matrix[i][1] = 0
            gene_loop_matrix[i][2] = 0
            gene_loop_matrix[i][3] = 0
            
gene_loop_matrix = stats.zscore(gene_loop_matrix, axis=1, ddof=1)
gene_loop_matrix[np.isnan(gene_loop_matrix)] = 0
        
        
        


#----------------------------------Merge_matrix---------------------------------------


matrix_all = np.hstack((matrix_b[:,:-1],matrix_peaks,gene_loop_matrix,matrix_b[:,-1].reshape(len(matrix_b),1)))
matrix_all[np.isnan(matrix_all)] = 0
#matrix_peaks_1 = []
#for i in matrix_peaks:
#    matrix_peaks_1.append(tuple(i))
#    
#matrix_peaks = np.array(matrix_peaks_1,dtype = datatype_2)

matrix_all_0 = np.array([a for a in matrix_all if a[-1] == 0.0])
matrix_all_1 = K_means_cluster(matrix_all_0,matrix_all_0[:,4:8],6)

for i in list(set(matrix_b[:,-1]))[1:]:
    matrix_all_0 = np.array([a for a in matrix_all if a[-1] == i])
    matrix_all_0 = K_means_cluster(matrix_all_0,matrix_all_0[:,4:8],6)
    matrix_all_1 = np.vstack((matrix_all_1,matrix_all_0))
        

#----------------------------------plot-----------------------------------------------
#n = range(16)
#pp = PdfPages('C:\\Users\\xxli\\Desktop\\gene_expression\\gene_expression_individual_classify_16.pdf')

left, bottom, width, height = 0.1, 0.1, 0.60, 0.80
size_axes = [left, bottom, width, height]
fig = plt.figure(figsize = (12, 12))
ax = fig.add_axes(size_axes)
#matrix_1 = np.hstack((matrix_b[:,:-1] , gene_loop_matrix_1))
matrix_1 = matrix_all_1[:,:-2]
vmax = np.percentile(matrix_1,95)
vmin = np.percentile(matrix_1,5)
im = ax.imshow(matrix_1,vmax=vmax,vmin=vmin,cmap='bwr',aspect = 'auto')
plt.colorbar(im)
x = ['CCS','NT5','NT6','fESC','CCS','NT5','NT6','fESC','CCS','NT5','NT6','fESC']
ax.set_xticks(np.arange(len(x)))
ax.set_xticklabels(x,fontsize = 15)
plt.title('Gene Expression,gene_numbers:' + str(len(matrix_1)),fontsize = 20)
#pp.savefig(fig)
#plt.close()

#-----------------------------------select_genes----------------------------------------

classify_gene = {}

n = 0
tmp = []
c = 3 
for i in matrix_all_1:
    if i[12] == c:
        if (max(i[4:8]) - min(i[4:8]) <= 0.5) & (min(i[8:12]) == i[8]):
            tmp.append(list(i) + [n])
    n += 1    
gene = {}
for i in tmp:
    values = list(i[:4])
    g = list(gene_fpkm.keys())[list(gene_fpkm.values()).index(values)]
    gene[g] = i
    
classify_gene[c] = {}
for i in gene:
    if gene[i][4:12] == [-1,-1,-1,-1,-1,-1,-1,-1]:
        continue
    else:
        classify_gene[c][i] = gene[i]
    
