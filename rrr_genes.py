# -*- coding: utf-8 -*-
"""
Created on Thu Oct 11 20:05:04 2018

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
from matplotlib.backends.backend_pdf import PdfPages

## Matplotlib Settings
matplotlib.rcParams['xtick.direction'] = 'out'
matplotlib.rcParams['ytick.direction'] = 'out'


# Our Own Color Map
my_cmap = plt.get_cmap('bwr')
my_cmap.set_bad('#D3D3D3')

cell = {'CCS':0 , 'NT5':1 , 'NT6':2 , 'fESC':3}

datatype = ({'names':['gene_id' , 'classify'],
             'formats':['S64' , 'S8' ]})
fpkm_type = ({'names':['gene_id' , 'gene_name' , 'ref' , 'strand' , 'start' , 'end' , 'CCS_1' , 'CCS_2' , 'CCS_3' , 'fESC_1' , 'fESC_2' , 'fESC_3' , 'NT5_1'	, 'NT5_2' ,	'NT5_3' , 'NT5_4' , 'NT6_1' , 'NT6_2' , 'NT6_3'],
             'formats':['S64' , 'S64' , 'S8' , 'S8' , np.int , np.int , np.float , np.float , np.float , np.float , np.float , np.float , np.float , np.float , np.float , np.float , np.float , np.float , np.float ]})
gtf_type_1 = ({'names':['gene_id' , 'gene_name' , 'chr' , 'start' , 'end'],
             'formats':['S64' , 'S64' , 'S8' , np.int , np.int]})
promoter_type = ({'names':['gene_name' , 'chr' , 'start' , 'end'],
             'formats':['S64' , 'S8' , np.int , np.int ]})
loop_type = ({'names':['chr' , 'start' , 'end'],
             'formats':['S8' , np.int , np.int]})
peak_type = ({'names':['chr' , 'start' , 'end' , 'score'],
             'formats':['<S8' , np.int , np.int , np.float]})

#-------------------------------------------------------Function--------------------------------------------------------------------

def Sort(a , s_1 , s_2):
    '''
    a: list of needed to be sorted
    '''
    a.sort(key = lambda x:(x[s_1],x[s_2]))
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
            peak_new.append((i[0] , x , y))
            a = list(a)
            for j in overlap:
                a.remove(j)
            a = np.array(a, dtype = peak_type)
        else:
            continue

    peak_new = Sort(peak_new , 0 ,1)
    peak_new = np.array(peak_new , dtype = peak_type)
    return peak_new


def union_loop(a,peak_type):
    '''
    a: np.array of peaks ,dtype = peak_type
    '''
    peak_new = []
    a = np.array(a , dtype = peak_type)
    for i in a:
        s_s = i['start'] - math.sqrt(2) * 40000
        s_e = i['start'] + math.sqrt(2) * 40000
        e_s = i['end'] - math.sqrt(2) * 40000
        e_e = i['end'] + math.sqrt(2) * 40000
        mask = (s_s <= a[a['chr'] == i['chr']]['start']) & (s_e >= a[a['chr'] == i['chr']]['start']) & \
               (e_s <= a[a['chr'] == i['chr']]['end']) & (e_e >= a[a['chr'] == i['chr']]['end'])
        overlap = a[a['chr'] == i['chr']][mask]
        if overlap.size != 0:
            x = overlap['start'].min()
            y = overlap['end'].max()
            peak_new.append((i[0] , x , y))
            a = list(a)
            for j in overlap:
                a.remove(j)
            a = np.array(a, dtype = peak_type)
        else:
            continue

    peak_new = Sort(peak_new , 0 ,1)
    peak_new = np.array(peak_new , dtype = peak_type)
    return peak_new




def fpkm_matrix(gene):
    matrix = np.zeros((len(gene) , 4))
    for i in range(len(gene)):
        if gene[i][4] != 0:
            ccs = math.log(gene[i][4] , 2)
        else:
            ccs = -12
        if gene[i][5] != 0:
            nt5 = math.log(gene[i][5] , 2)
        else:
            nt5 = -12
        if gene[i][6] != 0:
            nt6 = math.log(gene[i][6] , 2)
        else:
            nt6 = -12
        if gene[i][7] != 0:
            fesc = math.log(gene[i][7] , 2)
        else:
            fesc = -12
        matrix[i , 0] = ccs
        matrix[i , 1] = nt5
        matrix[i , 2] = nt6
        matrix[i , 3] = fesc
        
    return matrix


def matrix_index(matrix):
    index = np.zeros((len(matrix),1))
    for i in range(len(matrix)):
        index[i][0] = i
    matrix_1 = np.hstack((matrix , index))
    return matrix_1

def K_means_cluster(matrix_0,matrix,n):
    kmeans = KMeans(n_clusters=n, random_state=0).fit(matrix)
    a = kmeans.labels_
    a = a.reshape(len(a),1)

    matrix_1 = np.array(np.hstack((matrix_0,a)))
    matrix_1 = np.array(sorted(matrix_1, key=lambda x:x[-1]))
    return matrix_1

#--------------------------------------------------rrr_gene--------------------------------------------------------------------------
f = np.loadtxt('D:\\workspace\\rrr\\rrr_gene.txt' , skiprows = 1 , dtype = 'S64')

gene_fpkm = np.loadtxt('F:\\xxli\\data_analysis\\BDF1\\RNA_seq\\gene_expression\\all_genes.txt' , skiprows = 1 , dtype = fpkm_type)

gtf_1 = []
for i in gene_fpkm:
    if i['strand'] == '+':
        gtf_1.append((i['gene_id'].split(".")[0] , i[1] , i[2] , i[3] , i[4] , i[5] , i[6] , i[7] , i[8] , i[9] , i[10] , i[11] , i[12] , i[13] , i[14] , i[15] , i[16] , i[17] , i[18]))
    elif i['strand'] == '-':
        gtf_1.append((i['gene_id'].split(".")[0] , i[1] , i[2] , i[3] , i[5] , i[4] , i[6] , i[7] , i[8] , i[9] , i[10] , i[11] , i[12] , i[13] , i[14] , i[15] , i[16] , i[17] , i[18]))
    else:
        print 1
        continue
    
gene_fpkm = np.array(gtf_1 , dtype = fpkm_type)


rrr = []
for i in f:
    mask = (i == gene_fpkm['gene_name'])
    overlap = gene_fpkm[mask]
    if overlap.size == 0:
        continue
    elif overlap.size >1:
        print overlap.size
    ccs = (overlap[0]['CCS_1'] + overlap[0]['CCS_2'] + overlap[0]['CCS_3']) / 3
    nt5 = (overlap[0]['NT5_1'] + overlap[0]['NT5_2'] + overlap[0]['NT5_3'] + overlap[0]['NT5_4']) / 4
    nt6 = (overlap[0]['NT6_1'] + overlap[0]['NT6_2'] + overlap[0]['NT6_3']) / 3
    fesc = (overlap[0]['fESC_1'] + overlap[0]['fESC_2'] + overlap[0]['fESC_3']) / 3
    rrr.append((overlap[0]['gene_name'] , overlap[0]['ref'] , overlap[0]['start'] , overlap[0]['end'] , ccs , nt5 , nt6 , fesc))

        
#-------------------------------------------------------diff_gene_fpkm---------------------------------------------------------------     
        
vs = ['CCS_vs_NT5' , 'CCS_vs_NT6' , 'fESC_vs_CCS' , 'fESC_vs_NT5' , 'fESC_vs_NT6' , 'NT5_vs_NT6']
diff_gene = []
for c in vs:
    f = csv.reader(open('F:\\xxli\\data_analysis\\BDF1\\RNA_seq\\difference_expression\\baseMean_500_fc_1.5\\Filtered_new_' + c + '.csv'))
    n = 0
    for i in f:
        n += 1
        if n == 1:
            continue
        diff_gene.append(i[7].strip().split("\'")[1])
        
diff_gene = set(diff_gene)

diff_gene_age = {'older':[] , 'younger':[]}
for i in gene_age:
    for j in gene_age[i]:
        if j[0] in diff_gene:
            diff_gene_age[i].append(j)
        else:
            continue
        
gene_age = diff_gene_age

#
        
#----------------------------------------------gene_expression matrix----------------------------------------------------------------



matrix_a = fpkm_matrix(rrr)


##-------------------------------------------add matrix_index----------------------------------------------------------------------


matrix_a = matrix_index(matrix_a)

##-------------------------------------------K-means cluster-----------------------------------------------------------------------



matrix_a = K_means_cluster(matrix_a , matrix_a[:,:4] , 6)

        


#----------------------------------------------promoter----------------------------------------------------------------------------------
    
    

union_promoters = []
for c in cell:
    f = np.loadtxt('D:\\workspace\\IDR\\selected_peaks\\idr_select_peaks_chr\\' + c + '_IDR_Selected_promoter.narrowPeak' , usecols = ( -1 , 0 , 1 , 2 ) ,
                   dtype = promoter_type)
    for i in f:
        union_promoters.append((i[0].strip().split("\"")[1] , i[1] , i[2] , i[3]))
        
union_promoters = Sort(union_promoters , 1 , 2)
union_promoters = union_interval(union_promoters , promoter_type)        
    
p_percent = {'CCS':0 , 'NT5':0 , 'NT6':0 , 'fESC':0}
for c in cell:
    f = np.loadtxt('D:\\workspace\\IDR\\selected_peaks\\idr_select_peaks_chr\\' + c + '_IDR_Selected_promoter.narrowPeak' , usecols = ( -1 , 0 , 1 , 2 ) ,
                   dtype = promoter_type)

    n = 0
    for i in rrr:
        gene_name = i[0]
        chro = i[1]
        s = min(i[2] , i[3])
        e = max(i[2] , i[3])
        mask = (gene_name == union_promoters['gene_name'])
        overlap = union_promoters[mask]
        if overlap.size == 0:
            continue
        u_s = min(overlap['start'])
        u_e = max(overlap['end'])
        mask1 = (u_s <= f[f['chr'] == chro]['end']) & (u_e >= f[f['chr'] == chro]['start'])
        overlap1 = f[f['chr'] == chro][mask1]
        if overlap1.size != 0:
            n += 1
        elif overlap1.size >1 :
            print overlap1 , overlap
    percent = n / len(rrr)
    p_percent[c] = percent    
    
    
#----------------------------------------------loop-------------------------------------------------------------------------


loop_count = {'CCS':[] , 'NT5':[] , 'NT6':[] , 'fESC':[]}
loop_count_1 = {'CCS':[] , 'NT5':[] , 'NT6':[] , 'fESC':[]}
loop_count_2 = {'CCS':0 , 'NT5':0 , 'NT6':0 , 'fESC':0}

for i in rrr:
    chro = i[1].lstrip('chr')
    s = min(i[2] , i[3])
    e = max(i[2] , i[3])
    for c in cell:
        f = np.loadtxt('F:\\xxli\\data_analysis\\BDF1\HiC\\Loop_newSNP\\Clustered_' + c + '_hiccups_loops.txt' , skiprows = 1,
                       usecols = (0 , 1 , 2) , dtype = loop_type)
        mask1 = (s <= f[f['chr'] == chro]['end']) & (e >= f[f['chr'] == chro]['start'])
        overlap1 = f[f['chr'] == chro][mask1]
        if overlap1.size == 0:
            continue
        for k in overlap1:
            if min(abs(s - k['start']),abs(e - k['end'])) > 100000:
                continue
            loop_count[c].append(k)
            
for i in loop_count:
    loop_count_1[i] = union_loop(loop_count[i] , loop_type)
        
for i in loop_count_1:
    loop_count_2[i] = len(loop_count_1[i])
    
#------------------------------------------------loop_enhancer---------------------------------------------------------------------------


union_loops = []
union_peaks = []
peaks = {}


for c in cell:
    f1 = np.loadtxt('F:\\xxli\\data_analysis\\BDF1\HiC\\Loop_newSNP\\Clustered_' + c + '_hiccups_loops.txt' , skiprows = 1,
                    usecols = (0 , 1 , 2) , dtype = loop_type)
    f2 = np.loadtxt('D:\\workspace\\IDR\\selected_peaks\\idr_select_peaks_chr\\' + c + '_IDR_Selected_enhancer.narrowPeak' ,
                    usecols = (0 , 1 , 2) , dtype = loop_type)

    
    for i in f1:
        union_loops.append(('chr' + i[0] , i[1] , i[2]))
    for i in f2:
        union_peaks.append(i)
    peaks[c] = f2
    
union_loops = Sort(union_loops , 0 , 1)
union_peaks = Sort(union_peaks , 0 , 1)   

union_loops = union_loop(union_loops , loop_type)
union_peaks = union_interval(union_peaks , loop_type)

e_percent = {'CCS':0 , 'NT5':0 , 'NT6':0 , 'fESC':0}
for i in rrr:
    chro = i[1]
    s = min(i[2] , i[3])
    e = max(i[2] , i[3])
    mask1 = (s <= union_loops[union_loops['chr'] == chro]['end']) & (e >= union_loops[union_loops['chr'] == chro]['start'])
    overlap1 = union_loops[union_loops['chr'] == chro][mask1]
    if overlap1.size == 0:
        continue
    loop_enhancer = {}
    for k in overlap1:
        if min(abs(s - k['start']),abs(e - k['end'])) > 100000:
            continue
        dis = {}
        mask2 = (k['start'] <= union_peaks[union_peaks['chr'] == chro]['end']) & (k['end'] >= union_peaks[union_peaks['chr'] == chro]['start'])
        overlap2 = union_peaks[union_peaks['chr'] == chro][mask2]
        if overlap2.size == 0:
            continue
        for l in overlap2:
            d = min(abs(l['start'] - k['start']),abs(l['end'] - k['end']))
            if d > 100000:
                continue
            dis[d] = l
        if dis != {}:
            d = min(dis.keys())
            loop_enhancer[d] = (k , dis[d])
        else:
            continue
    if loop_enhancer != {}:
        keys = min(loop_enhancer.keys())
        enhancer = loop_enhancer[keys][1]
        for c in cell:
            mask = (enhancer[1] <= peaks[c][peaks[c]['chr'] == chro]['end']) & (enhancer[2] >= peaks[c][peaks[c]['chr'] == chro]['start'])
            overlap = peaks[c][peaks[c]['chr'] == chro][mask]
            if overlap.size != 0:
                e_percent[c] += 1
            
for i in e_percent:
    e_percent[i] = e_percent[i]/len(rrr)


#--------------------------------------------Draw heatmap------------------------------------------------------------------------
        
left, bottom, width, height = 0.1, 0.1, 0.60, 0.80
size_axes = [left, bottom, width, height]
fig = plt.figure(figsize = (12, 12))
ax = fig.add_axes(size_axes)
#matrix_1 = np.hstack((matrix_b[:,:-1] , gene_loop_matrix_1))
matrix = matrix_a[:,:-2]
vmax = np.percentile(matrix,95)
vmin = np.percentile(matrix,5)
im = ax.imshow(matrix,vmax=vmax,vmin=vmin,cmap=my_cmap,aspect = 'auto')
plt.colorbar(im)
x = ['CCS','NT5','NT6','fESC']
ax.set_xticks(np.arange(len(x)))
ax.set_xticklabels(x,fontsize = 15)
plt.title('Reprogramming Resistant Regions FPKM(log2),gene_numbers:' + str(len(matrix)),fontsize = 20)


#pp.savefig(fig)
            


#-------------------------------------------Draw bar_plot--------------------------------------------------------------------------
def Plot_bar(p , title , ylabel , classify , photo_name):
    left, bottom, width, height = 0.1, 0.1, 0.80, 0.80
    size_axes = [left, bottom, width, height]
    fig = plt.figure(figsize = (12, 12))
    ax = fig.add_axes(size_axes)
    #matrix_1 = np.hstack((matrix_b[:,:-1] , gene_loop_matrix_1))
    x = range(4)
    y = [p[title]['CCS'] , p[title]['NT5'] , p[title]['NT6'] , p[title]['fESC']]
    plt.bar(x,y)
    plt.xlabel('Cells' , fontsize = 20)
    plt.ylabel(ylabel, fontsize = 20)
    plt.title(title[0].upper() + title[1:] + ' Genes '+ classify , fontsize = 20)
    plt.ylim(0 , 0.8)
    xticks = ['CCS','NT5','NT6','fESC']
    ax.set_xticks(np.arange(len(x)))
    ax.set_xticklabels(xticks,fontsize = 20)
    plt.savefig('D:\\workspace\\gene_age\\' + photo_name)
    return fig
classify = 'classify_5'
fig1 = Plot_bar(p_percent , 'older' , 'Proportion of gained open promoter' ,  classify , classify + '\\older_genes_promoter_' + classify + '.png')
fig2 = Plot_bar(p_percent , 'younger' , 'Proportion of gained open promoter' , classify , classify + '\\younger_genes_promoter_' + classify + '.png')
fig3 = Plot_bar(e_percent , 'older' , 'Proportion of gained loop enhancer' , classify , classify + '\\older_genes_loop_enhancer' + classify + '.png')
fig4 = Plot_bar(e_percent , 'younger' , 'Proportion of gained loop enhancer' , classify , classify + '\\younger_genes_loop_enhancer' + classify + '.png')
        



def Plot_bar(p , title , ylabel , photo_name):
    left, bottom, width, height = 0.1, 0.1, 0.80, 0.80
    size_axes = [left, bottom, width, height]
    fig = plt.figure(figsize = (12, 12))
    ax = fig.add_axes(size_axes)
    #matrix_1 = np.hstack((matrix_b[:,:-1] , gene_loop_matrix_1))
    x = range(4)
    y = [p[title]['CCS'] , p[title]['NT5'] , p[title]['NT6'] , p[title]['fESC']]
    plt.bar(x,y)
    plt.xlabel('Cells' , fontsize = 20)
    plt.ylabel(ylabel, fontsize = 20)
    plt.title(title[0].upper() + title[1:] + ' Genes', fontsize = 20)
    plt.ylim(0 , 9000)
    xticks = ['CCS','NT5','NT6','fESC']
    ax.set_xticks(np.arange(len(x)))
    ax.set_xticklabels(xticks,fontsize = 20)
    plt.savefig('D:\\workspace\\gene_age\\' + photo_name)
    return fig

fig1 = Plot_bar(p_percent , 'older' , 'Proportion of gained open promoter' , 'older_genes_promoter.png')
fig2 = Plot_bar(p_percent , 'younger' , 'Proportion of gained open promoter' ,'younger_genes_promoter.png')
fig3 = Plot_bar(e_percent , 'older' , 'Proportion of gained loop enhancer' , 'older_genes_loop_enhancer.png')
fig4 = Plot_bar(e_percent , 'younger' , 'Proportion of gained loop enhancer' ,'younger_genes_loop_enhancer.png')
fig5 = Plot_bar(loop_count_2 , 'older' , 'Proportion of gained open Loop' , 'older_genes_Loop.png')        
fig6 = Plot_bar(loop_count_2 , 'younger' , 'Proportion of gained open Loop' , 'younger_genes_Loop.png')



def Plot_bar(p , title , ylabel , photo_name):
    left, bottom, width, height = 0.1, 0.1, 0.80, 0.80
    size_axes = [left, bottom, width, height]
    fig = plt.figure(figsize = (12, 12))
    ax = fig.add_axes(size_axes)
    #matrix_1 = np.hstack((matrix_b[:,:-1] , gene_loop_matrix_1))
    x = range(4)
    y = [p['CCS'] , p['NT5'] , p['NT6'] , p['fESC']]
    plt.bar(x,y)
    plt.xlabel('Cells' , fontsize = 20)
    plt.ylabel(ylabel, fontsize = 20)
    plt.title(title[0].upper() + title[1:] + ' Genes', fontsize = 20)
    plt.ylim(0 , 100)
    xticks = ['CCS','NT5','NT6','fESC']
    ax.set_xticks(np.arange(len(x)))
    ax.set_xticklabels(xticks,fontsize = 20)
    plt.savefig('D:\\workspace\\rrr' + photo_name)
    return fig

fig1 = Plot_bar(p_percent , 'Reprogramming resistant regions genes' , 'Proportion of gained open promoter' , '\\rrr_genes_promoter.png')
fig2 = Plot_bar(e_percent , 'Reprogramming resistant regions genes' , 'Proportion of gained loop enhancer' , '\\rrr_genes_loop_enhancer.png')
fig3 = Plot_bar(loop_count_2, 'Reprogramming resistant regions genes' , 'Proportion of gained Loop' , '\\rrr_genes_Loop.png')
