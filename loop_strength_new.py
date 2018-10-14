from __future__ import division
import numpy as np
import sys, cPickle
import os
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

def Load_loop(peak_file):
    """
    """
    loopType = np.dtype({'names':['chr','S','E'],
                     'formats':['S2', np.int, np.int]})

    loop = np.loadtxt(peak_file,dtype = loopType,skiprows = 1, usecols = (0,1,2))
    
    return loop
def Load_gene(gene_file):
    """
    """
    geneType = np.dtype({'names':['chr','S','E','FPKM','genename'],
                     'formats':['S8', np.int, np.int, np.float,'S30']})
    gene = np.loadtxt(gene_file , dtype = geneType, skiprows = 1, usecols =(2,4,5,7,1))

    return gene

def Get_Loop_Strength(inter,start,end,Len):
    """
    """
    maskA = (inter['bin1'] == start) & (inter['bin2'] == end)
    maskB = (inter['bin1'] >= start - Len) & (inter['bin1'] <= start + Len) & \
            (inter['bin2'] >= end - Len) & (inter['bin2'] <= end + Len)
#    maskC = (inter['bin1'] >= start - 4 * Len) & (inter['bin1'] <= start - 2 * Len) & \
#            (inter['bin2'] >= end + 2 * Len) & (inter['bin2'] <= end + 4 * Len)
    maskD = (inter['bin1'] >= start + 2 * Len) & (inter['bin1'] <= start + 4 * Len) & \
            (inter['bin2'] >= end - 4 * Len) & (inter['bin2'] <= end - 2 * Len)

    A = inter[maskA]
    B = inter[maskB]
#    C = inter[maskC]
    D = inter[maskD]
    n_B = len(np.nonzero(B)[0])
    n_D = len(np.nonzero(D)[0])
    A_S = A['IF']
    try:
        B_S = B['IF'].sum() / n_B
    except:
        B_S = 0
#        C_S = C['IF'].sum() / ((Len + 2) **2)
    try:
        D_S = D['IF'].sum() / n_D
    except:
        D_S = 1
    return A_S /D_S ,B_S /D_S


def read_npz(chro,HiCData):
    """
    HiCData in every chro
    """
    data = {}
    inter = HiCData[chro]
    max_bin = max(inter['bin1'].max(), inter['bin2'].max())
    for i in range(max_bin+1):
        data[i] = []
    for i in inter:
        data[i['bin1']].append((i['bin1'],i['bin2'],i['IF']))

    dtype = np.dtype({'names':['bin1','bin2','IF'],
                      'formats':[np.int, np.int, np.float]})
    for k,v in data.items():
        v = np.array(v,dtype = dtype)
        data[k] = v

    return data

def Get_inter(Index_npz,start,end):
    """
    """
#    Index_npz = read_npz(lib)
    inter = Index_npz[start]
    for i in range(end-start):
        inter = np.hstack((inter,Index_npz[start+i+1]))
    
    return inter

def Get_loops(cell):
    """
    cell : String ,The name of cells  e.g.'fESC','ESC' , 'NT5' , 'NT6' , 'CCS'
    get the loops' location as a dictionary
    """
    loops = {}
    fileSource = os.path.join(loopFolder, cell + '_cluster_filter-P_' + res + '.txt')
    Loop = np.loadtxt(fileSource, usecols = (0,1,2) , dtype = loopType, skiprows = 0)
    for i in Loop:
        if i['E'] - i['S'] >= 150000:
            if i['chr'] not in loops.keys():
                loops[i['chr']] = []
                loops[i['chr']].append(i)
            else:
                loops[i['chr']].append(i)
        else:
            continue
    for k,v in loops.items():
        v = np.array(v,loopType)
        loops[k] = v
    return loops

def get_union_loop(cell,Len):
    """
    cell : List ,The List of the cells'name , e.g.['ESC' , 'NT5' , 'NT6' , 'CCS']
    Len : The window of loop_overlap
    """
    loops = {}
    union = Get_loops(cell[0])
    loops[0] = Get_loops(cell[0])
    for i in range(1,len(cell)):
        loops[i] = Get_loops(cell[i])           
        for g in loops[i]:
            for j in loops[i][g]:
                mask = (j['S'] >=  union[j['chr']]['S'] - Len) & (j['S'] <=  union[j['chr']]['S'] + Len) & \
                       (j['E'] >=  union[j['chr']]['E'] - Len) & (j['E'] <=  union[j['chr']]['E'] + Len)
                tmp = union[j['chr']][mask]
                if tmp.size == 0:
                    t = list(union[j[0]])
                    t.append(j)
                    union[j[0]] = np.array(t , dtype = loopType)
                else:
                    continue
    return union

def get_common_loop(cell,Len):
    common = Get_loops(cell[0])
    common_tmp = []
    for i in range(1,len(cell)):
        loops = Get_loops(cell[i])
        for g in loops:
            for j in loops[g]:
                mask = (j['S'] >=  common[j['chr']]['S'] - Len) & (j['S'] <=  common[j['chr']]['S'] + Len) & \
                       (j['E'] >=  common[j['chr']]['E'] - Len) & (j['E'] <=  common[j['chr']]['E'] + Len)
                tmp = common[j['chr']][mask]
                if tmp.size == 0:
                   continue
                else:
                   common_tmp.append(tmp[0])
            common[g] = np.array(common_tmp , dtype = loopType)
            common_tmp = []
    return common


HiCFolder = '/public/home/xxli/data/BDF1/HiC/runHiC/workspace/Raw-GRCm38_68_chr'
loopFolder = '/public/home/xxli/data/BDF1/HiC/loop/cluster_filter_loop'
geneFolder = '/public/home/xxli/data/BDF1/HiC/loop/RNA'
cell = ['fESC','ESC','NT5','NT6','CCS']
enzyme = 'MboI'
res = '20K'
ResHiC = 20000
chrom = ['1', '2', '3', '4', '5', '6', '7', '8', '9','10','11', '12', '13', '14', '15', '16', '17', '18', '19', 'X']
loopType = np.dtype({'names':['chr','S','E'],
                     'formats':['S2', np.int, np.int]})
uniontype = np.dtype({'names':['chr','S','E','IF_fESC','IF_ESC','IF_NT5','IF_NT6','IF_CCS'],
                          'formats':['S2' , np.int , np.int, np.float,np.float,np.float,np.float,np.float]})
#---------------------------------union_loop--------------------------------
union = get_union_loop(cell , 0)
#--------------------------------Get_loop_strength--------------------------
final_loops = {}
final_loops_1 = {}
for i in range(5):
    ResLabel = str(ResHiC//1000) + 'K'
    Pre = '-'.join([cell[i], enzyme, 'allReps-filtered', ResLabel])
    HiCFil = Pre + '-sparse.npz'
    HiCSource = os.path.join(HiCFolder, HiCFil)
    lib = np.load(HiCSource)
    Out = open('/public/home/xxli/data/BDF1/HiC/loop/' + cell[i] + '_loop_point_strength.npz','wb')
    Out1 = open('/public/home/xxli/data/BDF1/HiC/loop/' + cell[i] + '_loop_ave_strength.npz','wb')
    final_loops[cell[i]] = {}
    final_loops_1[cell[i]] = {}
    for g in chrom:
        Index_npz = read_npz(g,lib)
        for l in union[g]:
            chro = l['chr']
            start = int(l['S'] / ResHiC)
            end = int(l['E'] / ResHiC)
            inter = Get_inter(Index_npz,start-10,start+10)
            (point,ave) = Get_Loop_Strength(inter,start,end,1)
            final_loops[cell[i]][tuple(l)] = point
            final_loops_1[cell[i]][tuple(l)] = ave
    cPickle.dump(final_loops[cell[i]],Out,2)
    cPickle.dump(final_loops_1[cell[i]],Out1,2)
    Out.close()
    Out1.close()