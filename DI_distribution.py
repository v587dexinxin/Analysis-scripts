# -*- coding: utf-8 -*-
"""
Created on Wed Jan 03 09:34:21 2018

@author: xxli
"""


from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
from scipy import sparse
import math
import matplotlib
from matplotlib.backends.backend_pdf import PdfPages
import pandas as pd
import seaborn as sns
from tadlib.hitad.chromLev import Chrom
from palettable.colorbrewer.qualitative import Dark2_8
from matplotlib.colors import LinearSegmentedColormap


datafolder = 'C:\\Users\\xxli\\Desktop\\scatter\\matrix_100K_NT5\\'
outfolder = 'C:\\Users\\xxli\\Desktop\\scatter\\Plot_heatmap_DI\\NT5_30K\\'
chrom=['1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19','X']
hexcolors = Dark2_8.hex_colors
interval = 300
def loadfile(filename):
    data = np.loadtxt(datafolder + filename,skiprows=1,dtype=np.str)
    matrix = np.array(data[:,2:],dtype=np.float)
    data = sparse.triu(sparse.coo_matrix(matrix))
    npz_data = np.array(zip(data.row,data.col,data.data) , dtype=[ ('bin1', '<i4'), ('bin2', '<i4'), ('IF', '<f8')])
    return matrix , npz_data
    
def properU(pos):
    """
    Express a genomic position in a proper unit (KB, MB, or both).
    
    """
    i_part = int(pos) // 1000000 # Integer Part
    d_part = (int(pos) % 1000000) // 1000 # Decimal Part
    
    if (i_part > 0) and (d_part > 0):
        return ''.join([str(i_part), 'M', str(d_part), 'K'])
    elif (i_part == 0):
        return ''.join([str(d_part), 'K'])
    else:
        return ''.join([str(i_part), 'M'])
    
def plot_matrix_DI(matrix , diff , ff , mm , al , filename , r):
    '''
    plot a compartment and heatmap mix fig
    matrix : the matrix of undistinguishable HiC data matrix
    diff: the difference value of DI between ff and mm
    ff: the value of ff DI
    mm: the value of mm DI
    al: the value of undistinguishable data DI
    filename: the name of pdfname
    r: Resolution
    '''
    ##ax size
    left, width = 0.15, 0.6
    bottom, height = 0.05, 0.6
    bottom_h = bottom + width
    size_heatmap = [left, bottom, width, height]
    size_signal = [left, bottom_h, width, 0.075]
    size_colorbar = [left + width + 0.08, bottom, width/20, height]
    

    pp = PdfPages(outfolder + filename)
    for i in range(int(len(matrix)/interval +1)):
        ##heatmap
        a = i*interval
        b = (i+1)*interval
        ticks = list(np.linspace(0, 30, 5).astype(float))
        pos = [((a + w*10) * r) for w in ticks]
        labels = [properU(p) for p in pos]
        fig = plt.figure(figsize=(12,12))
        ax = fig.add_axes(size_heatmap)
        nonzero = matrix[np.nonzero(matrix[a:b,a:b])]
        vmax = np.percentile(nonzero,90)
        lenth = len(matrix[a:b,a:b]) * r / 1000000
        sc = ax.imshow(matrix[a:b,a:b], cmap = LinearSegmentedColormap.from_list('interaction',['#FFFFFF','#CD0000']), vmax = vmax, vmin = 0,aspect = 'auto', interpolation = 'none', origin = 'lower',extent = (0, lenth, 0, lenth))
        ax.set_ylabel(filename,fontsize=20)
        ax.set_xticks(ticks)        
        ax.set_xticklabels(labels)
        ax.set_yticks(ticks)
        ax.set_yticklabels(labels)        
        ax = fig.add_axes(size_colorbar)
        fig.colorbar(sc,cax = ax)
    
        ##DI
        DI = {'diff':diff , 'mm':ff , 'pp':mm , 'all':al}
        keys = ['diff','mm','pp','all']
        t=1
        size_signal = [left, bottom_h, width, 0.075]
        if int(len(matrix)/interval) < i+1:
            for k in keys:
#                if k == 'diff':
                    ax = fig.add_axes(size_signal)
                    y = np.zeros(300,dtype = np.int)
                    y[0:(len(DI[k])%300)] = DI[k][a:]
                    ax.fill_between(np.arange(len(y)), y , color = hexcolors[t])
                    ax.set_xticks([])
                    ax.set_yticks([DI[k].min()/2,0,DI[k].max()/2])
                    ax.set_ylabel(k, fontsize=10)
                    ax.spines['top'].set_visible(False)
                    ax.spines['right'].set_visible(False)
                    ax.set_xlim((0,len(y)))
                    ax.set_ylim((DI[k].min()-0.5,DI[k].max()+0.5))
                    t += 1
                    size_signal[1] += 0.075

        else:
            for k in keys:
                    ax = fig.add_axes(size_signal)
                    ax.fill_between(np.arange(len(DI[k][a:b])), DI[k][a:b] , color = hexcolors[t])
                    ax.set_xticks([])
                    ax.set_yticks([DI[k].min()/2,0,DI[k].max()/2])
                    ax.set_ylabel(k, fontsize=10)
                    ax.spines['top'].set_visible(False)
                    ax.spines['right'].set_visible(False)
                    ax.set_xlim((0,len(DI[k][a:b])))
                    ax.set_ylim((DI[k].min()-0.5,DI[k].max()+0.5))
                    t += 1
                    size_signal[1] += 0.075
        pp.savefig(fig)
        plt.close()
    pp.close()
    return pp 
    


for g in chrom:
    filename = 'NT5_chr' + g + '_30K.pdf'
    ff = 'ff_chr' + g + '_30k.txt'
    mm = 'mm_chr' + g + '_30k.txt'
    al = 'all_chr' + g + '_30k.txt'
    matrix_ff , npz_data_ff = loadfile(ff)
    matrix_mm , npz_data_mm = loadfile(mm)
    matrix_all , npz_data_all = loadfile(al)
    output_ff = Chrom(g,30000,npz_data_ff,'ff',10000000)
    output_mm = Chrom(g,30000,npz_data_mm,'mm',10000000)
    output_all = Chrom(g,300000,npz_data_all,'mm',10000000)
    output_ff.windows = np.ones(output_ff.chromLen, dtype = np.int32) * 20
    output_mm.windows = np.ones(output_mm.chromLen, dtype = np.int32) * 20
    output_all.windows = np.ones(output_all.chromLen, dtype = np.int32) * 20
    output_ff.calDI(output_ff.windows, 0)
    output_mm.calDI(output_mm.windows, 0)
    output_all.calDI(output_all.windows, 0)
    lena = len(output_ff.DIs)
    lenb = len(output_mm.DIs)    
    diff = np.zeros(lena , dtype = np.float)

    if lena == lenb:
        for i in range(lena):
            diff[i] = output_ff.DIs[i] - output_mm.DIs[i]
    elif lena < lenb:
        for i in range(lena):
            diff[i] = output_ff.DIs[i] - output_mm.DIs[i]
        print 'ff<mm:',lenb-lena
    else:
        for i in range(lenb):
            diff[i] = output_ff.DIs[i] - output_mm.DIs[i]
        print 'ff>mm:',lena-lenb
    pp = plot_matrix_DI(matrix_ff , diff , output_ff.DIs , output_mm.DIs , output_all.DIs , filename , 30000 )
    
        
        
nonzero = a[np.nonzero(a)]
y = np.zeros((2479,2479),dtype = np.int)
y[nonzero] = a        
        