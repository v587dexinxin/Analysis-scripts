# -*- coding: utf-8 -*-
"""
Created on Fri Jan 12 20:05:13 2018

@author: DELL
"""

from __future__ import division
import numpy as np
from scipy import sparse
import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt

class Dis_distribution:
    '''
    Use .npz data to do Principal Component Analysis(PCA) and to selec suitable PCs as compartment.
    
    Parameters
    ----------
    data : npz
         numpy produced npz data. keys = chrs; dtype = [ ('bin1', '<i4'), ('bin2', '<i4'), ('IF', '<f8')]
    
    res : int
         Kb resolution for hic data. Defult: 200.
    
    '''
    chrs = ['1' ,'2' ,'3' ,'4' ,'5' ,'6' ,'7' ,'8' ,'9' ,'10' ,'11' ,'12' 
        ,'13' ,'14' ,'15' ,'16' ,'17' ,'18' ,'19' ,'X']
    
    def __init__(self, data, gap, res = 1000):
        self.data = data
        self.res = res
        self.gap = gap.copy()
        self.get_distribution()
#        self.get_pcs()
        
    def get_matrix(self): 
        '''
        from sparse to matrix
        '''
        self.matrix = {}
        self.lenth = {}
        for k in self.chrs:
            tuple_matrix = self.data[k]
            sizes = tuple_matrix['bin2'].max() + 1
            P1 = sparse.coo_matrix((tuple_matrix['IF'],(tuple_matrix['bin1'],tuple_matrix['bin2'])),shape=(sizes,sizes)).toarray()
            P2 = sparse.coo_matrix((tuple_matrix['IF'],(tuple_matrix['bin2'],tuple_matrix['bin1'])),shape=(sizes,sizes)).toarray()
            mat = np.zeros((sizes, sizes))
            mat[np.nonzero(P1)] = P1[np.nonzero(P1)]
            mat[np.nonzero(P2)] = P2[np.nonzero(P2)]
            self.matrix[k] = mat
            self.lenth[k] = sizes
            
    def get_distribution(self):
        '''
        get distance decay line
        '''
        self.get_matrix()
        dict_distribution = {}
        self.gap['lenth'] = self.gap['lenth']/self.res/1000
        self.gap['S'] = self.gap['S']/self.res/1000
        self.gap['E'] = self.gap['E']/self.res/1000
        gaps = self.gap[self.gap['lenth'] > 0 ]
        for k in self.chrs:
            lenth = self.lenth[k]
            matrix_lenth = np.zeros((lenth, lenth))
            for i in range(lenth):
                for j in range(lenth):
                    if i < j :
                        matrix_lenth[i,j] = abs( i - j )
    
            c_gap = gaps[gaps['chr'] == 'chr'+k]
            for g in c_gap:
                gap_site = range(g[1], g[2])
                matrix_lenth[gap_site,:] = -1
                matrix_lenth[:,gap_site] = -1
                
            max_lenth = matrix_lenth.max() + 1
            distribution = []
            for i in np.arange(1,max_lenth):
                distribution.append(self.matrix[k][matrix_lenth == i])
            dict_distribution[k] = np.array(distribution)
        self.distribution = dict_distribution
            
if __name__ == '__main__':
    data_fold = 'E:\\data\\nuclear_transfer\\Raw_matrix\\matrix_40K\\'
    R_1 = np.load(data_fold + 'NT5-MboI-R1-filtered-40K-sparse.npz')
    R_2 = np.load(data_fold + 'NT5-MboI-R2-filtered-40K-sparse.npz')
    R_c1 = np.load(data_fold + 'NT6-MboI-R1-filtered-40K-sparse.npz')
    R_c2 = np.load(data_fold + 'NT6-MboI-R2-filtered-40K-sparse.npz')
    gap_file = np.loadtxt('E:\\data\\nuclear_transfer\\gap.txt', usecols=(1,2,3,6), dtype=[('chr', '<S6'), ('S', '<i4'), ('E', '<i4'), ('lenth', '<i4')])
    TR_1 = Dis_distribution(R_1,gap_file,1000)
    TR_2 = Dis_distribution(R_2,gap_file,1000)
    TR_c1 = Dis_distribution(R_c1,gap_file,1000)
    TR_c2 = Dis_distribution(R_c2,gap_file,1000)
    dis_1 = TR_1.distribution['1'][:240].copy()
    dis_2 = TR_2.distribution['1'][:240].copy()
    dis_c1 = TR_c1.distribution['1'][:240].copy()
    dis_c2 = TR_c2.distribution['1'][:240].copy() 
    for i in TR_1.chrs[1:]:
        for j in range(35):
            lens = min(TR_1.distribution[i][j].size, TR_2.distribution[i][j].size,TR_c1.distribution[i][j].size,TR_c2.distribution[i][j].size)
            dis_1[j] = np.hstack((dis_1[j], TR_1.distribution[i][j][:lens]))
            dis_2[j] = np.hstack((dis_2[j], TR_2.distribution[i][j][:lens]))
            dis_c1[j] = np.hstack((dis_c1[j], TR_c1.distribution[i][j][:lens]))
            dis_c2[j] = np.hstack((dis_c2[j], TR_c2.distribution[i][j][:lens]))
    corr_12 = []
    corr_34 = []
    corr_13 = []
    corr_24 = []
    for i in range(35):
        corr_12.append(np.corrcoef(dis_1[i],dis_2[i])[0,1])
        corr_34.append(np.corrcoef(dis_c1[i],dis_c2[i])[0,1])
        corr_13.append(np.corrcoef(dis_1[i],dis_c1[i])[0,1])
        corr_24.append(np.corrcoef(dis_2[i],dis_c2[i])[0,1])
    
    plt.plot(range(1,36), corr_12[:35], label = 'NT5')
    plt.plot(range(1,36), corr_34[:35], label = 'NT6')
#    plt.plot(range(1,36), corr_13[:35], label = 'NT5_NT6_R1')
#    plt.plot(range(1,36), corr_24[:35], label = 'NT5_NT6_R2')
    plt.ylim([0.9,1.1])
    plt.xlim([0,50])
    plt.legend()

    fig = plt.figure( figsize=(12,10) )
    grid = gridspec.GridSpec(6, 4, top=0.95, bottom=0.05, left=0.1, right=0.9,
                             hspace=0.25, wspace=0.25)
    chrs = TR_1.chrs
    for i in range(20):
        x = int(i / 4)
        y = i % 4
        c = chrs[i]
        ax = plt.subplot(grid[x,y])
        corr_12 = []
        corr_13 = []
        corr_14 = []
        corr_23 = []
        corr_24 = []
        corr_34 = []
        for j in range(30):
            lens = min(TR_1.distribution[c][j].size, TR_2.distribution[c][j].size,TR_c1.distribution[c][j].size, TR_c2.distribution[c][j].size)
            corr_12.append(np.corrcoef(TR_1.distribution[c][j][:lens],TR_2.distribution[c][j][:lens])[0,1])
            corr_13.append(np.corrcoef(TR_1.distribution[c][j][:lens],TR_c1.distribution[c][j][:lens])[0,1])
            corr_14.append(np.corrcoef(TR_1.distribution[c][j][:lens],TR_c2.distribution[c][j][:lens])[0,1])
            corr_23.append(np.corrcoef(TR_2.distribution[c][j][:lens],TR_c1.distribution[c][j][:lens])[0,1])
            corr_24.append(np.corrcoef(TR_2.distribution[c][j][:lens],TR_c2.distribution[c][j][:lens])[0,1])
            corr_34.append(np.corrcoef(TR_c1.distribution[c][j][:lens],TR_c2.distribution[c][j][:lens])[0,1])
        l12 = plt.plot(range(1,31), corr_12[:30], label = '12')
        l13 = plt.plot(range(1,31), corr_13[:30], label = '13')
        l14 = plt.plot(range(1,31), corr_14[:30], label = '14')
        l23 = plt.plot(range(1,31), corr_23[:30], label = '23')
        l24 = plt.plot(range(1,31), corr_24[:30], label = '24')
        l34 = plt.plot(range(1,31), corr_34[:30], label = '34')
        if y != 0 :
            plt.yticks([])
        plt.xticks([])
        plt.title('chr'+c)
        for spine in ['right', 'top', 'left', 'bottom']:
            ax.spines[spine].set_visible(False)
        ax.set_ylim((0.9,1))
    ax = plt.subplot(grid[5,3])
    plt.xticks([])
    plt.yticks([])
    ax.legend(l12+l13+l14+l23+l24+l34,('12','13','14','23','24','34'),loc = 3,fontsize=8)
    for spine in ['right', 'top', 'left', 'bottom']:
            ax.spines[spine].set_visible(False)



#    data_fold = 'E:\\hi-c\\gm12878\\'
#    R_1 = np.load(data_fold + 'gm12878-MboI-R1-filtered-200K-sparse.npz')
#    gap_file = np.loadtxt('C:\\Users\\DELL\\Desktop\\gap.txt', usecols=(1,2,3,6), dtype=[('chr', '<S6'), ('S', '<i4'), ('E', '<i4'), ('lenth', '<i4')])
#    R_1 = Dis_distribution(R_1,gap_file,200)