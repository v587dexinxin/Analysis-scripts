# -*- coding: utf-8 -*-
"""
Created on Tue Jul 17 15:17:54 2018

@author: xxli
"""

import numpy as np


cell = 'NT6'
itype = np.dtype({'names':['gene_id','reads_count'],
                      'formats':['<S32' , np.int]})
f1 = np.loadtxt('C:\\Users\\xxli\\Desktop\\allele-specific\\RNA_seq\\RNA_newSNP_9\\Specific_' + cell + '_M_count.txt',dtype = itype)
f2 = np.loadtxt('C:\\Users\\xxli\\Desktop\\allele-specific\\RNA_seq\\RNA_newSNP_9\\Specific_' + cell + '_P_count.txt',dtype = itype)

f3 = np.loadtxt('C:\\Users\\xxli\\Desktop\\allele-specific\\RNA_seq\\Specific_' + cell + '_M_count.txt',dtype = itype)
f4 = np.loadtxt('C:\\Users\\xxli\\Desktop\\allele-specific\\RNA_seq\\Specific_' + cell + '_P_count.txt',dtype = itype)

a1 = list(f1['reads_count'])
a2 = list(f2['reads_count'])

a3 = list(f3['reads_count'])
a4 = list(f4['reads_count'])

del a1[-5:]
del a2[-5:]
del a3[-5:]
del a4[-5:]


a1 = np.array(a1)
a2 = np.array(a2)
a3 = np.array(a3)
a4 = np.array(a4)

M_new = a1.sum()
P_new = a2.sum()

M_old = a3.sum()
P_old = a4.sum()


print M_new/(M_new + P_new)
print M_old/(M_old + P_old)