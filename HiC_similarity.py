# -*- coding: utf-8 -*-
"""
Created on Wed Sep 13 09:26:28 2017

@author: han-luo
"""

from __future__ import division
import numpy as np
import sys,cPickle

    #Sparse Matrix To 2D Matrix for a chrmosome region
def getmatrix(inter,l_bin,r_bin):
    inter_matrix = np.zeros((r_bin - l_bin, r_bin - l_bin),dtype = float )
    #Extract the regional data
    mask = (inter['bin1'] >= l_bin) & (inter['bin1'] < r_bin) & \
           (inter['bin2'] >= l_bin) & (inter['bin2'] < r_bin)
    inter_extract = inter[mask]
    
    #Fill the matrix:
    for i in inter_extract:
        #Off-diagnoal parts
        if i['bin1'] != i['bin2']:
            inter_matrix[i['bin1'] - l_bin][i['bin2'] - l_bin] += i['IF']
            inter_matrix[i['bin2'] - l_bin][i['bin1'] - l_bin] += i['IF']
        else:
            #Diagonal part
            inter_matrix[i['bin1'] - l_bin][i['bin2'] - l_bin] += i['IF']
            
    return inter_matrix

    # load HiC data
HiCSource = sys.argv[1]
HiClib = np.load(HiCSource)
    #resolution
Res = int(HiClib['resolution'])
    #chromosome
#chrmo = ['1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19','X']
chrmo = [i for i in HiClib.keys() if i != 'resolution' and i != 'genomeInformation']

distance_contact = {}
    #max distance 10Mb
max_distance = 50000000
    #distance array step : Res
distance = np.arange(1,max_distance//Res,1)

for i in distance:
    distance_contact[i] = []
    
for c in chrmo:
    #Contact Matrix
    Matrix = getmatrix(HiClib[c],100,200)      
    for i in distance:
#        mask = (Matrix['bin2'] - Matrix['bin1']) == i
#        IF = Matrix[mask]['IF'].mean()
#        distance_contact[i] += [IF]
        tmp = []
        for m in range(50):
            tmp.append(Matrix[m][m+i])
        distance_contact[i] += tmp
    print "chromosome %s is done" % c

for dis,value in distance_contact.items():
    distance_contact[dis] = np.array(value)

out = sys.argv[2]    
f = open(out,'wb')
cPickle.dump(distance_contact,f,2)
f.close()