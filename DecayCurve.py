# -*- coding: utf-8 -*-
"""
Created on Tue Sep 05 15:35:21 2017

@author: han-luo
"""

#      HiC Contact decay curve with distance
from __future__ import division
import numpy as np
import sys, cPickle

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
max_distance = 10000000
    #distance array step : Res
distance = np.arange(1,max_distance//Res,1)

for i in distance:
    distance_contact[i] = []

for c in chrmo:
    # contact Matrix
    Matrix = HiClib[c]
    for i in distance:
        mask = (Matrix['bin2'] - Matrix['bin1']) == i
        IF = Matrix[mask]['IF'].mean()
        distance_contact[i] += [IF]
    print "chromosome %s is done" % c
for dis,value in distance_contact.items():
    distance_contact[dis] = np.array(value).mean()

    #output pickle file
out = sys.argv[2]    
f = open(out,'wb')
cPickle.dump(distance_contact,f,2)
f.close()
