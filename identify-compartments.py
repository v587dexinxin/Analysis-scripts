# -*- coding: utf-8 -*-
"""
Created on Fri Mar 10 15:47:48 2017

@author: wxt
"""

import sys, os
import numpy as np
from runHiC.chiclib import myGenome
from hiclib import binnedData
from mirnylib.h5dict import h5dict

def compartment_core(pos, res):
    
    gapsizes = pos[1:] - pos[:-1]
    endsIdx = np.where(gapsizes>1)[0]
    startsIdx = endsIdx + 1
    
    As = []
    for i in range(startsIdx.size - 1):
        start = pos[startsIdx[i]]
        end = pos[endsIdx[i+1]] + 1
        if end - start >= 1:
            As.append([start*res, end*res])
    
    if startsIdx.size > 0:
        start = pos[startsIdx[-1]]
        end = pos[-1] + 1
        if end - start >= 1:
            As.append([start*res, end*res])
        start = pos[0]
        end = pos[endsIdx[0]] + 1
        if end - start >= 1:
            As.append([start*res, end*res])

    if not startsIdx.size:
        start = pos[0]
        end = pos[-1]
        if end - start >= 1:
            As.append([start*res, end*res])
    
    return As

def compartmentFromPCA(pca, res, chrom, minsize=3):
    
    A_pos = np.where(pca>0)[0]
    As = compartment_core(A_pos,res)
    B_pos = np.where(pca<=0)[0]
    Bs = compartment_core(B_pos,res)
    compartments = []
    for a in As:
        if a[1] - a[0] >= res*minsize:
            compartments.append([chrom]+a+['A'])
        else:
            compartments.append([chrom]+a+['N'])
    for b in Bs:
        if b[1] - b[0] >= res*minsize:
            compartments.append([chrom]+b+['B'])
        else:
            compartments.append([chrom]+b+['N'])
    compartments.sort()
    
    return compartments

HiCFolder = '/public/home/hluo/Data/Mouse/HiC/compartments'

source = os.path.join(HiCFolder, sys.argv[1]) # file after calling runHiC binning
outfile = sys.argv[2]

# params from source
lib = h5dict(source, 'r')
res = lib['resolution']
gInfo = lib['genomeInformation']
genomeFolder = os.path.join(gInfo['dataFolder'], gInfo['genomeName'])

# Load binned data and perform PCA ...
genome_db = myGenome(genomeFolder, readChrms=gInfo['chroms'],
                     chrmFileTemplate=gInfo['template'],
                     gapFile=gInfo['gapFile'])
BD = binnedData.binnedData(res, genome_db)
name = os.path.split(source)[1].split('-')[0]
BD.simpleLoad(source, name)
BD.doCisPCADomains(3)
# Identify compartments ...
pcas = BD.PCDict[name].T
idx2label = gInfo['idx2label']
for i in idx2label:
    label = idx2label[i]
    mask = BD.chromosomeIndex==i
    tmp = pcas[mask]
    compartments = compartmentFromPCA(tmp[:,0], res, label)
    with open(outfile, 'a') as output:
        for c in compartments:
            line = '\t'.join(map(str,c))+'\n'
            output.write(line)
