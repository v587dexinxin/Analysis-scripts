# -*- coding: utf-8 -*-
"""
Created on Thu Nov 23 20:53:41 2017

@author: xxli
"""

f1=open("C:/Users/xxli/Desktop/BDF1_narrowpeaks/ES-1_peaks.narrowPeak","rt")
f2=open("C:/Users/xxli/Desktop/BDF1_narrowpeaks/NT5-1_peaks.narrowPeak","rt")
f3=open("C:/Users/xxli/Desktop/BDF1_narrowpeaks/NT6-2_peaks.narrowPeak","rt")
f4=open("C:/Users/xxli/Desktop/BDF1_narrowpeaks/ESC-ESC_NT5-chrX.txt","wt")
f5=open("C:/Users/xxli/Desktop/BDF1_narrowpeaks/ESC-ESC_NT6-chrX.txt","wt")
f6=open("C:/Users/xxli/Desktop/BDF1_narrowpeaks/NT5-NT5_ESC-chrX.txt","wt")
f7=open("C:/Users/xxli/Desktop/BDF1_narrowpeaks/NT5-NT5_NT6-chrX.txt","wt")
f8=open("C:/Users/xxli/Desktop/BDF1_narrowpeaks/NT6-NT6_ESC-chrX.txt","wt")
f9=open("C:/Users/xxli/Desktop/BDF1_narrowpeaks/NT6-NT6_NT5-chrX.txt","wt")
a=f1.readlines()
b=f2.readlines()
c=f3.readlines()
chrs=('chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9','chr10','chr11','chr12','chr13','chr14','chr15','chr16','chr17','chr18','chr19','chr20','chr21','chr22','chrX')
euchr=('chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9','chr10','chr11','chr12','chr13','chr14','chr15','chr16','chr17','chr18','chr19','chr20','chr21','chr22')
chrX=(['chrX'])
ES={}
NT5={}
NT6={}
ES_NT5={}
ES_NT6={}
NT5_NT6={}
NT5_ES={}
NT6_ES={}
NT6_NT5={}
es_nt5=0
nt5_es=0
es_nt6=0
nt6_es=0
nt5_nt6=0
nt6_nt5=0
es_nt5_nt6=0
for g in chrs:
    ES[g]=[]
    NT5[g]=[]
    NT6[g]=[]
    ES_NT5[g]=[]
    ES_NT6[g]=[]
    NT5_NT6[g]=[]
    NT5_ES[g]=[]
    NT6_ES[g]=[]
    NT6_NT5[g]=[]
#数据分类
for x in a:
    x=x.split()
    if ES.has_key(x[0]):
        ES[x[0]].append([x[1],x[2],x[6]])
for y in b:
    y=y.split()
    if NT5.has_key(y[0]):
        NT5[y[0]].append([y[1],y[2],y[6]])
for z in c:
    z=z.split()
    if NT6.has_key(z[0]):
        NT6[z[0]].append([z[1],z[2],z[6]])
#数据比较
for g in chrX:
    for i in ES[g]:
        for j in NT5[g]:
            if not (i[1]<=j[0] or i[0]>=j[1]):
                ES_NT5[g].append([i[0],i[1],i[2]])
                es_nt5 +=1
                break
for g in chrX:
    for i in NT5[g]:
        for j in ES[g]:
            if not (i[1]<=j[0] or i[0]>=j[1]):
                NT5_ES[g].append([i[0],i[1],i[2]])
                nt5_es +=1
                break
for g in chrX:
    for i in ES[g]:
        for j in NT6[g]:
            if not (i[1]<=j[0] or i[0]>=j[1]):
                ES_NT6[g].append([i[0],i[1],i[2]])
                es_nt6 +=1
                break
for g in chrX:
    for i in NT6[g]:
        for j in ES[g]:
            if not (i[1]<=j[0] or i[0]>=j[1]):
                NT6_ES[g].append([i[0],i[1],i[2]])
                nt6_es +=1
                break
for g in chrX:
    for i in NT5[g]:
        for j in NT6[g]:
            if not (i[1]<=j[0] or i[0]>=j[1]):
                NT5_NT6[g].append([i[0],i[1],i[2]])
                nt5_nt6 +=1
                break
for g in chrX:
    for i in NT6[g]:
        for j in NT5[g]:
            if not (i[1]<=j[0] or i[0]>=j[1]):
                NT6_NT5[g].append([i[0],i[1],i[2]])
                nt6_nt5 +=1
                break
for g in chrX:
    for i in ES_NT6[g]:
        for j in NT5[g]:
            if not (i[1]<=j[0] or i[0]>=j[1]):
                es_nt5_nt6 +=1
                break
#输出细胞特异性peaks

for g in chrX:
    for i in ES[g]:
        if i not in ES_NT5[g]:
            f4.writelines(g+'\t'+i[0]+'\t'+i[1]+'\t'+i[2]+'\n')
for g in chrX:
    for i in ES[g]:
        if i not in ES_NT6[g]:
            f5.writelines(g+'\t'+i[0]+'\t'+i[1]+'\t'+i[2]+'\n')
for g in chrX:
    for i in NT5[g]:
        if i not in NT5_ES[g]:
            f6.writelines(g+'\t'+i[0]+'\t'+i[1]+'\t'+i[2]+'\n')
for g in chrX:
    for i in NT5[g]:
        if i not in NT5_NT6[g]:
            f7.writelines(g+'\t'+i[0]+'\t'+i[1]+'\t'+i[2]+'\n')
for g in chrX:
    for i in NT6[g]:
        if i not in NT6_ES[g]:
            f8.writelines(g+'\t'+i[0]+'\t'+i[1]+'\t'+i[2]+'\n')
for g in chrX:
    for i in NT6[g]:
        if i not in NT6_NT5[g]:
            f9.writelines(g+'\t'+i[0]+'\t'+i[1]+'\t'+i[2]+'\n')

print 'es_nt5:'+str(es_nt5)+'\n'+'nt5_es:'+str(nt5_es)
print 'es_nt6:'+str(es_nt6)+'\n'+'nt6_es:'+str(nt6_es)
print 'nt5_nt6:'+str(nt5_nt6)+'\n'+'nt6_nt5:'+str(nt6_nt5)
print 'es_nt5_nt6:'+str(es_nt5_nt6)
f1.close()
f2.close()
f3.close()
f4.close()
f5.close()
f6.close()
f7.close()
f8.close()
f9.close()