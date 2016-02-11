#!/usr/bin/env python

import sys,os

from optparse import OptionParser
parser = OptionParser()

#popfile = sys.argv[1];      del sys.argv[1]
listdir = sys.argv[1]; del sys.argv[1]
outfile= sys.argv[1];      del sys.argv[1]
p = int(sys.argv[1]);      del sys.argv[1]
d1= int(sys.argv[1]);      del sys.argv[1]
d2= int(sys.argv[1]);      del sys.argv[1]
 
from logging import FileHandler, info
from logging import warning, basicConfig, DEBUG


x = outfile.split('.')
logfile = x[0] + '.log'

basicConfig(level=DEBUG,
    format='%(asctime)s %(levelname)s %(message)s',
    filename=logfile,
    filemode='w')

#parser.add_option('-d')
#parser.add_option('-p',type = 'int')
#parser.add_option('-t')
parser.add_option('-m')

(options,args) = parser.parse_args()

#x = options.d.split(':')
#d1,d2 = int(x[0]),int(x[1])

#p = options.p

#taxon = options.t.split(':')

AvTD = {}
VarTD = {}
AvTDm = {}
VarTDm = {}
AvTDmax = {}
VarTDmin = {}

for d in range(d1,d2+1):
    AvTD[d] = []
    VarTD[d] = []
    AvTDm[d] = []
    VarTDm[d] = []
    AvTDmax[d] = []
    VarTDmin[d] = []

if options.m:
    missing = options.m
else:
    missing = 'n'
"""program listdir outfile options

-d dimensions  -d d1:d2 no space between d1 and d2
-p permutations
-t taxon  taxon1:taxon2:taxon3, etc. no spaces
-m 'y' if missing data
"""




sample = {}; population = {}

#out = samplefile.split('.')
#output = out[0] + '.out'
#o = open(output,'a')

#saveout = sys.stdout
#sys.stdout = open(output, 'w')

from re import *

Taxon = {}; coef = {}; Taxon = {}; 

pathLengths= {}

###--------------functions-----------------------------###

def PathLength(population):
    taxonN = {}
    for t in taxon:
        Taxon[t] = {}
        x = [population[i][t] for i in sample]

        for i in set(x):
            Taxon[t][i] = x.count(i)
            
        taxonN[t] = len(Taxon[t])

    n = [float(len(Taxon[t])) for t in taxon]
    n.insert(0,1.0)
        
    #s = 100/float(N)
    raw = []
    for i in range((len(n)-1)):
        j = i + 1
        
        if n[i] > n[j]:
            c = 1
        else:
            c = (1 - n[i]/n[j])

        raw.append(c)

    s = sum(raw)

    adjco = [i*100/s for i in raw]
    coef = {}; pathLengths = {}
    for i in range(len(taxon)):
        t = taxon[i]
        coef[t] = sum(adjco[i:])
        pathLengths[t] = adjco[i]

    return coef, taxonN, pathLengths

def ATDmean(data,sample,coef):
    #[sample = data.keys()
    N = len(sample)

    Taxon = {}; taxonN = {}; AvTD = 0; n = 0
    #Taxon are counts of taxa at each level, taxonN are numbers of pairwise differences
    #at each level, with n being the accumlation of pairwise differences at that level. the difference
    #between n and TaxonN is the number of species that are in different taxa in that level
    #but not in upper levels

    for t in taxon:
        Taxon[t] = {}
        x = [data[i][t] for i in sample]
        for i in set(x):
            Taxon[t][i] = x.count(i)

    for t in taxon:
        taxonN[t] = sum([Taxon[t][i] * Taxon[t][j] for i in Taxon[t] for j in Taxon[t] if i != j])
        n = taxonN[t] - n
        AvTD = AvTD + (n * coef[t]) 
        n = taxonN[t]

    AvTD /= (N * (N - 1))
    
    return AvTD,taxonN, Taxon

def ATDvariance(taxonN, sample, atd, coef):
    vtd = []
    
    #N = sum(taxon)

    vtd = 0; N = 0; n = 0

    for t in taxon:
        n = taxonN[t] - n
        vtd = vtd + n * coef[t]**2 
        n = taxonN[t]

    N = len(sample)
    n = N * (N - 1)

    vtd = (vtd - ((atd*n)**2)/n)/n

    #vtd = (sum([tax1,tax2,tax3,tax4]) - (((atd * n)**2)/n))/n

    return vtd

def Sample(samplefile):
    sample = {}
    for i in open(samplefile):
        if match('Taxon:', i): continue
        elif match('Coefficients:', i): continue

        x = i.split()

        species = x[0]
        #sample[species] = {}
    
    #sample.append(species)
    """
        if missing == 'y':
            mtax = ''
            for t in taxon:
                
                if x[Index[t]] == '/':
                    #sample[species][t] = sample[species][t]
                    sample[species][t] = mtax
                else:
                    sample[species][t] = x[Index[t]]
                    mtax = x[Index[t]]

        else:
            for t in taxon:
                #y = Taxon[t]
                sample[species][t] = x[Index[t]]
    """
    sample[species] = population[species]
    
    return sample

def Funnel(p,d1,d2,population,coef):
        from random import sample
        pop = population.keys()
    
        dims = []; up = []; lo = []; means = []

        for d in range(d1, d2 + 1):
            perm = 'dimension: ' + str(d)
            info(perm)
            x.append(d)
            AvTDci = []; VarTDci = []
            for j in range(p):
                perm = 'permutation: ' + str(j + 1)
                info(perm)
                rsamp = sample(pop,d)
                atd,taxonN, Taxon = ATDmean(population,rsamp,coef); AvTDci.append(atd)
                vtd = ATDvariance(taxonN,rsamp,atd,coef); VarTDci.append(vtd)

            AvTDci.sort()
            VarTDci.sort()
        
            #AvTD = AvTDci[int(.05 * p)], sum(AvTDci)/p, AvTDci[int(.95 * p)], max(AvTDci)
            #VarTD = min(VarTDci), VarTDci[int(.05 * p)],sum(VarTDci)/p,VarTDci[int(.95 * p)] 
    
            AvTD[d].append(AvTDci[int(.05 * p)])
            VarTD[d].append(VarTDci[int(.95 * p)])
            AvTDm[d].append(sum(AvTDci)/p)
            VarTDm[d].append(sum(VarTDci)/p)
            AvTDmax[d].append(max(AvTDci))
            VarTDmin[d].append(min(VarTDci))

            #dims.append(d)
            #ciarray.append(AvTD[0])
            #carray.append(AvTD[1])

###------------------end functions------------------------###
n = len(os.listdir(listdir))

info('begin analysis of master lists')

for f in os.listdir(listdir):
    sample = {}; population = {}
    duplicates = []
    #l = '%s/%s' %(listdir,f)
    
    mlist = 'master list: ' + str(f)
    info(mlist)
    
    try:
        l = '%s/%s' %(listdir,f)
    except:
        l = '%s\%s' %(listdir,f)
    Index = {}
  
    
    for i in open(l):
        if match('Taxon:', i):
            x = i.split()
            taxon = x[1:]
            for t in taxon:
                j = taxon.index(t)
                Index[t] = j + 1
            continue

        elif match('Coefficients:', i):
            continue
        
        i.strip()
        x = i.split()    
        
        species = x[0]
    
        if species in sample.keys():
            duplicates.append(species)
        else:
            sample[species] = {}

        if missing == 'y':
            mtax = ''
            for t in taxon:
                if x[Index[t]] == '/':
                    sample[species][t] = mtax
                    
                else:
                    sample[species][t] = x[Index[t]]
                    mtax = x[Index[t]]
        else:
            for t in taxon:
                 sample[species][t] = x[Index[t]]
        
            #population[species] = {}
            #for t in taxon:
            #   population[species][t] = x[Index[t]]
        

    coef, popN, pathLengths = PathLength(sample)
    Funnel(p,d1,d2, sample,coef)
    #a, v = Funnel(p,d1,d2)


o  =open (outfile + '.psa','w')

l = 'd AvTD05%_2.5% AvTD05%_97.5% AvTDmean_2.5% AvTDmean_97.5% AvTDmax VarTD95%_2.5% VarTD95%_97.5% VarTDmean_2.5% VarTDmean_97.5% VarTDmin\n'
o.write(l)

n = len(os.listdir(listdir))

for d in range(d1, d2+1):
    AvTD[d].sort()
    VarTD[d].sort()
    AvTDm[d].sort()
    VarTDm[d].sort()

    print len(AvTD[d])
    print len(VarTD[d])
    print len(AvTDm[d])
    print len(VarTDm[d])
    l = '%i        %6.4f   %6.4f   %6.4f   %6.4f   %6.4f   %6.4f   %6.4f   %6.4f   %6.4f   %6.4f\n' \
    % (d, AvTD[d][int(n*0.025)], AvTD[d][int(n*0.975)],AvTDm[d][int(n*0.025)],AvTDm[d][int(n*0.975)], max(AvTDmax[d]), VarTD[d][int(n*0.025)],VarTD[d][int(n*0.975)],VarTDm[d][int(n*0.025)],VarTDm[d][int(n*0.975)], min(VarTDmin[d]))

    o.write(l)

