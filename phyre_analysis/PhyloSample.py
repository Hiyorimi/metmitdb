#!/usr/bin/env python

import sys
from random import *
popfile = sys.argv[1];      del sys.argv[1]
outfile= sys.argv[1];      del sys.argv[1]


mode = sys.argv[1];      del sys.argv[1]
Pdict = {}

if mode == 'a':
	p = sys.argv[1:];
	#p = p[0].split()
	Pdict['n'] = int(p[0])
	Pdict['t'] = p[1]
	Pdict['s'] = int(p[2])
	Pdict['e'] = int(p[3])
	Pdict['m'] = int(p[4])

#mode: 'a' automated -> a n t s e m
#n = number of new files; s = split; e = merge; m = move (number of each); t = taxon level
#mode: 'i' interactive


TA = {}; ta = []
master = {}; Tgroups = {}
#pop
changes = {}
ReverseTaxon = {}

KeepinOrder = []

for i in open(popfile):
	x = i.split(); s = x[0];
	if x[0] == 'Taxon:':
		for j in x[1:]:
			TA[j] = {}
			ta.append(j)
			changes[j] = {}
			Tgroups[j] = {}
			ReverseTaxon[j] = {}
		continue
	x = x[1:]
	master[s] = [i for i in x]
	x.append(s)
	KeepinOrder.append(s)
	for j in ta:
		if x[ta.index(j)] == '/': continue
		if x[ta.index(j)] not in TA[j]:
			ReverseTaxon[j][x[ta.index(j)]] = x[:ta.index(j)]
			Tgroups[j][x[ta.index(j)]] = [s]
			if x[ta.index(j) + 1] != '/':
				TA[j][x[ta.index(j)]] = [x[ta.index(j) + 1]]
			#else: TA[j][x[ta.index(j)]] = [x[ta.index(j) + 2]]
		else:
			if x[ta.index(j) + 1] != '/':
				if x[ta.index(j) + 1] not in TA[j][x[ta.index(j)]]:
					Tgroups[j][x[ta.index(j)]].append(s)
					TA[j][x[ta.index(j)]].append(x[ta.index(j) + 1])
			#else:
			#	if x[ta.index(j) + 2] not in TA[j][x[ta.index(j)]]:
			#		TA[j][x[ta.index(j)]].append(x[ta.index(j) + 2])


def MergeTaxa(ch, new,t, k=2):
	a = 0

	while a == 0:
		if ta.index(t) == 0:
			s = sample(TA[t],k)
			if s[0] in new or s[1] in new: continue
			if s[0] not in ch[t] or s[1] not in ch[t]: a = 1
		else:
			S = sample(TA[ta[ta.index(t)-1]],1)[0]
			#if S[0] in new or S[1] in new: continue
			if len(TA[ta[ta.index(t)-1]][S]) >= k:
				s = sample(TA[ta[ta.index(t)-1]][S],k)
				#if s[0] not in ch[t] or s[1] not in ch[t]: a = 1
				b = set(s)
				if len(b.intersection(ch[t].keys())) == 0:
					a = 1


	#while a == 0:
	#s = sample(TA[t],2)
	#if s[0]  s[1]
	for i in s:
		ch[t][s[s.index(i)]] = ['merge', ''.join(s)] 
	#ch[t][s[0]] = ['merge', s[0] + s[1]]
	#ch[t][s[1]] = ['merge', s[0] + s[1]]
	
	#C[t][s[0] + s[1]] = C[t][s[0]] + C[t][s[1]]
	#new.append(C[t][s[0] + s[1]])
	#C[t][s[0] + s[1]]
	return ch,new

def SplitTaxa(S, ch, new,t):
	a = 0

	while a == 0:
		s = sample(TA[t],1)[0]
		if len(TA[t][s]) >= 2: 
			if s not in ch[t].keys() and s not in new: a = 1
	
	S[t][s + '_1'] = []
	S[t][s + '_2'] = []

	ch[t][s] = ['split', s + '_1', s + '_2']

	for j in range(0,len(TA[t][s]),2):
		S[t][s + '_1'].append(TA[t][s][j])
	
	for j in range(1,len(TA[t][s]),2):
		S[t][s + '_2'].append(TA[t][s][j])

	new.append(S[t][s + '_2']); new.append(S[t][s + '_1'])
	#del C[t][s]
	return ch,new,S

def MoveTaxa(ch, new,t):
	a = 0

	if ta.index(t) == 0: print "cannot move taxa within uppermost level"
	
	elif ta.index(t) == 1:
		tt = ta[ta.index(t)-1]
		while a == 0:
			S = sample(TA[tt],2)
			#s = sample(T[t],1)[0]
			if len(TA[tt][S[0]]) >= 2: 
				s = sample(TA[tt][S[0]],1)[0]
				if s not in ch[t].keys(): a = 1
				
	elif ta.index(t) > 1: 
		tt = ta[ta.index(t)-2]

		while a == 0:
			S = sample(TA[tt],1)[0]
			if len(TA[tt][S]) >= 2:			
				S = sample(TA[tt][S],2)
				#s = sample(T[t],1)[0]
				try:
					if len(TA[ta[ta.index(tt)+1]][S[0]]) >= 2: 
						s = sample(TA[ta[ta.index(tt)+1]][S[0]],1)[0]
						if s not in ch[t].keys(): a = 1
				except:
					if len(TA[ta[ta.index(tt)+1]][S[1]]) >= 2: 
						s = sample(TA[ta[ta.index(tt)+1]][S[0]],1)[0]
						if s not in ch[t].keys(): a = 1
	ch[t][s] = ['move', S[1]]
	
	new.append(s)
	#C[tt][S[1]].append(s)
	#C[tt][S[0]].remove(s)
	return ch,new

if mode == 'i':
	z = 0
	S = {}
	new = []
	for i in ta:
		S[i] = {}
		
	while z == 0:
		prompt = 'what function would you like to perform? input m for move, e for merge, or s for split '
		i = str(raw_input(prompt))

		l = ''
		for j in ta[1:len(ta)]:
			l += '%d for %s ' %(ta.index(j), j)

		prompt = 'please input taxon level to manipulate: %s ' %l
		t = ta[input(prompt)]

		prompt = 'please input number of repetitions ' 
		n = input(prompt)

		if i == 'm':
			for r in range(n):
				changes,new = MoveTaxa(changes, new,t)
		elif i == 'e':
			prompt = 'how many taxon would you like to merge '
			K = int(input(prompt))
			for r in range(n):
				changes,new = MergeTaxa(changes,new,t, k=K)
		elif i == 's':
			for r in range(n):
				ch,new,S = SplitTaxa(S, changes, new,t)
		else:
			print 'that is not a legal option, please input again m for move, e for merge, or s for split'

		prompt = 'would you like to continue? please input y for yes or n for no'
		i = (raw_input(prompt))
		
		if i == 'n':
			print "printing new master list and changes file" 
			z = 1


	for t in ta:
		for i in changes[t]:
				for j in Tgroups[t][i]:
					if changes[t][i][0] == 'move':
						tt = ta[ta.index(t)-1]
						master[j][ta.index(tt)] = changes[t][i][1]
						for r in range(len(ReverseTaxon[tt][changes[t][i][1]])):
							master[j][r] = ReverseTaxon[tt][changes[t][i][1]][r]
					elif changes[t][i][0] == 'merge':
						master[j][ta.index(t)] = changes[t][i][1]

					elif changes[t][i][0] == 'split':
						tt = ta[ta.index(t)+1]
						for c in changes[t][i][1:]:
							for s in S[t][c]:
								if master[j][ta.index(tt)] == s:
									master[j][ta.index(t)] = c

	out = outfile + '.sim'
	o = open(out,'w')

	t = 'Taxon:\t%s\n' % '\t'.join(ta)
	o.write(t)

	for i in KeepinOrder:
	#for i in master:
		l = '%s\t%s\n' %(i,'\t'.join(master[i]))
		o.write(l)

	o.close()

	out = outfile + '.changes'
	o = open(out,'w')

	for t in changes:
		for i in changes[t]:
			if changes[t][i][0] == 'split':
				l = '%s %s split to %s and %s\n' %(t,i,changes[t][i][1],changes[t][i][2])
				o.write(l)
	for t in changes:
		for i in changes[t]:
			if changes[t][i][0] == 'merge':
				l = '%s %s merged to %s\n' %(t,i,changes[t][i][1])
				o.write(l)
	for t in changes:
		for i in changes[t]:
			if changes[t][i][0] == 'move':
				l = '%s %s moved to %s\n' %(t,i,changes[t][i][1])
				o.write(l)
	o.close()

####
import os
if mode == 'a':
	f = '%sMasterFiles' %outfile
	fc = '%sMasterChanges' %outfile
	try:
		os.mkdir(f)
		os.mkdir(fc)
	except:
		print 'directory exists'

	for n in range(Pdict['n']):
		ch = {}
		new = []
		S = {}
		
		M = {}
		
		for i in master:
			M[i] = master[i]
		
		for c in changes:
			ch[c] = {}
		
		for i in ta:
			S[i] = {}
		
		for m in range(Pdict['m']):
			ch,new = MoveTaxa(ch, new,Pdict['t'])
		for e in range(Pdict['e']):
			ch,new = MergeTaxa(ch, new,Pdict['t'], k=2)
		for s in range(Pdict['s']):
			ch,new,S = SplitTaxa(S, ch, new,Pdict['t'])
		
		#out = '%s_'

		G = {}

		for t in ta:
			for i in ch[t]:
					for j in Tgroups[t][i]:
						if ch[t][i][0] == 'merge':
							M[j][ta.index(t)] = ch[t][i][1]
						elif ch[t][i][0] == 'move':
							tt = ta[ta.index(t)-1]
							M[j][ta.index(tt)] = ch[t][i][1]
							for r in range(len(ReverseTaxon[tt][ch[t][i][1]])):
								M[j][r] = ReverseTaxon[tt][ch[t][i][1]][r]
						elif ch[t][i][0] == 'split':
							tt = ta[ta.index(t)+1]
							for c in ch[t][i][1:]:
								for s in S[t][c]:
									if M[j][ta.index(tt)] == s:
										M[j][ta.index(t)] = c

		
	
		#out = '%s/%s_%d.sim' %(f,outfile,n+1)
		out = '%s\%s_%d.sim' %(f,outfile,n+1)
		o = open(out,'w')

		t = 'Taxon:\t%s\n' % '\t'.join(ta)
		o.write(t)

		for i in KeepinOrder:
		#for i in master:
			l = '%s\t%s\n' %(i,'\t'.join(M[i]))
			o.write(l)

		o.close()

		#out = '%s/%s_%d.changes' %(fc,outfile,n+1)
		out = '%s\%s_%d.changes' %(fc,outfile,n+1)
		o = open(out,'w')

		for t in ch:
			for i in ch[t]:
				if ch[t][i][0] == 'split':
					l = '%s %s split to %s and %s\n' %(t,i,ch[t][i][1],ch[t][i][2])
					o.write(l)
		for t in ch:
			for i in ch[t]:
				if ch[t][i][0] == 'merge':
					l = '%s %s merged to %s\n' %(t,i,ch[t][i][1])
					o.write(l)
		for t in ch:
			for i in ch[t]:
				if ch[t][i][0] == 'move':
					l = '%s %s moved to %s\n' %(t,i,ch[t][i][1])
					o.write(l)
		o.close()

