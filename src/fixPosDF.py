# -*- coding: utf-8 -*-
"""
Created on Thu Nov  2 08:42:08 2017

@author: gritti
"""

import pickle
import numpy as np
import glob
import pandas as pd

paths = [ 'Y:\\Simone\\161030_lag2YFP+histmCherry' ]

worms1 = []
for path in paths:
	worms1.append( glob.glob( os.path.join( path, '*cellPos.pickle' ) ) )
worms = []
for exp in worms1:
	for w in exp:
		worms.append(w)
worms.sort()

for path in paths:
	N=5
	for worm in worms[N:N+1]:
		wName = worm.split('\\')[-1].split('_')[0]
		print(wName)
		
		df = pickle.load(open(worm,'rb'))
		
		#fix 1. cell position df
		t1p = df['1.p'].ix[pd.notnull(df['1.p'].X),'tidx'].values
		t1pp = df['1.pp'].ix[pd.notnull(df['1.pp'].X),'tidx'].values
		t1ppp = df['1.ppp'].ix[pd.notnull(df['1.ppp'].X),'tidx'].values
		t1ppa = df['1.ppa'].ix[pd.notnull(df['1.ppa'].X),'tidx'].values
		tb1 = df['b_1'].ix[pd.notnull(df['b_1'].X),'tidx'].values
		
		t4a = df['4.a'].ix[pd.notnull(df['4.a'].X),'tidx'].values
		t4aa = df['4.aa'].ix[pd.notnull(df['4.aa'].X),'tidx'].values
		t4aaa = df['4.aaa'].ix[pd.notnull(df['4.aaa'].X),'tidx'].values
		t4aap = df['4.aap'].ix[pd.notnull(df['4.aap'].X),'tidx'].values
		tb4 = df['b_4'].ix[pd.notnull(df['b_4'].X),'tidx'].values

		print(t1p)
		print(t1pp)
		print(t1ppp)
		print(t1ppa)
		print(tb1)
		print('\n')
		print(t4a)
		print(t4aa)
		print(t4aaa)
		print(t4aap)
		print(tb4)
		'''
		UNCOMMENT THE NEXT PART ONLY IF ENTIRELY SURE OF WHAT YOU ARE DOING!!!
		'''
		'''
#		for key in ['1.ppp','1.ppa','b_1','4.aaa','4.aap','b_4']:
#			df[key].ix[df[key].tidx==t1ppp[0],'X']=np.nan
#			df[key].ix[df[key].tidx==t1ppp[0],'Y']=np.nan
#			df[key].ix[df[key].tidx==t1ppp[0],'Z']=np.nan
#		
#		pickle.dump( df, open( worm, 'wb' ), protocol=2 )
		'''
		
		### check that nans have been properly inserted
#		print('\n\n\n')
#		t1p = df['1.p'].ix[pd.notnull(df['1.p'].X),'tidx'].values
#		t1pp = df['1.pp'].ix[pd.notnull(df['1.pp'].X),'tidx'].values
#		t1ppp = df['1.ppp'].ix[pd.notnull(df['1.ppp'].X),'tidx'].values
#		t1ppa = df['1.ppa'].ix[pd.notnull(df['1.ppa'].X),'tidx'].values
#		tb1 = df['b_1'].ix[pd.notnull(df['b_1'].X),'tidx'].values
#		
#		print(t1p)
#		print(t1pp)
#		print(t1ppp)
#		print(t1ppa)
#		print(tb1)
