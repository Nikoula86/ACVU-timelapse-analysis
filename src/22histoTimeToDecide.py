# -*- coding: utf-8 -*-
"""
Created on Tue Dec  1 10:31:00 2015

@author: kienle
"""

import numpy as np
import matplotlib.pyplot as plt
import pickle
from matplotlib import cm
from generalFunctions import *
from skimage import filters


import matplotlib as mpl
mpl.rcParams['pdf.fonttype'] = 42

def findTDecision( path, worm, ax1, filt = False, sigma = 20/60 ):

	timesDF = load_data_frame( path, worm + '_01times.pickle' )
	gonadPosDF = load_data_frame( path, worm + '_02gonadPos.pickle' )
	cellPosDF = load_data_frame( path, worm + '_04cellPos.pickle' )
	cellFluoDF = load_data_frame( path, worm + '_15cellFluo_q&d.pickle' )
	
	### extract division times
	tdiv1 = np.min( timesDF.ix[ pd.notnull( cellPosDF[ '1.ppp' ].X ) ].timesRel )
	tdiv2 = np.min( timesDF.ix[ pd.notnull( cellPosDF[ '4.aaa' ].X ) ].timesRel )

	### find ecdysis timepoint
	ecd = np.loadtxt( open( os.path.join( path, 'skin.txt'), 'rb' ) )
	# load ecdysis data
	index = np.where( ecd[:,0] == float(worm[1:]) )
	mintp = np.min( [ i for i in ecd[index, 1:6][0][0] if i >= 0 ] )
	lethtidx = ecd[ index, 1:6 ][0][0]
	lethtidx = lethtidx[ lethtidx>=0 ]
	tpL1 = timesDF.ix[ timesDF.tidxRel == 0, 'timesRel' ].values[0]
	tpL2 = timesDF.ix[ timesDF.tidxRel == ( lethtidx[1] - mintp ), 'timesRel' ].values[0]
	if len(lethtidx) > 2:
		tpL3 = timesDF.ix[ timesDF.tidxRel == ( lethtidx[2] - mintp ), 'timesRel' ].values[0]	

	# # relative to L2 start
	tpRel = tpL2
	# relative to first division
	tpRel = np.min( [ tdiv1, tdiv2 ] )

	# extract data from dataframe
	cellData = pd.DataFrame({})
	index = 0
	for idx, tRow in timesDF.iterrows():
		# print(tRow.tidxRel)

		#extract current cells
		currentCells = extract_current_cell_pos( cellPosDF, tRow.tidxRel )

		if len(currentCells) > 0:
			(c1,c4,b1,b4) = find_interesting_cells(currentCells)

			# print(c1.cname,c4.cname)
			newRow = pd.DataFrame( {'cell1':cellFluoDF[c1.cname].ix[ cellFluoDF[c1.cname].times == tRow.timesRel, '_488nm' ].values[0]/1000.,
				'cell4':cellFluoDF[c4.cname].ix[ cellFluoDF[c4.cname].times == tRow.timesRel, '_488nm' ].values[0]/1000.,
				'times':tRow.timesRel}, index= [index] )
			cellData = pd.concat( [cellData,newRow] )
			index += 1


	cellData = cellData.ix[ cellData.times >= np.max( [ tdiv1, tdiv2 ] ) ]
	# perform gaussian filter
	if filt:

		for val in zip(['cell1','cell4'],[tdiv1,tdiv2]):
			key = val[0]
			tdiv = val[1]

			### perform a gaussian filter on mother
			filtData = np.nan
			x_eval = cellData.ix[cellData.times >= tdiv, 'times'].values

			delta_x = x_eval[:,None] - cellData.ix[cellData.times >= tdiv, 'times'].values
			weights = np.exp(-delta_x*delta_x / (2*sigma*sigma)) / (np.sqrt(2*np.pi) * sigma)
			weights /= np.sum(weights, axis=1, keepdims=True)
			val = np.dot(weights, cellData.ix[cellData.times >= tdiv,key] )

			cellData.ix[cellData.times >= tdiv,key] = val

	cellData.ratio = np.abs( ( cellData.cell1 - cellData.cell4 ) / ( cellData.cell1 + cellData.cell4 ) )

	### find time to decide
	# print(tdiv1,tdiv2)
	thr = 1/3
	tDec = cellData.ix[ cellData.ratio < thr, 'times' ]
	if len( tDec ) > 0:
		tDec = tDec.values[-1]
	else:
		tDec = np.max( [ tdiv1,tdiv2 ] )
	print( np.abs( tdiv1 - tdiv2 ), tDec - np.max( [ tdiv1,tdiv2 ] ) )

	return ( np.abs( tdiv1 - tdiv2 ), tDec - np.max( [ tdiv1,tdiv2 ] ) )

if __name__ == '__main__':

	### setup figure for the phasePlot
	fig1 = plt.figure(figsize=(6.8,4.8))
	ax1 = fig1.add_subplot(111)
	fig1.subplots_adjust(left=0.15, right=.95, top=.95, bottom=0.15)
	for tl in ax1.get_xticklabels():
		tl.set_fontsize(18)
	for tl in ax1.get_yticklabels():
		tl.set_fontsize(18)
	# ax1.set_ylim((0,3.5))

	ax1.set_xlim(0,120/60.)
	ax1.set_ylim(0,10)

	# ax1.plot( [-2,10], [0,0], '--', color = 'black', lw=1)
	# ax1.plot( [-2,10], [1/3,1/3], '--', color = 'black', lw=1)
	# ax1.plot( [-2,10], [-1/3,-1/3], '--', color = 'black', lw=1)

	ax1.set_xlabel('dt Div', fontsize = 18)
	ax1.set_ylabel('t Decision', fontsize = 18)

	'''
	lag2YFP data
	'''
	### plot one worm
	#worms = [ 'W:\\Simone\\161111_lag2YFP+mCherry\\C17' ]

	### plot all available worms
	paths = [ 'Y:\\Simone\\161030_lag2YFP+histmCherry', 'Y:\\Simone\\161108_lag2YFP+mCherry','Y:\\Simone\\161111_lag2YFP+mCherry', 'Y:\\Simone\\170312_lag2YFP+histmCherry' ]
	worms1 = []
	for path in paths:
		worms1.append( glob.glob( os.path.join( path, '*cellFluo_q&d.pickle' ) ) )
	worms = []
	for exp in worms1:
		for w in exp:
			worms.append(w)
	worms.sort()

	# # ### plot every worm in a separate figure

	# figs = []
	# for idx, worm in enumerate( worms ):
	# 	plt.close(fig1)
	
	# 	figs.append(plt.figure(figsize=(7.8,7.8)))
	# 	ax1 = figs[-1].add_subplot(111)
	# 	ax1.set_xlabel('Z1 fluorescence',fontsize = 18)
	# 	ax1.set_ylabel('Z4 fluorescence',fontsize = 18)
	# 	ax1.plot( [0,5], [0,5], '--', color = 'black', lw=1)
	# 	ax1.set_xlim(0,4)
	# 	ax1.set_ylim(0,4)
	# 	ax1.set_title( 'worm number: ' + str(idx) )
	# 	figs[-1].subplots_adjust(left=0.15, right=.95, top=.95, bottom=0.15)
	# 	for tl in ax1.get_xticklabels():
	# 		tl.set_fontsize(18)
	# 	for tl in ax1.get_yticklabels():
	# 		tl.set_fontsize(18)

	# 	w = worm.split('\\')[-1].split('_')[0]
	# 	path = worm[:-len(worm.split('\\')[-1])]
				
	# 	plotPhase( path, w, ax1, filt = True, sigma = 15/60 )
	# 	plt.savefig('W:\\Nicola\\timelapse_data\\phasePlots\\'+str(idx)+'.pdf')
	# 	ax1.set_title( 'worm number: ' + str(idx) )
	# 	# plt.show()

	'''
	lag2GFP data
	'''
	### plot one worm
	# worms = ['X:\\Simone\\161002_lag2GFP+mCherry\\C15']
	### outlayers
	# worms = ['W:\\Simone\\160930_lag2GFP+mCherry\\C09','W:\\Simone\\160930_lag2GFP+mCherry\\C10','W:\\Simone\\160930_lag2GFP+mCherry\\C11','W:\\Simone\\160930_lag2GFP+mCherry\\C15','W:\\Simone\\160930_lag2GFP+mCherry\\C19','W:\\Simone\\161002_lag2GFP+mCherry\\C01']
	### regular worms
	# worms = ['W:\\Simone\\160930_lag2GFP+mCherry\\C01','W:\\Simone\\160930_lag2GFP+mCherry\\C02','W:\\Simone\\160930_lag2GFP+mCherry\\C03','W:\\Simone\\160930_lag2GFP+mCherry\\C06','W:\\Simone\\160930_lag2GFP+mCherry\\C08','W:\\Simone\\161002_lag2GFP+mCherry\\C04']
	### plot all
	# paths = [ 'X:\\Simone\\160930_lag2GFP+mCherry', 'X:\\Simone\\161002_lag2GFP+mCherry' ]
	# worms1 = []
	# for path in paths:
	# 	worms1.append( glob.glob( os.path.join( path, '*cellFluo.pickle' ) ) )
	# worms = []
	# for exp in worms1:
	# 	for w in exp:
	# 		worms.append(w)
	# worms.sort()

	'''
	HLH2GFP data
	'''
	### compute one worm
	# worms = ['X:\\Simone\\160407_HLH2_GFP_hist_mCherry\\C10']

	'''
	PLOT THE RATIO
	'''

	dtDiv = []
	tDec = []
	for idx, worm in enumerate( worms ):
		print(idx, worm)

		w = worm.split('\\')[-1].split('_')[0]
		path = worm[:-len(worm.split('\\')[-1])]
				
		dt,tdec = findTDecision( path, w, ax1, filt = True, sigma = 15/60 )
		dtDiv.append( dt )
		tDec.append( tdec )

		# ax1.text( dt*60., tdec, str(idx) )

	df = pd.DataFrame( { 'dt': np.array(dtDiv)*60., 'tDec': tDec } )
	bins = np.linspace(0,130,14)
	groups = df.groupby(pd.cut(df.dt,bins - 5))
	_mean = groups.mean().tDec
	_std = groups.std().tDec
	print(bins,_mean,_std)

	ax1.scatter( df.dt/60., df.tDec, color = 'black' , s = 50, lw=0,alpha = 1. )
	# ax1.errorbar( bins[:-1], groups.mean().tDec, yerr = groups.std().tDec, fmt='o', color = 'black' )
		
	plt.show()
