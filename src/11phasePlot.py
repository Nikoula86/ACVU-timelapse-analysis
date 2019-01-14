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

def plotPhase( path, worm, ax1, filt = False, sigma = 20/60 ):
	print( path, worm )

	timesDF = load_data_frame( path, worm + '_01times.pickle' )
	gonadPosDF = load_data_frame( path, worm + '_02gonadPos.pickle' )
	cellPosDF = load_data_frame( path, worm + '_04cellPos.pickle' )
	cellFluoDF = load_data_frame( path, worm + '_15cellFluo_q&d.pickle' )
	
	### extract division times
	tdiv1 = np.min( timesDF.ix[ pd.notnull( cellPosDF[ '1.ppp' ].X ) ].timesRel )
	tdiv2 = np.min( timesDF.ix[ pd.notnull( cellPosDF[ '4.aaa' ].X ) ].timesRel )
	tlast = np.max( timesDF.ix[ pd.notnull( cellPosDF[ '4.aaa' ].X ) ].timesRel )

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


	# perform gaussian filter
	if filt:

		for val in zip(['cell1','cell4'],[tdiv1,tdiv2]):
			key = val[0]
			tdiv = val[1]

			### perform a gaussian filter on mother
			filtData = np.nan
			x_eval = cellData.ix[cellData.times < tdiv, 'times'].values

			delta_x = x_eval[:,None] - cellData.ix[cellData.times < tdiv, 'times'].values
			weights = np.exp(-delta_x*delta_x / (2*sigma*sigma)) / (np.sqrt(2*np.pi) * sigma)
			weights /= np.sum(weights, axis=1, keepdims=True)
			val = np.dot(weights, cellData.ix[cellData.times < tdiv,key] )

			cellData.ix[cellData.times < tdiv,key] = val

			### perform a gaussian filter on mother
			filtData = np.nan
			x_eval = cellData.ix[cellData.times >= tdiv, 'times'].values

			delta_x = x_eval[:,None] - cellData.ix[cellData.times >= tdiv, 'times'].values
			weights = np.exp(-delta_x*delta_x / (2*sigma*sigma)) / (np.sqrt(2*np.pi) * sigma)
			weights /= np.sum(weights, axis=1, keepdims=True)
			val = np.dot(weights, cellData.ix[cellData.times >= tdiv,key] )

			cellData.ix[cellData.times >= tdiv,key] = val


	### plot data
	# print(tdiv1,tdiv2)

	# plot mothers
	ax1.plot( cellData.ix[ cellData.times<=np.min([tdiv1,tdiv2]), 'cell1' ], 
				cellData.ix[ cellData.times<=np.min([tdiv1,tdiv2]), 'cell4' ], 
				'-', lw = 2, color = 'black' )

	# plot intermediates
	ax1.plot( cellData.ix[ (cellData.times>=np.min([tdiv1,tdiv2])) & (cellData.times<=np.max([tdiv1,tdiv2])), 'cell1' ],
				cellData.ix[ (cellData.times>=np.min([tdiv1,tdiv2])) & (cellData.times<=np.max([tdiv1,tdiv2])), 'cell4' ],
				'-', lw = 2, color = 'blue' )

	# plot daughters
	ax1.plot( cellData.ix[ cellData.times>=np.max([tdiv1,tdiv2]), 'cell1' ],
				cellData.ix[ cellData.times>=np.max([tdiv1,tdiv2]), 'cell4' ],
				'-', lw = 2, color = 'magenta' )

	### make clear which one divided first
	val1 = cellData.ix[cellData.times == np.min([tdiv1,tdiv2]),'cell1'].values[0]
	val4 = cellData.ix[cellData.times == np.min([tdiv1,tdiv2]),'cell4'].values[0]
	if tdiv1 < tdiv2:
		ax1.plot([val1,val1],[0,val4],'--k')
	elif tdiv2 < tdiv1:
		ax1.plot([0,val1],[val4,val4],'--k')
	else:
		ax1.plot([val1,val1],[0,val4],'--k')
		ax1.plot([0,val1],[val4,val4],'--k')





if __name__ == '__main__':

	### setup figure for the phasePlot
	fig1 = plt.figure(figsize=(5.8,5.8))
	ax1 = fig1.add_subplot(111)
	fig1.subplots_adjust(left=0.15, right=.95, top=.95, bottom=0.15)
	for tl in ax1.get_xticklabels():
		tl.set_fontsize(18)
	for tl in ax1.get_yticklabels():
		tl.set_fontsize(18)
	# ax1.set_ylim((0,3.5))

	ax1.set_xlim(0,3)
	ax1.set_ylim(0,3)

	ax1.plot( [0,5], [0,5], '--', color = 'black', lw=1)
	ax1.plot( [0,5], [0,10], '--', color = 'black', lw=1)
	ax1.plot( [0,10], [0,5], '--', color = 'black', lw=1)

	ax1.set_xlabel('Z1 fluorescence', fontsize = 18)
	ax1.set_ylabel('Z4 fluorescence', fontsize = 18)

	'''
	lag2YFP data
	'''
	### plot one worm
	worms = [ 'W:\\Simone\\161111_lag2YFP+mCherry\\C14' ]

	### plot all available worms
	paths = [ 'W:\\Simone\\161030_lag2YFP+mCherry', 'W:\\Simone\\161108_lag2YFP+mCherry','W:\\Simone\\161111_lag2YFP+mCherry' ]
	worms1 = []
	for path in paths:
		worms1.append( glob.glob( os.path.join( path, '*cellFluo_q&d.pickle' ) ) )
	worms = []
	for exp in worms1:
		for w in exp:
			worms.append(w)
	worms.sort()

	# ### plot every worm in a separate figure

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

	for idx, worm in enumerate( worms ):

		w = worm.split('\\')[-1].split('_')[0]
		path = worm[:-len(worm.split('\\')[-1])]
				
		plotPhase( path, w, ax1, filt = True, sigma = 15/60 )

	plt.show()
