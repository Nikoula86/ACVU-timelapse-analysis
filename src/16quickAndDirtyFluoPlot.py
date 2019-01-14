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

def plotFluorescence( path, worm, ax1, color = 'k', channel = '488nm', lineages = [['1.p','1.pp','1.ppp'],['4.a','4.aa','4.aaa']], sisters = True, filt = False, raw = False, sigma = 20/60 ):
	print( path, worm )

	timesDF = load_data_frame( path, worm + '_01times.pickle' )
	gonadPosDF = load_data_frame( path, worm + '_02gonadPos.pickle' )
	cellPosDF = load_data_frame( path, worm + '_04cellPos.pickle' )
	cellFluoDF = load_data_frame( path, worm + '_15cellFluo_q&d.pickle' )
	
	colors = np.array( [ np.array( [cm.Blues(int(i))[:3] for i in (np.arange(len(lineages[0]))*254./(len(lineages[0])-1) / 254. ) * 127 + 127] ),
				np.array( [cm.Reds(int(i))[:3] for i in (np.arange(len(lineages[0]))*254./(len(lineages[0])-1) / 254. ) * 127 + 127 ] ) ] )

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

	### set relative time

	# # relative to L2 start
	tpRel = tpL2

	# relative to hatching time
	# tpRel=tpL1

	# relative to first cell division
	# tdiv1 = np.min( timesDF.ix[ pd.notnull( cellPosDF[ '1.ppp' ].X ) ].timesRel )
	# tdiv2 = np.min( timesDF.ix[ pd.notnull( cellPosDF[ '4.aaa' ].X ) ].timesRel )
	# tpRel = np.min([tdiv1,tdiv2])

	# print(tpRel)

	for idx, lin in enumerate( lineages ):
		for jdx, key in enumerate( lin ):
			# extract cell info and bckg info
			bckgData = cellFluoDF['b_'+key[0]].ix[pd.notnull(cellFluoDF[key]['_'+channel])]
			cellData = cellFluoDF[key].ix[pd.notnull(cellFluoDF[key]['_'+channel])]
			# correct fluorescence for background
			cellData['_'+channel] -= bckgData['_'+channel]
			# calculate time relative to whatever relative time is
			cellData.times -= tpRel

			# plot raw Data
			if raw:
				ax1.plot( cellData.times, 
							cellData['_'+channel],
							'o', color = colors[idx,jdx], mew=0, alpha = .6, ms = 4 )

			# plot gaussian filtered data
			if filt:
				cellData.expFilt = np.nan
				x_eval = cellData.times.values

				delta_x = x_eval[:,None] - cellData.times.values
				weights = np.exp(-delta_x*delta_x / (2*sigma*sigma)) / (np.sqrt(2*np.pi) * sigma)
				weights /= np.sum(weights, axis=1, keepdims=True)
				val = np.dot(weights, cellData['_'+channel].values )

				cellData.expFilt = val

				ax1.plot( cellData.times, 
							cellData.expFilt,
							'-', color = colors[idx,jdx], lw=2 )

			# plot sister data
			if sisters and len(key) == 5:
				sisName = '1.ppa'
				if key == '4.aaa': sisName = '4.aap'

				bckgData = cellFluoDF['b'+sisName[0]].ix[pd.notnull(cellFluoDF[sisName]['_'+channel])]
				sisData = cellFluoDF[sisName].ix[pd.notnull(cellFluoDF[sisName]['_'+channel])]
				sisData['_'+channel] -= bckgData['_'+channel]
				sisData['times'] -= tpRel

				# plot raw Data
				if raw:
					ax1.plot( sisData.times, 
								sisData['_'+channel],
								'o', color = colors[idx,jdx], mew=0, alpha = .6 )

				# plot gaussian filtered data
				if filt:
					### perform a gaussian filter
					sisData.expFilt = np.nan
					x_eval = sisData.times.values

					delta_x = x_eval[:,None] - sisData.times.values
					weights = np.exp(-delta_x*delta_x / (2*sigma*sigma)) / (np.sqrt(2*np.pi) * sigma)
					weights /= np.sum(weights, axis=1, keepdims=True)
					val = np.dot(weights, sisData['_'+channel].values )

					sisData.expFilt = val

					ax1.plot( sisData.times, 
								sisData.expFilt,
								'--', dashes = [3,2], color = colors[idx,jdx], lw=2 )

	## plot ecdysis
	ax1.plot( [tpL2 - tpRel, tpL2 - tpRel], [0,4000], '--', color = 'black', lw=1)
	if len(lethtidx) > 2:
		ax1.plot( [tpL3 - tpRel, tpL3 - tpRel], [0,4000], '--', color = 'black', lw=1)




if __name__ == '__main__':

	### setup figure for the timeseries
	fig1 = plt.figure(figsize=(5.8,2.9))
	ax1 = fig1.add_subplot(111)
	fig1.subplots_adjust(left=0.15, right=.95, top=.95, bottom=0.15)
	for tl in ax1.get_xticklabels():
		tl.set_fontsize(18)
	for tl in ax1.get_yticklabels():
		tl.set_fontsize(18)
	ax1.set_ylim((0,3500))

	# ax1.set_xlim((9,25))
	ax1.set_xlim(-2,10)
	ax1.set_xlabel('Time after L2 [hours]')
	ax1.set_ylabel('Fluorescence')

	'''
	lag2YFP data
	'''
	### plot one worm
	worms = [ 'Y:\\Simone\\161030_lag2YFP+histmCherry\\C23' ]
	
#	### plot all available worms
#	paths = [ 
#			'Y:\\Simone\\161030_lag2YFP+histmCherry',
#			'Y:\\Simone\\161108_lag2YFP+histmCherry', 
#			'Y:\\Simone\\161111_lag2YFP+histmCherry', 
#			'Y:\\Simone\\170312_lag2YFP+histmCherry'
#			]
#	worms1 = []
#	for path in paths:
#		worms1.append( glob.glob( os.path.join( path, '*cellFluo_q&d.pickle' ) ) )
#	worms = []
#	for exp in worms1:
#		for w in exp:
#			worms.append(w)
#	worms.sort()

	# ### plot every worm in a separate figure

	# figs = []
	# for idx, worm in enumerate( worms ):
	# 	plt.close(fig1)
	
	# 	figs.append(plt.figure(figsize=(7.8,5.8)))
	# 	ax1 = figs[-1].add_subplot(111)
	# 	ax1.set_xlabel('Time after L2 [hours]',fontsize = 18)
	# 	ax1.set_ylabel('Fluorescence',fontsize = 18)
	# 	figs[-1].subplots_adjust(left=0.15, right=.95, top=.95, bottom=0.15)
	# 	for tl in ax1.get_xticklabels():
	# 		tl.set_fontsize(18)
	# 	for tl in ax1.get_yticklabels():
	# 		tl.set_fontsize(18)

	# 	w = worm.split('\\')[-1].split('_')[0]
	# 	path = worm[:-len(worm.split('\\')[-1])]
				
	# 	plotFluorescence( path, w, ax1, channel = '488nm', filt = False, sigma = 15/60, lineages = [['1.p','1.pp','1.ppp'],['4.a','4.aa','4.aaa']] )
	# 	ax1.set_title( 'worm number: ' + str(idx) )
	# 	plt.savefig('W:\\Nicola\\timelapse_data\\rawPlots\\'+str(idx)+'.pdf')
		# plt.show()
	
	# plt.show()

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
	lin12nullUnc data
	'''
	### plot one worm
#	worms = [ 'Y:\\Simone\\170907_LIN12unc+histmCherry\\C06' ]

	'''
	HLH2GFP data
	'''
	### compute one worm
	# worms = ['X:\\Simone\\160407_HLH2_GFP_hist_mCherry\\C10']

	'''
	lag2 integration line data
	'''
#	### plot one worm
#	worms = [ 'Y:\\Simone\\171028_lag2int+JVZ32\\C01' ]
	### plot all available worms
#	paths = [ 
#			'Y:\\Simone\\171028_lag2int+JVZ32',
#			]
#	worms1 = []
#	for path in paths:
#		worms1.append( glob.glob( os.path.join( path, '*cellFluo_q&d.pickle' ) ) )
#	worms = []
#	for exp in worms1:
#		for w in exp:
#			worms.append(w)
#	worms.sort()

	'''
	PLOT THE FLUORESCENCE
	'''

	for idx, worm in enumerate( worms ):

		w = worm.split('\\')[-1].split('_')[0]
		path = worm[:-len(worm.split('\\')[-1])]
				
		plotFluorescence( path, w, ax1, channel = '488nm', filt = True, raw = True, sigma = 15/60, lineages = [['1.p','1.pp','1.ppp'],['4.a','4.aa','4.aaa']], sisters = False )
#		plotFluorescence( path, w, ax1, channel = '561nm', filt = True, raw = True, sigma = 15/60, lineages = [['1.p','1.pp','1.ppp'],['4.a','4.aa','4.aaa']], sisters = False )
		ax1.set_title( 'worm number: ' + str(idx) )

	plt.show()
