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

def plotFluorescence( path, worm, ax1, color = 'k', channel = '488nm', lineages = [['1.p','1.pp','1.ppp','1.ppa'],['4.a','4.aa','4.aaa','4.aap']], filt = False, sigma = 20/60 ):
	print( path, worm )

	timesDF = load_data_frame( path, worm + '_01times.pickle' )
	gonadPosDF = load_data_frame( path, worm + '_02gonadPos.pickle' )
	cellPosDF = load_data_frame( path, worm + '_04cellPos.pickle' )
	cellOutDF = load_data_frame( path, worm + '_05cellOut.pickle' )
	cellFluoDF = load_data_frame( path, worm + '_06cellFluo.pickle' )

	colors = np.array( [ np.array( [cm.Blues(int(i))[:3] for i in (np.arange(len(lineages[0]))*254./(len(lineages[0])-1) / 254. ) * 127 + 127] ),
				np.array( [cm.Reds(int(i))[:3] for i in (np.arange(len(lineages[0]))*254./(len(lineages[0])-1) / 254. ) * 127 + 127 ] ) ] )

	### find ecdysis timepoint

	ecd = np.loadtxt( open( os.path.join( path, 'skin.txt'), 'rb' ) )
	# load ecdysis data
	index = np.where( ecd[:,0] == float(worm[1:]) )
	mintp = np.min( [ i for i in ecd[index, 1:6][0][0] if i >= 0 ] )
	lethtidx = ecd[ index, 2:6 ][0][0] - 1
	tpL2 = timesDF.ix[ timesDF.tidxRel == lethtidx[0], 'timesRel' ].values[0]
	tpL1 = timesDF.ix[ timesDF.tidxRel == mintp, 'timesRel' ].values[0]

	darkField = load_stack( 'W:\\Orca_calibration\\AVG_darkField.tif' )
	flatField = load_stack( os.path.join( path, 'AVG_flatField_'+channel+'.tif' ) )
	### set relative time

	# # relative to L2 start
	# tpRel = tpL2

	# relative to hatching time
	# tpRel=tpL1

	# relative to first cell division
	tdiv1 = np.min( timesDF.ix[ pd.notnull( cellPosDF[ '1.ppp' ].X ) ].timesRel )
	tdiv2 = np.min( timesDF.ix[ pd.notnull( cellPosDF[ '4.aaa' ].X ) ].timesRel )
	tpRel = np.min([tdiv1,tdiv2])

	print(tpRel)

	### plot the timeseries
	for idx, lin in enumerate( lineages ):

		for jdx, key in enumerate( lin ):
			print(key)

			if channel == '488nm':
				cell = cellFluoDF[ key ].ix[ pd.notnull( cellFluoDF[ key ]._488nm ) ]
				bckg = cellFluoDF[ 'b_'+key[0] ].ix[ pd.notnull( cellFluoDF[ key ]._488nm ) ]
				times = timesDF.ix[ pd.notnull( cellFluoDF[ key ]._488nm ) ]

				if len(key) == 5:
					if key == '1.ppp':
						cell = cellFluoDF[ key ].ix[ pd.notnull( cellFluoDF[ key ]._488nm ) & pd.notnull( cellFluoDF[ key[:-1] + 'a' ]._488nm ) ]
						bckg = cellFluoDF[ 'b_'+key[0] ].ix[ pd.notnull( cellFluoDF[ key ]._488nm ) & pd.notnull( cellFluoDF[ key[:-1] + 'a' ]._488nm ) ]
						times = timesDF.ix[ pd.notnull( cellFluoDF[ key ]._488nm ) & pd.notnull( cellFluoDF[ key[:-1] + 'a' ]._488nm ) ]	
						sis = cellFluoDF[ key[:-1] + 'a' ].ix[ pd.notnull( cellFluoDF[ key ]._488nm ) & pd.notnull( cellFluoDF[ key[:-1] + 'a' ]._488nm ) ]
						cell._488nm -= ( sis._488nm )
					if key == '4.aaa':
						cell = cellFluoDF[ key ].ix[ pd.notnull( cellFluoDF[ key ]._488nm ) & pd.notnull( cellFluoDF[ key[:-1] + 'p' ]._488nm ) ]
						bckg = cellFluoDF[ 'b_'+key[0] ].ix[ pd.notnull( cellFluoDF[ key ]._488nm ) & pd.notnull( cellFluoDF[ key[:-1] + 'p' ]._488nm ) ]
						times = timesDF.ix[ pd.notnull( cellFluoDF[ key ]._488nm ) & pd.notnull( cellFluoDF[ key[:-1] + 'p' ]._488nm ) ]	
						sis = cellFluoDF[ key[:-1] + 'p' ].ix[ pd.notnull( cellFluoDF[ key ]._488nm ) & pd.notnull( cellFluoDF[ key[:-1] + 'p' ]._488nm ) ]
						cell._488nm -= ( sis._488nm )

				if filt:
					### perform a gaussian filter
					cell._488nmFilt = np.nan
					tidx = times.tidxRel
					x_eval = times.timesRel

					delta_x = x_eval[:,None] - times.timesRel.values
					weights = np.exp(-delta_x*delta_x / (2*sigma*sigma)) / (np.sqrt(2*np.pi) * sigma)
					weights /= np.sum(weights, axis=1, keepdims=True)
					val = np.dot(weights, cell._488nm )

					cell._488nmFilt = val
				else:
					cell._488nmFilt = cell._488nm

				# with background correction
				ax1.plot( times.timesRel-tpRel, 
							( cell._488nmFilt ),# / bckg._488nm, 
							'-', color = colors[idx,jdx], lw=2 )

if __name__ == '__main__':

	### setup figure for the timeseries
	fig1 = plt.figure(figsize=(5.8,3.8))
	ax1 = fig1.add_subplot(111)
	fig1.subplots_adjust(left=0.15, right=.95, top=.95, bottom=0.15)
	for tl in ax1.get_xticklabels():
		tl.set_fontsize(18)
	for tl in ax1.get_yticklabels():
		tl.set_fontsize(18)
	# ax1.set_ylim((0,3.5))

	# ax1.set_xlim((9,25))
	# ax1.set_xlim(-5,15)

	ax1.plot( [0,0], [0,4000], '--', color = 'black', lw=1)

	'''
	lag2YFP data
	'''
	# ### plot one worm
	# worms = [ 'W:\\Simone\\160516_lag2_YFP_hist_mCherry\\C21' ]
	# ### plot all available worms
	# paths = [ 'W:\\Simone\\160516_lag2_YFP_hist_mCherry', 'W:\\Simone\\160226_lag2_YFP_hist_mCherry' ]
	# worms1 = [ glob.glob( os.path.join( paths[0], '*cellFluo.pickle' ) ), glob.glob( os.path.join( paths[1], '*cellFluo.pickle' ) ) ]
	# worms = []
	# for exp in worms1:
	# 	for w in exp:
	# 		worms.append(w)
	# worms.sort()

	'''
	lag2GFP data
	'''
	### plot one worm
	worms = ['X:\\Simone\\161002_lag2GFP+mCherry\\C25']
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

	for idx, worm in enumerate( worms ):

		w = worm.split('\\')[-1].split('_')[0]
		path = worm[:-len(worm.split('\\')[-1])]
				
		plotFluorescence( path, w, ax1, channel = '488nm', filt = True, sigma = 20/60, lineages = [['1.pp','1.ppp','1.ppa'],['4.aa','4.aaa','4.aap']] )

		
	plt.show()
