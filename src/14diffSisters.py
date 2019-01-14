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

def plotFluorescence( path, worm, ax1, color = 'k', channel = '488nm' ):
	print( path, worm )

	timesDF = load_data_frame( path, worm + '_01times.pickle' )
	gonadPosDF = load_data_frame( path, worm + '_02gonadPos.pickle' )
	cellPosDF = load_data_frame( path, worm + '_04cellPos.pickle' )
	cellOutDF = load_data_frame( path, worm + '_05cellOut.pickle' )
	cellFluoDF = load_data_frame( path, worm + '_06cellFluo.pickle' )

	### find ecdysis timepoint

	ecd = np.loadtxt( open( os.path.join( path, 'skin.txt'), 'rb' ) )
	# load ecdysis data
	index = np.where( ecd[:,0] == float(worm[1:]) )
	mintp = np.min( [ i for i in ecd[index, 1:6][0][0] if i >= 0 ] )
	lethtidx = ecd[ index, 2:6 ][0][0] - 1
	tpL2 = timesDF.ix[ timesDF.tidxRel == lethtidx[0], 'timesRel' ].values[0]
	tpL1 = timesDF.ix[ timesDF.tidxRel == mintp, 'timesRel' ].values[0]

	### plot the timeseries
	for idx, key in enumerate( ['1.pp','4.aa'] ):

		cell1 = cellFluoDF[ key+'p' ].ix[ pd.notnull( cellFluoDF[ key+'p' ]._488nm ) & pd.notnull( cellFluoDF[ key+'a' ]._488nm ) ]
		cell2 = cellFluoDF[ key+'a' ].ix[ pd.notnull( cellFluoDF[ key+'a' ]._488nm ) & pd.notnull( cellFluoDF[ key+'p' ]._488nm ) ]
		times = timesDF.ix[ pd.notnull( cellFluoDF[ key+'p' ]._488nm ) & pd.notnull( cellFluoDF[ key+'a' ]._488nm ) ]
		tpRel = np.max( timesDF.ix[ pd.notnull( cellPosDF[ key ].X ) ].timesRel )

		# ### perform a gaussian filter
		# tidx = times.tidxRel
		# x_eval = times.timesRel
		# sigma = 20/60
		# delta_x = x_eval[:,None] - times.timesRel.values
		# weights = np.exp(-delta_x*delta_x / (2*sigma*sigma)) / (np.sqrt(2*np.pi) * sigma)
		# weights /= np.sum(weights, axis=1, keepdims=True)
		# val = np.dot(weights, cell1._488nm )
		# cell1._488nm = val

		# tidx = times.tidxRel
		# x_eval = times.timesRel
		# sigma = 20/60
		# delta_x = x_eval[:,None] - times.timesRel.values
		# weights = np.exp(-delta_x*delta_x / (2*sigma*sigma)) / (np.sqrt(2*np.pi) * sigma)
		# weights /= np.sum(weights, axis=1, keepdims=True)
		# val = np.dot(weights, cell2._488nm )
		# cell2._488nm = val


		# with background correction
		if key == '1.pp':
			ax1.plot( times.timesRel-tpRel, 
					( cell1._488nm - cell2._488nm ),
					'-k', ms = 2 )
		if key == '4.aa':
			ax1.plot( times.timesRel-tpRel, 
					( cell2._488nm - cell1._488nm ),
					'-k', ms = 2 )

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
	worms = ['X:\\Simone\\160930_lag2GFP+mCherry\\C04']
	### outlayers
	# worms = ['W:\\Simone\\160930_lag2GFP+mCherry\\C09','W:\\Simone\\160930_lag2GFP+mCherry\\C10','W:\\Simone\\160930_lag2GFP+mCherry\\C11','W:\\Simone\\160930_lag2GFP+mCherry\\C15','W:\\Simone\\160930_lag2GFP+mCherry\\C19','W:\\Simone\\161002_lag2GFP+mCherry\\C01']
	### regular worms
	# worms = ['W:\\Simone\\160930_lag2GFP+mCherry\\C01','W:\\Simone\\160930_lag2GFP+mCherry\\C02','W:\\Simone\\160930_lag2GFP+mCherry\\C03','W:\\Simone\\160930_lag2GFP+mCherry\\C06','W:\\Simone\\160930_lag2GFP+mCherry\\C08','W:\\Simone\\161002_lag2GFP+mCherry\\C04']
	### plot all
	# paths = [ 'X:\\Simone\\160930_lag2GFP+mCherry', 'X:\\Simone\\161002_lag2GFP+mCherry' ]
	# worms1 = [ glob.glob( os.path.join( paths[0], '*cellFluo.pickle' ) ), glob.glob( os.path.join( paths[1], '*cellFluo.pickle' ) ) ]
	# worms = []
	# for exp in worms1:
	# 	for w in exp:
	# 		worms.append(w)
	# worms.sort()

	for idx, worm in enumerate( worms ):

		w = worm.split('\\')[-1].split('_')[0]
		path = worm[:-len(worm.split('\\')[-1])]
				
		plotFluorescence( path, w, ax1, channel = '488nm' )

	plt.show()
