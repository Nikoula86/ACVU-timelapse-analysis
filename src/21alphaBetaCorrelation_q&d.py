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


def performGaussianFilter( data, sigma ):
	filtData = np.nan
	x_eval = data.times.values

	delta_x = x_eval[:,None] - data.times.values
	weights = np.exp(-delta_x*delta_x / (2*sigma*sigma)) / (np.sqrt(2*np.pi) * sigma)
	weights /= np.sum(weights, axis=1, keepdims=True)
	val = np.dot(weights, data._488nm.values )

	return val


def plotCorrelation( worms, ax1, ax2, filt = False, sigma = 20/60 ):

	tStep = 20
	tBins = np.arange(10,601,tStep)
	data = [ [[],[]] for i in tBins ]
	print(tBins)

	for idx, worm in enumerate( worms ):

		w = worm.split('\\')[-1].split('_')[0]
		path = worm[:-len(worm.split('\\')[-1])]
			
		print( path, w )

		timesDF = load_data_frame( path, w + '_01times.pickle' )
		cellPosDF = load_data_frame( path, w + '_04cellPos.pickle' )
		cellFluoDF = load_data_frame( path, w + '_15cellFluo_q&d.pickle' )
	
		### extract division times
		tdiv = [ np.min( timesDF.ix[ pd.notnull( cellPosDF[ '1.ppp' ].X ) ].timesRel ),
				np.min( timesDF.ix[ pd.notnull( cellPosDF[ '4.aaa' ].X ) ].timesRel ) ]

		for jdx, cnames in enumerate( [['1.ppp','1.ppa'],['4.aaa','4.aap']] ):
			cname = cnames[0]
			sname = cnames[1]

			cell = cellFluoDF[ cname ].ix[ pd.notnull( cellFluoDF[ cname ]._488nm ) & pd.notnull( cellFluoDF[ sname ]._488nm ) ]
			sis = cellFluoDF[ sname ].ix[ pd.notnull( cellFluoDF[ sname ]._488nm ) & pd.notnull( cellFluoDF[ cname ]._488nm ) ]

			cell.times = ( cell.times - tdiv[jdx] ) * 60
			sis.times = ( sis.times - tdiv[jdx] ) * 60

			# perform gaussian filter
			if filt:
				cell._488nm = performGaussianFilter(cell,sigma)
				sis._488nm = performGaussianFilter(sis,sigma)

			for kdx, tcenter in enumerate( tBins ):
				cDataInBin = cell.ix[ ( cell.times >= ( tcenter - tStep / 2. ) ) & ( cell.times < ( tcenter + tStep / 2. ) ) ]
				sDataInBin = sis.ix[ ( sis.times >= ( tcenter - tStep / 2. ) ) & ( sis.times < ( tcenter + tStep / 2. ) ) ]

				if len(cDataInBin) > 0:
					[ data[ kdx ][0].append( val ) for val in cDataInBin._488nm ]
					[ data[ kdx ][1].append( val ) for val in sDataInBin._488nm ]

	data = np.array( data )
	# print(data)

	# plot mothers
	ax1.plot( data[5,0], 
				data[5,1], 
				'o', color = 'black' )

	corr = np.nan * np.zeros(len(tBins))
	for idx, d in enumerate( data ):
		corr[idx] = np.polyfit(d[0],d[1],1)[0]#np.corrcoef( d[0],d[1] )[1,0]
	# print(tBins/60,corr)
	ax2.plot(tBins/60.,corr)

	corr = np.nan * np.zeros(len(tBins))
	for idx, d in enumerate( data ):
		corr[idx] = np.corrcoef( d[0],d[1] )[1,0]
	# print(tBins/60,corr)
	ax2.plot(tBins/60.,corr)






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

	ax1.set_xlim(0,3000)
	ax1.set_ylim(0,3000)

	ax1.plot( [0,3000], [0,3000], '--', color = 'black', lw=1)

	ax1.set_xlabel('alpha fluorescence', fontsize = 18)
	ax1.set_ylabel('beta fluorescence', fontsize = 18)

	### setup figure for the phasePlot
	fig2 = plt.figure(figsize=(5.8,5.8))
	ax2 = fig2.add_subplot(111)
	fig2.subplots_adjust(left=0.15, right=.95, top=.95, bottom=0.15)
	for tl in ax2.get_xticklabels():
		tl.set_fontsize(18)
	for tl in ax2.get_yticklabels():
		tl.set_fontsize(18)
	# ax1.set_ylim((0,3.5))

	ax2.set_xlim(0,10)
	ax2.set_ylim(-1,1)

	ax2.plot( [0,10], [0,0], '--', color = 'black', lw=1)

	ax2.set_xlabel('Time [hours after division]', fontsize = 18)
	ax2.set_ylabel('Correlation', fontsize = 18)

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

	plotCorrelation( worms, ax1, ax2, filt = True, sigma = 15/60 )

	plt.show()
