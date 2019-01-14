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
from matplotlib.colors import LinearSegmentedColormap


import matplotlib as mpl
mpl.rcParams['pdf.fonttype'] = 42

def computeData( path, worm, dT = 10, channel = '488nm', color = 'black' ):
	print( path, worm )

	timesDF = load_data_frame( path, worm + '_01times.pickle' )
	gonadPosDF = load_data_frame( path, worm + '_02gonadPos.pickle' )
	cellPosDF = load_data_frame( path, worm + '_04cellPos.pickle' )
	cellFluoDF = load_data_frame( path, worm + '_15cellFluo_q&d.pickle' )

	### check if they divided at similar times
	tDiv1 = np.min( timesDF.ix[ pd.notnull( cellPosDF[ '1.ppp' ].X ) ].timesRel )
	tDiv4 = np.min( timesDF.ix[ pd.notnull( cellPosDF[ '4.aaa' ].X ) ].timesRel )
	tRel = np.max( [ tDiv1, tDiv4 ] )
	# if np.abs( tDiv1 - tDiv4 ) >= 1:
	# 	return

	XYdata = pd.DataFrame({ 'dTdiv':tDiv1-tDiv4  ,'1.ppp':0, '1.ppa':0, '4.aaa':0, '4.aap':0 }, index=[0])
	for idx, key in enumerate( [['1.ppp','1.ppa'],['4.aaa','4.aap']] ):

		### extract all info
		cell = cellFluoDF[ key[0] ].ix[ pd.notnull( cellFluoDF[ key[0] ]._488nm ) & pd.notnull( cellFluoDF[ key[0] ]._488nm ) ]
		sis = cellFluoDF[ key[1] ].ix[ pd.notnull( cellFluoDF[ key[1] ]._488nm ) & pd.notnull( cellFluoDF[ key[1] ]._488nm ) ]

		### perform a gaussian filter
		# sigma = 20/60
		# delta_x = sis.times[:,None] - sis.times.values
		# weights = np.exp(-delta_x*delta_x / (2*sigma*sigma)) / (np.sqrt(2*np.pi) * sigma)
		# weights /= np.sum(weights, axis=1, keepdims=True)
		# sis._488nm = np.dot(weights, sis._488nm )

		# sigma = 20/60
		# delta_x = cell.times[:,None] - cell.times.values
		# weights = np.exp(-delta_x*delta_x / (2*sigma*sigma)) / (np.sqrt(2*np.pi) * sigma)
		# weights /= np.sum(weights, axis=1, keepdims=True)
		# cell._488nm = np.dot(weights, cell._488nm )

		### find slope
		tDiv = np.min( timesDF.ix[ pd.notnull( cellPosDF[ key[0] ].X ) ].timesRel )
		tIdxMin = np.min( timesDF.ix[ pd.notnull( cellPosDF[ key[0] ].X ) ].tidxRel )
		tIdxMax = np.max( timesDF.ix[ timesDF.timesRel <= ( tDiv + 1.5 ) ].tidxRel )

		sis = sis.ix[ ( sis.tidx >= tIdxMin ) & ( sis.tidx <= tIdxMax ) ]

		slope = np.polyfit(sis.times,sis._488nm,1)[0]

		### find outcome after N hours
		if ( tRel + dT + 1 ) < np.max( cell.times ):
			tMask = ( cell.times >= ( tRel + dT ) ) & ( cell.times < ( tRel + dT + 1 ) )
			exp = np.mean( ( cell.ix[ tMask, '_488nm' ] ).values[0] )
		else:
			tMask = ( cell.times == np.max(cell.times) )
			exp = np.mean( ( cell.ix[ tMask, '_488nm' ] ).values[0] )


		### store data
		XYdata[key[0]] = exp
		XYdata[key[1]] = slope

	print(XYdata)

	return XYdata


if __name__ == '__main__':

	### setup figure for the timeseries
	fig1 = plt.figure(figsize=(5.8,5.8))
	ax1 = fig1.add_subplot(111)
	fig1.subplots_adjust(left=0.15, right=.95, top=.95, bottom=0.15)
	for tl in ax1.get_xticklabels():
		tl.set_fontsize(18)
	for tl in ax1.get_yticklabels():
		tl.set_fontsize(18)

	ax1.set_ylim((-750,750))
	ax1.set_xlim((-2,2))
	ax1.plot([0,0],[-1000,1000],'--k')
	ax1.plot([-3,3],[0,0],'--k')

	ax1.set_ylabel( 'Sister Slope Difference', fontsize =18 )
	ax1.set_xlabel( 'Division Time diff (TdivZ1 - TdivZ4)', fontsize =18 )
	ax1.text( 2.5,400,'Outcome (Z1-Z4)/(Z1+Z4)', fontsize = 18, rotation = -90)

	'''
	lag2YFP data
	'''
	# ### plot one worm
	worms = [ 'W:\\Simone\\160516_lag2_YFP_hist_mCherry\\C21' ]
	### plot all available worms
	paths = [ 'W:\\Simone\\161030_lag2YFP+mCherry', 'W:\\Simone\\161108_lag2YFP+mCherry', 'W:\\Simone\\161111_lag2YFP+mCherry' ]
	worms1 = []
	for path in paths:
		worms1.append( glob.glob( os.path.join( path, '*cellFluo*.pickle' ) ) )
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
	# paths = [ 'X:\\Simone\\160930_lag2GFP+mCherry', 'X:\\Simone\\161002_lag2GFP+mCherry', 'X:\\Simone\\161004_lag2GFP+mCherry' ]
	# worms1 = []
	# for path in paths:
	# 	worms1.append( glob.glob( os.path.join( path, '*cellFluo.pickle' ) ) )
	# worms = []
	# for exp in worms1:
	# 	for w in exp:
	# 		worms.append(w)
	# worms.sort()

	cdict1 = {'red':   ((0.0, 1.0, 1.0),
				(0.4, 0.8, 0.8),
				(0.5, 0.0, 0.0),
				(1.0, 0.0, 0.0)),

	     'green': ((0.0, 0.0, 0.0),
				(1.0, 0.0, 0.0)),

	     'blue':  ((0.0, 0.0, 0.0),
				(0.5, 0.0, 0.0),
				(0.6, 0.8, 0.8),
				(1.0, 1.0, 1.0))
	    }
	blue2red = LinearSegmentedColormap('mycolormap', cdict1)

	data = pd.DataFrame({})
	for idx, worm in enumerate( worms ):

		w = worm.split('\\')[-1].split('_')[0]
		path = worm[:-len(worm.split('\\')[-1])]
		
		wormData = computeData( path, w, channel = '488nm', dT=10)
		data = pd.concat( [ data, wormData ] )
		ax1.text( wormData.dTdiv, wormData['1.ppa']-wormData['4.aap'], str(idx), fontsize = 15 )

	img = ax1.scatter( data.dTdiv, data['1.ppa']-data['4.aap'], c = (data['1.ppp']-data['4.aaa'])/(data['1.ppp']+data['4.aaa']), vmin = -1, vmax = 1, cmap = blue2red, s = 100, lw=0 )
	print('\n')
	print(((data['1.ppp']-data['4.aaa'])/(data['1.ppp']+data['4.aaa'])).values)
	plt.colorbar(img)
	plt.show()
