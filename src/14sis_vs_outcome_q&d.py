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

def plotFluorescence( path, worm, ax1, channel = '488nm', color = 'black', widx = 0 ):
	print( path, worm )

	timesDF = load_data_frame( path, worm + '_01times.pickle' )
	gonadPosDF = load_data_frame( path, worm + '_02gonadPos.pickle' )
	cellPosDF = load_data_frame( path, worm + '_04cellPos.pickle' )
	cellFluoDF = load_data_frame( path, worm + '_15cellFluo_q&d.pickle' )

	### check if they divided at similar times
	tDiv1 = np.min( timesDF.ix[ pd.notnull( cellPosDF[ '1.ppp' ].X ) ].timesRel )
	tDiv4 = np.min( timesDF.ix[ pd.notnull( cellPosDF[ '4.aaa' ].X ) ].timesRel )
	# if np.abs( tDiv1 - tDiv4 ) >= 1:
	# 	return

	### plot the timeseries
	XYdata = [ [], [] ]
	for idx, key in enumerate( ['1.pp','4.aa'] ):
		print(key)

		### extract all info
		if key == '1.pp':
			cell = cellFluoDF[ key+'p' ].ix[ pd.notnull( cellFluoDF[ key+'p' ]._488nm ) & pd.notnull( cellFluoDF[ key+'a' ]._488nm ) ]
			cell['times'] = timesDF.timesRel
			sis = cellFluoDF[ key+'a' ].ix[ pd.notnull( cellFluoDF[ key+'a' ]._488nm ) & pd.notnull( cellFluoDF[ key+'p' ]._488nm ) ]
			sis['times'] = timesDF.timesRel
		else:
			sis = cellFluoDF[ key+'p' ].ix[ pd.notnull( cellFluoDF[ key+'p' ]._488nm ) & pd.notnull( cellFluoDF[ key+'a' ]._488nm ) ]
			sis['times'] = timesDF.timesRel
			cell = cellFluoDF[ key+'a' ].ix[ pd.notnull( cellFluoDF[ key+'a' ]._488nm ) & pd.notnull( cellFluoDF[ key+'p' ]._488nm ) ]
			cell['times'] = timesDF.timesRel
		# bckg = cellFluoDF[ 'b_' + key[0] ].ix[ pd.notnull( cellFluoDF[ key+'a' ]._488nm ) & pd.notnull( cellFluoDF[ key+'p' ]._488nm ) ]
		# bckg['times'] = timesDF.timesRel

		### find slope
		tDiv = np.min( timesDF.ix[ pd.notnull( cellPosDF[ key + 'p' ].X ) ].timesRel )
		tIdxMin = np.min( timesDF.ix[ pd.notnull( cellPosDF[ key + 'p' ].X ) ].tidxRel )
		tIdxMax = np.max( timesDF.ix[ timesDF.timesRel <= ( tDiv + 1.5 ) ].tidxRel )

		sis.times -= tDiv
		# print(exp)
		# print(times)

		### perform a gaussian filter
		# sigma = 20/60
		# delta_x = sis.times[:,None] - sis.times.values
		# weights = np.exp(-delta_x*delta_x / (2*sigma*sigma)) / (np.sqrt(2*np.pi) * sigma)
		# weights /= np.sum(weights, axis=1, keepdims=True)
		# sis._488nm = np.dot(weights, sis._488nm )

		sis = sis.ix[ ( sis.tidx >= tIdxMin ) & ( sis.tidx <= tIdxMax ) ]

		# compute slope
		slope = np.polyfit(sis.times,sis._488nm,1)[0]
		# slope = np.mean(sis._488nm)

		### find outcome
		tIdxFinal = np.max( cell.tidx )
		exp = ( cell.ix[ cell.tidx == tIdxFinal, '_488nm' ] ).values[0]

		### store data
		XYdata[0].append( slope )
		XYdata[1].append( exp )
		print(slope,exp)

	ax1.plot( -np.diff(XYdata[0]), -np.diff(XYdata[1])/np.sum(XYdata[1]), 'o', color = color, alpha=.6,ms= 10 )
	# if path == 'X:\\Simone\\160930_lag2GFP+mCherry\\':
	# 	txt = 'E1'+worm
	# if path == 'X:\\Simone\\161002_lag2GFP+mCherry\\':
	# 	txt = 'E2'+worm
	ax1.text( -np.diff(XYdata[0]), -np.diff(XYdata[1])/np.sum(XYdata[1]), str(widx), fontsize = 15 )


if __name__ == '__main__':

	### setup figure for the timeseries
	fig1 = plt.figure(figsize=(5.8,5.8))
	ax1 = fig1.add_subplot(111)
	fig1.subplots_adjust(left=0.15, right=.95, top=.95, bottom=0.15)
	for tl in ax1.get_xticklabels():
		tl.set_fontsize(18)
	for tl in ax1.get_yticklabels():
		tl.set_fontsize(18)

	# ax1.set_ylim((-1,1))
	# ax1.set_xlim((-1000,1000))
	ax1.plot([-1000,1000],[0,0],'--k')
	ax1.plot([0,0],[-1,1],'--k')

	ax1.set_ylabel( 'Outcome (Z1-Z4)/(Z1+Z4)', fontsize = 18 )
	ax1.set_xlabel( 'Sister Slope Difference (Z1.ppa-Z4.aap)', fontsize =18 )

	'''
	lag2YFP data
	'''
	# ### plot one worm
	# worms = [ 'W:\\Simone\\160516_lag2_YFP_hist_mCherry\\C21' ]
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
	# paths = [ 'X:\\Simone\\160930_lag2GFP+mCherry', 'X:\\Simone\\161002_lag2GFP+mCherry' ]
	# worms1 = []
	# for path in paths:
	# 	worms1.append( glob.glob( os.path.join( path, '*cellFluo.pickle' ) ) )
	# worms = []
	# for exp in worms1:
	# 	for w in exp:
	# 		worms.append(w)
	# worms.sort()

	if len(worms) > 1:
		colors = np.array( [cm.jet(int(i))[:3] for i in ( np.arange(len(worms))*254./(len(worms)-1) / 254. ) * 254 ] )

	for idx, worm in enumerate( worms ):

		w = worm.split('\\')[-1].split('_')[0]
		path = worm[:-len(worm.split('\\')[-1])]
				
		plotFluorescence( path, w, ax1, channel = '488nm', widx = idx)#, color = colors[idx] )

	plt.show()
