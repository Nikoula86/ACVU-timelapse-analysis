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

def calcDivTime( path, worm ):

	timesDF = load_data_frame( path, worm + '_01times.pickle' )
	cellPosDF = load_data_frame( path, worm + '_04cellPos.pickle' )
	
	### find ecdysis timepoint

	ecd = np.loadtxt( open( os.path.join( path, 'skin.txt'), 'rb' ) )
	# load ecdysis data
	index = np.where( ecd[:,0] == float(worm[1:]) )
	mintp = np.min( [ i for i in ecd[index, 1:6][0][0] if i >= 0 ] )
	lethtidx = ecd[ index, 2:6 ][0][0] - 1
	tpL2 = timesDF.ix[ timesDF.tidxRel == lethtidx[0], 'timesRel' ].values[0]
	tpL1 = timesDF.ix[ timesDF.tidxRel == mintp, 'timesRel' ].values[0]

	### find division times
	tdiv1 = np.min( timesDF.ix[ pd.notnull( cellPosDF[ '1.ppp' ].X ) ].timesRel ) - tpL2
	tdiv2 = np.min( timesDF.ix[ pd.notnull( cellPosDF[ '4.aaa' ].X ) ].timesRel ) - tpL2

	### plot data
	return (tdiv1,tdiv2)




if __name__ == '__main__':

	### setup figure for the timeseries
	fig1 = plt.figure(figsize=(6.8,6.8))
	ax1 = fig1.add_subplot(111)
	fig1.subplots_adjust(left=0.15, right=.95, top=.95, bottom=0.15)
	for tl in ax1.get_xticklabels():
		tl.set_fontsize(18)
	for tl in ax1.get_yticklabels():
		tl.set_fontsize(18)

	ax1.plot( [-10,10], [0,0], '--', color = 'black', lw=1)
	ax1.plot( [0,0], [-10,10], '--', color = 'black', lw=1)
	ax1.plot( [-10,10], [-10,10], '--', color = 'black', lw=1)

	ax1.set_xlim(-2,5)
	ax1.set_ylim(-2,5)

	'''
	lag2YFP data
	'''
	### plot one worm
	# worms = [ 'W:\\Simone\\161030_lag2YFP+mCherry\\C23' ]
	### plot all available worms
	#paths = [ 'Y:\\Simone\\161030_lag2YFP+histmCherry', 'Y:\\Simone\\161108_lag2YFP+histmCherry', 'Y:\\Simone\\161111_lag2YFP+histmCherry', 'Y:\\Simone\\170312_lag2YFP+histmCherry' ]
	
	#paths = [ 'Y:\\Simone\\170405_LIN12BAL+histmCherry', 'Y:\\Simone\\170407_LIN12BAL+histmCherry', 'Y:\\Simone\\170419_LIN12BAL+histmCherry', 'Y:\\Simone\\170421_LIN12BAL+histmCherry' ]
	
	#worms1 = []
	#for path in paths:
	#	worms1.append( glob.glob( os.path.join( path, '*cellPos.pickle' ) ) )
	#worms = []
	#for exp in worms1:
	#	for w in exp:
	#		worms.append(w)
	#worms.sort()

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
	# 	worms1.append( glob.glob( os.path.join( path, '*cellPos.pickle' ) ) )
	# worms = []
	# for exp in worms1:
	# 	for w in exp:
	# 		worms.append(w)
	# worms.sort()

	'''
	lin12nullunc data
	'''
	### plot one worm
	worms = [ 'T:\\Simone\\170714_LIN12unc+histmCherry\\C24' ]
	### plot all available worms
	#paths = [ 'Y:\\Simone\\161030_lag2YFP+histmCherry', 'Y:\\Simone\\161108_lag2YFP+histmCherry', 'Y:\\Simone\\161111_lag2YFP+histmCherry', 'Y:\\Simone\\170312_lag2YFP+histmCherry' ]
	
	#paths = [ 'Y:\\Simone\\170405_LIN12BAL+histmCherry', 'Y:\\Simone\\170407_LIN12BAL+histmCherry', 'Y:\\Simone\\170419_LIN12BAL+histmCherry', 'Y:\\Simone\\170421_LIN12BAL+histmCherry' ]
	
	#worms1 = []
	# for path in paths:
		# worms1.append( glob.glob( os.path.join( path, '*cellPos.pickle' ) ) )
	# worms = []
	# for exp in worms1:
		# for w in exp:
			# worms.append(w)
	# worms.sort()

	'''
	START
	'''
	tDiv = []
	for idx, worm in enumerate( worms ):

		w = worm.split('\\')[-1].split('_')[0]
		path = worm[:-len(worm.split('\\')[-1])]
		print( idx, path, w )
				
		tdiv = calcDivTime( path, w )
		ax1.text(tdiv[0],tdiv[1], str(idx) )
		tDiv.append(tdiv)

	tDiv = np.array( tDiv )
	print(tDiv)
	ax1.plot( tDiv[:,0], tDiv[:,1], 'ok', alpha = .5 )

	ax1.set_xlim(np.floor(np.min(tDiv)),np.floor(np.max(tDiv)+1))
	ax1.set_ylim(np.floor(np.min(tDiv)),np.floor(np.max(tDiv)+1))

	ax1.set_xlim(-2,4)
	ax1.set_ylim(-2,4)

	ax1.set_xlabel('Z1.pp division',fontsize = 18)
	ax1.set_ylabel('Z4.aa division',fontsize = 18)
	if 'lag2GFP' in worms[0]:
		ax1.set_title('lag2GFP',fontsize = 18)
	if 'lag2YFP' in worms[0]:
		ax1.set_title('lag2YFP',fontsize = 18)

	#print('CORRELATION COEFFICIENT:', np.corrcoeff(tDiv[:,0], tDiv[:,1]))
	plt.show()
