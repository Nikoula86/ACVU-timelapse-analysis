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

def plotPos( path, worm, ax1, filt = False, sigma = 20/60 ):
	print( path, worm )

	paramsDF = load_data_frame( path, worm + '_01params.pickle' )
	timesDF = load_data_frame( path, worm + '_01times.pickle' )
	gonadPosDF = load_data_frame( path, worm + '_02gonadPos.pickle' )
	cellPosDF = load_data_frame( path, worm + '_04cellPos.pickle' )
	cellFluoDF = load_data_frame( path, worm + '_15cellFluo_q&d.pickle' )
	# print(paramsDF)
	# return
	cellFluoDF.pop('b_1',None)
	cellFluoDF.pop('b_4',None)
	cellFluoDF.pop('4.a',None)
	cellFluoDF.pop('1.p',None)

	### extract division times
	tfirst = np.min( timesDF.ix[ pd.notnull( cellPosDF[ '1.pp' ].X ) ].timesRel )
	tdiv1 = np.min( timesDF.ix[ pd.notnull( cellPosDF[ '1.ppp' ].X ) ].timesRel )
	tdiv4 = np.min( timesDF.ix[ pd.notnull( cellPosDF[ '4.aaa' ].X ) ].timesRel )
	tlast = np.max( timesDF.ix[ pd.notnull( cellPosDF[ '4.aaa' ].X ) ].timesRel )

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

	# extract data from dataframe
	timesDF = timesDF.ix[ pd.notnull( cellPosDF[ '1.pp' ].X ) | pd.notnull( cellPosDF[ '1.ppp' ].X ) ]
	vmin = 10000
	vmax = 0
	for idx, tRow in timesDF.iterrows():
		# print(tRow.timesRel)

		#extract current cells
		currentCells = extract_current_cell_fluo_qd( cellPosDF, cellFluoDF, tRow.tidxRel )
		currentCells = currentCells.ix[ (currentCells.cname!='b_1')&(currentCells.cname!='b_4') ]
		cell1 = '1.ppp'
		if cell1 not in list(currentCells.cname):
			cell1 = '1.pp'
		refPos = np.array( [ currentCells.ix[currentCells.cname==cell1,'X'].values[0] * paramsDF.pxlSize, 
					currentCells.ix[currentCells.cname==cell1,'Y'].values[0] * paramsDF.pxlSize ])#, 
#					currentCells.ix[currentCells.cname==cell1,'Z'].values[0] ] )
		cell4 = '4.aaa'
		if cell4 not in list(currentCells.cname):
			cell4 = '4.aa'
		refPos += np.array( [ currentCells.ix[currentCells.cname==cell4,'X'].values[0] * paramsDF.pxlSize, 
					currentCells.ix[currentCells.cname==cell4,'Y'].values[0] * paramsDF.pxlSize])# , 
#					currentCells.ix[currentCells.cname==cell4,'Z'].values[0] ] )
		refPos /= 2.

		# print(currentCells)

		for cname in currentCells.cname:

			# print(tRow.timesRel, cname, leftcell, rightcell)
			posCell = np.array( [ currentCells.ix[currentCells.cname==cname,'X'].values[0] * paramsDF.pxlSize, 
						currentCells.ix[currentCells.cname==cname,'Y'].values[0] * paramsDF.pxlSize])#, 
#						currentCells.ix[currentCells.cname==cname,'Z'].values[0] ] )

			cellFluoDF[ cname ].ix[ cellFluoDF[cname].times == tRow.timesRel, 'dist' ] = ( -1 * (cname[0]=='1') + 1 * (cname[0]=='4') ) * np.sqrt( np.sum( ( posCell - refPos )**2 ) )
			vmin = np.min( [ vmin, np.min( cellFluoDF[cname]._488nm ) ] )
			vmax = np.max( [ vmax, np.max( cellFluoDF[cname]._488nm ) ] )

	# print(data)
	for cname in list( cellFluoDF.keys() ):
		celldf = cellFluoDF[ cname ].ix[ pd.notnull( cellFluoDF[cname]._488nm ) ]
		ax1.plot( celldf.dist, celldf.times - tpRel, '-k', lw=.5 )
		im = ax1.scatter( celldf.dist, celldf.times - tpRel, c=celldf._488nm,s=100,vmin = vmin,vmax = vmax/1.2, lw = 0.5 )

	# plt.colorbar(im)









if __name__ == '__main__':

	### setup figure for the positions Plot
	fig1 = plt.figure(figsize=(5.8,3.8))
	ax1 = fig1.add_subplot(111)
	fig1.subplots_adjust(left=0.15, right=.95, top=.95, bottom=0.15)
	for tl in ax1.get_xticklabels():
		tl.set_fontsize(18)
	for tl in ax1.get_yticklabels():
		tl.set_fontsize(18)

	ax1.set_xlim(-15,15)
	ax1.plot( [-15,15], [0,0], '--', color = 'black', lw=1)

	'''
	lag2YFP data
	'''
	### plot one worm
	worms = [ 'W:\\Simone\\161111_lag2YFP+mCherry\\C17' ]

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

	figs = []
	for idx, worm in enumerate( worms ):
		plt.close(fig1)
	
		figs.append(plt.figure(figsize=(5.8,10.8)))
		ax1 = figs[-1].add_subplot(111)
		ax1.set_xlabel('Z1.ppa<-Distance [um]->Z4.aap',fontsize = 18)
		ax1.set_ylabel('Time after L2 [hours]',fontsize = 18)
		# ax1.plot( [0,5], [0,5], '--', color = 'black', lw=1)
		# ax1.set_xlim(0,4)
		# ax1.set_ylim(0,4)
		ax1.set_xlim(-15,15)
		ax1.set_ylim(-3,10)
		ax1.plot( [-15,15], [0,0], '--', color = 'black', lw=1)
		ax1.set_title( 'worm number: ' + str(idx) )
		figs[-1].subplots_adjust(left=0.15, right=.95, top=.95, bottom=0.15)
		for tl in ax1.get_xticklabels():
			tl.set_fontsize(18)
		for tl in ax1.get_yticklabels():
			tl.set_fontsize(18)

		w = worm.split('\\')[-1].split('_')[0]
		path = worm[:-len(worm.split('\\')[-1])]
		ax1.invert_yaxis()
				
		plotPos( path, w, ax1, filt = True, sigma = 15/60 )
		ax1.set_title( 'worm number: ' + str(idx) )
		plt.savefig('W:\\Nicola\\timelapse_data\\posPlots\\'+str(idx)+'.pdf')
		plt.show()

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

	# for idx, worm in enumerate( worms ):

	# 	w = worm.split('\\')[-1].split('_')[0]
	# 	path = worm[:-len(worm.split('\\')[-1])]
				
	# 	plotPos( path, w, ax1, filt = True, sigma = 15/60 )

	ax1.invert_yaxis()
	plt.show()
