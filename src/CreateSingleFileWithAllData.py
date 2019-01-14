# -*- coding: utf-8 -*-
"""
Created on Fri Jan  5 14:22:09 2018

@author: gritti
"""

import numpy as np
import matplotlib.pyplot as plt
#import pickle
from matplotlib import cm
import os
import generalFunctions as gf
import pandas as pd
import glob
#from skimage import filters


import matplotlib as mpl
mpl.rcParams['pdf.fonttype'] = 42

def plotRatio( path, worm, sigma = 20/60 ):
	wormID = os.path.join( path.split('\\')[-2], worm )
	print( path, worm )

	timesDF = gf.load_data_frame( path, worm + '_01times.pickle' )
#	gonadPosDF = gf.load_data_frame( path, worm + '_02gonadPos.pickle' )
	cellPosDF = gf.load_data_frame( path, worm + '_04cellPos.pickle' )
	cellFluoDF = gf.load_data_frame( path, worm + '_15cellFluo_q&d.pickle' )
	
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
	tpL2 = timesDF.ix[ timesDF.tidxRel == ( lethtidx[1] - mintp ), 'timesRel' ].values[0]

	# relative to L2 start
	tpRel = tpL2

	# extract data from dataframe
	cellData = pd.DataFrame({})
	index = 0
	for idx, tRow in timesDF.iterrows():
		# print(tRow.tidxRel)
		

		#extract current cells
		currentCells = gf.extract_current_cell_pos( cellPosDF, tRow.tidxRel )

		if len(currentCells) > 0:
			(c1,c4,b1,b4) = gf.find_interesting_cells(currentCells)

			# print(c1.cname,c4.cname)
			val1 = cellFluoDF[c1.cname].ix[ cellFluoDF[c1.cname].times == tRow.timesRel, '_488nm' ].values[0]/1000.
			valb1 = cellFluoDF[b1.cname].ix[ cellFluoDF[b1.cname].times == tRow.timesRel, '_488nm' ].values[0]/1000.
			val4 = cellFluoDF[c4.cname].ix[ cellFluoDF[c4.cname].times == tRow.timesRel, '_488nm' ].values[0]/1000.
			valb4 = cellFluoDF[b4.cname].ix[ cellFluoDF[b4.cname].times == tRow.timesRel, '_488nm' ].values[0]/1000.
			val1corr = val1 - valb1
			val4corr = val4 - valb4
			newRow = pd.DataFrame( {'cell1':val1,
			   'b1':valb1,
			   'cell1corr':val1corr,
				'cell4':val4,
			   'b4':valb4,
			   'cell4corr':val4corr,
				'times':tRow.timesRel}, index= [index] )
			cellData = pd.concat( [cellData,newRow] )
			index += 1

	# print(tdiv1,tdiv2)
	mothersData = cellData.ix[ cellData.times < np.min([tdiv1,tdiv2]) ]
	intermediateData = cellData.ix[ (cellData.times>=np.min([tdiv1,tdiv2])) & (cellData.times<np.max([tdiv1,tdiv2])) ]
	daughtersData = cellData.ix[ cellData.times>=np.max([tdiv1,tdiv2]) ]
		
	absZ1Raw = np.array( [ mothersData.cell1.values, intermediateData.cell1.values, daughtersData.cell1.values ] )
	absZ4Raw = np.array( [ mothersData.cell4.values, intermediateData.cell4.values, daughtersData.cell4.values ] )
	absb1Raw = np.array( [ mothersData.b1.values, intermediateData.b1.values, daughtersData.b1.values ] )
	absb4Raw = np.array( [ mothersData.b4.values, intermediateData.b4.values, daughtersData.b4.values ] )

	ratioRaw = np.array( [ ( ( mothersData.cell1 - mothersData.cell4 ) / ( mothersData.cell1 + mothersData.cell4 ) ).values,
			 ( ( intermediateData.cell1 - intermediateData.cell4 ) / ( intermediateData.cell1 + intermediateData.cell4 ) ).values,
			 ( ( daughtersData.cell1 - daughtersData.cell4 ) / ( daughtersData.cell1 + daughtersData.cell4 ) ).values, ] )

	ratioCorr = np.array( [ ( ( mothersData.cell1corr - mothersData.cell4corr ) / ( mothersData.cell1corr + mothersData.cell4corr ) ).values,
			 ( ( intermediateData.cell1corr - intermediateData.cell4corr ) / ( intermediateData.cell1corr + intermediateData.cell4corr ) ).values,
			 ( ( daughtersData.cell1corr - daughtersData.cell4corr ) / ( daughtersData.cell1corr + daughtersData.cell4corr ) ).values, ] )
	
	# perform gaussian filter on raw  and bckg-corrected data
	for val in zip(['cell1','cell4'],[tdiv1,tdiv2]):
		key = val[0]
		tdiv = val[1]

		### perform a gaussian filter on mother
		x_eval = cellData.ix[cellData.times < tdiv, 'times'].values

		delta_x = x_eval[:,None] - cellData.ix[cellData.times < tdiv, 'times'].values
		weights = np.exp(-delta_x*delta_x / (2*sigma*sigma)) / (np.sqrt(2*np.pi) * sigma)
		weights /= np.sum(weights, axis=1, keepdims=True)

		val = np.dot(weights, cellData.ix[cellData.times < tdiv,key] )
		cellData.ix[cellData.times < tdiv,key+'Filt'] = val

		val = np.dot(weights, cellData.ix[cellData.times < tdiv,key+'corr'] )
		cellData.ix[cellData.times < tdiv,key+'corrFilt'] = val

		### perform a gaussian filter on mother
		x_eval = cellData.ix[cellData.times >= tdiv, 'times'].values

		delta_x = x_eval[:,None] - cellData.ix[cellData.times >= tdiv, 'times'].values
		weights = np.exp(-delta_x*delta_x / (2*sigma*sigma)) / (np.sqrt(2*np.pi) * sigma)
		weights /= np.sum(weights, axis=1, keepdims=True)

		val = np.dot(weights, cellData.ix[cellData.times >= tdiv,key] )
		cellData.ix[cellData.times >= tdiv,key+'Filt'] = val

		val = np.dot(weights, cellData.ix[cellData.times >= tdiv,key+'corr'] )
		cellData.ix[cellData.times >= tdiv,key+'corrFilt'] = val

	mothersData = cellData.ix[ cellData.times < np.min([tdiv1,tdiv2]) ]
	intermediateData = cellData.ix[ (cellData.times>=np.min([tdiv1,tdiv2])) & (cellData.times<np.max([tdiv1,tdiv2])) ]
	daughtersData = cellData.ix[ cellData.times>=np.max([tdiv1,tdiv2]) ]

	ratioFiltRaw = np.array( [ ( ( mothersData.cell1Filt - mothersData.cell4Filt ) / ( mothersData.cell1Filt + mothersData.cell4Filt ) ).values,
			 ( ( intermediateData.cell1Filt - intermediateData.cell4Filt ) / ( intermediateData.cell1Filt + intermediateData.cell4Filt ) ).values,
			 ( ( daughtersData.cell1Filt - daughtersData.cell4Filt ) / ( daughtersData.cell1Filt + daughtersData.cell4Filt ) ).values, ] )

	ratioFiltCorr = np.array( [ ( ( mothersData.cell1corrFilt - mothersData.cell4corrFilt ) / ( mothersData.cell1corrFilt + mothersData.cell4corrFilt ) ).values,
			 ( ( intermediateData.cell1corrFilt - intermediateData.cell4corrFilt ) / ( intermediateData.cell1corrFilt + intermediateData.cell4corrFilt ) ).values,
			 ( ( daughtersData.cell1corrFilt - daughtersData.cell4corrFilt ) / ( daughtersData.cell1corrFilt + daughtersData.cell4corrFilt ) ).values, ] )
	
	times = np.array( [ mothersData.times.values, intermediateData.times.values, daughtersData.times.values ] ) - tpRel

	# convert the ratio from (1-4)/(1+4) to (1st-2nd)/(1st+2nd)
	if tdiv1>=tdiv2:
		ratioRaw*=-1
		ratioFiltRaw*=-1
		ratioCorr*=-1
		ratioFiltCorr*=-1
	
	# save the ACid
	ACid = 4
	if ratioFiltCorr[2][-1]	> 0:
		ACid = 1
	
	# explanation:
	exp = ['folder', 'ACid', 'timesToL2', 'absValZ1', 'absValb1', 'absValZ4', 'absValb4', 'ratioRaw', 'ratioFiltRaw', 'ratioCorr', 'ratioFiltCorr']
	return [ exp, wormID, ACid, times, absZ1Raw, absb1Raw, absZ4Raw, absb4Raw, ratioRaw, ratioFiltRaw, ratioCorr, ratioFiltCorr ]



if __name__ == '__main__':


#%%
	cwd = os.getcwd()
	parentDir = os.path.join(os.path.dirname(os.path.dirname(cwd)),'ACVU_data','timelapse')

	'''
	lag2YFP data
	'''
	### plot one worm
	worms = [ os.path.join( parentDir, '161111_lag2YFP+histmCherry', 'C13' ) ]

	### plot all available worms
	dataFolders = [ '161030_lag2YFP+histmCherry', '161108_lag2YFP+histmCherry', '161111_lag2YFP+histmCherry', '170312_lag2YFP+histmCherry' ]
	paths = [os.path.join(parentDir,dataFolder) for dataFolder in dataFolders]
	worms1 = []
	for path in paths:
		worms1.append( glob.glob( os.path.join( path, '*cellFluo_q&d.pickle' ) ) )
	worms = []
	for exp in worms1:
		for w in exp:
			worms.append(w)
	worms.sort()

	'''
	lag2multi data
#	'''
#	### plot one worm
#	worms = [ os.path.join( parentDir, '180106_lag2multi+JVZ32', 'C66' ) ]
#
#	### plot all available worms
#	dataFolders = [ '180106_lag2multi+JVZ32' ]
#	paths = [os.path.join(parentDir,dataFolder) for dataFolder in dataFolders]
#	worms1 = []
#	for path in paths:
#		worms1.append( glob.glob( os.path.join( path, '*cellFluo_q&d.pickle' ) ) )
#	worms = []
#	for exp in worms1:
#		for w in exp:
#			worms.append(w)
#	worms.sort()
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
#	dataFolders = [ '160930_lag2GFP+mCherry', '161002_lag2GFP+mCherry' ]
#	paths = [os.path.join(parentDir,dataFolder) for dataFolder in dataFolders]
#	worms1 = []
#	for path in paths:
#		worms1.append( glob.glob( os.path.join( path, '*cellFluo_q&d.pickle' ) ) )
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
#%%
	'''
	PLOT THE RATIO
	'''

	allData = []
	for idx, wID in enumerate( worms ):

		worm = wID.split('\\')[-1].split('_')[0]
		path = wID[:-len(wID.split('\\')[-1])]
				
		data = plotRatio( path, worm, sigma = 15/60 )
		allData.append(data)
		
#	pickle.dump(allData,open(os.path.join(parentDir,'allData_lag2multi.pickle'),'wb'))
	