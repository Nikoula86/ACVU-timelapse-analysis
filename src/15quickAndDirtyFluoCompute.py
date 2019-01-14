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

#def computeFluorescence( path, worm, channel = '488nm', lineages = [['1.p','1.pp','1.ppp'],['4.a','4.aa','4.aaa']] ):
	#print( path, worm )
def computeFluorescence( path, worm, channel = '488nm', lineages = [['1.p','1.pp','1.ppp'],['4.a','4.aa','4.aaa']] ):
	print( path, worm )

	if not os.path.isfile( path + worm + '_analyzedImages/combined_withSisters_'+channel+'.tif' ):
		print('cellMovie not found!')
		return

	timesDF = load_data_frame_pandas( path, worm + '_01times.pickle' )
	gonadPosDF = load_data_frame_pandas( path, worm + '_02gonadPos.pickle' )
	cellPosDF = load_data_frame_pandas( path, worm + '_04cellPos.pickle' )
	
	# create new dataframe or load old one
	if os.path.isfile( path + worm + '_15cellFluo_q&d.pickle' ):
		cellFluo = load_data_frame_pandas( path, worm + '_15cellFluo_q&d.pickle' )
	else:
		cellFluo = cellPosDF.copy()
		for key in cellFluo.keys():
			cellFluo[key]['_488nm'] = np.nan
			cellFluo[key]['_561nm'] = np.nan
			cellFluo[key]['times'] = timesDF.timesRel
			cellFluo[key]['X488nm'] = cellPosDF[key].X
			cellFluo[key]['Y488nm'] = cellPosDF[key].Y
			cellFluo[key]['X561nm'] = cellPosDF[key].X
			cellFluo[key]['Y561nm'] = cellPosDF[key].Y
			cellFluo[key] = cellFluo[key].drop('X', 1)
			cellFluo[key] = cellFluo[key].drop('Y', 1)
			cellFluo[key] = cellFluo[key].drop('Z', 1)
		
			
	### compute cells
	for idx, lin in enumerate( lineages ):
		count = 0
		movieCell = load_stack(path + worm + '_analyzedImages/cell_'+channel+'_'+lin[-1].replace('.','')+'.tif' )
		sisName = lin[-1][:-1] + ('a'*(lin[-1][-1]=='p')+'p'*(lin[-1][-1]=='a'))
		movieSis = load_stack(path + worm + '_analyzedImages/cell_'+channel+'_'+sisName.replace('.','')+'.tif' )

		center = int( movieCell[0].shape[0]/2 )

		expCell = np.array( [ np.mean( frame[center-5:center+5,center-5:center+5] ) for frame in movieCell ] )
		expCell = expCell[expCell!=0]

		expSis = np.array( [ np.mean( frame[center-5:center+5,center-5:center+5] ) for frame in movieSis ] )
		expSis = expSis[expSis!=0]

		# print(len(movieCell))
		for jdx, key in enumerate( lin ):
			print('\nComputing ' +channel + ' fluorescence for ' + key )
			cellData = cellPosDF[key].ix[pd.notnull(cellPosDF[key].X)]

			if len(cellData) == 0: continue

			start = count
			count += len(cellData)

			exp = expCell[start:count]
			print('length of the full ' + key[:2] + ' nonNull fluo movie data: ' + str( len(expCell) ) )
			print('length of ' + key + ' pos data: ' + str( len(cellData) ) )
			print('tidx in which ' + key + ' is found: ', list(cellData.tidx) )

			cellData['exp'] = exp
			cellData = cellData.ix[ pd.notnull(cellData.exp) ]

			# save data in the new dataframe
			cellFluo[key].ix[pd.notnull(cellPosDF[key].X), '_'+channel ] = cellData.exp

		# plot sister data
		print('\nComputing ' +channel + ' fluorescence for ' + sisName )
		sisData = cellPosDF[sisName].ix[pd.notnull(cellPosDF[sisName].X)]
		print('length of ' + sisName + ' nonNull fluo movie data: ' + str( len(expSis) ) )
		print('length of ' + sisName + ' pos data: ' + str( len(sisData) ) )
		print('tidx in which ' + sisName + ' is found: ', list(sisData.tidx) )

		sisData['exp'] = expSis
		sisData = sisData.ix[ pd.notnull(sisData.exp) ]

		# save data in the new dataframe
		cellFluo[sisName].ix[pd.notnull(cellPosDF[sisName].X), '_'+channel ] = sisData.exp

	### compute background
	for idx, lin in enumerate( ['b_1','b_4'] ):
		count = 0
		movieCell = load_stack(path + worm + '_analyzedImages/bckg_'+channel+'_'+lin[-1]+'.tif' )
		# print(len(movieCell))

		key = lin
		print('\nComputing ' +channel + ' fluorescence for bckg: ' + key)
		cellData = cellPosDF[key].ix[pd.notnull(cellPosDF[key].X)]

		start = count
		count += len(cellData)

		center = int( movieCell[0].shape[0]/2 )
		exp = np.array( [ np.mean( frame[center-5:center+5,center-5:center+5] ) for frame in movieCell[start:count] ] )
		exp[exp==0] = np.nan
		# print(len(exp),len(cellData))

		cellData['exp'] = exp
		cellData = cellData.ix[ pd.notnull(cellData.exp) ]

		# save data in the new dataframe
		cellFluo[key].ix[pd.notnull(cellPosDF[key].X), '_'+channel ] = cellData.exp


	save_data_frame( cellFluo, path, worm + '_15cellFluo_q&d.pickle'  )



if __name__ == '__main__':
	
	'''
	lag2YFP data
	'''
	# ### compute one worm
	absPath = '/Users/ngritti/Dropbox/ACVU_project_private/ACVU_data/timelapse/pickledFiles/'
	worms = [ absPath + '180524_lag2ts+JVZ32/C07' ]
	### plot all available worms
#	paths = [ 
#			'Y:\\Simone\\161030_lag2YFP+histmCherry',
#			'Y:\\Simone\\161108_lag2YFP+histmCherry', 
#			'Y:\\Simone\\161111_lag2YFP+histmCherry', 
#			'Y:\\Simone\\170312_lag2YFP+histmCherry'
#			]
#	worms1 = []
#	for path in paths:
#		worms1.append( glob.glob( os.path.join( path, '*cellPos.pickle' ) ) )
#	worms = []
#	for exp in worms1:
#		for w in exp:
#			worms.append(w)
#	worms.sort()
#	print(worms)
	
	'''
	lag2GFP data
	'''
	### compute one worm
	# worms = ['X:\\Simone\\161002_lag2GFP+mCherry\\C15']
	### outlayers
	# worms = ['W:\\Simone\\160930_lag2GFP+mCherry\\C09','W:\\Simone\\160930_lag2GFP+mCherry\\C10','W:\\Simone\\160930_lag2GFP+mCherry\\C11','W:\\Simone\\160930_lag2GFP+mCherry\\C15','W:\\Simone\\160930_lag2GFP+mCherry\\C19','W:\\Simone\\161002_lag2GFP+mCherry\\C01']
	### regular worms
	# worms = ['W:\\Simone\\160930_lag2GFP+mCherry\\C01','W:\\Simone\\160930_lag2GFP+mCherry\\C02','W:\\Simone\\160930_lag2GFP+mCherry\\C03','W:\\Simone\\160930_lag2GFP+mCherry\\C06','W:\\Simone\\160930_lag2GFP+mCherry\\C08','W:\\Simone\\161002_lag2GFP+mCherry\\C04']
	### compute all
	# paths = [ 'X:\\Simone\\160930_lag2GFP+mCherry', 'X:\\Simone\\161002_lag2GFP+mCherry', 'X:\\Simone\\161004_lag2GFP+mCherry' ]
	# paths = [ 'X:\\Simone\\161002_lag2GFP+mCherry' ]
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
	#worms = ['X:\\Simone\\160407_HLH2_GFP_hist_mCherry\\C11']
	
	'''
	LIN12nullunc data
	'''
	# ### compute one worm
#	worms = ['Y:\\Simone\\170907_LIN12unc+histmCherry\\C09']
	
	### plot all available worms
	# paths = [ 'T:\\Simone\\170312_lag2YFP+histmCherry' ]#'T:\\Simone\\161030_lag2YFP+mCherry', 'T:\\Simone\\161108_lag2YFP+mCherry' ,'T:\\Simone\\161111_lag2YFP+mCherry' ]
	# worms1 = []
	# for path in paths:
	# 	worms1.append( glob.glob( os.path.join( path, '*cellPos.pickle' ) ) )
	# worms = []
	# for exp in worms1:
	# 	for w in exp:
	# 		worms.append(w)
	# worms.sort()
	# print(worms)
	
	'''
	lag2 integration line data
	'''
	# ### compute one worm
#	worms = [ 'Y:\\Simone\\171028_lag2int+JVZ32\\C08' ]
	### plot all available worms
	# paths = [ 'W:\\Simone\\170312_lag2YFP+histmCherry' ]#'W:\\Simone\\161030_lag2YFP+mCherry', 'W:\\Simone\\161108_lag2YFP+mCherry' ,'W:\\Simone\\161111_lag2YFP+mCherry' ]
	# worms1 = []
	# # for path in paths:
	# 	# worms1.append( glob.glob( os.path.join( path, '*cellPos.pickle' ) ) )
	# worms = []
	# for exp in worms1:
	# 	for w in exp:
	# 		worms.append(w)
	# worms.sort()
	# print(worms)
	
	'''
	COMPUTE THE FLUORESCENCE
	'''
	
	channels = ['488nm','561nm']
	for idx, worm in enumerate( worms ):
		
		w = worm.split('/')[-1].split('_')[0]
		path = worm[:-len(worm.split('/')[-1])]
		
		for channel in channels:
			computeFluorescence( path, w, channel = channel )
