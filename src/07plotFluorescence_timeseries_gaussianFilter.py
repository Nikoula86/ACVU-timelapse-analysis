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

def plotFluorescence( path, worm, ax1, ax2, color = 'k', channel = '488nm', lineages = [['1.p','1.pp','1.ppp'],['4.a','4.aa','4.aaa']], filt = False, sigma = 20/60 ):
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

			if channel == '488nm':
				cell = cellFluoDF[ key ].ix[ pd.notnull( cellFluoDF[ key ]._488nm ) ]
				bckg = cellFluoDF[ 'b_'+key[0] ].ix[ pd.notnull( cellFluoDF[ key ]._488nm ) ]
				times = timesDF.ix[ pd.notnull( cellFluoDF[ key ]._488nm ) ]

				if filt:
					### perform a gaussian filter
					cell._488nmFilt = np.nan
					tidx = times.tidxRel
					x_eval = times.timesRel

					delta_x = x_eval[:,None] - times.timesRel.values
					weights = np.exp(-delta_x*delta_x / (2*sigma*sigma)) / (np.sqrt(2*np.pi) * sigma)
					weights /= np.sum(weights, axis=1, keepdims=True)
					val = np.dot(weights, cell._488nm - bckg._488nm )

					cell._488nmFilt = val
				else:
					cell._488nmFilt = cell._488nm

				# with background correction
				if (key == '1.ppa') or (key == '4.aap'):
					ax1.plot( times.timesRel-tpRel, 
								( cell._488nmFilt ),# / bckg._488nm, 
								'--', dashes = [3,2],
								color = colors[idx,jdx], lw=2 )
				else:
					ax1.plot( times.timesRel-tpRel, 
							( cell._488nmFilt ),# / bckg._488nm, 
							'-', color = colors[idx,jdx], lw=2 )


			elif channel == '561nm':
				cell = cellFluoDF[ key ].ix[ pd.notnull( cellFluoDF[ key ]._561nm ) ]
				if len(cell.tidx>0):
					# print(cell)
					bckg = cellFluoDF[ 'b_'+key[0] ].ix[ pd.notnull( cellFluoDF[ key ]._561nm ) ]
					times = timesDF.ix[ pd.notnull( cellFluoDF[ key ]._561nm ) ]

					### calculate signal with otsu
					area = []
					for tidx in cell.tidx:
						currentCells = extract_current_cell_fluo( cellPosDF, cellOutDF, cellFluoDF, tidx )
						currentCell = currentCells.ix[ currentCells.cname == key ].squeeze()
						gonadPos = extract_pos( gonadPosDF.ix[ gonadPosDF.tidx == tidx ].squeeze() )

						imgpxl = currentCell.imgPxl
						medianCorrection = np.median( ( flatField - darkField ) )
						darkF = crop_image( darkField, gonadPos, 512 )
						flatF = crop_image( flatField, gonadPos, 512 )
						den = np.clip( flatF.astype(np.float) - darkF.astype(np.float), 0, None )

						cellPos = extract_3Dpos( currentCell )
						cellOut = extract_out( currentCell ) * imgpxl / 1000.
						drift = np.array([currentCell.ix['X'+channel]-currentCell.X,currentCell.ix['Y'+channel]-currentCell.Y])
				
						#decorrect image	
						fileName = os.path.join( path, w + '_analyzedImages', timesDF.ix[timesDF.tidxRel==tidx,'fName'].values[0] + channel + '.tif')
						stack = load_stack( fileName )
						img = stack[cellPos[2]]

						imgCorr = np.clip( img.astype(np.float) - darkF.astype(np.float), 0, None ) * medianCorrection / den
						imgCorrSmall = imgCorr[cellPos[1]-imgpxl/2:cellPos[1]+imgpxl/2+1,cellPos[0]-imgpxl/2:cellPos[0]+imgpxl/2+1]

						thr = filters.threshold_otsu(imgCorrSmall)
						global_otsu = imgCorrSmall >= thr
						# plt.imshow(global_otsu,interpolation='nearest')
						# plt.show()
						# plt.imshow(imgCorrSmall,interpolation='nearest')
						# plt.show()
						signal = np.sum(imgCorrSmall*global_otsu) / np.sum(global_otsu)
						cell.ix[cell.tidx==tidx,'_561nm'] = signal
						area = np.sum(global_otsu)

						# ### find the new signal
						# if is_outline_cell(currentCell):
						# 	area.append( calculate_area( imgCorrSmall, drift, cellOut ) )
						# else:
						# 	area.append(36)
						print(fileName, signal, area)
					area = np.array(area)



					# without background correction
					ax1.plot( times.timesRel-tpRel, 
								( cell._561nm - bckg._561nm ),
								'-', color = 'b')#colors[idx,jdx], lw=2 )
					ax1.plot( times.timesRel-tpRel, 
								( bckg._561nm ),
								'-', color = 'b')# colors[idx,jdx], lw=2 )
					# ax1.plot( times.timesRel-tpRel, 
					# 			( cell._561nm ),
					# 			'-', color = colors[idx,jdx], lw=2 )

	### plot the ratio
	data = pd.DataFrame( { 'tidx': timesDF.tidxRel,
							'times': timesDF.timesRel,
							'ratio': np.nan,
							'lineage': np.nan } )
	# print(data)
	colorsRatio = ['black','blue','magenta']

	### perform a gaussian filter
	if filt:
		for key in cellFluoDF.keys():
			if key[0]!='b':
				cell = cellFluoDF[ key ].ix[ pd.notnull( cellFluoDF[ key ]._488nm ) ]
				bckg = cellFluoDF[ 'b_'+key[0] ].ix[ pd.notnull( cellFluoDF[ key ]._488nm ) ]
				times = timesDF.ix[ pd.notnull( cellFluoDF[ key ]._488nm ) ]

				x_eval = times.timesRel

				delta_x = x_eval[:,None] - times.timesRel.values
				weights = np.exp(-delta_x*delta_x / (2*sigma*sigma)) / (np.sqrt(2*np.pi) * sigma)
				weights /= np.sum(weights, axis=1, keepdims=True)
				cell._488nm = np.dot(weights, cell._488nm - bckg._488nm )
				for idx, row in cell.iterrows():
					cellFluoDF[key].ix[ cellFluoDF[key].tidx == row.tidx, '_488nm' ] = row._488nm

	for idx, trow in timesDF.iterrows():

		currentCells = extract_current_cell_fluo( cellPosDF, cellOutDF, cellFluoDF, trow.tidxRel )
		# print('\n',trow,currentCells)

		if len(currentCells) > 2:

			(cell1, cell2, b1, b2) = find_interesting_cells( currentCells )

			if channel == '488nm':
				# with bckg correction
				c1signal = ( cell1._488nm - b1._488nm )# / b1._488nm
				c2signal = ( cell2._488nm - b2._488nm )# / b2._488nm

				# # without bckg correction
				# c1signal = cell1.fluo488nm
				# c2signal = cell2.fluo488nm
			elif channel == '561nm':
				# with bckg correction
				c1signal = ( cell1._561nm - b1._561nm )# / b1._561nm
				c2signal = ( cell2._561nm - b2._561nm )# / b2._561nm

				# # without bckg correction
				# c1signal = cell1.fluo561nm
				# c2signal = cell2.fluo561nm

			# print(cell1.cname,cell1.fluo488nm,cell2.cname,cell2.fluo488nm)
			data.ix[ data.tidx == trow.tidxRel, 'ratio' ] = ( c1signal - c2signal ) / ( c1signal + c2signal )
			data.ix[ data.tidx == trow.tidxRel, 'lineage' ] = int(np.min([len(cell1.cname),len(cell2.cname)])-3)
			# print(cell1.cname,cell2.cname,int(np.min([len(cell1.cname),len(cell2.cname)])-3),colorsRatio[int(np.min([len(cell1.cname),len(cell2.cname)])-3)])

	data = data.ix[ pd.notnull( data.ratio ) ].reset_index()
	# print(data)

	for idx, d in data.iterrows():
		if d.lineage >= 0:
			# print(idx,idx+1)
			ax2.plot( [ d.times-tpRel, data.times.values[np.clip(idx+1,0,len(data)-1)]-tpRel ], [ d.ratio, data.ratio.values[np.clip(idx+1,0,len(data)-1)] ], 'o', color = colorsRatio[int(d.lineage)], lw=2 )

	# ### plot ecdysis
	# for idx, tidx in enumerate( lethtidx ):
	# 	if tidx >= 0:
	# 		tp = timesDF.ix[ timesDF.tidxRel == tidx, 'timesRel' ].values[0]
	# 		ax1.plot([tp-tpRel,tp-tpRel],[0,2000],'--', color = color, lw=1)
	# 		ax2.plot([tp-tpRel,tp-tpRel],[-1,1],'--', color = color, lw=1)
			# ax2.plot([tp,tp],[-1,1],'--', color = color, lw=1)

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

	### setup figure for the ratio
	fig2 = plt.figure(figsize=(7,2.75))
	ax2 = fig2.add_subplot(111)
	fig2.subplots_adjust(left=0.15, right=.95, top=.95, bottom=0.15)
	for tl in ax2.get_xticklabels():
		tl.set_fontsize(18)
	for tl in ax2.get_yticklabels():
		tl.set_fontsize(18)
	ax2.set_ylim(-1,1)
	
	# ax2.set_xlim((9,25))
	# ax2.set_xlim(-5,15)

	ax1.plot( [0,0], [0,4000], '--', color = 'black', lw=1)
	ax2.plot( [0,0], [-6,6], '--', color = 'black', lw=1)

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
#	worms = ['X:\\Simone\\161002_lag2GFP+mCherry\\C15']
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
	lag2 integration line data
	'''
	# ### plot one worm
	worms = [ 'Y:\\Simone\\171028_lag2int+JVZ32\\C01' ]
	# ### plot all available worms
	# paths = [ 'W:\\Simone\\160516_lag2_YFP_hist_mCherry', 'W:\\Simone\\160226_lag2_YFP_hist_mCherry' ]
	# worms1 = [ glob.glob( os.path.join( paths[0], '*cellFluo.pickle' ) ), glob.glob( os.path.join( paths[1], '*cellFluo.pickle' ) ) ]
	# worms = []
	# for exp in worms1:
	# 	for w in exp:
	# 		worms.append(w)
	# worms.sort()

	for idx, worm in enumerate( worms ):

		w = worm.split('\\')[-1].split('_')[0]
		path = worm[:-len(worm.split('\\')[-1])]
				
		plotFluorescence( path, w, ax1, ax2, channel = '488nm', filt = True, sigma = 15/60, lineages = [['1.pp','1.ppp','1.ppa'],['4.aa','4.aaa','4.aap']] )

	ax2.plot([-6,8],[0,0],'--k', lw=1)
	ax2.set_ylim(-1,1)
	ax2.set_xlim(-5,10)
		
	plt.show()
