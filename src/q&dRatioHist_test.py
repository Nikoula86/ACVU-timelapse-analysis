# -*- coding: utf-8 -*-
"""
Created on Thu Jan  4 13:48:24 2018

@author: gritti
"""

# -*- coding: utf-8 -*-
"""
Created on Wed Jan  3 11:45:30 2018

@author: gritti
"""

import numpy as np
import matplotlib.pyplot as plt
import pickle
#from matplotlib import cm
import os
import generalFunctions as gf
import pandas as pd
import glob
#from skimage import filters


import matplotlib as mpl
mpl.rcParams['pdf.fonttype'] = 42

def plotRatio( data, path, worm ):
	print( path, worm )
	
	timesDF = gf.load_data_frame( path, worm + '_01times.pickle' )
	cellPosDF = gf.load_data_frame( path, worm + '_04cellPos.pickle' )
	
	### extract division times
	tdiv1 = np.min( timesDF.ix[ pd.notnull( cellPosDF[ '1.ppp' ].X ) ].timesRel )
	tdiv4 = np.min( timesDF.ix[ pd.notnull( cellPosDF[ '4.aaa' ].X ) ].timesRel )
	
	times = data[0]
	# set tpRel as 1st dividing cell
	if len(times[1])>0:
		times = times - times[1][0]

	ratio = data[1]
	# adjust so to compute (1st-2nd)/(1st+2nd) instead of (Z1-Z4)/(Z1+Z4)
	if tdiv1>=tdiv4:
		ratio*=-1
	# make lines connected
	if len( times[1] ) > 0:
		times[0] = np.append(times[0],times[1][0])
		times[1] = np.append(times[1],times[2][0])
		ratio[0] = np.append(ratio[0],ratio[1][0])
		ratio[1] = np.append(ratio[1],ratio[2][0])
	else:
		times[0] = np.append(times[0],times[2][0])
		ratio[0] = np.append(ratio[0],ratio[2][0])
		
	### plot data
	startOtherDir = ( (ratio[2][0]*ratio[2][-1]) < 0 )
	dtOtherDir = np.sum( np.diff( times[2][(ratio[2]*ratio[2][-1])<0] ) )
	print(times[2][ratio[2]>0])
	
	colors = ['k','b','m']
	
#	if ( startOtherDir ) & ( dtOtherDir > 2.5 ):# & ( np.abs(tdiv1-tdiv4) < 0.6 ):
		
#	for idx, val in enumerate( zip(times,ratio) ):
##		if (tdiv1!=tdiv4) and (np.sign(ratio[2][-1])==-1):
#		# plot mothers
#		ax1.plot( val[0], val[1],
#					'-', lw = .5, color = colors[idx], alpha = alphaFilt )
#		if (tdiv1<tdiv4) and (np.sign(ratio[2][-1])==1):
#			# plot mothers
#			ax1.plot( val[0], val[1],
#						'-', lw = .5, color = colors[idx], alpha = alphaFilt )

#	print(tdiv1,tdiv4,ratio[2][-1])
#	print((tdiv1>=tdiv4)==(ratio[2][-1]>0))
	## plot ecdysis
#	if plotEcd:
#		ax1.plot( [tpL2 - tpRel, tpL2 - tpRel], [-1,1], '--', color = 'black', lw=.25, dashes = (2,2))
#		if len(lethtidx) > 2 and plotL3:
#			ax1.plot( [tpL3 - tpRel, tpL3 - tpRel], [-1,1], '--', color = 'black', lw=.25, dashes = (2,2))


	return (ratio[0][0], ratio[2][0])


if __name__ == '__main__':


	#%% setup figure for the timeseries
	fig1 = plt.figure(figsize=(1.8,1.1))
	ax1 = fig1.add_subplot(111)
#	fig1.subplots_adjust(left=0.25, right=.99, top=.99, bottom=0.25)
#	for tl in ax1.get_xticklabels():
#		tl.set_fontsize(8)
#	for tl in ax1.get_yticklabels():
#		tl.set_fontsize(8)
#	ax1.set_ylim((-1.1,1.1))
#
#	# ax1.set_xlim((9,25))
#	ax1.set_xlim(-2.5,9.5)
#	ax1.set_xlabel('Time after L2 (h)', fontsize = 8, fontname = 'Calibri',labelpad=2)
#	ax1.set_ylabel('Norm fluor\n difference (a.u.)', fontsize = 8, fontname = 'Calibri',labelpad=-0)
#	
#	# remove label axis and so on
#	ax1.tick_params(
#		axis='y',          # changes apply to the x-axis
#		which='both',      # both major and minor ticks are affected
#		left='off',      # ticks along the bottom edge are off
#		right='off',         # ticks along the top edge are off
#		labelleft='off', # labels along the bottom edge are off
#		length = 1.5,
#		pad = 1 )
#	ax1.tick_params(
#		axis='x',          # changes apply to the x-axis
#		which='both',      # both major and minor ticks are affected
#		bottom='off',      # ticks along the bottom edge are off
#		top='off',         # ticks along the top edge are off
#		labelbottom='off', # labels along the bottom edge are off
#		length = 1.5,
#		pad = 1 )
#	
#	for axis in ['top','bottom','left','right']:
#		ax1.spines[axis].set_linewidth(0.25)
#	
#	ax1.xaxis.set_tick_params(width=.25)
#	ax1.yaxis.set_tick_params(width=.25)
#
#	ax1.xaxis.set_ticks_position('bottom')
#	ax1.yaxis.set_ticks_position('left')
#
#	ax1.set_xticks([-2,0,2,4,6,8])
#	ax1.set_yticks([-1,-.5,0,.5,1])
#	
#	for label in (ax1.get_xticklabels() + ax1.get_yticklabels()):
#		label.set_fontname('Calibri')
#		label.set_fontsize(8)
#	
#	plt.text(-1.5,.7,'L1',fontname='Calibri',fontsize = 8)
#	plt.text(3,.7,'L2',fontname='Calibri',fontsize = 8)
#	plt.text(8.2,.7,'L3',fontname='Calibri',fontsize = 8)
#	
#	plt.plot([-2.5,9.5],[-1,-1],'--k',dashes = (2,2),lw=.25)
#	plt.plot([-2.5,9.5],[0,0],'--k',dashes = (2,2),lw=.25)
#	plt.plot([-2.5,9.5],[1,1],'--k',dashes = (2,2),lw=.25)
#%%
	cwd = os.getcwd()
	parentDir = os.path.join(os.path.dirname(os.path.dirname(cwd)),'ACVU_data','timelapse')

	'''
	lag2YFP data
	'''
	### plot one worm
	worms = [ os.path.join( parentDir, '161111_lag2YFP+histmCherry', 'C07' ) ]

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

	allFiltData = pickle.load( open( os.path.join( parentDir, 'allDataFilt.pickle' ), 'rb' ) )
	vals = []
	for idx, wormName in enumerate( worms ):

		worm = wormName.split('\\')[-1].split('_')[0]
		path = wormName[:-len(wormName.split('\\')[-1])]
				
		data = allFiltData[idx]
		val = plotRatio( data, path, worm )
		vals.append(val)
		
	vals = np.transpose(vals)

	ax1.hist(vals[0],bins = 20, range = (-1,1), alpha=.3, width = 0.06)
	ax1.hist(vals[1],bins = 20, range = (-1,1), alpha=.3, width = 0.06)

	ax1.set_xlim([-1,1])
	plt.show();
