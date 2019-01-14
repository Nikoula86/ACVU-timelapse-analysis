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
#import generalFunctions as gf
#import pandas as pd
import glob
#from skimage import filters


import matplotlib as mpl
mpl.rcParams['pdf.fonttype'] = 42

def plotRatio( data, count = 0, alpha = 1. ):
	
	times = data[ 1 + data[0].index('timesToL2') ]
	# set tpRel as 2nd dividing cell
#	times = times - times[2][0]
	# set tpRel as 1st dividing cell
#	if len(times[1])>0
#		times = times - times[1][0]

	ratio = data[ 1 + data[0].index('ratioFiltCorr') ]
#	_sum = data[ 1 + data[0].index('absValZ1') ] + data[ 1 + data[0].index('absValZ4') ]
#	ratio *= _sum
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
	dtOtherDirThr = np.sum( np.diff( times[2][(ratio[2]*ratio[2][-1])<0] ) ) >= 1.5
#	print(times[2][ratio[2]>0])
	
	colors = ['k','b','m']
	
#	if ( startOtherDir ) & ( dtOtherDirThr ):
#		count+=1
	for idx, val in enumerate( zip(times,ratio) ):
#		if (tdiv1!=tdiv4) and (np.sign(ratio[2][-1])==-1):
		# plot mothers
		ax1.plot( val[0], val[1],
					'-o', lw = .5, color = colors[idx], alpha = alpha )
	return count

if __name__ == '__main__':


	#%% setup figure for the timeseries
	fig1 = plt.figure(figsize=(5.8,3.1))
	ax1 = fig1.add_subplot(111)
	fig1.subplots_adjust(left=0.15, right=.99, top=.99, bottom=0.15)
	for tl in ax1.get_xticklabels():
		tl.set_fontsize(8)
	for tl in ax1.get_yticklabels():
		tl.set_fontsize(8)
	ax1.set_ylim((-1.1,1.1))

	# ax1.set_xlim((9,25))
	ax1.set_xlim(-2.5,9.5)
	ax1.set_xlabel('Time after L2 (h)', fontsize = 8, fontname = 'Calibri',labelpad=2)
	ax1.set_ylabel('Norm fluor\n difference (a.u.)', fontsize = 8, fontname = 'Calibri',labelpad=-0)
	
	# remove label axis and so on
	ax1.tick_params(
		axis='y',          # changes apply to the x-axis
		which='both',      # both major and minor ticks are affected
		left='off',      # ticks along the bottom edge are off
		right='off',         # ticks along the top edge are off
		labelleft='off', # labels along the bottom edge are off
		length = 1.5,
		pad = 1 )
	ax1.tick_params(
		axis='x',          # changes apply to the x-axis
		which='both',      # both major and minor ticks are affected
		bottom='off',      # ticks along the bottom edge are off
		top='off',         # ticks along the top edge are off
		labelbottom='off', # labels along the bottom edge are off
		length = 1.5,
		pad = 1 )
	
	for axis in ['top','bottom','left','right']:
		ax1.spines[axis].set_linewidth(0.25)
	
	ax1.xaxis.set_tick_params(width=.25)
	ax1.yaxis.set_tick_params(width=.25)

	ax1.xaxis.set_ticks_position('bottom')
	ax1.yaxis.set_ticks_position('left')

#	ax1.set_xticks([-2,0,2,4,6,8])
#	ax1.set_yticks([-1,-.5,0,.5,1])
	
	for label in (ax1.get_xticklabels() + ax1.get_yticklabels()):
		label.set_fontname('Calibri')
		label.set_fontsize(18)
	
	plt.text(-1.5,.7,'L1',fontname='Calibri',fontsize = 8)
	plt.text(3,.7,'L2',fontname='Calibri',fontsize = 8)
	plt.text(8.2,.7,'L3',fontname='Calibri',fontsize = 8)
	
	plt.plot([-2.5,9.5],[-1,-1],'--k',dashes = (2,2),lw=.25)
	plt.plot([-2.5,9.5],[0,0],'--k',dashes = (2,2),lw=.25)
	plt.plot([-2.5,9.5],[1,1],'--k',dashes = (2,2),lw=.25)
	plt.plot([-2.5,9.5],[1/3,1/3],'--k',dashes = (2,2),lw=.25)
	plt.plot([-2.5,9.5],[-1/3,-1/3],'--k',dashes = (2,2),lw=.25)
#%%
	cwd = os.getcwd()
	parentDir = os.path.join(os.path.dirname(os.path.dirname(cwd)),'ACVU_data','timelapse')

#%%
	'''
	PLOT THE RATIO
	'''

#	allData = pickle.load( open( os.path.join( parentDir, 'allData_lag2multi.pickle' ), 'rb' ) )
	allData = pickle.load( open( os.path.join( parentDir, 'allData.pickle' ), 'rb' ) )
	vals = []
	counts = 0
	for data in allData:
		if '161111_lag2YFP+histmCherry\\C05' in data[1]:
#		if '180106_lag2multi+JVZ32\\C69' in data[1]:
			print(data[1])
	
			c = plotRatio( data )
			counts += c

	print(counts)
	plt.show();
