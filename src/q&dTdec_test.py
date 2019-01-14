# -*- coding: utf-8 -*-
"""
Created on Tue Dec  1 10:31:00 2015

@author: kienle
"""

import numpy as np
import matplotlib.pyplot as plt
import pickle
import os
#from matplotlib import cm
from generalFunctions import *
#from skimage import filters


import matplotlib as mpl
mpl.rcParams['pdf.fonttype'] = 42

def findTDecision( data ):

	### find dTdiv
	times = data[data[0].index('timesToL2')+1]
	# set tpRel as 2nd dividing cell
	times = times - times[2][0]
	# set tpRel as 1st dividing cell
#	if len(times[1])>0:
#		times = times - times[1][0]

	if len( times[1] ) == 0:
		dTdiv = 0.
	else:
		dTdiv = times[2][0] - times[1][0]
	
	ratio = data[data[0].index('ratioFiltCorr')+1]
	### find time to decide
	# print(tdiv1,tdiv2)
	thr = 1/3
	tDec = times[2][ np.abs(ratio[2]) < thr ]
	if len( tDec ) > 0:
		tDec = tDec[-1]
	else:
		tDec = times[2][-1]

	return ( dTdiv, tDec )

def plotTdec(allData,color):

	dtDiv = []
	tDec = []
	
	for idx, data in enumerate( allData ):

		dt,tdec = findTDecision( data )
		dtDiv.append( dt )
		tDec.append( tdec )

		print(idx,data[1], dt, tdec)
#		ax1.text( dt, tdec )

	ax1.scatter( dtDiv, tDec, color = color , s = 100, lw=0, alpha = .7 )

	bins = np.linspace(-5,305,32) / 60.
	digitized = np.digitize(dtDiv, bins)
	bin_means = [np.array(tDec)[digitized == i].mean() for i in range(1, len(bins))]
	bin_std = [np.array(tDec)[digitized == i].std() for i in range(1, len(bins))]
	
	dtDivMean = np.array(bins)[np.isfinite(bin_means)] + 5/60.
	tDecMean = np.array(bin_means)[np.isfinite(bin_means)]
	tDecStd = np.array(bin_std)[np.isfinite(bin_means)]
	ax1.plot(dtDivMean, tDecMean, '-', color = color, lw=2)
#	ax1.fill_between(dtDivMean, tDecMean-tDecStd, tDecMean+tDecStd, facecolor = color, alpha = .2, lw=0. )
	
	print('\n tDec for dtDiv < 30min:') # 5.51 +- 0.77
	print(np.mean(tDecMean[:4]), '+-', np.std(tDecMean[:4]))
	print('\n tDec for dtDiv > 30min:') # 3.60 +- 0.84
	print(np.mean(tDecMean[5:]), '+-', np.std(tDecMean[5:]))
	print('\n tDec total:') # 3.60 +- 0.84
	print(np.mean(tDec), '+-', np.std(tDec),'\n\n')

if __name__ == '__main__':

	### setup figure for the phasePlot
	fig1 = plt.figure(figsize=(6.8,4.8))
	ax1 = fig1.add_subplot(111)
	fig1.subplots_adjust(left=0.15, right=.95, top=.95, bottom=0.15)
	for tl in ax1.get_xticklabels():
		tl.set_fontsize(18)
	for tl in ax1.get_yticklabels():
		tl.set_fontsize(18)
	# ax1.set_ylim((0,3.5))

	ax1.set_xlim(0,3.1)
	ax1.set_ylim(0,8)

	# ax1.plot( [-2,10], [0,0], '--', color = 'black', lw=1)
	# ax1.plot( [-2,10], [1/3,1/3], '--', color = 'black', lw=1)
	# ax1.plot( [-2,10], [-1/3,-1/3], '--', color = 'black', lw=1)

	ax1.set_xlabel('dt Div', fontsize = 18)
	ax1.set_ylabel('t Decision', fontsize = 18)

	cwd = os.getcwd()
	parentDir = os.path.join(os.path.dirname(os.path.dirname(cwd)),'ACVU_data','timelapse')




	'''
	plot WT data
	'''
	allData = pickle.load( open( os.path.join( parentDir, 'allData.pickle' ), 'rb' ) )
	plotTdec(allData,'grey')

	'''
	plot lag2multi copies data
	'''
	allData = pickle.load( open( os.path.join( parentDir, 'allData_lag2multi.pickle' ), 'rb' ) )
	plotTdec(allData,'red')


		
	plt.show()
