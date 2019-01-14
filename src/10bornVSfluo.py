import glob
from tifffile import *
from generalFunctions import *
import numpy as np
import PIL
from PIL import Image, ImageDraw, ImageFont
import os.path
import matplotlib.pyplot as plt
import pickle
import scipy.interpolate as ip
from scipy.ndimage import interpolation
import matplotlib as mpl
from matplotlib.colors import LinearSegmentedColormap
mpl.rcParams['pdf.fonttype'] = 42

'''
lag2GFP strain
'''
# paths = [ 'W:\\Simone\\160930_lag2GFP+mCherry', 'W:\\Simone\\161002_lag2GFP+mCherry', 'W:\\Simone\\161004_lag2GFP+mCherry' ]

'''
lag2YFP strain
'''
paths = [ 'W:\\Simone\\161030_lag2YFP+mCherry', 'W:\\Simone\\161108_lag2YFP+mCherry', 'W:\\Simone\\161111_lag2YFP+mCherry' ]

'''
START
'''
worms1 = []
for path in paths:
	worms1.append( glob.glob( os.path.join( path, '*cellFluo_q&d.pickle' ) ) )
worms = []
for exp in worms1:
	for w in exp:
		worms.append(w)
worms.sort()
print(worms)

def tBorn( cname, tDF, cPosDF ):
	cDF = cPosDF[ cname ].ix[ pd.notnull( cPosDF[ cname ].X ) ]
	tidxborn = np.min( cDF.tidx )
	tborn = tDF.ix[ tDF.tidxRel == tidxborn, 'timesRel' ].values[0]
	return ( tborn, tidxborn )

def findData( cname, tLastBorn, tDF, cFluoDF, dt ):
	cDF = cFluoDF[ cname ].ix[ pd.notnull( cFluoDF[ cname ]._488nm ) ]
	# print(cDF)
	# print(tLastBorn,dt)

	if dt != -1:
		tidxmin = tDF.ix[ tDF.timesRel > ( tLastBorn + dt ), 'tidxRel' ].values[0]
		print(tidxmin)
		tidx = cDF.ix[ cDF.tidx >= tidxmin, 'tidx' ].values[0]
	else:	# if dt = -1, find the last tp
		tidx = np.max(cDF.tidx)

	times = tDF.ix[ tDF.tidxRel == tidx, 'timesRel' ].values[0]
	cellFluo = cDF.ix[ cDF.tidx == tidx, '_488nm' ].values[0]

	# bckg = cFluoDF[ 'b_' + cname[0] ].ix[ pd.notnull( cFluoDF[ 'b_' + cname[0] ]._488nm ) ]
	# bckgFluo = bckg.ix[ bckg.tidx == tidx, '_488nm' ].values[0]

	return ( cellFluo, times )

data = []

dts = [-1]#np.linspace(0,5,11)
print(dts)

for dtidx, dt in enumerate(dts):
	print(dt)

	fig1 = plt.figure(figsize=(6.8,6.8))
	ax1 = fig1.add_subplot(111)
	fig1.subplots_adjust(left=0.1, right=.9, top=.9, bottom=0.1)
	for tl in ax1.get_xticklabels():
		tl.set_fontsize(18)
	for tl in ax1.get_yticklabels():
		tl.set_fontsize(18)
	# ax1.set_ylim((-6,6))

	ax1.set_xlim((-2,2))
	# ax1.set_xlim(-.6,1.1)
	ax1.plot([0,0],[-1,1],'--k')
	ax1.plot([-2,2],[0,0],'--k')
	ax1.set_xlabel( 'Division Time diff (TdivZ1 - TdivZ4)', fontsize = 18 )
	ax1.set_ylabel( 'Outcome (Z1-Z4)/(Z1+Z4)', fontsize = 18 )
	
	for widx, worm in enumerate(worms):
		print(worm)

		### load parameters and times dataframes
		w = worm.split('\\')[-1].split('_')[0]
		path = worm[:-len(worm.split('\\')[-1])]
		# print(path,w)
		timesDF = load_data_frame( path, w + '_01times.pickle' )
		cellPosDF = load_data_frame( path, w + '_04cellPos.pickle' )
		cellFluoDF = load_data_frame( path, w + '_15cellFluo_q&d.pickle' )

		# background correction
		for key in cellFluoDF.keys():
			if key[0]!='b':
				cellFluoDF[key]._488nm -= cellFluoDF['b_'+key[0]]._488nm

		# print(cellPosDF['1.ppp'].keys())
		# print(cellFluoDF['1.ppp'].ix[ cellFluoDF['1.ppp'].tidx == 100 ])
		# print(timesDF.keys())

		### CREATE CELL DICTIONARY
		cell1 = { 'cname': '1.ppp', 'tidxborn': np.nan, 'tborn': np.nan, 'exp': np.nan, 'texp': np.nan }
		cell4 = { 'cname': '4.aaa', 'tidxborn': np.nan, 'tborn': np.nan, 'exp': np.nan, 'texp': np.nan }

		### FIND DIVISION TIMES
		cell1['tborn'], cell1['tidxborn'] = tBorn(cell1['cname'],timesDF,cellPosDF)
		cell4['tborn'], cell4['tidxborn'] = tBorn(cell4['cname'],timesDF,cellPosDF)

		### FIND FLUORESCENCE
		cell1['exp'], cell1['texp'] = findData( cell1['cname'], np.max( [ cell1['tborn'], cell4['tborn'] ] ), timesDF, cellFluoDF, dt=dt )
		cell4['exp'], cell4['texp'] = findData( cell4['cname'], np.max( [ cell1['tborn'], cell4['tborn'] ] ), timesDF, cellFluoDF, dt=dt )

		ax1.plot( cell1['tborn'] - cell4['tborn'], np.min( [ 1, (cell1['exp'] - cell4['exp'])/(cell1['exp'] + cell4['exp']) ] ), 'ok', alpha = .5 )
		ax1.text( cell1['tborn'] - cell4['tborn'], np.min( [ 1, (cell1['exp'] - cell4['exp'])/(cell1['exp'] + cell4['exp']) ] ), str(widx), fontsize = 15 )
		print(np.max( [ cell1['tborn'], cell4['tborn'] ] ), cell1['texp'],cell4['texp'])
		print(cell1)
		print(cell4)

	plt.show()
	# plt.savefig('X:\\Simone\\160930_lag2GFP+mCherry\\' + 'tp%s.png'%dtidx)
	# plt.close(fig1)
