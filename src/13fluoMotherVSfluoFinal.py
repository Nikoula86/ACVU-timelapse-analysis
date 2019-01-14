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
### plot one worm
# worms = ['X:\\Simone\\160930_lag2GFP+mCherry\\C09','X:\\Simone\\160930_lag2GFP+mCherry\\C10','X:\\Simone\\160930_lag2GFP+mCherry\\C11','X:\\Simone\\160930_lag2GFP+mCherry\\C15','X:\\Simone\\160930_lag2GFP+mCherry\\C19','X:\\Simone\\161002_lag2GFP+mCherry\\C01']
# worms = ['X:\\Simone\\160930_lag2GFP+mCherry\\C01','X:\\Simone\\160930_lag2GFP+mCherry\\C02','X:\\Simone\\160930_lag2GFP+mCherry\\C03','X:\\Simone\\160930_lag2GFP+mCherry\\C06','X:\\Simone\\160930_lag2GFP+mCherry\\C08','X:\\Simone\\161002_lag2GFP+mCherry\\C04']
### plot all
paths = [ 'W:\\Simone\\160930_lag2GFP+mCherry', 'W:\\Simone\\161002_lag2GFP+mCherry' ]
worms1 = []
for path in paths:
	worms1.append( glob.glob( os.path.join( path, '*cellFluo.pickle' ) ) )
worms = []
for exp in worms1:
	for w in exp:
		worms.append(w)
print(worms)

'''
lag2YFP strain
'''
# paths = [ 'W:\\Simone\\160516_lag2_YFP_hist_mCherry', 'W:\\Simone\\160226_lag2_YFP_hist_mCherry' ]
# worms1 = []
# for path in paths:
# 	worms1.append( glob.glob( os.path.join( path, '*cellFluo.pickle' ) ) )
# worms = []
# for exp in worms1:
# 	for w in exp:
# 		worms.append(w)
# print(worms)

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

	bckg = cFluoDF[ 'b_' + cname[0] ].ix[ pd.notnull( cFluoDF[ 'b_' + cname[0] ]._488nm ) ]
	bckgFluo = bckg.ix[ bckg.tidx == tidx, '_488nm' ].values[0]

	return ( cellFluo - bckgFluo, times )

data = []
dts = []

dt = -1#np.linspace(0,5,11)
print(dt)

fig1 = plt.figure(figsize=(5.8,5.8))
ax1 = fig1.add_subplot(111)
fig1.subplots_adjust(left=0.1, right=.9, top=.9, bottom=0.1)
for tl in ax1.get_xticklabels():
	tl.set_fontsize(18)
for tl in ax1.get_yticklabels():
	tl.set_fontsize(18)
# ax1.set_ylim((-6,6))

ax1.set_xlim((-1,1))
# ax1.set_xlim(-.6,1.1)
ax1.plot([0,0],[-1,1],'--k')
ax1.plot([-1,1],[0,0],'--k')

xValues = []
yValues = []
for widx, worm in enumerate(worms):
	print(worm)

	### load parameters and times dataframes
	w = worm.split('\\')[-1].split('_')[0]
	path = worm[:-len(worm.split('\\')[-1])]
	# print(path,w)
	paramsDF = load_data_frame( path, w + '_01params.pickle' )
	timesDF = load_data_frame( path, w + '_01times.pickle' )
	gpDF = load_data_frame( path, w + '_02gonadPos.pickle' )
	cellPosDF = load_data_frame( path, w + '_04cellPos.pickle' )
	cellOutDF = load_data_frame( path, w + '_05cellOut.pickle' )
	cellFluoDF = load_data_frame( path, w + '_06cellFluo.pickle' )
	# apdvPosDF = load_data_frame( path, w + '_08apdvPos.pickle' )

	# print(cellPosDF['1.ppp'].keys())
	# print(cellFluoDF['1.ppp'].ix[ cellFluoDF['1.ppp'].tidx == 100 ])
	# print(timesDF.keys())

	# find time of divisions
	tdiv1 = np.min( timesDF.ix[ pd.notnull( cellPosDF[ '1.ppp' ].X ) ].timesRel )
	tdiv4 = np.min( timesDF.ix[ pd.notnull( cellPosDF[ '4.aaa' ].X ) ].timesRel )
	tDiv = [ np.min([tdiv1,tdiv4]) , np.max([tdiv1,tdiv4]) ]
	# find tidx of divisions
	tDivIdx = [ timesDF.ix[ timesDF.timesRel == tDiv[0], 'tidxRel' ].values[0]-1, timesDF.ix[ timesDF.timesRel == tDiv[1], 'tidxRel' ].values[0]-1 ]
	# find tMax
	tMax = np.max( timesDF.ix[ pd.notnull( cellPosDF[ '1.ppp' ].X ) ].timesRel )
	tMaxIdx = timesDF.ix[ timesDF.timesRel == tMax, 'tidxRel' ].values[0]

	# extract mother cells info
	# print(cellFluoDF['1.ppp'].keys())
	mother1 = cellFluoDF['1.pp'].ix[ ( pd.notnull( cellFluoDF['1.pp']._488nm ) ) & ( cellFluoDF['1.pp'].tidx<tDivIdx[0] ) & ( cellFluoDF['1.pp'].tidx>(tDivIdx[0]-3) ), ['_488nm','tidx'] ]
	bckg1 = cellFluoDF[ 'b_1' ].ix[ ( pd.notnull( cellFluoDF['1.pp']._488nm ) ) & ( cellFluoDF['1.pp'].tidx<tDivIdx[0] ) & ( cellFluoDF['1.pp'].tidx>(tDivIdx[0]-3) ), ['_488nm','tidx'] ]
	mother1._488nm = ( mother1._488nm - bckg1._488nm ) / ( bckg1._488nm )
	# mother1['times'] = timesDF.ix[ pd.notnull( cellFluoDF['1.pp']._488nm ), 'timesRel' ].values - tDiv[0]

	mother4 = cellFluoDF['4.aa'].ix[ ( pd.notnull( cellFluoDF['4.aa']._488nm ) ) & ( cellFluoDF['4.aa'].tidx<tDivIdx[0] ) & ( cellFluoDF['1.pp'].tidx>(tDivIdx[0]-3) ), ['_488nm','tidx'] ]
	bckg4 = cellFluoDF[ 'b_4' ].ix[ ( pd.notnull( cellFluoDF['4.aa']._488nm ) ) & ( cellFluoDF['4.aa'].tidx<tDivIdx[0] ) & ( cellFluoDF['1.pp'].tidx>(tDivIdx[0]-3) ), ['_488nm','tidx'] ]
	mother4._488nm = ( mother4._488nm - bckg4._488nm ) / ( bckg4._488nm )
	# mother4['times'] = timesDF.ix[ pd.notnull( cellFluoDF['4.aa']._488nm ), 'timesRel' ].values - tDiv[0]

	xValues.append( np.mean((mother1._488nm - mother4._488nm)/(mother1._488nm + mother4._488nm)) )

	# extract cells info
	# print(cellFluoDF['1.ppp'].keys())
	cell1 = cellFluoDF['1.ppp'].ix[ cellFluoDF['1.ppp'].tidx==tMaxIdx, ['_488nm','tidx'] ]
	bckg1 = cellFluoDF[ 'b_1' ].ix[ cellFluoDF['1.ppp'].tidx==tMaxIdx, ['_488nm','tidx'] ]
	fluo1 = ( ( cell1._488nm - bckg1._488nm ) / ( bckg1._488nm ) ).values[0]

	cell4 = cellFluoDF['4.aaa'].ix[ cellFluoDF['4.aaa'].tidx==tMaxIdx, ['_488nm','tidx'] ]
	bckg4 = cellFluoDF[ 'b_4' ].ix[ cellFluoDF['4.aaa'].tidx==tMaxIdx, ['_488nm','tidx'] ]
	fluo4 = ( ( cell4._488nm - bckg4._488nm ) / ( bckg4._488nm ) ).values[0]

	yValues.append( (fluo1 - fluo4)/(fluo1 + fluo4) )

print(xValues,yValues)
ax1.plot(xValues,yValues,'o')
for idx,w in enumerate(worms):
	ax1.text(xValues[idx],yValues[idx],w.split('\\')[-1].split('_')[0])
plt.show()
