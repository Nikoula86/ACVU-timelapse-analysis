import glob
from tifffile import *
from generalFunctions import *
import numpy as np
import PIL
from PIL import Image, ImageDraw, ImageFont
import os.path
import matplotlib.pyplot as plt
from matplotlib import cm
import pickle
import scipy.interpolate as ip
from scipy.ndimage import interpolation
import matplotlib as mpl
from matplotlib.colors import LinearSegmentedColormap
mpl.rcParams['pdf.fonttype'] = 42


def plotPhaseFluorescence( path, worm, ax1, color = 'k', channel = '488nm', filt = False, sigma = 15/60 ):

	timesDF = load_data_frame( path, worm + '_01times.pickle' )
	cellPosDF = load_data_frame( path, worm + '_04cellPos.pickle' )
	cellOutDF = pd.DataFrame({})#load_data_frame( path, worm + '_05cellOut.pickle' )
	cellFluoDF = load_data_frame( path, worm + '_15cellFluo_q&d.pickle' )

	lineages = [['1.p','1.pp','1.ppp'],['4.a','4.aa','4.aaa']]
	colors = np.array( [ np.array( [cm.Blues(int(i))[:3] for i in (np.arange(len(lineages[0]))*254./(len(lineages[0])-1) / 254. ) * 127 + 127] ),
				np.array( [cm.Reds(int(i))[:3] for i in (np.arange(len(lineages[0]))*254./(len(lineages[0])-1) / 254. ) * 127 + 127 ] ) ] )

	### find ecdysis timepoint

	ecd = np.loadtxt( open( os.path.join( path, 'skin.txt'), 'rb' ) )
	# load ecdysis data
	index = np.where( ecd[:,0] == float(worm[1:]) )
	mintp = np.min( [ i for i in ecd[index, 1:6][0][0] if i >= 0 ] )
	lethtidx = ecd[ index, 2:6 ][0][0] - mintp
	tpL2 = timesDF.ix[ timesDF.tidxRel == lethtidx[0], 'timesRel' ].values[0]

	### set relative time

	# # relative to L2 start
	# tpRel = tpL2

	# # relative to hatching time
	# tpRel=0

	# relative to first cell division
	tdiv1 = np.min( timesDF.ix[ pd.notnull( cellPosDF[ '1.ppp' ].X ) ].timesRel )
	tdiv2 = np.min( timesDF.ix[ pd.notnull( cellPosDF[ '4.aaa' ].X ) ].timesRel )
	tpRel = np.min([tdiv1,tdiv2])

	print(tpRel)

	### phase plot
	dataMother = pd.DataFrame( { 'tidx': timesDF.tidxRel,
							'times': timesDF.timesRel,
							'signal1': np.nan,
							'signal4': np.nan,
							'lineage': np.nan } )
	# print(data)
	# colorsRatio = ['black','blue','magenta']

	for idx, trow in timesDF.iterrows():

		currentCells = extract_current_cell_pos( cellPosDF, trow.tidxRel )
		for cname in currentCells.cname:
			# print(cellFluoDF[cname].ix[ cellFluoDF[cname].tidx == trow.tidxRel, '_488nm' ].values[0])
			currentCells.ix[ currentCells.cname == cname, '_488nm' ] = cellFluoDF[cname].ix[ cellFluoDF[cname].tidx == trow.tidxRel, '_488nm' ].values
		print(currentCells)

		# print(cell1.cname,cell1.fluo488nm,cell2.cname,cell2.fluo488nm)
		dataMother.ix[ data.tidx == trow.tidxRel, 'signal1' ] = cell1._488nm
		dataMother.ix[ data.tidx == trow.tidxRel, 'signal4' ] = cell4._488nm
		# print(cell1.cname,cell2.cname,int(np.min([len(cell1.cname),len(cell2.cname)])-3),colorsRatio[int(np.min([len(cell1.cname),len(cell2.cname)])-3)])

	data = data.ix[ pd.notnull( data.signal1 ) ].reset_index()
	# print(data)

	ax1.plot(data.signal1,data.signal4)
	# for idx, d in data.iterrows():
	# 	# print(idx,idx+1)
	# 	ax1.plot( [ d.times-tpRel, data.times.values[np.clip(idx+1,0,len(data)-1)]-tpRel ], [ d.ratio, data.ratio.values[np.clip(idx+1,0,len(data)-1)] ], '-', color = colorsRatio[int(d.lineage)], lw=2 )


'''
lag2YFP data
'''
# ### plot one worm
worms = [ 'W:\\Simone\\161030_lag2YFP+mCherry\\C23' ]
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

### setup figure for the timeseries
fig1 = plt.figure(figsize=(5.8,5.8))
ax1 = fig1.add_subplot(111)
fig1.subplots_adjust(left=0.15, right=.95, top=.95, bottom=0.15)
for tl in ax1.get_xticklabels():
	tl.set_fontsize(18)
for tl in ax1.get_yticklabels():
	tl.set_fontsize(18)

for idx, worm in enumerate( worms ):

	w = worm.split('\\')[-1].split('_')[0]
	path = worm[:-len(worm.split('\\')[-1])]
			
	plotPhaseFluorescence( path, w, ax1, channel = '488nm', filt = True, sigma = 15/60 )


plt.show()