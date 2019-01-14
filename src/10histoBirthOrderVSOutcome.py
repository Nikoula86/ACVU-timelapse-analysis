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
paths = [ 'Y:\\Simone\\161030_lag2YFP+histmCherry',
			'Y:\\Simone\\161108_lag2YFP+histmCherry',
			'Y:\\Simone\\161111_lag2YFP+histmCherry',
			'Y:\\Simone\\170312_lag2YFP+histmCherry' ]

'''
LIN12BALhistmCherry
'''
# paths = [ 'T:\\Simone\\170405_LIN12BAL+histmCherry', 'T:\\Simone\\170407_LIN12BAL+histmCherry', 'T:\\Simone\\170419_LIN12BAL+histmCherry' ] 

'''
LIN12GFP strain
'''
# paths = [ 'W:\\Simone\\160420_LIN12GFP+histmCherry', 'W:\\Simone\\161212_LIN12GFP+histmCherry', 'W:\\Simone\\170220_LIN12GFP+histmCherry', 'W:\\Simone\\170221_LIN12GFP+histmCherry' ]

'''
HLH2GFP strain
'''
# paths = [ 'W:\\Simone\\160407_HLH2GFP+histmCherry', 'W:\\Simone\\160502_HLH2GFP+histmCherry', 'W:\\Simone\\170211_HLH2GFP+histmCherry' ]

'''
LIN12nullunc strain
'''
#paths = [ 'Y:\\Simone\\170711_LIN12unc+histmCherry',
#			'Y:\\Simone\\170714_LIN12unc+histmCherry',
#			'Y:\\Simone\\170722_LIN12unc+histmCherry',
#			'Y:\\Simone\\170727_LIN12unc+histmCherry',
#			'Y:\\Simone\\170730_LIN12unc+histmCherry'] 
			
'''
LAGintegration strain
'''
#paths = [ 
#		'Y:\\Simone\\180106_lag2multi+JVZ32',
#		'Y:\\Simone\\180117_lag2multi+JVZ32'
#		] 

'''
START
'''
worms1 = []
for path in paths:
	worms1.append( glob.glob( os.path.join( path, '*01times.pickle' ) ) )
worms = []
for exp in worms1:
	for w in exp:
		worms.append(w)
worms.sort()
print(worms)
print('\nN:',len(worms),'\n')

binDT = 10 #binning in minutes

tMax = 3.
data = pd.DataFrame({'absdt':np.linspace(0,tMax*60./binDT,tMax*60./binDT+1),
					'_second': 0*np.linspace(0,tMax*60./binDT,tMax*60./binDT+1),
					'_first': 0*np.linspace(0,tMax*60./binDT,tMax*60./binDT+1)})
# print(data)
for widx, worm in enumerate(worms):
	# print(worm)

	### load parameters and times dataframes
	w = worm.split('\\')[-1].split('_')[0]
	path = worm[:-len(worm.split('\\')[-1])]
#	print(path,w)
	# print(path,w)
	timesDF = load_data_frame( path, w + '_01times.pickle' )

	# print(cellPosDF['1.ppp'].keys())
	# print(cellFluoDF['1.ppp'].ix[ cellFluoDF['1.ppp'].tidx == 100 ])
	# print(timesDF.keys())

	### load txt file
	fname = os.path.join( path, 'birthOrderToOutcome.txt')

	# Using the newer with construct to close the file automatically.
	with open(fname, 'r') as f:
	    cdata = f.readlines()
	try:
		cdata = np.array([[int(i) for i in x.strip().split('\t')] for x in cdata[1:]])
		cdata = cdata[cdata[:,0]==int(w[1:])][0,1:]
	except:
		continue
	# print(cdata)

	tb1 = timesDF.ix[ timesDF.tidxRel == cdata[0], 'timesRel' ].values[0] 
	tb4 = timesDF.ix[ timesDF.tidxRel == cdata[1], 'timesRel' ].values[0] 

	### CALCULATE DATA
	dT = round( np.abs(tb1-tb4) * 60. / 10. )
	# print('ciao',dT,tb1,tb4)
	if tb1 > tb4:
		if cdata[2] == 1:
			data.ix[ data.absdt == dT, '_second' ] += 1
		if cdata[2] == 4:
			data.ix[ data.absdt == dT, '_first' ] += 1
	elif tb4 > tb1:
		if cdata[2] == 4:
			data.ix[ data.absdt == dT, '_second' ] += 1
		if cdata[2] == 1:
			data.ix[ data.absdt == dT, '_first' ] += 1

data.absdt /= ( 60. / 10 )
data = data.ix[ (data.absdt<=tMax) & (data.absdt>0) ]

'''
PLOT ABSOLUTE VALUES
'''

# print(data)
fig1 = plt.figure(figsize=(3.5,3.5))
ax1 = fig1.add_subplot(111)
fig1.subplots_adjust(left=0.20, right=.93, top=.93, bottom=0.20)
for tl in ax1.get_xticklabels():
	tl.set_fontsize(15)
for tl in ax1.get_yticklabels():
	tl.set_fontsize(15)

# newdata = 

ax1.set_xlim(0,tMax)
ax1.set_xlabel('Time between births', fontsize = 15)
ax1.set_ylabel('Number of animals', fontsize = 15)
# regular histograms
# ax1.set_ylim(0,10)
# width = 0.35
# print((data.dt.values - width/2.), data._second.values, width)
# tot = data._first + data._second

# ax1.bar( (data.dt )*10, data._second, width*10, color = 'blue')
# ax1.bar( (data.dt - width)*10, data._first, width*10, color = 'red')
# ax1.set_xticks(np.arange(0, 121, 30))

# histograms Stacked
#ax1.set_ylim(0,1.)
width = 1.
print(data)
tot = data._first + data._second
f = np.sum( data.ix[data.absdt<=0.5,'_first'] )
s = np.sum( data.ix[data.absdt<=0.5,'_second'] )
print( f, s, f/(f+s), np.sum(tot) )

ax1.bar( (data.absdt-width/(60./binDT)), data._second, width/(60./binDT), color = '#67a9cf')
ax1.bar( (data.absdt-width/(60./binDT)), data._first, width/(60./binDT), color = '#ef8a62',  bottom=data._second)
ax1.set_xticks(np.arange(0, 2.01, 0.5))
ax1.set_title( 'Nt=' + str(int(np.sum(tot))) + ', f=' + str(int(f)) )

'''
PLOT FRACTIONS
'''

# print(data)
fig2 = plt.figure(figsize=(3.5,3.5))
ax2 = fig2.add_subplot(111)
fig2.subplots_adjust(left=0.20, right=.93, top=.93, bottom=0.20)
for tl in ax2.get_xticklabels():
	tl.set_fontsize(15)
for tl in ax2.get_yticklabels():
	tl.set_fontsize(15)
tMax = 3.

# newdata = 

ax2.set_xlim(0,tMax)
ax2.set_xlabel('Time between births', fontsize = 15)
ax2.set_ylabel('Fraction of animals', fontsize = 15)
# regular histograms
# ax1.set_ylim(0,10)
# width = 0.35
# print((data.dt.values - width/2.), data._second.values, width)
# tot = data._first + data._second

# ax1.bar( (data.dt )*10, data._second, width*10, color = 'blue')
# ax1.bar( (data.dt - width)*10, data._first, width*10, color = 'red')
# ax1.set_xticks(np.arange(0, 121, 30))

# histograms Stacked
ax2.set_ylim(0,1.)
width = 1.

tot = data._first + data._second
f = np.sum( data.ix[data.absdt<=0.5,'_first'] )
s = np.sum( data.ix[data.absdt<=0.5,'_second'] )
print( f, s, f/(f+s), np.sum(tot) )

data._first = data._first / tot
data = data.fillna(0.)
data._second = 1. - data._first
print(data)
# data._second = data._second / tot

ax2.bar( (data.absdt-width/(60./binDT)), data._second, width/(60./binDT), color = '#67a9cf')
ax2.bar( (data.absdt-width/(60./binDT)), data._first, width/(60./binDT), color = '#ef8a62',  bottom=data._second)
ax2.set_xticks(np.arange(0, 2.01, 0.5))
ax2.set_title( 'Nt=' + str(int(np.sum(tot))) + ', f=' + str(int(f)) )


plt.show()


