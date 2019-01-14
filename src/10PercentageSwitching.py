# -*- coding: utf-8 -*-
"""
Created on Thu Jan  4 16:52:25 2018

@author: gritti
"""

import glob
import pandas as pd
#from tifffile import *
import generalFunctions as gf
import numpy as np
#import PIL
#from PIL import Image, ImageDraw, ImageFont
import os.path
import os
import matplotlib.pyplot as plt
import pickle
#import scipy.interpolate as ip
#from scipy.ndimage import interpolation
import matplotlib as mpl
#from matplotlib.colors import LinearSegmentedColormap
mpl.rcParams['pdf.fonttype'] = 42

#%%
cwd = os.getcwd()
parentDir = os.path.join(os.path.dirname(os.path.dirname(cwd)),'ACVU_data','timelapse')
'''
lag2GFP strain
'''
# paths = [ '160930_lag2GFP+mCherry', '161002_lag2GFP+mCherry', '161004_lag2GFP+mCherry' ]

'''
lag2YFP strain
'''
dataFolders = [ '161030_lag2YFP+histmCherry', '161108_lag2YFP+histmCherry', '161111_lag2YFP+histmCherry', '170312_lag2YFP+histmCherry' ]

'''
LIN12GFP strain
'''
# paths = [ '160420_LIN12GFP+histmCherry', '161212_LIN12GFP+histmCherry', '170220_LIN12GFP+histmCherry', '170221_LIN12GFP+histmCherry' ]

'''
HLH2GFP strain
'''
# paths = [ '160407_HLH2GFP+histmCherry', '160502_HLH2GFP+histmCherry', '170211_HLH2GFP+histmCherry' ]


paths = [os.path.join(parentDir,dataFolder) for dataFolder in dataFolders]

#%%
'''
LOAD DATA
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

#%%
binDT = 10 #binning in minutes

data = pd.DataFrame({'absdt':np.linspace(0,12,13),
					'_switch': 0*np.linspace(0,12,13),
					'_noswitch': 0*np.linspace(0,12,13)})

# print(data)
allFiltData = pickle.load( open( os.path.join( parentDir, 'allData.pickle' ), 'rb' ) )
count = 0

for widx, worm in enumerate(worms):
	# print(worm)

	### load parameters and times dataframes
	w = worm.split('\\')[-1].split('_')[0]
	path = worm[:-len(worm.split('\\')[-1])]
	# print(path,w)

	### calculate if switched or not
	d = allFiltData[widx]

	times = d[ 1 + d[0].index('timesToL2') ]
	# set tpRel as 1st dividing cell
	if len(times[1])>0:
		times = times - times[1][0]

	ratio = d[ 1 + d[0].index('ratioFilt') ]
		
	### plot data
	startOtherDir = ( (ratio[2][0]*ratio[2][-1]) < 0 )
	dtOtherDirThr = ( np.sum( np.diff( times[2][(ratio[2]*ratio[2][-1])<0] ) ) >= 1.5 )
	
	### CALCULATE DATA
#	dT = round( np.abs(tb1-tb4) * 60. / 10. )
	if len(times[1])>0:
		dT = round( 60 * ( times[2][0] - times[1][0] ) / 10 )
	else:
		dT = 0
#	print('ciao',dT)

	if startOtherDir and dtOtherDirThr:
		count += 1
#	if dtOtherDirThr:
		data.ix[ data.absdt == dT, '_switch' ] += 1
	else:
		data.ix[ data.absdt == dT, '_noswitch' ] += 1
print(count)
data.absdt /= ( 60. / 10 )
#data = data.ix[ (data.absdt<=2) & (data.absdt>0) ]

#%%
### setup figure
fig1 = plt.figure(figsize=(3.7,3.7))
ax1 = fig1.add_subplot(111)
fig1.subplots_adjust(left=0.25, right=.95, top=.95, bottom=0.25)
for tl in ax1.get_xticklabels():
	tl.set_fontsize(8)
for tl in ax1.get_yticklabels():
	tl.set_fontsize(8)

ax_min = 0
ax_max = 2.

#ax1.set_ylim(0.,1.)
ax1.set_xlim(ax_min,ax_max)

ax1.set_xlabel('Time between births (h)', fontsize = 8, fontname = 'Calibri',labelpad=4)
ax1.set_ylabel('Fraction of animals', fontsize = 8, fontname = 'Calibri',labelpad=4)

# remove label axis and so on
ax1.tick_params(
	axis='y',          # changes apply to the x-axis
	which='both',      # both major and minor ticks are affected
	left='off',      # ticks along the bottom edge are off
	right='off',         # ticks along the top edge are off
	labelleft='off', # labels along the bottom edge are off
	length = 1.5,
	pad = 2 )
ax1.tick_params(
	axis='x',          # changes apply to the x-axis
	which='both',      # both major and minor ticks are affected
	bottom='off',      # ticks along the bottom edge are off
	top='off',         # ticks along the top edge are off
	labelbottom='off', # labels along the bottom edge are off
	length = 1.5,
	pad = 4 )

for axis in ['top','bottom','left','right']:
	ax1.spines[axis].set_linewidth(0.25)

ax1.xaxis.set_tick_params(width=.25)
ax1.yaxis.set_tick_params(width=.25)

ax1.xaxis.set_ticks_position('bottom')
ax1.yaxis.set_ticks_position('left')

ax1.set_xticks([0,.5,1,1.5,2])
#ax1.set_yticks([0,.5,1])

for label in (ax1.get_xticklabels() + ax1.get_yticklabels()):
	label.set_fontname('Calibri')
	label.set_fontsize(8)

#%%
# regular histograms
# ax1.set_ylim(0,10)
# width = 0.35
# print((data.dt.values - width/2.), data._second.values, width)
# tot = data._first + data._second

# ax1.bar( (data.dt )*10, data._second, width*10, color = 'blue')
# ax1.bar( (data.dt - width)*10, data._first, width*10, color = 'red')
# ax1.set_xticks(np.arange(0, 121, 30))

# histograms Stacked
width = 1.
print(data)
tot = data._switch + data._noswitch
f = np.sum( data.ix[data.absdt<=0.5,'_switch'] )
s = np.sum( data.ix[data.absdt<=0.5,'_noswitch'] )
print( f, s, f/(f+s), np.sum(tot) )

data._switch = data._switch / tot
data = data.fillna(0.)
data._noswitch = 1. - data._switch
print(data)

ax1.bar( (data.absdt-width/(60./binDT)), data._noswitch, width/(60./binDT), color = '#67a9cf',lw=.25)
ax1.bar( (data.absdt-width/(60./binDT)), data._switch, width/(60./binDT), color = '#ef8a62',  bottom=data._noswitch,lw=.25)

plt.show()


