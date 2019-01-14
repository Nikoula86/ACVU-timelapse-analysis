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

pathsAll = []
strains = []
'''
lag2GFP strain
'''
# pathsAll.append( [ 'W:\\Simone\\160930_lag2GFP+mCherry', 'W:\\Simone\\161002_lag2GFP+mCherry', 'W:\\Simone\\161004_lag2GFP+mCherry' ] )
# strains.append('lag2GFP')
'''
lag2YFP strain
'''
pathsAll.append( [ 'W:\\Simone\\161030_lag2YFP+histmCherry', 'W:\\Simone\\161108_lag2YFP+histmCherry', 'W:\\Simone\\161111_lag2YFP+histmCherry', 'W:\\Simone\\170312_lag2YFP+histmCherry' ] )
strains.append('lag2YFP')
'''
LIN12GFP strain
'''
pathsAll.append( [ 'W:\\Simone\\160420_LIN12GFP+histmCherry', 'W:\\Simone\\161212_LIN12GFP+histmCherry', 'W:\\Simone\\170220_LIN12GFP+histmCherry', 'W:\\Simone\\170221_LIN12GFP+histmCherry' ] )
strains.append('LIN12GFP')
'''
HLH2GFP strain
'''
pathsAll.append( [ 'W:\\Simone\\160407_HLH2GFP+histmCherry', 'W:\\Simone\\160502_HLH2GFP+histmCherry', 'W:\\Simone\\170211_HLH2GFP+histmCherry' ] )
strains.append('HLH2GFP')
'''
START
'''
allData = []
for idx, paths in enumerate( pathsAll ):
	worms1 = []
	for path in paths:
		worms1.append( glob.glob( os.path.join( path, '*01times.pickle' ) ) )
	worms = []
	for exp in worms1:
		for w in exp:
			worms.append(w)
	worms.sort()
	# print(worms)
	print(strains[idx]+'\nN:',len(worms),'\n')

	binDT = 10 #binning in minutes

	data = pd.DataFrame({'absdt':np.linspace(0,12,13),
						'_second': 0*np.linspace(0,12,13),
						'_first': 0*np.linspace(0,12,13)})
	# print(data)
	for widx, worm in enumerate(worms):
		# print(worm)

		### load parameters and times dataframes
		w = worm.split('\\')[-1].split('_')[0]
		path = worm[:-len(worm.split('\\')[-1])]
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
		cdata = np.array([[int(i) for i in x.strip().split('\t')] for x in cdata[1:]])
		cdata = cdata[cdata[:,0]==int(w[1:])][0,1:]
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

	tMax = 2.
	data.absdt /= ( 60. / 10 )
	data = data.ix[ (data.absdt<=tMax) & (data.absdt>0) ]

	# histograms Stacked
	ftot = np.sum( data.ix[:,'_first'] )
	stot = np.sum( data.ix[:,'_second'] )
	f = np.sum( data.ix[data.absdt<=0.5,'_first'] )
	s = np.sum( data.ix[data.absdt<=0.5,'_second'] )


	allData.append( [ftot+stot, ftot, 100*ftot/(ftot+stot), f+s, f, 100*f/(f+s)] )

allData = np.array(allData)


'''
PLOT
'''
fig = plt.figure(figsize=(5,5))
ax1 = fig.add_subplot(111)
fig.subplots_adjust(left=0.20, right=.93, top=.93, bottom=0.20)
for tl in ax1.get_xticklabels():
	tl.set_fontsize(15)
for tl in ax1.get_yticklabels():
	tl.set_fontsize(15)


width = 0.5

ax1.set_xticks(np.arange(len(allData))+1)
ax1.set_xticklabels(strains,rotation=45)
ax1.set_ylabel('Fraction of animals',fontsize = 15)

### histogram of <30 minutes
# ax1.bar(np.arange(len(allData)) + 1 - width/2., allData[:,5], width,color = '#ef8a62',edgecolor = '#ef8a62')
# for i in np.arange(len(allData)):
# 	ax1.text(i+1-width/3.,allData[i,5]+0.5, '%d/%d'%(int(allData[i,4]),int(allData[i,3])),fontsize = 15)
	
### histogram of all worms
ax1.bar(np.arange(len(allData)) + 1 - width/2., allData[:,2], width,color = '#ef8a62',edgecolor = '#ef8a62')
for i in np.arange(len(allData)):
	ax1.text(i+1-width/3.,allData[i,2]+0.5, '%d/%d'%(int(allData[i,1]),int(allData[i,0])),fontsize = 15)

plt.show()
