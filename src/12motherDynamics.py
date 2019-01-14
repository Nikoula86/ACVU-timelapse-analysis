import numpy as np
import matplotlib.pyplot as plt
import pickle
from matplotlib import cm
from generalFunctions import *
from skimage import filters


import matplotlib as mpl
mpl.rcParams['pdf.fonttype'] = 42

fig = plt.figure(figsize=(5.8,11.6))
ax1 = fig.add_subplot(211)
ax2 = fig.add_subplot(212)
fig.subplots_adjust(left=0.15, right=.95, top=.95, bottom=0.15)
for ax in [ax1,ax2]:
	for tl in ax.get_xticklabels():
		tl.set_fontsize(18)
	for tl in ax.get_yticklabels():
		tl.set_fontsize(18)

'''
lag2YFP
'''
### plot single worms
# worms = ['W:\\Simone\\160930_lag2GFP+mCherry\\C09','W:\\Simone\\160930_lag2GFP+mCherry\\C10','W:\\Simone\\160930_lag2GFP+mCherry\\C11','W:\\Simone\\160930_lag2GFP+mCherry\\C15','W:\\Simone\\160930_lag2GFP+mCherry\\C19','W:\\Simone\\161002_lag2GFP+mCherry\\C01']
# worms = ['W:\\Simone\\160930_lag2GFP+mCherry\\C01','W:\\Simone\\160930_lag2GFP+mCherry\\C02','W:\\Simone\\160930_lag2GFP+mCherry\\C03','W:\\Simone\\160930_lag2GFP+mCherry\\C06','W:\\Simone\\160930_lag2GFP+mCherry\\C08','W:\\Simone\\161002_lag2GFP+mCherry\\C04']
### plot all worms
# paths = [ 'W:\\Simone\\161030_lag2YFP+mCherry' ]
# worms1 = []
# for path in paths:
# 	worms1.append( glob.glob( os.path.join( path, '*cellFluo_q&d.pickle' ) ) )
# worms = []
# for exp in worms1:
# 	for w in exp:
# 		worms.append(w.split('_06')[0])

'''
lag2GFP
'''
### plot single worms
# worms = ['W:\\Simone\\160930_lag2GFP+mCherry\\C09','W:\\Simone\\160930_lag2GFP+mCherry\\C10','W:\\Simone\\160930_lag2GFP+mCherry\\C11','W:\\Simone\\160930_lag2GFP+mCherry\\C15','W:\\Simone\\160930_lag2GFP+mCherry\\C19','W:\\Simone\\161002_lag2GFP+mCherry\\C01']
# worms = ['W:\\Simone\\160930_lag2GFP+mCherry\\C01','W:\\Simone\\160930_lag2GFP+mCherry\\C02','W:\\Simone\\160930_lag2GFP+mCherry\\C03','W:\\Simone\\160930_lag2GFP+mCherry\\C06','W:\\Simone\\160930_lag2GFP+mCherry\\C08','W:\\Simone\\161002_lag2GFP+mCherry\\C04']
### plot all worms
paths = [ 'W:\\Simone\\160930_lag2GFP+mCherry', 'W:\\Simone\\161002_lag2GFP+mCherry' ]
worms1 = []
for path in paths:
	worms1.append( glob.glob( os.path.join( path, '*cellFluo.pickle' ) ) )
worms = []
for exp in worms1:
	for w in exp:
		worms.append(w.split('_06')[0])

'''
HLH2::GFP
'''
### plot single worms
worms = [ 'X:\\Simone\\160407_HLH2_GFP_hist_mCherry\\C02' ]


print(worms)

dfFirst = pd.DataFrame({})
dfSecond = pd.DataFrame({})
dfCell = pd.DataFrame({})
for wormFull in worms:
	print(wormFull)

	worm = wormFull.split('\\')[-1].split('_')[0]
	path = wormFull[:-len(wormFull.split('\\')[-1])]

	timesDF = load_data_frame( path, worm + '_01times.pickle' )
	gonadPosDF = load_data_frame( path, worm + '_02gonadPos.pickle' )
	cellPosDF = load_data_frame( path, worm + '_04cellPos.pickle' )
	# cellOutDF = load_data_frame( path, worm + '_05cellOut.pickle' )
	cellFluoDF = load_data_frame( path, worm + '_06cellFluo.pickle' )

	# find time of first division and name of first dividing cell
	tdiv1 = np.min( timesDF.ix[ pd.notnull( cellPosDF[ '1.ppp' ].X ) ].timesRel )
	tdiv2 = np.min( timesDF.ix[ pd.notnull( cellPosDF[ '4.aaa' ].X ) ].timesRel )
	tpRel = np.min([tdiv1,tdiv2])
	if tdiv1 == tdiv2:
		firstDividing = 'same'
		secondDividing = 'same'
	elif tpRel == tdiv1:
		firstDividing = '1.pp'
		secondDividing = '4.aa'
	elif tpRel == tdiv2:
		firstDividing = '4.aa'
		secondDividing = '1.pp'

	# find tidx of first division
	tDivIdx = timesDF.ix[ timesDF.timesRel == tpRel, 'tidxRel' ].values[0]

	# find tidx of second division
	tDiv2 = np.max([tdiv1,tdiv2])
	tDiv2Idx = timesDF.ix[ timesDF.timesRel == tDiv2, 'tidxRel' ].values[0]

	# extract mother cells info
	# print(cellFluoDF['1.ppp'].keys())
	mother1 = cellFluoDF['1.pp'].ix[ pd.notnull( cellFluoDF['1.pp']._488nm ), ['_488nm','tidx'] ]
	mother1._488nm -= cellFluoDF[ 'b_1' ].ix[ pd.notnull( cellFluoDF['1.pp']._488nm ) ]._488nm
	mother1._488nm /= cellFluoDF[ 'b_1' ].ix[ pd.notnull( cellFluoDF['1.pp']._488nm ) ]._488nm
	mother1['_488nmAbs'] = mother1._488nm
	mother1['times'] = timesDF.ix[ pd.notnull( cellFluoDF['1.pp']._488nm ), 'timesRel' ].values - tpRel

	mother4 = cellFluoDF['4.aa'].ix[ pd.notnull( cellFluoDF['4.aa']._488nm ), ['_488nm','tidx'] ]
	mother4._488nm -= cellFluoDF[ 'b_4' ].ix[ pd.notnull( cellFluoDF['4.aa']._488nm ) ]._488nm
	mother4._488nm /= cellFluoDF[ 'b_4' ].ix[ pd.notnull( cellFluoDF['4.aa']._488nm ) ]._488nm
	mother4['_488nmAbs'] = mother4._488nm
	mother4['times'] = timesDF.ix[ pd.notnull( cellFluoDF['4.aa']._488nm ), 'timesRel' ].values - tpRel

	# extract cell info
	if firstDividing == 'same':
		cell = pd.DataFrame({})
	elif firstDividing == '1.pp':
		mother1._488nm /= np.mean(mother4.ix[ mother4.times<0, '_488nmAbs'])
		mother4._488nm /= np.mean(mother4.ix[ mother4.times<0, '_488nmAbs'])
		cell = cellFluoDF['1.ppp'].ix[ pd.notnull( cellFluoDF['1.ppp']._488nm ), ['_488nm','tidx'] ]
		cell._488nm -= cellFluoDF[ 'b_1' ].ix[ pd.notnull( cellFluoDF['1.ppp']._488nm ) ]._488nm
		cell._488nm /= cellFluoDF[ 'b_1' ].ix[ pd.notnull( cellFluoDF['1.ppp']._488nm ) ]._488nm
		cell['times'] = timesDF.ix[ pd.notnull( cellFluoDF['1.ppp']._488nm ), 'timesRel' ].values - tpRel
		cell._488nm /= np.mean(mother4.ix[ mother4.times<0, '_488nmAbs'])
		cell = cell.ix[  (cell.tidx < tDiv2Idx) ]
	elif firstDividing == '4.aa':
		mother1._488nm /= np.mean(mother1.ix[ mother1.times<0, '_488nmAbs'])
		mother4._488nm /= np.mean(mother1.ix[ mother1.times<0, '_488nmAbs'])
		cell = cellFluoDF['4.aaa'].ix[ pd.notnull( cellFluoDF['4.aaa']._488nm ), ['_488nm','tidx'] ]
		cell._488nm -= cellFluoDF[ 'b_4' ].ix[ pd.notnull( cellFluoDF['4.aaa']._488nm ) ]._488nm
		cell._488nm /= cellFluoDF[ 'b_4' ].ix[ pd.notnull( cellFluoDF['4.aaa']._488nm ) ]._488nm
		cell['times'] = timesDF.ix[ pd.notnull( cellFluoDF['4.aaa']._488nm ), 'timesRel' ].values - tpRel
		cell._488nm /= np.mean(mother1.ix[ mother1.times<0, '_488nmAbs'])
		cell = cell.ix[  (cell.tidx < tDiv2Idx) ]

	# print(firstDividing)
	# print(mother1)
	# print(mother4)
	# print(cell)

	if firstDividing == 'same':
		dfFirst = pd.concat([dfFirst,mother1])
		dfFirst = pd.concat([dfFirst,mother4])
		ax1.plot( mother1.times, mother1._488nm, color = 'red', alpha = .5 )
		ax1.plot( mother4.times, mother4._488nm, color = 'red', alpha = .5 )
	elif firstDividing == '1.pp':
		dfFirst = pd.concat([dfFirst,mother1])
		dfSecond = pd.concat([dfSecond,mother4])
		dfCell = pd.concat([dfCell,cell])
		ax1.plot( mother1.times, mother1._488nm, color = 'red', alpha = .5 )
		ax2.plot( mother4.times, mother4._488nm, color = 'blue', alpha = .5 )
		ax1.plot( cell.times, cell._488nm, color = 'green', alpha = .5 )
	elif firstDividing == '4.aa':
		dfFirst = pd.concat([dfFirst,mother4])
		dfSecond = pd.concat([dfSecond,mother1])
		dfCell = pd.concat([dfCell,cell])
		ax1.plot( mother4.times, mother4._488nm, color = 'red', alpha = .5, lw=.5 )
		ax2.plot( mother1.times, mother1._488nm, color = 'blue', alpha = .5, lw=.5 )
		ax1.plot( cell.times, cell._488nm, color = 'green', alpha = .5 )

print(dfFirst)
print(dfSecond)

colors = ['red','blue','green']
ax = [ax1,ax2,ax1]
for idx, df in enumerate( [dfFirst, dfSecond, dfCell] ):
	# group data by time and remove missing datapoint
	bins = np.linspace(df.times.min(), df.times.max(), 20)
	print(bins)
	groups = df.groupby(np.digitize(df.times, bins))

	# Get the mean of each bin:
	_mean = groups.mean()
	print(_mean)
	_std = groups.std()
	_std = _std._488nm.fillna(value=0)
	print(_std)

	ax[idx].plot(_mean.times,_mean._488nm, lw = 2, color = colors[idx])

for ax in [ax1,ax2]:
	ax.plot([-5,3], [1,1], '--k')
	ax.plot([0,0], [0,4], '--k')
	ax.set_xlim((-5,3))
	ax.set_ylim((0,4))
plt.show()
