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

ax1.set_ylabel( 'Normalized fluorescence' )
ax2.set_ylabel( 'Normalized fluorescence' )
ax2.set_xlabel( 'Time after division [hours]' )

'''
lag2YFP
'''
### plot single worms
# worms = ['W:\\Simone\\160930_lag2GFP+mCherry\\C09','W:\\Simone\\160930_lag2GFP+mCherry\\C10','W:\\Simone\\160930_lag2GFP+mCherry\\C11','W:\\Simone\\160930_lag2GFP+mCherry\\C15','W:\\Simone\\160930_lag2GFP+mCherry\\C19','W:\\Simone\\161002_lag2GFP+mCherry\\C01']
# worms = ['W:\\Simone\\160930_lag2GFP+mCherry\\C01','W:\\Simone\\160930_lag2GFP+mCherry\\C02','W:\\Simone\\160930_lag2GFP+mCherry\\C03','W:\\Simone\\160930_lag2GFP+mCherry\\C06','W:\\Simone\\160930_lag2GFP+mCherry\\C08','W:\\Simone\\161002_lag2GFP+mCherry\\C04']
### plot all worms
paths = [ 'W:\\Simone\\161030_lag2YFP+mCherry', 'W:\\Simone\\161108_lag2YFP+mCherry', 'W:\\Simone\\161111_lag2YFP+mCherry' ]
worms1 = []
for path in paths:
	worms1.append( glob.glob( os.path.join( path, '*cellFluo_q&d.pickle' ) ) )
worms = []
for exp in worms1:
	for w in exp:
		worms.append(w)
worms.sort()

'''
lag2GFP
'''
### plot single worms
# worms = ['W:\\Simone\\160930_lag2GFP+mCherry\\C09','W:\\Simone\\160930_lag2GFP+mCherry\\C10','W:\\Simone\\160930_lag2GFP+mCherry\\C11','W:\\Simone\\160930_lag2GFP+mCherry\\C15','W:\\Simone\\160930_lag2GFP+mCherry\\C19','W:\\Simone\\161002_lag2GFP+mCherry\\C01']
# worms = ['W:\\Simone\\160930_lag2GFP+mCherry\\C01','W:\\Simone\\160930_lag2GFP+mCherry\\C02','W:\\Simone\\160930_lag2GFP+mCherry\\C03','W:\\Simone\\160930_lag2GFP+mCherry\\C06','W:\\Simone\\160930_lag2GFP+mCherry\\C08','W:\\Simone\\161002_lag2GFP+mCherry\\C04']
### plot all worms
# paths = [ 'W:\\Simone\\160930_lag2GFP+mCherry', 'W:\\Simone\\161002_lag2GFP+mCherry', 'W:\\Simone\\161004_lag2GFP+mCherry' ]
# worms1 = []
# for path in paths:
# 	worms1.append( glob.glob( os.path.join( path, '*cellFluo_q&d.pickle' ) ) )
# worms = []
# for exp in worms1:
# 	for w in exp:
# 		worms.append(w)
# worms.sort()

'''
HLH2::GFP
'''
### plot single worms
# worms = [ 'X:\\Simone\\160407_HLH2_GFP_hist_mCherry\\C02' ]

'''
START SCRIPT
'''

print(worms)

dfMotherFirst = pd.DataFrame({})
dfMotherSecond = pd.DataFrame({})
dfCellFirst = pd.DataFrame({})
dfCellSecond = pd.DataFrame({})
dfSisFirst = pd.DataFrame({})
dfSisSecond = pd.DataFrame({})
for wormFull in worms:
	print(wormFull)

	worm = wormFull.split('\\')[-1].split('_')[0]
	path = wormFull[:-len(wormFull.split('\\')[-1])]

	timesDF = load_data_frame( path, worm + '_01times.pickle' )
	gonadPosDF = load_data_frame( path, worm + '_02gonadPos.pickle' )
	cellPosDF = load_data_frame( path, worm + '_04cellPos.pickle' )
	# cellOutDF = load_data_frame( path, worm + '_05cellOut.pickle' )
	cellFluoDF = load_data_frame( path, worm + '_15cellFluo_q&d.pickle' )

	# find time of first division and name of first dividing cell
	tdiv1 = np.min( timesDF.ix[ pd.notnull( cellPosDF[ '1.ppp' ].X ) ].timesRel )
	tdiv4 = np.min( timesDF.ix[ pd.notnull( cellPosDF[ '4.aaa' ].X ) ].timesRel )
	tpRel = np.min([tdiv1,tdiv4])
	tpRel1 = tdiv1
	tpRel4 = tdiv4
	if tdiv1 == tdiv4:
		firstDividing = 'same'
		secondDividing = 'same'
	elif tdiv1 < tdiv4:
		firstDividing = '1.pp'
		secondDividing = '4.aa'
	elif tdiv4 < tdiv1:
		firstDividing = '4.aa'
		secondDividing = '1.pp'

	### extract mother cells info
	# print(cellFluoDF['1.ppp'].keys())
	mother1 = cellFluoDF['1.pp'].ix[ pd.notnull( cellFluoDF['1.pp']._488nm ), ['_488nm','tidx','times'] ]
	mother1.times -= tpRel1

	mother4 = cellFluoDF['4.aa'].ix[ pd.notnull( cellFluoDF['4.aa']._488nm ), ['_488nm','tidx','times'] ]
	mother4.times -= tpRel4

	### extract cell info
	cell1 = cellFluoDF['1.ppp'].ix[ pd.notnull( cellFluoDF['1.ppp']._488nm ), ['_488nm','tidx','times'] ]
	cell1.times -= tpRel1

	cell4 = cellFluoDF['4.aaa'].ix[ pd.notnull( cellFluoDF['4.aaa']._488nm ), ['_488nm','tidx','times'] ]
	cell4.times -= tpRel4

	### extract cell info
	sis1 = cellFluoDF['1.ppa'].ix[ pd.notnull( cellFluoDF['1.ppa']._488nm ), ['_488nm','tidx','times'] ]
	sis1.times -= tpRel1

	sis4 = cellFluoDF['4.aap'].ix[ pd.notnull( cellFluoDF['4.aap']._488nm ), ['_488nm','tidx','times'] ]
	sis4.times -= tpRel4

	### renormalize mother data
	ren = np.mean( pd.concat([mother1.ix[mother1.times <= tpRel],mother4.ix[mother4.times <= tpRel]])._488nm )
	mother1._488nm /= ren
	mother4._488nm /= ren
	cell1._488nm /= ren
	cell4._488nm /= ren
	sis1._488nm /= ren
	sis4._488nm /= ren

	# print(firstDividing)
	# print(mother1)
	# print(mother4)
	# print(cell)

	if firstDividing == 'same':
		print('ciao')
		# dfMotherFirst = pd.concat([dfMotherFirst,mother1])
		# dfMotherFirst = pd.concat([dfMotherFirst,mother4])
		# dfCellFirst = pd.concat([dfCellFirst,cell1])
		# dfCellFirst = pd.concat([dfCellFirst,cell4])
		# ax1.plot( mother1.times, mother1._488nm, color = 'red', alpha = .5 )
		# ax1.plot( mother4.times, mother4._488nm, color = 'red', alpha = .5 )
		# ax1.plot( cell1.times, cell1._488nm, color = 'blue', alpha = .5 )
		# ax1.plot( cell4.times, cell4._488nm, color = 'blue', alpha = .5 )
	elif firstDividing == '1.pp':
		dfMotherFirst = pd.concat([dfMotherFirst,mother1])
		dfMotherSecond = pd.concat([dfMotherSecond,mother4])
		dfCellFirst = pd.concat([dfCellFirst,cell1])
		dfCellSecond = pd.concat([dfCellSecond,cell4])
		dfSisFirst = pd.concat([dfSisFirst,sis1])
		dfSisSecond = pd.concat([dfSisSecond,sis4])
		# ax1.plot( mother1.times, mother1._488nm, color = 'red', alpha = .5 )
		# ax2.plot( mother4.times, mother4._488nm, color = 'red', alpha = .5 )
		# ax1.plot( cell1.times, cell1._488nm, color = 'blue', alpha = .5 )
		# ax2.plot( cell4.times, cell4._488nm, color = 'blue', alpha = .5 )
		# ax1.plot( sis1.times, sis1._488nm, color = 'yellow', alpha = .5 )
		# ax2.plot( sis4.times, sis4._488nm, color = 'yellow', alpha = .5 )
	elif firstDividing == '4.aa':
		dfMotherFirst = pd.concat([dfMotherFirst,mother4])
		dfMotherSecond = pd.concat([dfMotherSecond,mother1])
		dfCellFirst = pd.concat([dfCellFirst,cell4])
		dfCellSecond = pd.concat([dfCellSecond,cell1])
		dfSisFirst = pd.concat([dfSisFirst,sis4])
		dfSisSecond = pd.concat([dfSisSecond,sis1])
		# ax1.plot( mother4.times, mother4._488nm, color = 'red', alpha = .5, lw=.5 )
		# ax2.plot( mother1.times, mother1._488nm, color = 'red', alpha = .5, lw=.5 )
		# ax1.plot( cell4.times, cell4._488nm, color = 'blue', alpha = .5 )
		# ax2.plot( cell1.times, cell1._488nm, color = 'blue', alpha = .5 )
		# ax1.plot( sis4.times, sis4._488nm, color = 'yellow', alpha = .5 )
		# ax2.plot( sis1.times, sis1._488nm, color = 'yellow', alpha = .5 )

# print(dfMotherFirst)
# print(dfMotherSecond)

colors = ['red','red','blue','blue','orange','orange']
ax = [ax1,ax2,ax1,ax2,ax1,ax2]
for idx, df in enumerate( [dfMotherFirst, dfMotherSecond, dfCellFirst, dfCellSecond, dfSisFirst, dfSisSecond] ):
	print(idx)

	# group data by time and remove missing datapoint
	bins = np.linspace(int(df.times.min()), int(df.times.max()+1), ( int(df.times.max()+1) - int(df.times.min()) ) * 2 + 1 )
	# print(bins)
	groups = df.groupby(np.digitize(df.times, bins))

	# Get the mean of each bin:
	_mean = groups.mean()
	# print('\n mean:\n',_mean)
	_std = groups.std()
	_std = _std._488nm.fillna(value=0)
	# print('\n std:\n',_std.values)

	# print(len(bins),len(_mean._488nm))
	ax[idx].errorbar( _mean.times, _mean._488nm, _std.values / 2, lw = 2, color = colors[idx] )

for ax in [ax1,ax2]:
	ax.plot([-5,10], [1,1], '--k')
	ax.plot([0,0], [0,4], '--k')
	# ax.set_xlim((-5,3))
	ax.set_ylim((0,4))
plt.show()
