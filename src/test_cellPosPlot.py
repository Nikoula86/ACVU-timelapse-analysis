# -*- coding: utf-8 -*-
"""
Created on Tue Jan 23 13:35:34 2018

@author: gritti
"""

# -*- coding: utf-8 -*-
"""
Created on Mon Jan 22 11:26:11 2018

@author: gritti
"""

import numpy as np
import matplotlib.pyplot as plt
import pickle
from matplotlib import cm
import generalFunctions as gf
import os
import pandas as pd
from scipy.ndimage import interpolation
import pandas as pd
import glob

#%%
def isDaughter(cname, cdata):
    return ( ( cdata.cname.str.startswith( cname ) ) & ( cdata.cname.str.len() == (len(cname)+1) ) ).any()

def coordTransformO2C(r,size):
    return np.array( [ size/2+r[0], size/2-r[1] ] )

def coordTransformC2O(r,size):
    return np.array( [ r[0]-size/2, size/2-r[1] ] )

def rotateCoords(r, alpha):
    return np.array( [ r[0]*np.cos(alpha) - r[1]*np.sin(alpha), r[0]*np.sin(alpha) + r[1]*np.cos(alpha) ] )

def Yreflect(r):
	return np.array( [ r[0], -r[1] ] )

def inverse(pos1,pos2,posn):
	sign = (pos2[0]-pos1[0])*(posn[1]-pos1[1])-(pos2[1]-pos1[0])*(posn[0]-pos1[0])
	return np.sign(sign) == 1

def extractCoords(cell,channel):
	# print(cell)
	# return np.array([cell.X,cell.Y])
	if channel == '488nm':
		if not 'X488nm' in cell.keys():
			cell.X488nm = cell.X
			cell.Y488nm = cell.Y
		if np.isnan( cell.X488nm ):
			cell.X488nm = cell.X
			cell.Y488nm = cell.Y
		return np.array([cell.X488nm,cell.Y488nm])
	if channel == '561nm':
		if not 'X561nm' in cell.keys():
			cell.X561nm = cell.X
			cell.Y561nm = cell.Y
		if np.isnan( cell.X561nm ):
			cell.X561nm = cell.X
			cell.Y561nm = cell.Y
		return np.array([cell.X561nm,cell.Y561nm])

def extractCoordsRot(cell):
	return np.array([cell.cXposRot,cell.cYposRot])

def extractCoordsO(celltp,cname):
	return np.array([celltp.ix[celltp.cname==cname,'cXposO'].values[0],celltp.ix[celltp.cname==cname,'cYposO'].values[0]])

def extractRefPointsCoordsO(postp,pname):
	return np.array([postp.ix[postp.pname==pname,'cXposO'].values[0],postp.ix[postp.pname==pname,'cYposO'].values[0]])

def extractCoordsORot(celltp,cname):
	return np.array([celltp.ix[celltp.cname==cname,'cXposORot'].values[0],celltp.ix[celltp.cname==cname,'cYposORot'].values[0]])

def extractRefPointsCoordsORot(postp,pname):
	return np.array([postp.ix[postp.pname==pname,'cXposORot'].values[0],postp.ix[postp.pname==pname,'cYposORot'].values[0]])

def flatFieldCorrection( imgs, darkField, flatField, bodyData, tidx ):
	size = imgs[0].shape[0]
	medianCorrection = np.median( ( flatField - darkField ) )

	gp = np.floor( bodyData.ix[ bodyData.tidx == tidx, 'gonadPos' ].values[0] )

	darkF = np.zeros( ( size, size ) )
	flatF = np.zeros( ( size, size ) )

	darkF[ -np.min( [ gp[1]-size/2, 0 ] ) : size-np.max( [ gp[1]+size/2-2047, 0 ] ) , 
							-np.min( [ gp[0]-size/2, 0 ] ) : size-np.max( [ gp[0]+size/2-2047, 0 ] ) ] = darkField[
							np.max( [ gp[1]-size/2, 0 ] ) : np.min( [ gp[1]+size/2, 2047 ] ) , 
							np.max( [ gp[0]-size/2, 0 ] ) : np.min( [ gp[0]+size/2, 2047 ] ) ]
	flatF[ -np.min( [ gp[1]-size/2, 0 ] ) : size-np.max( [ gp[1]+size/2-2047, 0 ] ) , 
							-np.min( [ gp[0]-size/2, 0 ] ) : size-np.max( [ gp[0]+size/2-2047, 0 ] ) ] = flatField[
							np.max( [ gp[1]-size/2, 0 ] ) : np.min( [ gp[1]+size/2, 2047 ] ) , 
							np.max( [ gp[0]-size/2, 0 ] ) : np.min( [ gp[0]+size/2, 2047 ] ) ]

	# print(flatF,medianCorrection)

	# plt.figure()
	# plt.imshow(darkF,cmap='gray')
	# plt.figure()
	# plt.imshow(flatF,cmap='gray')

	return ((imgs - darkF)/((flatF-darkF)/medianCorrection)).astype(np.uint16)

def rotateImage(imgs, alpha):
	return np.array( [ interpolation.rotate(img,alpha*180/np.pi,reshape=False) for img in imgs ] )

def flipImage(refPoints):
	if refPoints.ix[refPoints.pname=='d','cYposRot'].values[0] > refPoints.ix[refPoints.pname=='a','cYposRot'].values[0]:
		return True
	else:
		return False

def find_interesting_cells( currentCells ):
	cnames = list( currentCells.cname )

	if '1.p' in cnames:
		cell1ppa = currentCells.ix[ currentCells.cname == '1.p' ]
		cell1ppp = currentCells.ix[ currentCells.cname == '1.p' ]
		if '4.a' in cnames:
			cell4aaa = currentCells.ix[ currentCells.cname == '4.a' ]
			cell4aap = currentCells.ix[ currentCells.cname == '4.a' ]
		if '4.aa' in cnames:
			cell4aaa = currentCells.ix[ currentCells.cname == '4.aa' ]
			cell4aap = currentCells.ix[ currentCells.cname == '4.aa' ]

	if '1.pp' in cnames:
		cell1ppa = currentCells.ix[ currentCells.cname == '1.pp' ]
		cell1ppp = currentCells.ix[ currentCells.cname == '1.pp' ]
		if '4.a' in cnames:
			cell4aaa = currentCells.ix[ currentCells.cname == '4.a' ]
			cell4aap = currentCells.ix[ currentCells.cname == '4.a' ]
		if '4.aa' in cnames:
			cell4aaa = currentCells.ix[ currentCells.cname == '4.aa' ]
			cell4aap = currentCells.ix[ currentCells.cname == '4.aa' ]
		if '4.aaa' in cnames:
			cell4aaa = currentCells.ix[ currentCells.cname == '4.aaa' ]
			if '4.aap' in cnames:
				cell4aap = currentCells.ix[ currentCells.cname == '4.aap' ]
			else:
				cell4aap = currentCells.ix[ currentCells.cname == '4.aaa' ]
				cell4aap.X = np.nan

	if '1.ppp' in cnames:
		cell1ppp = currentCells.ix[ currentCells.cname == '1.ppp' ]
		if '1.ppa' in cnames:
			cell1ppa = currentCells.ix[ currentCells.cname == '1.ppa' ]
		else:
			cell1ppa = currentCells.ix[ currentCells.cname == '1.ppp' ]
			cell1ppa.X = np.nan
			
		if '4.aa' in cnames:
			cell4aaa = currentCells.ix[ currentCells.cname == '4.aa' ]
			cell4aap = currentCells.ix[ currentCells.cname == '4.aa' ]
		if '4.aaa' in cnames:
			cell4aaa = currentCells.ix[ currentCells.cname == '4.aaa' ]
			if '4.aap' in cnames:
				cell4aap = currentCells.ix[ currentCells.cname == '4.aap' ]
			else:
				cell4aap = currentCells.ix[ currentCells.cname == '4.aaa' ]
				cell4aap.X = np.nan

	cell1ppa.cname = '1.ppa'
	cell1ppa = cell1ppa.squeeze()
	cell1ppp.cname = '1.ppp'
	cell1ppp = cell1ppp.squeeze()
	cell4aaa.cname = '4.aaa'
	cell4aaa = cell4aaa.squeeze()
	cell4aap.cname = '4.aap'
	cell4aap = cell4aap.squeeze()
	# print(cell1,cell2)
	return ( cell1ppa, cell1ppp, cell4aaa, cell4aap )

#%%
import matplotlib as mpl
mpl.rcParams['pdf.fonttype'] = 42

cwd = os.getcwd()
parentDir = os.path.join(os.path.dirname(os.path.dirname(cwd)),'ACVU_data','timelapse')
#print(worms)

#%%
colors = ['b','k','r','m']

fig = plt.figure(figsize=(4,7))
ax = fig.add_subplot(111)


allData = pickle.load(open(os.path.join(parentDir,'allData_cellPos.pickle'),'rb'))
worms = allData[0]

for idx, wormPath in enumerate( worms ):
	print(wormPath)
	
	dataList = allData[1][idx]
	for i in np.arange(4):
		ax.plot(dataList[i+1][np.isfinite(dataList[i+1])],dataList[0][np.isfinite(dataList[i+1])],color = colors[i],alpha = .25)
	


ax.xaxis.tick_top()
ax.set_ylim(-1, 10)
plt.gca().invert_yaxis()
ax.set_xlim(-25, 25)	
ax.plot([-25,25],[0,0],'--k',lw=.5)
ax.plot([0,0],[-5,10],'--k',lw=.5)

#%% compute mean and std
colors = ['b','k','r','m']

fig = plt.figure(figsize=(4,7))
ax = fig.add_subplot(111)

newdata = [ 
				np.linspace(0,10,61),
			   [[] for i in np.linspace(0,10,61)],
			   [[] for i in np.linspace(0,10,61)],
			   [[] for i in np.linspace(0,10,61)],
			   [[] for i in np.linspace(0,10,61)]
		   ]

for idx,wormPath in enumerate(worms):
	
	dataList = allData[1][idx]
	dataListPositive = [ i[dataList[0]>=0] for i in dataList ]
	newjdxs = np.digitize(dataListPositive[0],newdata[0]) - 1
	
	for jdx, newjdx in enumerate(newjdxs):
		for i in [1,2,3,4]:
			if np.isfinite( dataListPositive[i][jdx] ):
				newdata[i][newjdx].append(dataListPositive[i][jdx])


times = newdata[0]
meanPos = [ [ np.mean(i) for i in j ] for j in newdata[1:]  ]
stdPos = [ [ np.std(i) for i in j ] for j in newdata[1:]  ]

for i in np.arange(4):
	meany = np.array(meanPos[i])
	stdy = np.array(stdPos[i])
	ax.plot(meany,times,color = colors[i],alpha = .25)
	ax.fill_betweenx(times,meany-stdy,meany+stdy,facecolor=colors[i],alpha=.2)

ax.xaxis.tick_top()
ax.set_ylim(-1, 10)
plt.gca().invert_yaxis()
ax.set_xlim(-25, 25)	
ax.plot([-25,25],[0,0],'--k',lw=.5)
ax.plot([0,0],[-5,10],'--k',lw=.5)

plt.show()