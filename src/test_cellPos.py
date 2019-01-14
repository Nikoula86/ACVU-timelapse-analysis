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
expFolders = [ 
	 		 '161030_lag2YFP+histmCherry',
			 '161108_lag2YFP+histmCherry',
			 '161111_lag2YFP+histmCherry',
			 '170312_lag2YFP+histmCherry',
			 ]
worms1 = []
for path in expFolders:
	worms1.append( glob.glob( os.path.join( parentDir, path, '*08apdvPos.pickle' ) ) )
worms = []
for exp in worms1:
	for w in exp:
		worms.append(w)
worms.sort()
#print(worms)

colors = ['b','k','r','m']
channel = '488nm'
imgSize = 512

fig = plt.figure(figsize=(4,7))
ax = fig.add_subplot(111)

allData = [[],[]]

#%%
for wormPath in worms:
	print(wormPath)
	
	worm = wormPath.split('\\')[-1].split('_')[0]
	path = wormPath.split(worm)[0]
	
	paramsDF = gf.load_data_frame( os.path.join( path ), worm + '_01params.pickle' )
	timesDF = gf.load_data_frame( os.path.join( path ), worm + '_01times.pickle' )
	gpDF = gf.load_data_frame( os.path.join( path ), worm + '_02gonadPos.pickle' )
	cellPosDF = gf.load_data_frame( os.path.join( path ), worm + '_04cellPos.pickle' )
	apdvPosDF = gf.load_data_frame( os.path.join( path ), worm + '_08apdvPos.pickle' )
	
	#%%	
	### extract division times
	tdiv1 = np.min( timesDF.ix[ pd.notnull( cellPosDF[ '1.ppp' ].X ) ].timesRel )
	tdiv2 = np.min( timesDF.ix[ pd.notnull( cellPosDF[ '4.aaa' ].X ) ].timesRel )
	
	timesDF.timesRel -= np.max([tdiv1,tdiv2])
	
	#%%

	dataList = [[],[],[],[],[]]
	
	for idx, trow in timesDF.iterrows():
	
		### skip problematic timepoints
		# if trow.tidxRel == 42:
		# 	continue
	
		currentCells = gf.extract_current_cell_pos( cellPosDF, trow.tidxRel )
		currentCellsOriginalPos = gf.extract_current_cell_pos( cellPosDF, trow.tidxRel )
	
		### if there are only 2 cells (one is the bckg) then skip the timepoint
		if len( currentCells.cname ) > 2:
			dataList[0].append(trow.timesRel)
	
			### find out which cells we want to show
			(c1ppa,c1ppp,c4aaa,c4aap) = find_interesting_cells( currentCells )
			interestingCells = pd.concat([c1ppa, c1ppp, c4aaa, c4aap],axis=1).transpose().reset_index()
			print(interestingCells)
	
			refPoints = gf.extract_current_apdv_pos( apdvPosDF, trow.tidxRel )
#			print(trow.tidxRel,list(interestingCells.cname))
	
			### correct for XY
			gonadPos = gf.extract_pos( gpDF.ix[ gpDF.tidx == trow.tidxRel ].squeeze() )
	
	
			### get coordinates with respect to center of the image (origin)
			for jdx, cell in interestingCells.iterrows():
	
				posC = extractCoords(cell,channel)
				posO = coordTransformC2O(posC, imgSize)
	
				interestingCells.ix[ interestingCells.cname == cell.cname, 'cXposO' ] = posO[0]
				interestingCells.ix[ interestingCells.cname == cell.cname, 'cYposO' ] = posO[1]
	
			for jdx, pos in refPoints.iterrows():
				# print(pos)
				posC = gf.extract_pos(pos) - gonadPos + 256
				posO = coordTransformC2O(posC, imgSize)
	
				refPoints.ix[ refPoints.pname == pos.pname, 'cXposO' ] = posO[0]
				refPoints.ix[ refPoints.pname == pos.pname, 'cYposO' ] = posO[1]
	
			vect = extractRefPointsCoordsO(refPoints,'p') - extractRefPointsCoordsO(refPoints,'a')
			alpha = -np.arctan2(vect[1],vect[0])
	
			### rotate coordinates
			for jdx, cell in interestingCells.iterrows():
	
				posO = extractCoordsO(interestingCells, cell.cname)
				posORot = rotateCoords(posO,alpha)
	
				interestingCells.ix[interestingCells.cname==cell.cname,'cXposORot'] = posORot[0]
				interestingCells.ix[interestingCells.cname==cell.cname,'cYposORot'] = posORot[1]
	
				interestingCells.ix[interestingCells.cname==cell.cname,'cXposRot'] = coordTransformO2C( extractCoordsORot(interestingCells,cell.cname), imgSize )[0]
				interestingCells.ix[interestingCells.cname==cell.cname,'cYposRot'] = coordTransformO2C( extractCoordsORot(interestingCells,cell.cname), imgSize )[1]
	
			for jdx, pos in refPoints.iterrows():
	
				posO = extractRefPointsCoordsO(refPoints,pos.pname)
				posORot = rotateCoords(posO,alpha)
	
				refPoints.ix[refPoints.pname==pos.pname,'cXposORot'] = posORot[0]
				refPoints.ix[refPoints.pname==pos.pname,'cYposORot'] = posORot[1]
	
				refPoints.ix[refPoints.pname==pos.pname,'cXposRot'] = coordTransformO2C( extractRefPointsCoordsORot(refPoints,pos.pname), imgSize )[0]
				refPoints.ix[refPoints.pname==pos.pname,'cYposRot'] = coordTransformO2C( extractRefPointsCoordsORot(refPoints,pos.pname), imgSize )[1]
	
				### flip image
				flipTF = flipImage( refPoints )
	
				### flip coordinates
				if flipTF:
					# print('flipping...')
					for jdx, cell in interestingCells.iterrows():
						posORot = extractCoordsORot(interestingCells, cell.cname)
	
						interestingCells.ix[interestingCells.cname==cell.cname,'cXposORot'] = posORot[0]
						interestingCells.ix[interestingCells.cname==cell.cname,'cYposORot'] = -posORot[1]
	
						interestingCells.ix[interestingCells.cname==cell.cname,'cXposRot'] = coordTransformO2C( extractCoordsORot(interestingCells,cell.cname), imgSize )[0]
						interestingCells.ix[interestingCells.cname==cell.cname,'cYposRot'] = coordTransformO2C( extractCoordsORot(interestingCells,cell.cname), imgSize )[1]
						
#			print(interestingCells)
			dataList[1].append(interestingCells.ix[interestingCells.cname=='1.ppa','cXposRot'].values[0])
			dataList[2].append(interestingCells.ix[interestingCells.cname=='1.ppp','cXposRot'].values[0])
			dataList[3].append(interestingCells.ix[interestingCells.cname=='4.aaa','cXposRot'].values[0])
			dataList[4].append(interestingCells.ix[interestingCells.cname=='4.aap','cXposRot'].values[0])
	
	dataList = np.array(dataList)
	centerMass = np.nanmean( dataList[1:], 0 )
	
	dataList[1:] -= centerMass
	dataList[1:] *= 0.108
		
	for i in np.arange(4):
		ax.plot(dataList[i+1][np.isfinite(dataList[i+1])],dataList[0][np.isfinite(dataList[i+1])],color = colors[i])
	allData[0].append(wormPath)
	allData[1].append(dataList)
ax.xaxis.tick_top()
ax.set_ylim(-3, 12)
plt.gca().invert_yaxis()
ax.set_xlim(-25, 25)	
plt.show()