# -*- coding: utf-8 -*-
"""
Created on Fri Apr 10 17:00:56 2015

@author: gritti


NB: CHANGE outpath ACCORDING TO THE COMPUTER YOU ARE USING!
"""

import glob
import pandas as pd
from tifffile import *
import generalFunctions as gf
import numpy as np
import PIL
from PIL import Image, ImageDraw, ImageFont
import os.path
import matplotlib.pyplot as plt
import pickle
import scipy.interpolate as ip
from scipy.ndimage import interpolation
import matplotlib as mpl
mpl.rcParams['pdf.fonttype'] = 42

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

	darkF[ -int(np.min( [ gp[1]-size/2, 0 ] )) : int(size-np.max( [ gp[1]+size/2-2047, 0 ] )) , 
							-int(np.min( [ gp[0]-size/2, 0 ] )) : int(size-np.max( [ gp[0]+size/2-2047, 0 ] )) ] = darkField[
							int(np.max( [ gp[1]-size/2, 0 ] )) : int(np.min( [ gp[1]+size/2, 2047 ] )) , 
							int(np.max( [ gp[0]-size/2, 0 ] )) : int(np.min( [ gp[0]+size/2, 2047 ] )) ]
	flatF[ -int(np.min( [ gp[1]-size/2, 0 ] )) : int(size-np.max( [ gp[1]+size/2-2047, 0 ] )) , 
							-int(np.min( [ gp[0]-size/2, 0 ] )) : int(size-np.max( [ gp[0]+size/2-2047, 0 ] )) ] = flatField[
							int(np.max( [ gp[1]-size/2, 0 ] )) : int(np.min( [ gp[1]+size/2, 2047 ] )) , 
							int(np.max( [ gp[0]-size/2, 0 ] )) : int(np.min( [ gp[0]+size/2, 2047 ] )) ]

	# print(flatF,medianCorrection)

	# plt.figure()
	# plt.imshow(darkF,cmap='gray')
	# plt.figure()
	# plt.imshow(flatF,cmap='gray')

	return ((imgs - darkF)/((flatF-darkF)/medianCorrection)).astype(np.uint16)

def rotateImage(imgs, alpha):
	return np.array( [ interpolation.rotate(img,alpha*180/np.pi,reshape=False) for img in imgs ] )

def flipImage(imgs,refPoints):
	if refPoints.ix[refPoints.pname=='d','cYposRot'].values[0] > refPoints.ix[refPoints.pname=='a','cYposRot'].values[0]:
		return np.array( [ img[::-1, :] for img in imgsRot ] ), True
	else:
		return imgs, False

'''
lag2YFP data
'''
#worms = [ 'Y:\\Simone\\161030_lag2YFP+histmCherry\\C23' ]
'''
lin12nullunc data
'''
#worms = [ 'Y:\\Simone\\170907_LIN12unc+histmCherry\\C14' ]
'''
lag2 integration line data
'''
#worms = [ 'Y:\\Simone\\171028_lag2int+JVZ32\\C08' ]
'''
lag2 ts line data
'''
absPath = '/Users/ngritti/Dropbox/ACVU_project_private/ACVU_data/timelapse/pickledFiles/'
worms = [ absPath + '180524_lag2ts+JVZ32/C07' ]

### plot all available worms
# paths = [ 'W:\\Simone\\161030_lag2YFP+mCherry', 'W:\\Simone\\161108_lag2YFP+mCherry', 'W:\\Simone\\161111_lag2YFP+mCherry' ]
# worms1 = []
# for path in paths:
# 	worms1.append( glob.glob( os.path.join( path, '*cellFluo_q&d.pickle' ) ) )
# worms = []
# for exp in worms1:
# 	for w in exp:
# 		worms.append(w)
# worms.sort()

channels = ['488nm','561nm']

size=40

### INCLUDE DRIFT CORRECTION

for worm in worms:
	# w = 'C01'
	w = worm.split('/')[-1].split('_')[0]
	path = worm[:-len(worm.split('/')[-1])]

	### Nicola computer
	outpath = os.path.join( path, w + '_analyzedImages' )
	# ### Simone computer
	# outpath = 'Y:\\Presentation\\Figures\\160211_ACVU\\movie_' + path.split('\\')[-2] + '_%s' %w

	print(outpath)

	### load parameters and times dataframes
	paramsDF = gf.load_data_frame_pandas( path, w + '_01params.pickle' )
	timesDF = gf.load_data_frame_pandas( path, w + '_01times.pickle' )
	gpDF = gf.load_data_frame_pandas( path, w + '_02gonadPos.pickle' )
	cellPosDF = gf.load_data_frame_pandas( path, w + '_04cellPos.pickle' )
	apdvPosDF = gf.load_data_frame_pandas( path, w + '_08apdvPos.pickle' )
	
	### load darkField
	#darkField = load_stack( 'X:\\Orca_calibration\\AVG_darkField.tif' )
	# darkField = gf.load_stack( 'Y:/Orca_calibration/AVG_darkField.tif' )
	darkField = gf.load_stack( '/Users/ngritti/Dropbox/ACVU_project_PRIVATE/ACVU_data/timelapse/pickledFiles/AVG_darkField.tif' )
	
	for channel in channels:

		### load flat field
		flatField = gf.load_stack( os.path.join( path, 'AVG_flatField_'+channel+'.tif' ) )
	
		### MAKE MOVIE FOR INTERESTING CELLS -> 1.p,1.pp,1.ppp and 4.a,4.aa,4.aaa
		### create the movies
		movie = [[],[],[],[],[],[]]

		# for idx, trow in timesDF.ix[timesDF.tidxRel==47].iterrows(): # test single timepoints
		for idx, trow in timesDF.iterrows():

			### skip problematic timepoints
			# if trow.tidxRel == 42:
			# 	continue

			currentCells = gf.extract_current_cell_pos( cellPosDF, trow.tidxRel )
			currentCellsOriginalPos = gf.extract_current_cell_pos( cellPosDF, trow.tidxRel )
			# print(currentCells)

			### if there are only 2 cells (one is the bckg) then skip the timepoint
			if len( currentCells.cname ) > 2:

				### find out which cells we want to show
				(c1,c2,b1,b2) = gf.find_interesting_cells( currentCells )
				interestingCells = pd.concat([c1, c2, b1, b2],axis=1).transpose().reset_index()

				refPoints = gf.extract_current_apdv_pos( apdvPosDF, trow.tidxRel )
				print(trow.tidxRel,list(interestingCells.cname))

				### load stack
				imgs = gf.load_stack( os.path.join( outpath, trow.fName + channel + '.tif') )
				# print(imgfile)
				imgSize = imgs[0].shape[0]

				### correct for XY
				gonadPos = gf.extract_pos( gpDF.ix[ gpDF.tidx == trow.tidxRel ].squeeze() ).astype(np.uint16)
				imgsXYCorr = gf.flat_field_correction( imgs, darkField, flatField, gonadPos )
				# print(imgsXYCorr[22])

				# ### plot raw data
				# plt.figure()
				# plt.xlim(0,imgs[21].shape[0])
				# plt.ylim(0,imgs[21].shape[1])
				# plt.gca().invert_yaxis()
				# plt.imshow(imgs[21],cmap='gray')
				# colors = ['r','g']
				# for kdx, cell in interestingCells.iterrows():
				# 	pos = extractCoords(cell)
				# 	plt.plot( pos[0], pos[1], 'o',color = colors[kdx])
				# for kdx, cell in refPoints.iterrows():
				# 	pos = extract_pos(cell) - gonadPos + 256
				# 	plt.plot( pos[0], pos[1], 'o',color = 'b')

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

				### rotate image

				vect = extractRefPointsCoordsO(refPoints,'p') - extractRefPointsCoordsO(refPoints,'a')
				alpha = -np.arctan2(vect[1],vect[0])
				# print(alpha)
				imgsRot = rotateImage(imgsXYCorr, alpha)

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
				imgsFinal, flipTF = flipImage( imgsRot, refPoints )

				### flip coordinates
				if flipTF:
					# print('flipping...')
					for jdx, cell in interestingCells.iterrows():
						posORot = extractCoordsORot(interestingCells, cell.cname)

						interestingCells.ix[interestingCells.cname==cell.cname,'cXposORot'] = posORot[0]
						interestingCells.ix[interestingCells.cname==cell.cname,'cYposORot'] = -posORot[1]

						interestingCells.ix[interestingCells.cname==cell.cname,'cXposRot'] = coordTransformO2C( extractCoordsORot(interestingCells,cell.cname), imgSize )[0]
						interestingCells.ix[interestingCells.cname==cell.cname,'cYposRot'] = coordTransformO2C( extractCoordsORot(interestingCells,cell.cname), imgSize )[1]

				# ### plot rotated and flipped image
				# plt.figure()
				# plt.xlim(0,imgsFinal[21].shape[0])
				# plt.ylim(0,imgsFinal[21].shape[1])
				# plt.gca().invert_yaxis()
				# plt.imshow(imgsFinal[21],cmap='gray')
				# colors = ['r','g']
				# for kdx, cell in interestingCells.iterrows():
				# 	pos = extractCoordsRot(cell)
				# 	plt.plot( pos[0], pos[1], 'o',color = colors[kdx])
				# for kdx, cell in refPoints.iterrows():
				# 	pos = extractCoordsRot(cell)
				# 	plt.plot( pos[0], pos[1], 'o',color = 'b')

				### crop images
				for jdx, cell in interestingCells.iterrows():

					pos = extractCoordsRot(cell)
					# small = imgsFinal[ currentCellsOriginalPos.ix[currentCellsOriginalPos.cname == cell.cname,'Z'], pos[1]-size:pos[1]+size, pos[0]-size:pos[0]+size ]
					small = gf.crop_image( imgsFinal[ int(currentCellsOriginalPos.ix[currentCellsOriginalPos.cname == cell.cname,'Z'].values[0]) ], np.array([pos[0],pos[1]]), 2*size )
					# print(cell)
					# print(small.shape)

					# print(cell, small.shape)
					# plt.figure()
					# plt.imshow(small,cmap='gray')

					movie[jdx].append(small)
					# print(np.array(movie[jdx]).shape)

				### APPEND FRAME FOR SISTER CELLS, IF ANY 
				for jdx, cname in enumerate(['1.ppa','4.aap']):
					isSis = False
					sis = None
					if cname in list( currentCells.cname ):
						sis = pd.DataFrame(currentCells.ix[ currentCells.cname == cname ])
						isSis = True

					if not isSis:
						movie[jdx+4].append( np.zeros((2*size,2*size)).astype(np.uint16) )

					else:
						### get coordinates with respect to center of the image (origin)
						posC = extractCoords(sis.squeeze(),channel)
						posO = coordTransformC2O(posC, imgSize)

						currentCells.ix[ currentCells.cname == cname, 'cXposO' ] = [ posO[0] ]
						currentCells.ix[ currentCells.cname == cname, 'cYposO' ] = [ posO[1] ]

						### rotate coordinates
						posO = extractCoordsO(currentCells, cname)
						posORot = rotateCoords(posO,alpha)

						currentCells.ix[ currentCells.cname == cname, 'cXposORot' ] = posORot[0]
						currentCells.ix[ currentCells.cname == cname, 'cYposORot' ] = posORot[1]

						currentCells.ix[ currentCells.cname == cname, 'cXposRot' ] = coordTransformO2C( extractCoordsORot(currentCells,cname), imgSize )[0]
						currentCells.ix[ currentCells.cname == cname, 'cYposRot' ] = coordTransformO2C( extractCoordsORot(currentCells,cname), imgSize )[1]

						### flip coordinates
						if flipTF:
							# print('flipping...')
							posORot = extractCoordsORot(currentCells, cname)

							currentCells.ix[ currentCells.cname == cname, 'cXposORot' ] = posORot[0]
							currentCells.ix[ currentCells.cname == cname, 'cYposORot' ] = -posORot[1]

							currentCells.ix[ currentCells.cname == cname, 'cXposRot' ] = coordTransformO2C( extractCoordsORot(currentCells,cname), imgSize )[0]
							currentCells.ix[ currentCells.cname == cname, 'cYposRot' ] = coordTransformO2C( extractCoordsORot(currentCells,cname), imgSize )[1]

						### crop images
						pos = extractCoordsRot(currentCells.ix[ currentCells.cname == cname].squeeze())
						# small = imgsFinal[ currentCellsOriginalPos.ix[currentCellsOriginalPos.cname == cname,'Z'], pos[1]-size:pos[1]+size, pos[0]-size:pos[0]+size ]
						small = gf.crop_image( imgsFinal[ int( currentCellsOriginalPos.ix[currentCellsOriginalPos.cname == cname,'Z'] ) ], np.array( [ int(pos[0]), int(pos[1]) ] ), 2*size )
						# print(cell)
						# print(small.shape)

						# print(cell, small.shape)
						# plt.figure()
						# plt.imshow(small,cmap='gray')

						movie[jdx+4].append(small)

		# print(np.array(movie[0]).shape,np.array(movie[0]).dtype)
		# print(np.array(movie[1]).shape,np.array(movie[1]).dtype)
		# print(np.array(movie[2]).shape,np.array(movie[2]).dtype)
		# print(np.array(movie[3]).shape,np.array(movie[3]).dtype)
		## save movies
		imsave( os.path.join( outpath, 'cell_'+channel+'_1ppp.tif' ), np.array(movie[0]).astype(np.uint16) )
		imsave( os.path.join( outpath, 'cell_'+channel+'_4aaa.tif' ), np.array(movie[1]).astype(np.uint16) )
		imsave( os.path.join( outpath, 'bckg_'+channel+'_1.tif' ), np.array(movie[2]).astype(np.uint16) )
		imsave( os.path.join( outpath, 'bckg_'+channel+'_4.tif' ), np.array(movie[3]).astype(np.uint16) )
		imsave( os.path.join( outpath, 'cell_'+channel+'_1ppa.tif' ), np.array(movie[4]).astype(np.uint16) )
		imsave( os.path.join( outpath, 'cell_'+channel+'_4aap.tif' ), np.array(movie[5]).astype(np.uint16) )
		imsave( os.path.join( outpath, 'combined_withSisters_'+channel+'.tif' ), np.concatenate( 
																				( np.array(movie[4]),
																				np.array(movie[0]),
																				np.array(movie[1]),
																				np.array(movie[5]) ), 2 ).astype(np.uint16) )

plt.show()

