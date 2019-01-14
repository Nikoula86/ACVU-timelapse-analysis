

"""
PyQt seam cells analysis GUI

NB: the python package tifffile (Gohlke) needs to be installed.

author: Nicola Gritti
last edited: June 2015
"""

import sys
from tifffile import *
from generalFunctions import *
import pickle
import os
from PyQt4 import QtGui, QtCore
import numpy as np
from matplotlib.figure import Figure
import matplotlib.pyplot as plt
from matplotlib import cm
import glob
import pandas as pd
import time
from matplotlib.colors import LinearSegmentedColormap
from skimage import morphology

# create colormaps
cdict1 = {'red':   ((0.0, 0.0, 0.0),
                    (1.0, 1.0, 1.0)),

         'green': ((0.0, 0.0, 0.0),
                   (1.0, 0.0, 0.0)),

         'blue':  ((0.0, 0.0, 0.0),
                   (1.0, 0.0, 0.0))
         }

black_red = LinearSegmentedColormap('BlueRed1', cdict1)

cdict2 = {'red':   ((0.0, 0.0, 0.0),
                    (1.0, 0.0, 0.0)),

         'green': ((0.0, 0.0, 0.0),
                   (1.0, 0.0, 0.0)),

         'blue':  ((0.0, 0.0, 0.0),
                   (1.0, 1.0, 1.0))
         }

black_blue = LinearSegmentedColormap('BlueRed1', cdict2)

#############################################



path = 'W:\\Nicola\\160701_LAG2gfp_1Xoutcrossed_line2_withBeads'
worm = 'C01'

pathDial = path + '\\' + worm + '_analyzedImages'

### load parameters and times dataframes
paramsDF = load_data_frame( path, worm + '_01params.pickle' )
timesDF = load_data_frame( path, worm + '_01times.pickle' )
gpDF = load_data_frame( path, worm + '_02gonadPos.pickle' )
cellPosDF = load_data_frame( path, worm + '_04cellPos.pickle' )
cellOutDF_inner = load_data_frame( path, worm + '_05cellOut_1.pickle' )
cellOutDF_outer = load_data_frame( path, worm + '_05cellOut_2.pickle' )

lineage = '1.ppp'

cellFluoDF = load_data_frame( path, worm + '_06cellFluo_1.pickle' )
key = lineage
cell = cellFluoDF[ key ].ix[ pd.notnull( cellFluoDF[ key ]._488nm ) ]
area_cell = cellOutDF_inner[ key ].ix[ pd.notnull( cellFluoDF[ key ]._488nm ), 'area' ]		

plt.plot(timesDF[ pd.notnull( cellOutDF_inner[ lineage ].area ) ].timesRel,cell._488nm)
plt.plot(timesDF[ pd.notnull( cellOutDF_inner[ lineage ].area ) ].timesRel,area_cell*10)

x_eval = timesDF[ pd.notnull( cellOutDF_inner[ lineage ].area ) ].timesRel.values#np.linspace(np.min(timesDF[ pd.notnull( cellOutDF_inner[ lineage ].area ) ].timesRel), np.max(timesDF[ pd.notnull( cellOutDF_inner[ lineage ].area ) ].timesRel), 101)
sigma = 0.2

delta_x = x_eval[:,None] - timesDF[ pd.notnull( cellOutDF_inner[ lineage ].area ) ].timesRel.values
weights = np.exp(-delta_x*delta_x / (2*sigma*sigma)) / (np.sqrt(2*np.pi) * sigma)
weights /= np.sum(weights, axis=1, keepdims=True)
y_eval = np.dot(weights, cell._488nm)

plt.plot(x_eval,y_eval)


# delta_x = x_eval[:, None] - timesDF[ pd.notnull( cellOutDF_inner[ lineage ].area ) ].timesRel

# print(timesDF[pd.notnull(cellOutDF_inner[lineage].area)])

# tseries = [[],[]]
# for idx, tRow in timesDF[ pd.notnull( cellOutDF_inner[ lineage ].area ) ].iterrows():
# 	print(tRow.tidxRel)
# 	tp = tRow.tidxRel
# 	tRow = timesDF.ix[ timesDF.tidxRel == tp ].squeeze()
# 	fileName = os.path.join( pathDial, tRow.fName + '488nm.tif')
# 	firststack = load_stack( fileName )

# 	# find cell image
# 	cell = extract_current_cell_out( cellPosDF, cellOutDF_inner, tp )
# 	cell = cell.ix[ cell.cname == '1.ppp' ].squeeze()
# 	cellimg = firststack[cell.Z, cell.Y-cell.imgPxl/2.:cell.Y+cell.imgPxl/2., cell.X-cell.imgPxl/2.:cell.X+cell.imgPxl/2.]

# 	# find inner mask
# 	cell = extract_current_cell_out( cellPosDF, cellOutDF_inner, tp )
# 	cell = cell.ix[ cell.cname == '1.ppp' ].squeeze()
# 	vertices = np.array( [ np.append(cell.Xout,cell.Xout[0]) * cell.imgPxl / 1000., np.append(cell.Yout[:],cell.Yout[0]) * cell.imgPxl / 1000. ] ).T
# 	p = Path(vertices)
# 	# create the mask (image full of 0 and 1, the 1s are wherethe cell is)
# 	points = [ (i,j) for i in np.arange(cell.imgPxl) for j in np.arange(cell.imgPxl) ]
# 	maskinner = p.contains_points(points).reshape(cell.imgPxl,cell.imgPxl).T

# 	# find outer mask
# 	cell = extract_current_cell_out( cellPosDF, cellOutDF_outer, tp )
# 	cell = cell.ix[ cell.cname == '1.ppp' ].squeeze()
# 	vertices = np.array( [ np.append(cell.Xout,cell.Xout[0]) * cell.imgPxl / 1000., np.append(cell.Yout[:],cell.Yout[0]) * cell.imgPxl / 1000. ] ).T
# 	p = Path(vertices)
# 	# create the mask (image full of 0 and 1, the 1s are wherethe cell is)
# 	points = [ (i,j) for i in np.arange(cell.imgPxl) for j in np.arange(cell.imgPxl) ]
# 	maskouter = p.contains_points(points).reshape(cell.imgPxl,cell.imgPxl).T

# 	# plt.figure()
# 	# img1 = plt.imshow(cellimg, cmap = 'gray', interpolation = 'nearest')
# 	# img2 = plt.imshow(maskinner, cmap = black_blue, interpolation = 'nearest', alpha = .5)
# 	# img3 = plt.imshow(maskouter, cmap = black_red, interpolation = 'nearest', alpha = .1)

# 	# find out how many layers there are in between
# 	newmask = np.copy(maskinner)
# 	layernumber = 0
# 	overlap = 1
# 	data = [[],[]]
# 	while overlap > 0:

# 		# plt.figure()
# 		# img1 = plt.imshow(cellimg, cmap = 'gray', interpolation = 'nearest')
# 		# img2 = plt.imshow(newmask, cmap = black_blue, interpolation = 'nearest', alpha = .1)
# 		# img3 = plt.imshow(maskouter, cmap = black_red, interpolation = 'nearest', alpha = .1)

# 		overlap = np.sum( maskouter - ( maskouter * newmask ) )
# 		layernumber += 1

# 		data[0].append( np.sum( newmask * maskouter ) )
# 		data[1].append( np.sum( newmask * maskouter * cellimg ) )

# 		# if np.mod( layernumber, 2 ) == 0:
# 		# 	newmask = morphology.binary_dilation( newmask, morphology.square(3) )
# 		# else:
# 		newmask = morphology.binary_dilation( newmask, morphology.disk(1) )
# 		# print(np.sum(maskinner), np.sum(maskouter),overlap, layernumber)

# 	# plt.figure()
# 	# img1 = plt.imshow(cellimg, cmap = 'gray', interpolation = 'nearest')
# 	# img2 = plt.imshow(newmask, cmap = black_blue, interpolation = 'nearest', alpha = .5)
# 	# img3 = plt.imshow(maskouter, cmap = black_red, interpolation = 'nearest', alpha = .1)

# 	data = np.array(data)

# 	weights = np.arange(len(data[0]))[::-1]+1
# 	# print(data[0],data[1],weights)
# 	# print( np.sum( data[1] * weights ) / ( np.sum(weights) ) )

# 	tseries[0].append(tRow.tidxRel)
# 	tseries[1].append( np.sum( data[1] / data[0] * weights ) / ( np.sum(weights) ) )


# plt.plot(tseries[0],tseries[1])

plt.show()
