import numpy as np
import matplotlib.pyplot as plt
import pickle
from matplotlib import cm
from generalFunctions import *
from skimage.filters import threshold_otsu
from skimage import morphology
from scipy.ndimage.morphology import binary_fill_holes
import cv2

import matplotlib as mpl
mpl.rcParams['pdf.fonttype'] = 42

path = 'W:\\Simone\\161002_lag2GFP+mCherry\\'
worm = 'C05'

paramsDF = load_data_frame( path, worm + '_01params.pickle' )
timesDF = load_data_frame( path, worm + '_01times.pickle' )
gpDF = load_data_frame( path, worm + '_02gonadPos.pickle' )
cellPosDF = load_data_frame( path, worm + '_04cellPos.pickle' )
cellOutDF = load_data_frame( path, worm + '_05cellOut.pickle' )
# cellFluoDF = load_data_frame( path, worm + '_06cellFluo.pickle' )

tidx = 97

fname = timesDF.ix[ timesDF.tidxRel == tidx, 'fName' ].values[0]
print(fname)
imgs = load_stack(path+worm+'_analyzedImages\\'+fname+'488nm.tif')

currentCellsData = extract_current_cell_out( cellPosDF, cellOutDF, tidx )
print('cells in this timepoint:\n',currentCellsData)

cname = '1.ppp'
cellData = currentCellsData.ix[currentCellsData.cname == cname]
print('\nCell of interest:\n',cellData)

cropimg = imgs[ int(cellData.Z),
	int(cellData.Y-cellData.imgPxl/2):int(cellData.Y+cellData.imgPxl/2),
	int(cellData.X-cellData.imgPxl/2):int(cellData.X+cellData.imgPxl/2)]

plt.figure()
plt.imshow( cropimg, 
	cmap = 'gray', interpolation = 'nearest' )


mask = cropimg >= threshold_otsu(cropimg)
plt.figure()
plt.imshow(mask, cmap = 'gray', interpolation = 'nearest' )
mask = binary_fill_holes( mask )
mask = morphology.remove_small_objects( mask, 100 )
mask = morphology.binary_dilation( cropimg >= threshold_otsu(cropimg), selem = np.ones((5,5)) )
lbl = morphology.label(mask)
nlbl = lbl[int(mask.shape[0]/2),int(mask.shape[1]/2)]
print(nlbl)
mask = lbl==1
print(mask.shape, np.sum(mask))

cont = cv2.findContours(mask.astype(np.uint8), cv2.RETR_TREE, cv2.CHAIN_APPROX_SIMPLE)[1][0]
cont = cv2.convexHull(cont)
cont = cont[:,0,:]

print(cont)
# print(cont[1])
# print(cont[1][0][:,0,:])
vertices = [list(cont[:,0]),list(cont[:,1])]
vertices[0].append(vertices[0][0])
vertices[1].append(vertices[1][0])
vertices = np.array(vertices)
plt.plot(vertices[0],vertices[1],'-y')

plt.show()



