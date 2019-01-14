# -*- coding: utf-8 -*-
"""
Created on Tue Dec  1 10:31:00 2015

@author: kienle
"""

import numpy as np
import matplotlib.pyplot as plt
import pickle
from matplotlib import cm
from mpl_toolkits.mplot3d import Axes3D
from generalFunctions import *
from skimage import filters, morphology
import scipy.spatial as spatial


import matplotlib as mpl
mpl.rcParams['pdf.fonttype'] = 42

def plotFluorescence(path,worm,ax1,ax2,channel='488nm',correction =True):

  # import worm data
  timesDF = load_data_frame( path, worm + '_01times.pickle' )
  gonadPosDF = load_data_frame( path, worm + '_02gonadPos.pickle' )
  cellPosDF = load_data_frame( path, worm + '_04cellPos.pickle' )
  cellOutDF = load_data_frame( path, worm + '_05cellOut.pickle' )

  # import darkField and flatField
  darkField = load_stack( 'W:\\Orca_calibration\\AVG_darkField.tif' )
  flatField = load_stack( os.path.join( path, 'AVG_flatField_' + channel + '.tif' ) )

  ### find ecdysis timepoint
  ecd = np.loadtxt( open( os.path.join( path, 'skin.txt'), 'rb' ) )
  # load ecdysis data
  index = np.where( ecd[:,0] == float(worm[1:]) )
  mintp = np.min( [ i for i in ecd[index, 1:6][0][0] if i >= 0 ] )
  lethtidx = ecd[ index, 2:6 ][0][0] - 1
  tpL2 = timesDF.ix[ timesDF.tidxRel == lethtidx[0], 'timesRel' ].values[0]
  tpL1 = timesDF.ix[ timesDF.tidxRel == mintp, 'timesRel' ].values[0]

  # relative to hatching time
  tpRel=tpL1

  ### check FF correction
  # mag = np.median( flatField - darkField ) / ( flatField - darkField )
  # data = [[],[],[],[]]
  # for idx, y in enumerate( np.arange(2048)[::64] ):
  #   for jdx, x in enumerate( np.arange(2048)[::64] ):
  #     print(x,y)
  #     data[0].append(x)
  #     data[1].append(y)
  #     data[2].append(y)
  #     data[3].append(mag[y,x])
  # return data

  ### create dataframe for fluorescence  
  cellList = list( cellPosDF.keys() )
  singleCell = pd.DataFrame( { 'tidx': timesDF.tidxRel,
          'times' : timesDF.timesRel,
          'fluo': np.nan,
          'area': np.nan,
          'bfluo': np.nan,
          'barea': np.nan,
          'X': np.nan,
          'Y': np.nan,
          'Z': np.nan } )
  # singleCell.fluo = singleCell.fluo.astype(floating)

  cellList = [ cn for cn in cellList if cn[0] != 'b' ]
  cellFluoDF = { cn: singleCell.copy() for cn in cellList }

  ### filter timesDF for tests (much faster!!!)
  # tFilter = ( timesDF.tidxRel > 75 ) & ( timesDF.tidxRel < 80 )
  # timesDF = timesDF[ tFilter ]

  for idx, tRow in timesDF.iterrows():
    
    print(tRow.fName)

    #extract position of the gonad in the raw image and information of all labeled cells
    gonadPos = gonadPosDF.ix[ gonadPosDF.tidx == tRow.tidxRel ].squeeze()
    currentCells = extract_current_cell_out( cellPosDF, cellOutDF, tRow.tidxRel )
    # print(currentCells)

    if len( currentCells ) > 0:
    
      ### for each timepoint, if there are cells labeled, load the cropped image
      stack = load_stack( os.path.join( path, worm + '_analyzedImages', tRow.fName + channel + '.tif') )
      if correction:
        stack = flat_field_correction( stack, darkField, flatField, np.array( gonadPos ) )
        # stack = np.mean( flatField - darkField ) * ( stack - darkField ) / ( flatField - darkField )

      ### check the stack
      # plt.figure()
      # plt.imshow( stack[10], interpolation = 'nearest', cmap = 'gray', vmin = 0, vmax = 3000 )
      # plt.plot( gonadPos.X, gonadPos.Y, 'or' )

      ### check the flat field correction - NB: it doeasn't work now!!!
      # plt.figure()
      # plt.imshow( np.mean( flatField - darkField ) * ( stack[10] - darkField ) / ( flatField - darkField ), interpolation = 'nearest', cmap = 'gray', vmin = 0, vmax = 3000 )
      # plt.plot( gonadPos.X, gonadPos.Y, 'or' )
      # plt.show()

      currentCellList = list( currentCells.cname )
      currentCellList = [ c for c in currentCellList if c[0]!='b' ]
      currentCellList.sort()
      for cname in currentCellList:
        # print(cname)

        ### crop the full image to the area of the cell
        cell = currentCells.ix[ currentCells.cname == cname ].squeeze()
        cropimg = stack[ int( cell.Z ), int( cell.Y - cell.imgPxl/2 ) : int( cell.Y + cell.imgPxl/2 ), int( cell.X - cell.imgPxl/2 ) : int( cell.X + cell.imgPxl/2 )  ]
        # cropimg1 = stack[ int( 15 ), int( cell.Y - cell.imgPxl/2 ) : int( cell.Y + cell.imgPxl/2 ), int( cell.X - cell.imgPxl/2 ) : int( cell.X + cell.imgPxl/2 )  ]

        ### check the cropimg
        # plt.figure()
        # plt.imshow( cropimg, interpolation = 'nearest', cmap = 'gray', vmin = 0, vmax = 3000 )
        # plt.show()

        ### create the mask
        cellOut = extract_out( cell ) * cell.imgPxl / 1000.
        vertices = np.array( [ np.append(cellOut[:,0],cellOut[0,0]), np.append(cellOut[:,1],cellOut[0,1]) ] ).T
        p = Path(vertices)
        points = [ (i,j) for i in np.arange(cell.imgPxl) for j in np.arange(cell.imgPxl) ]
        mask = p.contains_points(points).reshape(int(cell.imgPxl),int(cell.imgPxl)).T

        ### create the background mask
        bmask = morphology.binary_dilation(mask,morphology.disk(10)) != mask

        ### check the mask
        # plt.figure()
        # plt.imshow(mask,interpolation = 'nearest',cmap = cm.OrRd, alpha = .5)
        # plt.figure()
        # plt.imshow(bmask,interpolation = 'nearest',cmap = cm.OrRd, alpha = .5)
        # plt.show()

        ### compute the signal
        signal = np.sum( mask * cropimg )# / np.sum( mask )
        bsignal = np.sum( bmask * cropimg )# / np.sum( mask )
        area = np.sum( mask )
        barea = np.sum( bmask )

        ### retrieve cell Pos relative to original image
        cellPos = np.array( [ gonadPos.X, gonadPos.Y ] ) + np.array( [ cell.X, cell.Y ] ) - 256

        ### save data
        cellFluoDF[ cname ].ix[ cellFluoDF[ cname ].tidx == tRow.tidxRel, 'X' ] = cellPos[0]
        cellFluoDF[ cname ].ix[ cellFluoDF[ cname ].tidx == tRow.tidxRel, 'Y' ] = cellPos[1]
        cellFluoDF[ cname ].ix[ cellFluoDF[ cname ].tidx == tRow.tidxRel, 'Z' ] = cell.Z
        cellFluoDF[ cname ].ix[ cellFluoDF[ cname ].tidx == tRow.tidxRel, 'area' ] = area
        cellFluoDF[ cname ].ix[ cellFluoDF[ cname ].tidx == tRow.tidxRel, 'fluo' ] = np.float( signal )
        cellFluoDF[ cname ].ix[ cellFluoDF[ cname ].tidx == tRow.tidxRel, 'barea' ] = barea
        cellFluoDF[ cname ].ix[ cellFluoDF[ cname ].tidx == tRow.tidxRel, 'bfluo' ] = np.float( bsignal )

        # return
  # print(cellList)

  ### plot timeseries
  for cname in cellList:
    area = cellFluoDF[ cname ].ix[ pd.notnull( cellFluoDF[ cname ].X ), 'area' ]
    barea = cellFluoDF[ cname ].ix[ pd.notnull( cellFluoDF[ cname ].X ), 'barea' ]
    fluo = cellFluoDF[ cname ].ix[ pd.notnull( cellFluoDF[ cname ].X ), 'fluo' ]
    bfluo = cellFluoDF[ cname ].ix[ pd.notnull( cellFluoDF[ cname ].X ), 'bfluo' ]
    times = cellFluoDF[ cname ].ix[ pd.notnull( cellFluoDF[ cname ].X ), 'times' ]

    ax1.plot( times - tpL1, (fluo / area - bfluo / barea ) * area )

  ### 3D plot
  #reorganize data
  data = [[],[],[],[]]
  for cname in cellList:
    cell = cellFluoDF[ cname ].ix[ pd.notnull( cellFluoDF[ cname ].X ) ]
    for x in cell.X:
      data[0].append(x)
    for y in cell.Y:
      data[1].append(y)
    for z in cell.Z:
      data[2].append(z)
    for f in cell.fluo:
      data[3].append(f)

  # ax2.scatter(data[0], data[1], data[2], c=data[3], marker='o', s = 75)
  return data

      


#################################################################################################
def fmt(x, y, z):
    return 'x: {x:0.2f}\ny: {y:0.2f}\nz: {z:0.2f}'.format(x=x, y=y, z=z)
      
class FollowDotCursor(object):
    """Display the x,y location of the nearest data point.
    http://stackoverflow.com/a/4674445/190597 (Joe Kington)
    http://stackoverflow.com/a/13306887/190597 (unutbu)
    http://stackoverflow.com/a/15454427/190597 (unutbu)
    """
    def __init__(self, ax, x, y, z, tolerance=5, formatter=fmt, offsets=(-20, 20)):
        try:
            x = np.asarray(x, dtype='float')
        except (TypeError, ValueError):
            x = np.asarray(mdates.date2num(x), dtype='float')
        y = np.asarray(y, dtype='float')
        mask = ~(np.isnan(x) | np.isnan(y))
        x = x[mask]
        y = y[mask]
        self.z = z
        self._points = np.column_stack((x, y))
        self.offsets = offsets
        y = y[np.abs(y-y.mean()) <= 3*y.std()]
        self.scale = x.ptp()
        self.scale = y.ptp() / self.scale if self.scale else 1
        self.tree = spatial.cKDTree(self.scaled(self._points))
        self.formatter = formatter
        self.tolerance = tolerance
        self.ax = ax
        self.fig = ax.figure
        self.ax.xaxis.set_label_position('top')
        self.dot = ax.scatter(
            [x.min()], [y.min()], s=130, color = 'green', alpha=0.7)
        self.annotation = self.setup_annotation()
        plt.connect('motion_notify_event', self)

    def scaled(self, points):
        points = np.asarray(points)
        return points * (self.scale, 1)

    def __call__(self, event):
        ax = self.ax
        # event.inaxes is always the current axis. If you use twinx, ax could be
        # a different axis.
        if event.inaxes == ax:
            x, y = event.xdata, event.ydata
        elif event.inaxes is None:
            return
        else:
            inv = ax.transData.inverted()
            x, y = inv.transform([(event.x, event.y)]).ravel()
        annotation = self.annotation
        (x, y), idx = self.snap(x, y)
        annotation.xy = x, y
        annotation.set_text(self.formatter(x, y, self.z[idx]))
        self.dot.set_offsets((x, y))
        bbox = ax.viewLim
        event.canvas.draw()

    def setup_annotation(self):
        """Draw and hide the annotation box."""
        annotation = self.ax.annotate(
            '', xy=(0, 0), ha = 'right',
            xytext = self.offsets, textcoords = 'offset points', va = 'bottom',
            bbox = dict(
                boxstyle='round,pad=0.5', fc='yellow', alpha=0.75),
            arrowprops = dict(
                arrowstyle='->', connectionstyle='arc3,rad=0'))
        return annotation

    def snap(self, x, y):
        """Return the value in self.tree closest to x, y."""
        dist, idx = self.tree.query(self.scaled((x, y)), k=1, p=1)
        try:
            return self._points[idx], idx
        except IndexError:
            # IndexError: index out of bounds
            return self._points[0], idx


####################################################################################################

if __name__ == '__main__':

  ### setup figure for the timeseries
  fig1 = plt.figure(figsize=(5.8,3.8))
  ax1 = fig1.add_subplot(111)
  fig1.subplots_adjust(left=0.15, right=.95, top=.95, bottom=0.15)
  for tl in ax1.get_xticklabels():
    tl.set_fontsize(18)
  for tl in ax1.get_yticklabels():
    tl.set_fontsize(18)
  # ax1.set_ylim((0,2048))
  # ax1.set_xlim((0,2048))

  ### setup figure for 3D plot
  fig2 = plt.figure(figsize=(5.8,3.8))
  ax2 = fig2.add_subplot(111, projection = '3d')
  fig1.subplots_adjust(left=0.15, right=.95, top=.95, bottom=0.15)
  ax2.set_ylim((0,2048))
  ax2.set_xlim((0,2048))
  # cursor = FollowDotCursor(ax2, x, y, tolerance=20)

  path = 'W:\\Nicola\\160701_LAG2gfp_1Xoutcrossed_line2_withBeads'
  worms = ['C01']#,'C02']
  lineages = ['1.ppp']

  for idx, w in enumerate( worms ):
    
    data = plotFluorescence( path, w, ax1, ax2, channel = '488nm', correction = True )

    ax2.scatter(data[0], data[1], data[2], c=data[3], marker='o', s = 75)
    # cursor = FollowDotCursor(ax2, data[0], data[1], data[3], tolerance=2)

    
  plt.show()

