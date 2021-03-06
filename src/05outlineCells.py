

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
from matplotlib.backend_bases import key_press_handler
from matplotlib.backends.backend_qt4agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends import qt_compat
import glob
import pandas as pd
from scipy.ndimage.morphology import binary_fill_holes
from skimage import morphology
from skimage.filters import threshold_otsu
import cv2
use_pyside = qt_compat.QT_API == qt_compat.QT_API_PYSIDE

class GUI(QtGui.QWidget):
    
	def __init__(self):

		super(GUI, self).__init__()

		self.setWindowTitle( 'Outline Cells' )
		self.cellNames = ['1.p','4.a','1.pp','4.aa','1.ppa','1.ppp','4.aaa','4.aap','b_1','b_4']
		self.initUI()
    
    #-----------------------------------------------------------------------------------------------
    # INITIALIZATION OF THE WINDOW - DEFINE AND PLACE ALL THE WIDGETS
    #-----------------------------------------------------------------------------------------------

	def initUI(self):
		# SET THE GEOMETRY

		mainWindow = QtGui.QVBoxLayout()
		mainWindow.setSpacing(15)

		fileBox = QtGui.QHBoxLayout()
		spaceBox1 = QtGui.QHBoxLayout()
		rawDataBox = QtGui.QHBoxLayout()

		mainWindow.addLayout(fileBox)
		mainWindow.addLayout(spaceBox1)
		mainWindow.addLayout(rawDataBox)

		Col1 = QtGui.QGridLayout()
		Col2 = QtGui.QHBoxLayout()
		Col3 = QtGui.QVBoxLayout()

		rawDataBox.addLayout(Col1)
		rawDataBox.addLayout(Col2)
		rawDataBox.addLayout(Col3)

		self.setLayout(mainWindow)

		# DEFINE ALL WIDGETS AND BUTTONS

		loadBtn = QtGui.QPushButton('Load DataSet')
		saveBtn = QtGui.QPushButton('Save data (F12)')

		tpLbl = QtGui.QLabel('Relative Tp:')
		slLbl = QtGui.QLabel('Slice:')
		fNameLbl = QtGui.QLabel('File name:')
		zoomLbl = QtGui.QLabel('Zoom (imgPxl):')

		self.tp = QtGui.QSpinBox(self)
		self.tp.setValue(-5)
		self.tp.setMaximum(100000)

		self.sl = QtGui.QSpinBox(self)
		self.sl.setValue(0)
		self.sl.setMaximum(100000)

		self.fName = QtGui.QLabel('')

		self.zoom = QtGui.QSpinBox(self)
		self.zoom.setValue(50)
		self.zoom.setMinimum(10)
		self.zoom.setMaximum(150)
		self.zoom.setSingleStep(10)

		self._488nmBtn = QtGui.QRadioButton('488nm')
		self._561nmBtn = QtGui.QRadioButton('561nm')
		self.CoolLEDBtn = QtGui.QRadioButton('CoolLED')

		outlineBtn = QtGui.QPushButton('Automatic outline')
		fullOutlineBtn = QtGui.QPushButton('Automatic outline for all images')

		self.sld1 = QtGui.QSlider(QtCore.Qt.Vertical, self)
		self.sld1.setMaximum(2**16-1)
		self.sld1.setValue(0)
		self.sld2 = QtGui.QSlider(QtCore.Qt.Vertical, self)
		self.sld2.setMaximum(2**16)
		self.sld2.setValue(2**16-1)

		self.fig1 = Figure((8.0, 8.0), dpi=100)
		self.fig1.subplots_adjust(left=0., right=1., top=1., bottom=0.)
		self.ax1 = self.fig1.add_subplot(111)
		self.canvas1 = FigureCanvas(self.fig1)
		self.canvas1.setFocusPolicy( QtCore.Qt.ClickFocus )
		self.canvas1.setFocus()
		self.canvas1.setFixedSize(QtCore.QSize(600,600))
		self.canvas1.setSizePolicy( QtGui.QSizePolicy.Expanding, QtGui.QSizePolicy.Expanding )

		self.fig2 = Figure((4.0, 4.0), dpi=100)
		self.fig2.subplots_adjust(left=0., right=1., top=1., bottom=0.)
		self.ax2 = self.fig2.add_subplot(111)
		self.canvas2 = FigureCanvas(self.fig2)
		self.canvas2.setFixedSize(QtCore.QSize(300,300))
		self.canvas2.setSizePolicy( QtGui.QSizePolicy.Expanding, QtGui.QSizePolicy.Expanding )

		# PLACE ALL THE WIDGET ACCORDING TO THE GRIDS

		fileBox.addWidget(loadBtn)
		fileBox.addWidget(saveBtn)

		spaceBox1.addWidget(self.HLine())

		Col1.addWidget(tpLbl, 0, 0)#, 1, 1, Qt.AlignTop)
		Col1.addWidget(self.tp, 0, 1)#, 1, 1, Qt.AlignTop)
		Col1.addWidget(slLbl, 1, 0)#, 1, 1, Qt.AlignTop)
		Col1.addWidget(self.sl, 1, 1)#, 1, 1, Qt.AlignTop)
		Col1.addWidget(fNameLbl, 2, 0)
		Col1.addWidget(self.fName, 2, 1)
		Col1.addWidget(zoomLbl, 3, 0)
		Col1.addWidget(self.zoom, 3, 1)
		Col1.addWidget(self._488nmBtn, 4, 0 )
		Col1.addWidget(self._561nmBtn, 5, 0 )
		Col1.addWidget(self.CoolLEDBtn, 6, 0 )
		Col1.addWidget(outlineBtn,7,0,1,2)
		Col1.addWidget(fullOutlineBtn,8,0,1,2)

		Col2.addWidget(self.sld1)
		Col2.addWidget(self.sld2)
		Col2.addWidget(self.canvas1)

		Col3.addWidget(self.canvas2)

		self.setFocus()
		self.show()

		# BIND BUTTONS TO FUNCTIONS

		loadBtn.clicked.connect(self.selectWorm)
		saveBtn.clicked.connect(self.saveData)
		outlineBtn.clicked.connect(self.automaticOutline)
		fullOutlineBtn.clicked.connect(self.fullAutomaticOutline)

		self.tp.valueChanged.connect(self.loadNewStack)
		self.sl.valueChanged.connect(self.updateAllCanvas)
		self.zoom.valueChanged.connect(self.changeZoom)
		self.sld1.valueChanged.connect(self.updateBC)
		self.sld2.valueChanged.connect(self.updateBC)

		self._488nmBtn.toggled.connect(self.radio488Clicked)
		self._561nmBtn.toggled.connect(self.radio561Clicked)
		self.CoolLEDBtn.toggled.connect(self.radioCoolLEDClicked)

		self.fig1.canvas.mpl_connect('button_press_event',self.onMouseClickOnCanvas1)        
		self.fig1.canvas.mpl_connect('scroll_event',self.wheelEvent)        

	#-----------------------------------------------------------------------------------------------
	# FORMATTING THE WINDOW
	#-----------------------------------------------------------------------------------------------

	def center(self):
        
		qr = self.frameGeometry()
		cp = QtGui.QDesktopWidget().availableGeometry().center()
		qr.moveCenter(cp)
		self.move(qr.topLeft())
        
	def HLine(self):

		toto = QtGui.QFrame()
		toto.setFrameShape(QtGui.QFrame.HLine)
		toto.setFrameShadow(QtGui.QFrame.Sunken)
		return toto

	def VLine(self):
	    
	    toto = QtGui.QFrame()
	    toto.setFrameShape(QtGui.QFrame.VLine)
	    toto.setFrameShadow(QtGui.QFrame.Sunken)
	    return toto

	def heightForWidth(self, width):

		return width

	#-----------------------------------------------------------------------------------------------
	# BUTTON FUNCTIONS
	#-----------------------------------------------------------------------------------------------

	def selectWorm(self):

		### store the folders
		self.pathDial = QtGui.QFileDialog.getExistingDirectory(self, 'Select a folder', 'X:\\Simone\\160129_MCHERRY_HLH2GFP_onHB101')#'Y:\\Images')
		self.worm = self.pathDial.split("\\")[-1].split('_')[0]
		self.path = os.path.dirname( self.pathDial )
		self.setWindowTitle('Outline - ' + self.pathDial.split("\\")[-2] + ' - ' + self.pathDial.split("\\")[-1][:3] )

		### give error message if there is no CoolLED movie in the selected folder
		flist = glob.glob( self.pathDial + '\\*_movie.tif' )
		if len(flist)==0:#not os.path.isfile( os.path.join( self.pathDial, '*_movie.tif' ) ):
			QtGui.QMessageBox.about(self,'Warning!','There is no movie in this folder! Create a movie first!')
			return

		### load parameters and times dataframes
		self.paramsDF = load_data_frame( self.path, self.worm + '_01params.pickle' )
		self.timesDF = load_data_frame( self.path, self.worm + '_01times.pickle' )
		self.gpDF = load_data_frame( self.path, self.worm + '_02gonadPos.pickle' )
		self.cellPosDF = load_data_frame( self.path, self.worm + '_04cellPos.pickle' )

		# extract some info
		self.compression = self.paramsDF.compression
		self.hatchingtidx = int( self.paramsDF.tidxHatch )

		### if the cellOutline pickle file already exists, load it, otherwise create a blank one
		if os.path.isfile( os.path.join(self.path, self.worm + '_05cellOut.pickle' ) ):
			self.cellOutDF = load_data_frame( self.path, self.worm + '_05cellOut.pickle' )

		else:
			self.cellOutDF = create_cell_out( self.timesDF, self.cellNames )

		# detect available channels
		self.channels = []
		chns = ['CoolLED','488nm','561nm']
		for c in chns:

			if os.path.isfile( os.path.join( self.pathDial, c + '_movie.tif' ) ):

				self.channels.append(c)
		print(self.channels)
		self.currentChannel = self.channels[0]

		### detect size of the cropped images
		tp = first_tidx_pos_all_cells( self.cellPosDF )
		self.prevtp = tp-1
		tRow = self.timesDF.ix[ self.timesDF.tidxRel == tp ].squeeze()
		fileName = os.path.join( self.pathDial, tRow.fName + self.currentChannel + '.tif')
		firststack = load_stack( fileName )
		self.cropsize = firststack.shape[1]
		self.nslices = firststack.shape[0]

		### intialize canvases
		self.initializeCanvas1()
		self.initializeCanvas2()

		### extract current cells already labeled
		self.currentCells = extract_current_cell_out( self.cellPosDF, self.cellOutDF, self.tp.value(), self.zoom.value() )
		self.analyzedCell = '---'

		### update the text of the fileName
		self.fName.setText( self.timesDF.ix[ self.timesDF.tidxRel == tp, 'fName' ].values[0])

		### set the timepoint to the hatching time
		self.tp.setMinimum(np.min(self.timesDF.tidxRel))
		self.tp.setMaximum(np.max(self.timesDF.tidxRel))

		### set the max slice number
		self.sl.setMaximum( self.nslices-1 )

		self.tp.setValue( tp )

		if self.currentChannel == 'CoolLED':
			self.CoolLEDBtn.setChecked(True)    # this uppdates the canvas1 once more
		elif self.currentChannel == '561nm':
			self._561nmBtn.setChecked(True)    # this uppdates the canvas1 once more
		elif self.currentChannel == '488nm':
			self._488nmBtn.setChecked(True)    # this uppdates the canvas1 once more

		# self.pathDial.show()
		self.setFocus()

	def loadNewStack(self):

		# update the cell outline data frame before updating the images and retrireving new cells
		newCellOutDF = update_cell_out_DF( self.currentCells, self.cellOutDF, self.prevtp )
		self.cellOutDF = newCellOutDF
		self.prevtp = self.tp.value()

		# before changing timepoint, print labeled cells and check if they are OK
		# print( 'cells labeled:\n ', self.currentCells )

		# print(self.fList['gfp'][self.tp.value()])
		tRow = self.timesDF.ix[ self.timesDF.tidxRel == self.tp.value() ].squeeze()

		### update the text of the fileName
		self.fName.setText( self.timesDF.ix[ self.timesDF.tidxRel == self.tp.value(), 'fName' ].values[0])

		print( 'Loading... ', self.pathDial, tRow.fName )

		# calculate the max value of the previous stack
		try:
			prevmax = np.max( [ np.max(self.stacks[ch]) for ch in self.channels ] )
		# if it's the first time a stack is to be loaded (so if there is no previous stack), set it to zero
		except:
			prevmax = 0

		# load all the available stacks - this is the slowest part of the code!!!
		self.stacks = {}
		for ch in self.channels:
			fileName = os.path.join( self.pathDial, tRow.fName + ch + '.tif')
			if os.path.isfile( fileName ):
				# print(MultiImage('X:\\Simone\\160129_MCHERRY_HLH2GFP_onHB101\\C02_analyzedImages\\Z003_488nm.tif'))
				# print(fileName, MultiImage( fileName ))
				# self.stacks[ch] = MultiImage( fileName )
				self.stacks[ch] = load_stack( fileName )
			# if there are no files for the timepoint, create a blank image
			else:
				self.stacks[ch] = prevmax * np.ones((self.nslices,self.cropsize,self.cropsize))

		### extract current cells already labeled
		if len(self.currentCells) > 0:
			print('previous cells:')
			print(self.currentCells)
			
		self.currentCells = extract_current_cell_out( self.cellPosDF, self.cellOutDF, self.tp.value(), self.zoom.value() )
		# print(self.currentCells)

		# if there are cells labeled and if the previously currently analyzed cell is present, set it as the currently labeled cell and select the right slice
		if len(self.currentCells) > 0:
			print('cells detected:')
			print(self.currentCells)

			### update currently analyzed cell
			if self.analyzedCell not in list( self.currentCells.cname ):
				self.analyzedCell = self.currentCells.cname[0]

			### update slice, if it's the same slice number, manually replot the images
			newslice = self.currentCells.ix[ self.currentCells.cname == self.analyzedCell, 'Z' ].values[0]
			if newslice != self.sl.value():
				self.sl.setValue( newslice )

			else:
				self.updateAllCanvas()
		# if no cells are found, manually plot the blank images
		elif len(self.currentCells) == 0:
			self.updateAllCanvas()

		# update the BC
		self.setBCslidersMinMax()            

	def saveData(self):
	    
		save_data_frame( self.cellOutDF, self.path, self.worm + '_05cellOut.pickle' )       
		self.setFocus()
	    
	def updateAllCanvas(self):
		self.updateCanvas1()
		self.updateCanvas2()
	    
	def radio488Clicked(self, enabled):
		# print('radio 488 clicked')

		if enabled:
			if '488nm' in self.channels:
				self.currentChannel = '488nm'
				self.setFocus()
				self.updateCanvas1()
			else:
				if self.currentChannel == 'CoolLED':
					self.CoolLEDBtn.setChecked(True)    # this uppdates the canvas1 once more
				elif self.currentChannel == '561nm':
					self._561nmBtn.setChecked(True)    # this uppdates the canvas1 once more
				QtGui.QMessageBox.about(self, 'Warning', 'No 488nm channel!')

	def radio561Clicked(self, enabled):
	    # print('radio 561 clicked')

		if enabled:
			if '561nm' in self.channels:
				self.currentChannel = '561nm'
				self.setFocus()
				self.updateCanvas1()
			else:
				if self.currentChannel == 'CoolLED':
					self.CoolLEDBtn.setChecked(True)    # this uppdates the canvas1 once more
				elif self.currentChannel == '488nm':
					self._488nmBtn.setChecked(True)    # this uppdates the canvas1 once more
				QtGui.QMessageBox.about(self, 'Warning', 'No 561nm channel!')

	def radioCoolLEDClicked(self, enabled):
	    # print('radio LED clicked')

		if enabled:
			if 'CoolLED' in self.channels:
				self.currentChannel = 'CoolLED'
				self.setFocus()
				self.updateCanvas1()
			else:
				if self.currentChannel == '561nm':
					self._561nmBtn.setChecked(True)    # this uppdates the canvas1 once more
				elif self.currentChannel == '488nm':
					self._488nmBtn.setChecked(True)    # this uppdates the canvas1 once more
				QtGui.QMessageBox.about(self, 'Warning', 'No CoolLED channel!')

	def updateBC(self):
		# change brightness and contrast
		self.imgplot1.set_clim( self.sld1.value(), self.sld2.value() )  
		self.imgplot2.set_clim( self.sld1.value(), self.sld2.value() )  
		self.canvas1.draw()
		self.canvas2.draw()

	#-----------------------------------------------------------------------------------------------
	# DEFAULT FUNCTION FOR KEY AND MOUSE PRESS ON WINDOW
	#-----------------------------------------------------------------------------------------------

	def keyPressEvent(self, event):

		# print(event.key())

		# change timepoint
		if event.key() == QtCore.Qt.Key_Right:
			cellsLabeled = 0
			nextTp = self.tp.value()
			while cellsLabeled == 0:
				nextTp += 1
				cellsLabeled = len(extract_current_cell_pos(self.cellPosDF,nextTp))
			dt = nextTp - self.tp.value()
			self.changeSpaceTime( 'time', dt )

		elif event.key() == QtCore.Qt.Key_Left:
			cellsLabeled = 0
			prevTp = self.tp.value()
			while cellsLabeled == 0:
				prevTp -= 1
				cellsLabeled = len(extract_current_cell_pos(self.cellPosDF,prevTp))
			dt = prevTp - self.tp.value()
			self.changeSpaceTime( 'time', dt )

		# change slice
		elif event.key() == QtCore.Qt.Key_Up:
			self.changeSpaceTime( 'space', +1 )
		    
		elif event.key() == QtCore.Qt.Key_Down:
			self.changeSpaceTime( 'space', -1 )

		elif event.key() == QtCore.Qt.Key_Space:

			if len( self.currentCells ) > 0:
				idx = np.mod( self.currentCells.ix[ self.currentCells.cname == self.analyzedCell ].index + 1, len(self.currentCells) )
				self.analyzedCell = self.currentCells.cname.values[idx][0]
				self.sl.setValue( self.currentCells.ix[ self.currentCells.cname == self.analyzedCell, 'Z' ] )

				self.updateAllCanvas()

		self.setFocus()

	def wheelEvent(self,event):
		if self.canvas1.underMouse():
			step = event.step
		else:          
			step = event.delta()/abs(event.delta())
		self.sl.setValue( self.sl.value() + step) 

	#-----------------------------------------------------------------------------------------------
	# ADDITIONAL FUNCTIONS FOR KEY AND MOUSE PRESS ON CANVASES AND BUTTONS
	#-----------------------------------------------------------------------------------------------

	def onMouseClickOnCanvas1(self, event):

		pos = np.array( [ int(np.round(event.xdata)), int(np.round(event.ydata)) ] )

		outline = extract_out( self.currentCells.ix[ self.currentCells.cname == self.analyzedCell ].squeeze() )

		if event.button == 1:

			if np.isnan( outline[0,0] ):
				outline = np.array( [ pos ] )
			else:
				outline = np.vstack( [ outline, pos ] )

		elif event.button == 3:

			if len( outline[ :, 0 ] ) == 1:
				outline = np.array( [ [ np.nan, np.nan ] ] )
			else:
				outline = outline[:-1,:]

		idx = self.currentCells.ix[ self.currentCells.cname == self.analyzedCell ].index

		self.currentCells.Xout.values[ idx ] = [ outline[:,0] ]
		self.currentCells.Yout.values[ idx ] = [ outline[:,1] ]

		self.updateCanvas1()
		self.setFocus()
		# print(event.button,event.xdata,event.ydata)

	def automaticOutline(self):

		imgpxl = self.zoom.value()
		# if not np.isnan( self.currentCells.ix[ self.currentCells.cname == self.analyzedCell ].imgPxl.values[0] ):
		# 	imgpxl = self.currentCells.ix[ self.currentCells.cname == self.analyzedCell ].imgPxl

		# extract current cell data
		pos = extract_3Dpos( self.currentCells.ix[ self.currentCells.cname == self.analyzedCell ].squeeze() )
		# crop image
		cropimg = self.stacks[self.currentChannel][self.sl.value(),pos[1]-imgpxl/2:pos[1]+imgpxl/2+1,pos[0]-imgpxl/2:pos[0]+imgpxl/2+1]

		# create first mask
		mask = cropimg >= threshold_otsu(cropimg)

		if mask[int(mask.shape[0]/2),int(mask.shape[1]/2)] != 0:
			# adjust the mask to be smooth
			mask = binary_fill_holes( mask )
			mask = morphology.remove_small_objects( mask, 100 )
			mask = morphology.binary_dilation( cropimg >= threshold_otsu(cropimg), selem = np.ones((3,3)) )
			lbl = morphology.label(mask)
			nlbl = lbl[int(mask.shape[0]/2),int(mask.shape[1]/2)]
			mask = lbl==nlbl

			# find vertices of the polygon
			cont = cv2.findContours(mask.astype(np.uint8), cv2.RETR_TREE, cv2.CHAIN_APPROX_SIMPLE)[1][0]
			cont = cv2.convexHull(cont)
			cont = cont[:,0,:]

			vertices = [list(cont[:,0]),list(cont[:,1])]
			vertices[0].append(vertices[0][0])
			vertices[1].append(vertices[1][0])
			outline = np.array(vertices)

			# save them in the right cell data
			idx = self.currentCells.ix[ self.currentCells.cname == self.analyzedCell ].index
			# print(outline.shape, imgpxl)
			self.currentCells.Xout.values[ idx ] = [ outline[0]/imgpxl*1000. ]
			self.currentCells.Yout.values[ idx ] = [ outline[1]/imgpxl*1000. ]

		self.updateCanvas1()
		self.setFocus()

	def fullAutomaticOutline( self ):

		if self.CoolLEDBtn.isChecked():
			QtGui.QMessageBox.about(self, 'Warning', 'Select a proper fluorescence channel!')
			return

		firsttp = first_tidx_pos_all_cells( self.cellPosDF )
		lasttp = last_tidx_pos_all_cells( self.cellPosDF )

		for idx, trow in self.timesDF.ix[ ( self.timesDF.tidxRel <= lasttp ) & ( self.timesDF.tidxRel >= firsttp ) ].iterrows():
			self.tp.setValue(trow.tidxRel)

			for jdx, cell in self.currentCells.iterrows():
				self.analyzedCell = cell.cname
				if cell.cname[0] != 'b':
					self.automaticOutline()

		self.tp.setValue( firsttp )

	#-----------------------------------------------------------------------------------------------
	# UTILS
	#-----------------------------------------------------------------------------------------------

	def changeZoom( self ):
		imgpxl = self.zoom.value()
		self.currentCells.ix[ self.currentCells.cname == self.analyzedCell, 'imgPxl' ] = imgpxl
		self.updateCanvas1()

	def setBCslidersMinMax(self):
		self.sld1.setMaximum( np.max( [ np.max(self.stacks[ch]) for ch in self.channels ] ) )
		self.sld1.setMinimum(0)
		self.sld2.setMaximum (np.max( [ np.max(self.stacks[ch]) for ch in self.channels ] ) )
		self.sld2.setMinimum(0)

	def initializeCanvas1(self):
		# print('initializing canvas1')

		self.fig1.clf()
		self.fig1.subplots_adjust(left=0., right=1., top=1., bottom=0.)
		self.ax1 = self.fig1.add_subplot(111)
		self.canvas1.draw()

		# plot the first blank image with the right size
		self.ax1.cla()
		self.imgplot1 = self.ax1.imshow( np.zeros((1000,1000)), cmap = 'gray', interpolation = 'nearest' )

		# remove the white borders
		self.ax1.autoscale(False)
		self.ax1.axis('Off')
		self.fig1.subplots_adjust(left=0., right=1., top=1., bottom=0.)

		# plot cell pos and name
		self.text_c1 = []
		self.plot_c1 = []

		# redraw the canvas
		self.canvas1.draw()
		self.setFocus()

	def updateCanvas1(self):
		# print('updating canvas1')

		# clear cell text and points
		# print(self.text1,self.points1)
		for text in self.text_c1:
			text.remove()
		self.text_c1 = []

		for points in self.plot_c1:
			self.ax1.lines.remove(points)
		self.plot_c1 = []


		# if no cells labeled, leave the image blank
		if len( self.currentCells ) == 0:
			self.imgplot1.set_data( np.ones((self.zoom.value(),self.zoom.value())) )
			self.canvas1.draw()
			return

		# detect and update zoom size
		imgpxl = self.zoom.value()
		if not np.isnan( self.currentCells.ix[ self.currentCells.cname == self.analyzedCell, 'imgPxl' ].values[0] ):
			imgpxl = self.currentCells.ix[ self.currentCells.cname == self.analyzedCell, 'imgPxl' ].values[0]
			self.zoom.setValue( imgpxl )
		else:
			self.currentCells.ix[ self.currentCells.cname == self.analyzedCell, 'imgPxl' ] = imgpxl

		# extract current cell data
		pos = extract_3Dpos( self.currentCells.ix[ self.currentCells.cname == self.analyzedCell ].squeeze() )

		# plot the image
		self.imgplot1.set_data( self.stacks[self.currentChannel][self.sl.value(),pos[1]-imgpxl/2:pos[1]+imgpxl/2+1,pos[0]-imgpxl/2:pos[0]+imgpxl/2+1] )

		# change brightness and contrast
		self.imgplot1.set_clim( self.sld1.value(), self.sld2.value() )  

		# print cell name
		if pos[2] == self.sl.value():
			self.text_c1.append( self.ax1.text( 10, 50, self.analyzedCell, color='yellow', size='medium', alpha=.8,
					rotation=0, fontsize = 20 ) )

		### draw outline
		outline = extract_out( self.currentCells.ix[ self.currentCells.cname == self.analyzedCell ].squeeze() )

		# print(self.currentCells.ix[ self.currentCells.cname == self.analyzedCell ].squeeze())

		if len( outline ) > 1:
			outline = np.vstack( [ outline, outline[0] ] )

		self.plot_c1.append( self.ax1.plot( outline[:,0], outline[:,1], '-x', color='yellow', ms=2, alpha=1., lw = .5 )[0] )

		# # redraw the canvas
		self.canvas1.draw()
		self.setFocus()

	def initializeCanvas2(self):
		# print('initializing canvas2')

		self.fig2.clf()
		self.fig2.subplots_adjust(left=0., right=1., top=1., bottom=0.)
		self.ax2 = self.fig2.add_subplot(111)
		self.canvas2.draw()

		# plot the first blank image with the right size
		self.ax2.cla()
		self.imgplot2 = self.ax2.imshow( np.zeros((self.cropsize,self.cropsize)), cmap = 'gray' )

		# remove the white borders
		self.ax2.autoscale(False)
		self.ax2.axis('Off')
		self.fig2.subplots_adjust(left=0., right=1., top=1., bottom=0.)

		# plot cell pos and name
		self.text_c2 = []
		self.plot_c2 = []

		# redraw the canvas
		self.canvas2.draw()
		self.setFocus()

	def updateCanvas2(self):
		# print('updating canvas2')

		# plot the image
		self.imgplot2.set_data( self.stacks[self.currentChannel][self.sl.value()] )

		# change brightness and contrast
		self.imgplot2.set_clim( self.sld1.value(), self.sld2.value() )  

		# clear cell text and points
		# print(self.text1,self.points1)
		for text in self.text_c2:
			text.remove()
		self.text_c2 = []

		for points in self.plot_c2:
			self.ax2.lines.remove(points)
		self.plot_c2 = []

		# extract current cell data
		for idx, cell in self.currentCells.iterrows():

			if cell.Z == self.sl.value():

				color = 'red'
				if cell.cname == self.analyzedCell:
					color = 'yellow'

				self.text_c2.append( self.ax2.text( cell.X, cell.Y + 18, cell.cname, color=color, size='medium', alpha=.8,
						rotation=0 ) )
				self.plot_c2.append( self.ax2.plot( cell.X, cell.Y, 'o', color=color, alpha = .8, mew = 0 )[0] )


		# redraw the canvas
		self.canvas2.draw()
		self.setFocus()

	def changeSpaceTime(self, whatToChange, increment):

		if whatToChange == 'time':
			# if they are OK (and not going to negative times), change timepoint
			self.tp.setValue( self.tp.value() + increment )

		if whatToChange == 'space':
			self.sl.setValue( self.sl.value() + increment )

if __name__ == '__main__':
    
    app = QtGui.QApplication.instance() # checks if QApplication already exists 
    if not app: # create QApplication if it doesnt exist 
        app = QtGui.QApplication(sys.argv)
    
    gui = GUI()
    app.setStyle("plastique")
    # app.installEventFilter(gui)
    sys.exit(app.exec_())
    


