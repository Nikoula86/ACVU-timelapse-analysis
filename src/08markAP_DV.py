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
use_pyside = qt_compat.QT_API == qt_compat.QT_API_PYSIDE

class GUI(QtGui.QWidget):
    
    def __init__(self):

        super(GUI, self).__init__()
        
        self.setWindowTitle( 'AP-DV Analysis' )
        self.lbltxt = '"wheel" press: change side, currently %s\n"i" or "u" press: change cell sides'
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
        
        rawDataBox.addLayout(Col1)
        rawDataBox.addLayout(Col2)
        
        self.setLayout(mainWindow)

        # DEFINE ALL WIDGETS AND BUTTONS
        
        loadBtn = QtGui.QPushButton('Load DataSet')
        saveBtn = QtGui.QPushButton('Save data (F12)')
        
        tpLbl = QtGui.QLabel('Relative Tp:')
        fNameLbl = QtGui.QLabel('File name:')
        
        self.tp = QtGui.QSpinBox(self)
        self.tp.setValue(0)
        self.tp.setMinimum(-100000)
        self.tp.setMaximum(100000)

        self.fName = QtGui.QLabel('')
        
        self._488nmBtn = QtGui.QRadioButton('488nm')
        self._561nmBtn = QtGui.QRadioButton('561nm')
        self.CoolLEDBtn = QtGui.QRadioButton('CoolLED')
        
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

        # PLACE ALL THE WIDGET ACCORDING TO THE GRIDS

        fileBox.addWidget(loadBtn)
        fileBox.addWidget(saveBtn)

        spaceBox1.addWidget(self.HLine())

        Col1.addWidget(tpLbl, 0, 0)#, 1, 1, Qt.AlignTop)
        Col1.addWidget(self.tp, 0, 1)#, 1, 1, Qt.AlignTop)
        Col1.addWidget(fNameLbl, 1, 0)#, 1, 1, Qt.AlignTop)
        Col1.addWidget(self.fName, 1, 1)#, 1, 1, Qt.AlignTop)
        Col1.addWidget(self._488nmBtn, 2, 0 )
        Col1.addWidget(self._561nmBtn, 3, 0 )
        Col1.addWidget(self.CoolLEDBtn, 4, 0 )
        
        Col2.addWidget(self.sld1)
        Col2.addWidget(self.sld2)
        Col2.addWidget(self.canvas1)
        
        self.setFocus()
        self.show()
        
        # BIND BUTTONS TO FUNCTIONS
        
        loadBtn.clicked.connect(self.selectWorm)
        saveBtn.clicked.connect(self.saveData)

        self.tp.valueChanged.connect(self.updateAllCanvas)
        self.sld1.valueChanged.connect(self.updateAllCanvas)
        self.sld2.valueChanged.connect(self.updateAllCanvas)

        self._488nmBtn.toggled.connect(self.radioClicked)
        self._561nmBtn.toggled.connect(self.radioClicked)
        self.CoolLEDBtn.toggled.connect(self.radioClicked)

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
        self.pathDial = QtGui.QFileDialog.getExistingDirectory(self, 'Select a folder', 'X:/Simone/160226_lag2YFP_histmCherry')
        self.worm = self.pathDial.split("/")[-1].split('_')[0]
        self.path = os.path.dirname( self.pathDial )
        self.setWindowTitle('Outline - ' + self.pathDial.split("/")[-2] + ' - ' + self.pathDial.split("/")[-1][:3] )

        ### give error message if there is no CoolLED movie in the selected folder
        if not os.path.isfile( os.path.join( self.pathDial, 'CoolLED_movie.tif' ) ):
            QtGui.QMessageBox.about(self,'Warning!','There is no movie in this folder! Create a movie first!')
            return

        ### load all movies (without timestamps, we will add it later on)
        self.channels = {}
        
        if os.path.isfile( os.path.join( self.pathDial, '488nm_movie.tif' ) ):
            self.channels['488nm'] = load_stack( os.path.join( self.pathDial, '488nm_movie.tif') )
        
        if os.path.isfile( os.path.join( self.pathDial, '561nm_movie.tif') ):
            self.channels['561nm'] = load_stack( os.path.join( self.pathDial, '561nm_movie.tif') )
        
        if os.path.isfile( os.path.join( self.pathDial, 'CoolLED_movie.tif' ) ):
            self.channels['CoolLED'] = load_stack( os.path.join( self.pathDial, 'CoolLED_movie.tif' ) )

        self.currentChannel = 'CoolLED'
        
        ### load parameters and times dataframes
        self.paramsDF = load_data_frame_pandas( self.path, self.worm + '_01params.pickle' )
        self.timesDF = load_data_frame_pandas( self.path, self.worm + '_01times.pickle' )
        self.gpDF = load_data_frame_pandas( self.path, self.worm + '_02gonadPos.pickle' )
        self.cellPosDF = load_data_frame_pandas( self.path, self.worm + '_04cellPos.pickle' )
        
        # extract some info
        self.compression = self.paramsDF.compression
        self.hatchingtidx = int( self.paramsDF.tidxHatch )

        ### if the AP pickle file already exists, load it, otherwise create a blank one
        if os.path.isfile( os.path.join(self.path, self.worm + '_08apdvPos.pickle' ) ):
            print('Found an old')
            self.apdvPosDF = load_data_frame_pandas( self.path, self.worm + '_08apdvPos.pickle' )
        
        else:
            print('Creating new')
            self.apdvPosDF = create_apdv_pos( self.timesDF )

        ### extract current cells already labeled
        self.currentPos = extract_current_apdv_pos( self.apdvPosDF, first_tidx_pos_all_cells( self.cellPosDF ) )
        print(self.currentPos)

        ### set the timepoint to the hatching time
        self.tp.setMinimum(np.min(self.timesDF.tidxRel))
        self.tp.setMaximum(np.max(self.timesDF.tidxRel))
        self.tp.setValue( first_tidx_pos_all_cells( self.cellPosDF ) )
        self.tp.setMinimum( first_tidx_pos_all_cells( self.cellPosDF ) )
        self.tp.setMaximum( last_tidx_pos_all_cells( self.cellPosDF ) + 1 )
        
        ### update the text of the fileName
        self.fName.setText(self.timesDF.ix[self.timesDF.tidxRel == self.tp.value(), 'fName'].values[0])

        self.setFocus()

    def saveData(self):

        save_data_frame( self.apdvPosDF, self.path, self.worm + '_08apdvPos.pickle' )
        
    def updateAllCanvas(self):

        self.updateRadioBtn()
        self.updateCanvas1()
        
    def radioClicked(self):

        if self._488nmBtn.isChecked():
            if '488nm' in self.channels.keys():
                self.currentChannel = '488nm'
            else:
                QtGui.QMessageBox.about(self, 'Warning', 'No 488nm channel!')

        elif self._561nmBtn.isChecked():
            if '561nm' in self.channels.keys():
                self.currentChannel = '561nm'
            else:
                QtGui.QMessageBox.about(self, 'Warning', 'No 561nm channel!')

        elif self.CoolLEDBtn.isChecked():
            if 'CoolLED' in self.channels.keys():
                self.currentChannel = 'CoolLED'
            else:
                QtGui.QMessageBox.about(self, 'Warning', 'No CoolLED channel!')

        self.setBCslidersMinMax()
        self.resetBC()
        self.setFocus()
        self.updateAllCanvas()

    def keyPressEvent(self, event):

        # print(event.key())

        # change timepoint
        if event.key() == QtCore.Qt.Key_Right:
            self.changeSpaceTime( 'time', +1 )

        elif event.key() == QtCore.Qt.Key_Left:
            self.changeSpaceTime( 'time', -1 )

        # change slice
        elif event.key() == QtCore.Qt.Key_Up:
            self.changeSpaceTime( 'space', +1 )
            
        elif event.key() == QtCore.Qt.Key_Down:
            self.changeSpaceTime( 'space', -1 )
        
        # key press on cropped image
        if self.canvas1.underMouse():
            self.onKeyPressOnCanvas1(event)
            
        self.setFocus()

    #-----------------------------------------------------------------------------------------------
    # DEFAULT FUNCTION FOR KEY AND MOUSE PRESS ON WINDOW
    #-----------------------------------------------------------------------------------------------

    def wheelEvent(self,event):

        if self.canvas1.underMouse():
            step = event.step
        
        else:          
            step = event.delta()/abs(event.delta())

        ### update daytaframe with previously labeled cells
        newapdvDF = update_apdv_pos_DF( self.currentPos, self.apdvPosDF, self.tp.value() )
        self.apdvPosDF = newapdvDF
        print(self.currentPos) # print previously labeled positions

        ### extract current cells already labeled
        self.currentPos = extract_current_apdv_pos( self.apdvPosDF, self.tp.value() + step )

        self.tp.setValue( self.tp.value() + step ) 

        self.fName.setText( self.timesDF.ix[ self.timesDF.tidxRel == self.tp.value(), 'fName' ].values[0] )

    def changeSpaceTime(self, whatToChange, increment):

        ### update daytaframe with previously labeled cells
        newapdvDF = update_apdv_pos_DF( self.currentPos, self.apdvPosDF, self.tp.value() )
        self.apdvPosDF = newapdvDF
        print(self.currentPos) # print previously labeled positions

        ### extract current cells already labeled
        self.currentPos = extract_current_apdv_pos( self.apdvPosDF, self.tp.value() + increment )

        if whatToChange == 'time':
            # if they are OK (and not going to negative times), change timepoint
            self.tp.setValue( self.tp.value() + increment )

        if whatToChange == 'space':
            self.sl.setValue( self.sl.value() + increment )

        self.fName.setText( self.timesDF.ix[ self.timesDF.tidxRel == self.tp.value(), 'fName' ].values[0] )

    #-----------------------------------------------------------------------------------------------
    # ADDITIONAL FUNCTIONS FOR KEY AND MOUSE PRESS ON CANVASES
    #-----------------------------------------------------------------------------------------------

    def onKeyPressOnCanvas1(self, event):
        
        posNameList = [ QtCore.Qt.Key_A, QtCore.Qt.Key_P, QtCore.Qt.Key_D ]

        # find the position of the cursor relative to the image in pixel
        imgshape = self.channels[self.currentChannel][0].shape
        canshape = self.canvas1.size()
        cf = imgshape[0]/canshape.width()
        refpos = self.canvas1.mapFromGlobal(QtGui.QCursor.pos())
        refpos = np.array([ int( refpos.x() * cf ), int( refpos.y() * cf )]) * self.compression

        ### find the closest cell to the cursor
        idx = closer_2Dpos( refpos.astype( np.uint16 ), self.currentPos )

        ### assign the name to the cell
        if any( [ event.key() == cn for cn in posNameList ] ):

            self.currentPos.ix[ idx, 'pname' ] = QtGui.QKeySequence(event.key()).toString().lower()

        self.updateCanvas1()
        self.setFocus()

    def onMouseClickOnCanvas1(self, event):
        
        refpos = np.array( [ event.xdata, event.ydata ] ) * self.compression 

        if event.button == 1:

            # create an empty cell in the currentCells df: the only entries are tidx, xyzpos and cname
            newpos = create_single_apdv_pos( refpos.astype(np.uint16), self.tp.value() )
            ### assign automatic name
            if len( self.currentPos ) == 0: newpos.pname = 'a'
            elif len( self.currentPos ) == 1: newpos.pname = 'p'
            elif len( self.currentPos ) == 2: newpos.pname = 'd'
            self.currentPos = pd.concat( [ self.currentPos, newpos ] )
            
        elif event.button == 3:

            # remove a cell (the closest to the cursor at the moment of right-click)
            idx = closer_2Dpos( refpos.astype(np.uint16), self.currentPos )
            self.currentPos = self.currentPos.drop( [ idx ] )

        self.currentPos = self.currentPos.reset_index(drop=True)
        
        self.updateCanvas1()
        self.setFocus()

    #-----------------------------------------------------------------------------------------------
    # UTILS
    #-----------------------------------------------------------------------------------------------

    def updateRadioBtn(self):
        if self.currentChannel == '488nm':
            self._488nmBtn.setChecked(True)
        elif self.currentChannel == '561nm':
            self._561nmBtn.setChecked(True)
        elif self.currentChannel == 'CoolLED':
            self.CoolLEDBtn.setChecked(True)
        self.setFocus()

    def setBCslidersMinMax(self):
        self.sld1.setMaximum(np.max(self.channels[self.currentChannel]))
        self.sld1.setMinimum(np.min(self.channels[self.currentChannel]))
        self.sld2.setMaximum(np.max(self.channels[self.currentChannel]))
        self.sld2.setMinimum(np.min(self.channels[self.currentChannel]))

    def resetBC(self):
        self.sld1.setValue(np.min(self.channels[self.currentChannel]))
        self.sld2.setValue(np.max(self.channels[self.currentChannel]))
        
    def updateCanvas1(self):
        
        # plot the image
        self.ax1.cla()
        imgplot = self.ax1.imshow( self.channels[self.currentChannel][self.tp.value() + self.hatchingtidx], cmap = 'gray' )
        
        # remove the white borders and plot outline and spline
        self.ax1.autoscale(False)
        self.ax1.axis('Off')
        self.fig1.subplots_adjust(left=0., right=1., top=1., bottom=0.)

        # change brightness and contrast
        self.sld1.setValue(np.min([self.sld1.value(),self.sld2.value()]))
        self.sld2.setValue(np.max([self.sld1.value(),self.sld2.value()]))
        imgplot.set_clim(self.sld1.value(), self.sld2.value())  

        # print gonad position
        gonadPos = extract_pos( self.gpDF.ix[ self.gpDF.tidx == self.tp.value() ].squeeze() ) / self.compression
        if len( gonadPos.shape ) > 0:
            self.ax1.plot( gonadPos[0], gonadPos[1], 'o', color='green', ms=10, mew=0, alpha=.7, lw = 0 ) 

        # print cell positions
        cells = extract_current_cell_pos(self.cellPosDF, self.tp.value())
        for idx, cell in cells.iterrows():
            if cell.cname[0] != 'b': 
                pos = np.array( [ gonadPos[0]*self.compression-256+cell.X, gonadPos[1]*self.compression-256+cell.Y ] ) / self.compression
                self.ax1.plot( pos[0], pos[1], 'o', color = 'blue', ms = 5, mew = 0, alpha = .5 )

        # pritn apdv pos
        for idx, pos in self.currentPos.iterrows():
            p = extract_pos( pos ) / self.compression

            self.ax1.plot( p[0], p[1], 'o', color='white', ms=10, mew=0, alpha=.8, lw = 0 ) 
            self.ax1.text( p[0], p[1] + 20, pos.pname.capitalize(), color='white', size='medium', alpha=.8,
                    rotation=0 )

        ### find ecdysis timepoint
        ecd = np.loadtxt( open( os.path.join( self.path, 'skin.txt'), 'rb' ) )
        # load ecdysis data
        index = np.where( ecd[:,0] == float(self.worm[1:]) )
        mintp = np.min( [ i for i in ecd[index, 1:6][0][0] if i >= 0 ] )
        lethtidx = ecd[ index, 1:6 ][0][0]
        lethtidx = lethtidx[ lethtidx>=0 ]
        tpL2 = self.timesDF.ix[ self.timesDF.tidxRel == ( lethtidx[1] - mintp ), 'timesRel' ].values[0]
        # print time
        # print(self.timesDF.ix[ self.timesDF.tidxRel == self.tp.value(), 'timesRel' ])
        self.ax1.text( 5, 15, '%.2f' % ( self.timesDF.ix[ self.timesDF.tidxRel == self.tp.value(), 'timesRel' ].values[0] - tpL2 ), color = 'white' )     

        # redraw the canvas
        self.canvas1.draw()
        self.setFocus()

if __name__ == '__main__':
    
    app = QtGui.QApplication.instance() # checks if QApplication already exists 
    if not app: # create QApplication if it doesnt exist 
        app = QtGui.QApplication( sys.argv )
    
    gui = GUI()
    app.setStyle( "plastique" )
    # app.installEventFilter(gui)
    sys.exit( app.exec_() )
