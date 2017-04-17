# apphotRect.py

# WIC 2017-04-16: use photutils to implement our rectangular aperture
# photometry for resolved objects

import numpy as np
import os, glob
import sys 
import time
import matplotlib.pylab as plt
from matplotlib.colors import LogNorm
plt.ion()

from astropy.table import Table, vstack
from astropy.io import fits

# photutils
from photutils import RectangularAperture
from photutils import aperture_photometry

class ApertureSet(object):

    """Set of apertures for the photometry"""

    # (I'm not sure if I want to stick with photutils, since that
    # currently enforces the same aperture size for all
    # apertures. Either way, we separate the apertures list out from
    # the individual image.

    def __init__(self):

        # Some parameters for aperture construction
        self.rectWidth = 220.
        self.rectHeight = 90.
        self.rectTheta = np.radians(325.)
        

        # sky rectangles
        self.skyWidth = np.copy(self.rectWidth)
        self.skyHeight = np.copy(self.rectHeight)

        # sky nudge
        self.skyNudge = np.array([0., 30.])

        # some default positions for apertures
        self.aperCens = np.array([])
        self.skyCens = np.array([])

        # INCLUDE ENTRIES FOR THE SKY REGIONS!!!

    def setDefaultCenters(self):

        """Sets up default aperture centers for testing on the
        2017-04-14 Jupiter data"""

        # the third one is a dummy to help me get the array dimensions
        # the right way round...
        vX = np.array([944., 892., 1105., 693., 1297.])
        vY = np.array([520., 446., 365., 592., 250.5])

        self.aperCens = np.vstack(( vX, vY ))

        # set up the sky apertures
        self.skyCens = np.copy(self.aperCens)

        # DO THE OFFSET HERE.
        self.skyCens[1] -= ( self.skyNudge[1] + self.skyHeight )
        self.skyCens[1,0] += 2.0 * (self.skyNudge[1] + self.skyHeight)
        
    def buildApertures(self):

        """Builds the apertures from the centers"""

        if np.size(self.aperCens) < 1:
            return

        self.apersObj = RectangularAperture(positions=self.aperCens, \
                                                w=self.rectWidth, \
                                                h=self.rectHeight, \
                                                theta=self.rectTheta)

        self.apersSky = RectangularAperture(positions=self.skyCens, \
                                                w=self.rectWidth, \
                                                h=self.rectHeight, \
                                                theta=self.rectTheta)

    def showApertures(self, figNum=1, inAx=None):

        """Uses photutils apertures built-in to show the apertures"""

        # This is a really stupid cargo-cult way to proceed... Come
        # back to this later!!

        if not inAx:
            fig = plt.figure(figNum)
            fig.clf()
            ax = fig.add_subplot(111)
            ax.plot([0., 1024], [0., 1024], 'w.', alpha=0.)
        else:
            ax = inAx

        self.apersObj.plot(ax=ax, color='g', alpha=0.7)
        self.apersSky.plot(ax=ax, color='b', ls='--', alpha=0.7)

class OneImage(object):

    """Object holding a single image, its header, and photometry on
    that image. Inherits apertures from the call. Optionally, can pass
    the image data as an array and the header as a list."""

    def __init__(self, inPath='', apertures=[], \
                     apPath = '', \
                     aImg=np.array([]), hdr=[], \
                     iExt=0, \
                     Verbose=True):

        # path to the image
        self.inPath = inPath[:]
        
        # path to the apertures if given
        self.apPath = apPath[:]

        # image data and header
        self.aImg = aImg
        self.hdr = hdr
        
        # Which FITS extension has the data we want?
        self.iExt = iExt

        # source and sky apertures
        self.apsSource = []
        self.apsSky = []

        # photometry table
        self.tPhotObj = Table()
        self.tPhotSky = Table()

        # some header keywords
        self.keyDateNum= 'JD'
        self.keyObsdate= 'DATE-OBS'
        self.keyTemp= 'CCD-TEMP'
        self.keyItime='EXPTIME'
        self.keyGain='EGAIN'
        self.keyFilter='FILTER'

        # meta information
        self.DMeta = {}

        # Control variables
        self.Verbose = Verbose

    def loadImage(self):

        """Loads the image to memory"""

        if not os.access(self.inPath, os.R_OK):
            if self.Verbose:
                print "OneImage.loadImage FATAL - cannot read path %s" \
                    % (self.inPath)
            return

        # if here, then we can proceed
        self.aImg, self.hdr = fits.getdata(self.inPath, self.iExt, header=True)

    def doPhotom(self):

        """Extracts the photometry on the sky and on the object"""

        # might put try/except clause in here
        self.tPhotObj = aperture_photometry(self.aImg, self.apsSource, \
                                                method='subpixel')
        self.tPhotSky = aperture_photometry(self.aImg, self.apsSky, \
                                                method='subpixel')

    def pollHeader(self):

        """Polls the header for useful information (like observation
        time, etc.)"""

        # This is going to be a bit tedious...
        self.DMeta = {}
        for sKey in [self.keyDateNum, self.keyObsdate, \
                         self.keyTemp, self.keyItime, \
                         self.keyGain, self.keyFilter]:
            try:
                self.DMeta[sKey] = self.hdr[sKey]
            except:
                badKey=True

    def attachMetaToPhot(self):

        """Attaches metadata to the photometry table"""

        for sKey in self.DMeta.keys():
            self.tPhotObj.meta[sKey] = self.DMeta[sKey]

    def arrangePhotAsTableRow(self):

        """The calling routine may want to treat this set as a table
        row. This method re-arranges the photom and sky tables into a
        table with one column per object."""

        self.tRow = Table()
        
        # pass the meta information
        for sMeta in [self.keyDateNum, self.keyItime, \
                          self.keyGain, self.keyTemp]:
            self.tRow[sMeta] = [np.float(self.DMeta[sMeta])]

        # now add the photometry
        for iRow in range(len(self.tPhotObj)):
            sCount = str(iRow).zfill(2)

            sPho = 'src_%s' % (sCount)
            sSky = 'sky_%s' % (sCount)
            sAreaSrc = 'srcArea_%s' % (sCount)
            sAreaSky = 'skyArea_%s' % (sCount)

            # now put the values in
            self.tRow[sPho] = [self.tPhotObj[iRow]['aperture_sum'] ]
            self.tRow[sSky] = [self.tPhotSky[iRow]['aperture_sum'] ]
            self.tRow[sAreaSrc] = [self.apsSource.area() ]
            self.tRow[sAreaSky] = [self.apsSky.area() ]

    def showImage(self, useLog=True, figNum=1, colormap='gray'):

        """Utility method - plots the image and shows the apertures"""

        # do nothing if the image is not defined
        if np.shape(self.aImg) < 1:
            return

        # ensure the colormap is set
        try:
            cmap = plt.cm.get_cmap(colormap)
        except:
            cmap = plt.cm.get_cmap('gray')
        
        plt.figure(figNum)
        plt.clf()
        if not useLog:
            plt.imshow(self.aImg, origin='lower', interpolation='nearest', \
                           cmap=cmap)
        else:
            plt.imshow(self.aImg, origin='lower', interpolation='nearest', \
                           cmap=cmap, norm=LogNorm(), vmin=900., vmax=3000.)

        plt.xlabel(r"$X_{\rm image}$")
        plt.ylabel(r"$Y_{\rm image}$")

        plt.xlim(500., 1400.)
        plt.ylim(50., 750.)

        # show colorbar
        #plt.colorbar()

def TestDefaultApertures(Verbose=True, filOut='testPhot.csv', \
                             writeInterval=50):

    """Tests setting up the source apertures and offsetting them for
    the sky"""

    # Testing on laptop: run from directory
    # /Users/clarkson/Data/UMD_Observatory/20170408/Aligned_R_Good

    AS = ApertureSet()
    AS.setDefaultCenters()
    AS.buildApertures()

    # construct table to hold all the photometry
    tMaster = Table()

    # loop through all the images
    # dirIms = '/Users/clarkson/Data/UMD_Observatory/20170408/Aligned_R_Good'
    dirIms = os.getcwd()
    LIms = glob.glob('%s/Aligned*fits' % (dirIms))

    if len(LIms) < 1:
        print "TestDefaultApertures WARN - no images in %s"
        print "TestDefaultApertures WARN - exitting."
        return

    # directory for plots
    dirPlots = './tmpPlots'
    if not os.access(dirPlots, os.R_OK):
        os.makedirs(dirPlots)

    # For testing the timing...
    t0 = time.time()

    for iImg in range(len(LIms)):
        pathImg = LIms[iImg]
        fileImg = os.path.split(pathImg)[-1]
        fileStem = os.path.splitext(fileImg)[0]

        IM = OneImage(pathImg)
        IM.loadImage()
        IM.pollHeader()

        # Attach the aperture to the image...
        IM.apsSource = AS.apersObj
        IM.apsSky = AS.apersSky

        # ... extract the photometry for this image...
        IM.doPhotom()
        # IM.attachMetaToPhot()

        # create a table out of the data we've just extracted
        IM.arrangePhotAsTableRow()

        # now accumulate onto the master photometry table
        if len(tMaster) < 1:
            tMaster = IM.tRow.copy()
        else:
            tMaster = vstack([tMaster, IM.tRow])

        # write every so often (we don't want to have to re-run a long
        # set of measurements twice).
        if iImg % writeInterval < 1 and iImg > 1:
            tMaster.write(filOut, overwrite=True)


        if Verbose:
            nImgs = np.float(len(LIms))
            dt = time.time() - t0
            
            timePerRow = dt / np.float(iImg+1.)
            timeLeft = timePerRow * np.float(nImgs - iImg)

            sys.stdout.write("\r %4i of %4i:: %5.2e per row, about %.2fs remain" % \
                                 (iImg, nImgs, timePerRow, timeLeft) )
            sys.stdout.flush()

        # plot the image, then plot the apertures over it
        IM.showImage()
        ax = plt.gca()
        AS.showApertures(inAx=ax)

        # Since we're doing this as a set, we know what the first time
        # will be for the plot...
        daysElapsed = IM.tRow['JD'] - tMaster['JD'][0]
        sTime = '%.3f' % (daysElapsed * 1440.)
        plt.title('Time elapsed: %s min' % (sTime))

        figName = 'APs_%s.jpg' % (fileStem)
        pathFig = '%s/%s' % (dirPlots, figName)
        plt.savefig(pathFig)

    tMaster.write(filOut, overwrite=True)

    print tMaster
