# apphotRect.py

# WIC 2017-04-16: use photutils to implement our rectangular aperture
# photometry for resolved objects

import numpy as np
import os
import time
import matplotlib.pylab as plt
from matplotlib.colors import LogNorm
plt.ion()

from astropy.table import Table
from astropy.io import fits

class ApertureSet(object):

    """Set of apertures for the photometry"""

    # (I'm not sure if I want to stick with photutils, since that
    # currently enforces the same aperture size for all
    # apertures. Either way, we separate the apertures list out from
    # the individual image.

    def __init__(self):

        # Some parameters for aperture construction
        self.rectWidth = 50.
        self.rectHeight = 20.

        # sky rectangles
        self.skyWidth = np.copy(self.rectWidth)
        self.skyHeight = np.copy(self.rectHeight)

        # some default positions for apertures
        self.aperCens = np.array([])
        self.skyCens = np.array([])

        # INCLUDE ENTRIES FOR THE SKY REGIONS!!!

    def setDefaultCenters(self):

        """Sets up default aperture centers for testing on the
        2017-04-14 Jupiter data"""

        # the third one is a dummy to help me get the array dimensions
        # the right way round...
        vX = np.array([400., 400., 250.]) 
        vY = np.array([50., 70., 70.])]

        self.aperCens = np.vstack(( vX, vY ))

        # set up the sky apertures
        self.skyCens = np.copy(self.aperCens)

        # DO THE OFFSET HERE.

        

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

        # photometry table
        self.tPhot = Table()

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

    def showApertures(self, useLog=True, figNum=1, colormap='viridis'):

        """Utility method - plots the image and shows the apertures"""

        # do nothing if the image is not defined
        if np.shape(self.aImg) < 1:
            return

        # ensure the colormap is set
        try:
            cmap = plt.cm.get_cmap('gray')
        except:
            cmap = plt.cm.get_cmap(colormap)
        
        plt.figure(figNum)
        plt.clf()
        if not useLog:
            plt.imshow(self.aImg, origin='lower', interpolation='nearest', \
                           cmap=cmap)
        else:
            plt.imshow(self.aImg, origin='lower', interpolation='nearest', \
                           cmap=cmap, norm=LogNorm())

        # show colorbar
        plt.colorbar()
