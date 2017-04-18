#
# CheckLinearity.py
#

# Given a set of flatfield frames, estimate the linearity for each
# pixel

# WIC 2015-05-05 - started updating this to take a SECTION of the
# image (may be necessary to spatially clump if the signal is weak)

import os, glob, time, sys
import numpy as N
import pyfits

# useful for smoothing
from scipy import ndimage
from scipy import signal

import pickle

import FitAgainstExptime

from pylab import cm as CM

import pylab as P
P.ion()

# WIC Apr 2 2014: I've written this as an object+methods rather than a
# sequence of methods because it's clearer to me... plus this way all
# the pieces are collected in the object that we can tweak, rather
# than sending them through from one to the other.

# Object that holds a "stack" of flatfields
class FlatStack(object):

    """Stack of image-planes and associated methods"""

    def __init__(self, DirImgs='', sStart='Flat_', sExt='fits', \
                 ThreshBright=40000.0, MinTExp=0.00):

        # Source directory for fits frames
        self.DirImgs=DirImgs[:]
        if len(DirImgs) < 1:
            self.DirImgs=os.getcwd()

        # Search terms for flatfield frames
        self.sStart='Flat_'
        self.sExt = 'fit'

        self.sStart='Linearity_Dome_060V_VBand-00'
        self.sExt='fit'

        # Take from the argument
        self.sStart = sStart[:]
        self.sExt = sExt[:]

        # List of flatfield frames
        self.LFrames = []

        # actual data info
        self.NFrames = 0
        self.LHeaders = []
        self.vTimes = N.array([])
        self.aStack = N.array([])

        # similar vector for times
        self.aTimes = N.array([])
        self.vTheseTimes = N.array([])

        # Threshold for what we call a "high" and "low" count
        # value. We leave it up to the observer to take enough points
        # above and below this range.
        self.ThreshBright= N.copy(ThreshBright) # was 30000.0

        # Min Exposure time to accept
        self.MinTExp = N.copy(MinTExp)

        # It's useful to designate pixel-counts vectors for the
        # current pixel-stack of interest since we're likely to be
        # doing multiple operations on it
        self.ThisCol = 350 # 50 #295 # column-index of this pixel
        self.ThisRow = 275 # 129 # 190 # row-index of this pixel

        # 2015-05-05 new feature: group regions! Default is length-1
        # region in each direction.
        self.RegLenCol = 1
        self.RegLenRow = 1

        # 2015-05-06 store the rowmax and colmax as separate entities
        # so that other methods can use them. For example, could call
        # a method that populates the max values based on the array
        # dimensions
        self.ThisColMax = self.ThisCol + 1  
        self.ThisRowMax = self.ThisRow + 1

        # 2015-05-06 should think about how best to encode the
        # indices... could use a boolean plane-array that we can then
        # just read out from all the other parts we might
        # need. (E.g. when writing coefficients to output frames)

        self.vThisPixel = N.array([])
        self.gThisHi = N.array([]) # Selection of points >= threshbright
        self.gThisLo = N.array([]) # selection of poitns < threshbright

        # Fitted parameters for the currently-selected pixel
        self.ThisLinr = N.array([])
        self.ThisQuad = N.array([])

        # Fitted parameters for the deltas
        self.vParsAgainstLinear = N.array([])
        self.vParsDeficit = N.array([])

        # Finely-sampled exptime grid for plotting
        self.tFine = N.array([])

        # Array of gradients
        self.aLinr = N.array([])
        self.aQuad = N.array([])

        # Array of linear fit to delta vs linear:
        self.aLinrVsLinr = N.array([])

        # Array of coefficients for the thing we actually want to
        # measure: (Linear minus Observed) as a function of (Observed)
        self.aExpon = N.array([])  # Expfit coefficients
        self.aPokes = N.array([])  # curve outside nSigma for the first time
        self.aThresh = N.array([]) # When deficit > MinFracNonlin times counts

        # Vector of delta-counts
        self.vThisDeficit = N.array([])
        
        # Some characterizing information about the pixel:
        self.nSigLinear = 6.0 # 5.0 # 7.0 # range to use for judging which piece is linear
        self.ThisStdLinear = 0.0 # stddev of the linear part of the curve

        self.ThisLowThresh = 0.0
        self.ThisLowPoke = 0.0

        # What threshold nonlinearity do we use? Let's start with 1 percent
        self.MinFracNonlin = 0.02

        self.FileLinr  = 'EstLinear.pickle'
        self.FileQuad  = 'EstQuad.pickle'
        self.FileExpon = 'LinearityExpon.pickle'
        self.FileLL = 'LinearityVsLinear.pickle'
        self.FileThresh = 'LinearityThresholds.pickle'

        self.FileStack = 'StackUsed.pickle'

        # closeness for discrepancies
        self.nClose = 20
        self.nDiscrep = 20

        # predicted deficit-image
        self.aImgPred = N.array([])
        self.aImgSmooth = N.array([])

    def AssembleImageSet(self):

        """Wrapper for methods defined below: load images from disk and stack
into an (NFrames x nY x nX) array that we'll use to look for linearity."""

        self.FindSourceFrames()
#        print self.LFrames
        for iPlane in range(len(self.LFrames)):
            self.LoadSinglePlane(iPlane)

            # NOT YET
            sys.stdout.write("\rAssemble INFO: %4i of %4i" % (iPlane+1, len(self.LFrames) ))
            sys.stdout.flush()

        # Flush the buffer
        print " "

        # Now lift the exptimes from the headers we've read in
        self.AssembleExptimes()
        
        # Produce vector of finely-sampled exposure times
        self.tFine = N.linspace(start=0, stop=N.max(self.vTimes), num=100)

    def WriteImageSet(self):

        """Writes assembled image-set to disk"""
        
        pickle.dump(N.asarray(self.aStack, 'int32'), open(self.FileStack, 'w'))
        
    def FindSourceFrames(self):

        """Finds frames that match our search-term filenames 'self.sStart' and
        'self.sExt' """

        sGlobString = "%s/%s*.%s" % (self.DirImgs, self.sStart, self.sExt)
        self.LFrames = glob.glob(sGlobString)

    def LoadSinglePlane(self, iPlane=0):

        """Load a single image and pile it on to the stack"""

        # Only do something if the iPlane'th entry in the file list exists!
        try:
            ThisFits = self.LFrames[iPlane]
        except:
            return

        # file must be readable
        try:
            ThisImg = pyfits.getdata(ThisFits)
        except:
            return

        # header must be readable to get the exptime
        try:
            ThisHeader = pyfits.getheader(ThisFits)
        except:
            return

        # WIC 2015-04-06 - probably don't need this for the Dobsonian
        # ST402. Also watchout for the filename search string.

        # print ThisHeader['INSTRUME']

        # WIC: Evil bodge for single-linear region...

        DoingDark=True
        if ThisHeader['INSTRUME'].find('402') < 0:
            if not DoingDark:
                if ThisHeader['EXPTIME'] > 32.0 or \
                   ThisHeader['EXPTIME'] < 3.0:
                    return

        else:

            # For the 402, 0.04s is about the minimum. Rerun.
            if ThisHeader['EXPTIME'] < self.MinTExp:  # was 0.04!!!
                return

                #        print iPlane
                
        # create plane of exptimes as well
        aTimes = ThisImg * 0.0 + ThisHeader['EXPTIME']

        # now do the appending
        self.LHeaders.append(ThisHeader)

        if N.size(N.shape(self.aStack)) < 2:
            self.aStack = ThisImg
            self.aTimes = N.copy(aTimes)
        else:
            self.aStack = N.dstack(( self.aStack, ThisImg ))
            self.aTimes = N.dstack(( self.aTimes, aTimes  ))

        # If we've got this far, we increment the number of planes
        # we've got so far.
        self.NFrames += 1

    def AssembleExptimes(self):

        """Get the exposure times from the headers"""

        # By doing a second loop for the times like this, we can be
        # flexible about how we lift the integration time.        
        self.vTimes = N.zeros(self.NFrames)
        for iTime in range(self.NFrames):
            ThisHeader = self.LHeaders[iTime]
            self.vTimes[iTime] = N.float(ThisHeader['EXPTIME'])

    def RunStatsOnAllPixels(self, UseSimple=True):

        """Fit the growth-curve to all the pixels. Might take a while..."""

        ShapeStack = N.shape(self.aStack)
        nRows = ShapeStack[0]
        nCols = ShapeStack[1]

        # Initialise the results arrays
        self.aLinr = N.zeros((nRows, nCols, 2 ))
        self.aQuad = N.zeros((nRows, nCols, 3 ))

        self.aExpon = N.zeros(( nRows, nCols, 2 ))
        self.aPokes = N.zeros(( nRows, nCols ))
        self.aThresh = N.zeros(( nRows, nCols ))

        self.aLinrVsLinr = N.zeros(( nRows, nCols, 2 ))

        tStart = time.time()

        print nRows, nCols
        print ShapeStack

        iDone = 0

        # Now fill this in...
        for jRow in range(nRows):
            for iCol in range(nCols):
                
                # select this pixel                
                self.ThisCol = iCol
                self.ThisRow = jRow                

                # Initialise the pixel
                self.InitSinglePixel()
                self.SelectPixelVector()

                if UseSimple:
                    
                    # Create the object and pass in the arrays
                    PT = FitAgainstExptime.PixTrend(Verbose=False)

                    # Pass in any hard limits  # WARN HARDCODE 2015-04-30
                    #                    PT.tLimLo = 2.5
                    #                    PT.tLimHi = 15.

                    PT.nClose = self.nClose
                    PT.nDiscrep = self.nDiscrep

                    PT.vExptimes = N.copy(self.vTimes)
                    PT.vCounts = N.copy(self.vThisPixel)
                    PT.vIDS = N.zeros(N.size(self.vThisPixel))
                    PT.sWhichCol = '%i' % (self.ThisCol)
                    PT.sWhichRow = '%i' % (self.ThisRow)

                    # Clump, perform stats
                    PT.ClumpByExptime()
                    PT.StatsPerExptime()    
                    try:
                        PT.WrapFindPointsOnLine()
                    except:
                        print "---"
                        print "PROBLEM with col, row (%4i, %4i)" % (jRow, iCol)
                        continue
                    PT.FindDiscrepantChunks()
                    PT.EstLoIntersection()

                    # Translate into pieces this routine expects
                    self.aLinr[jRow, iCol] = N.copy(PT.parsRough)
                    self.aQuad[jRow, iCol] = N.copy(PT.parsRef)
#                    self.aQuad[jRow, iCol][0] = PT.covaRef[0][0]
#                    self.aQuad[jRow, iCol][1] = PT.covaRef[1][1]
#                    self.aQuad[jRow, iCol][2] = PT.xcmRef
                    self.aPokes[jRow, iCol] = N.copy(PT.valTimeCross)
                    self.aThresh[jRow, iCol] = N.copy(PT.xcmRef)

                    # All of the above is if we're using the simple method. 
                    # If not...
                else:

                    # Fit the curve of growth
                    self.FitGrowthSelected()
                    self.GetDeltaCounts()
                    self.FitDeltaCounts()
                    self.GetDeltaLimits()

                    # Fit against *incident* flux
                    self.FitDeltasAgainstLinear()

                    # Now store in the master-arrays
                    self.aLinr[jRow, iCol] = N.copy(self.ThisLinr)
                    self.aQuad[jRow, iCol] = N.copy(self.ThisQuad)
                    self.aExpon[jRow, iCol] = N.copy(self.vParsDeficit)
                    self.aPokes[jRow, iCol] = N.copy(self.ThisLowPoke)
                    self.aThresh[jRow, iCol] = N.copy(self.ThisLowThresh)

                    self.aLinrVsLinr[jRow, iCol] = \
                                                   N.copy(self.vParsAgainstLinear)

                # Whichever method we are using, will want to write
                # every so often.
                if iDone % 50000 == 20:
                    ThatTime = time.time()
                    print "%i of %i: %i %i: %.2f" % (iDone, nRows*nCols, jRow, iCol, ThatTime - tStart)
                    # dump here
                    self.WriteFittedStats()
                
                    # If we're using the simple method, plot graph of
                    # the pixels too.
                    if UseSimple:
                        PT.ShowPixel()

                iDone += 1

        # OK well that takes about a minute to go through, that's not
        # terrible...

        self.WriteFittedStats()

    def WriteFittedStats(self):

        # dump them both to disk:
        pickle.dump( self.aLinr, open(self.FileLinr, 'w' ))
        pickle.dump( self.aQuad, open(self.FileQuad, 'w' ))
        pickle.dump( self.aExpon, open(self.FileExpon, 'w'))
        pickle.dump( N.dstack((self.aPokes, self.aThresh)), open(self.FileThresh, 'w') )
        pickle.dump( self.aLinrVsLinr, open(self.FileLL, 'w'))

    def InitSinglePixel(self):

        """Re-initialises the single-pixel quantities."""

        self.vThisDeficit = N.array([])
        self.ThisStdLinear = 0.0
        self.ThisLowThresh = 0.0
        self.ThisLowPoke = 0.0
        self.vThisPixel = N.array([])
        self.gThisHi = N.array([]) 
        self.gThisLo = N.array([]) 

        # Fitted parameters
        self.vParsDeficit = N.array([])
        self.vParsAgainstLinear = N.array([])
        self.ThisLinr = N.array([])
        self.ThisQuad = N.array([])

    def SelectPixelVector(self):

        """Lifts the 1D counts vector for the [ThisCol'th, self.ThisRow] pixel
for all exposure-times."""

        # 2015-05-05 New feature: take self.RegLenCol x self.RegLenRow
        # region
        RowMin = self.ThisRow
        ColMin = self.ThisCol
        
        # Don't try to read past the image bounaries!
        shStack = N.shape(self.aStack)
        RowMax = N.min( [RowMin + self.RegLenRow, shStack[0]] )
        ColMax = N.min( [ColMin + self.RegLenCol, shStack[1]] )

        # ensure all are integers
        RowMax = int(RowMax)
        ColMax = int(ColMax)

        # Trust the user to have selected a sensible pixel... we'll
        # deal with that later.
        #self.vThisPixel = self.aStack[self.ThisRow, self.ThisCol]

        ArrayOfRegion = self.aStack[RowMin:RowMax, ColMin:ColMax]
        self.vThisPixel = N.ravel(ArrayOfRegion)

        # same for times
        self.vTheseTimes = N.ravel(self.aTimes[RowMin:RowMax, ColMin:ColMax])

        print "DEBUG INFO:", N.shape(self.vThisPixel)

        # Select "high" and "low" values (WIC - this might be obsolete)
        self.gThisHi = N.where( (self.vThisPixel >= self.ThreshBright) & \
                                (self.vThisPixel < 62000.0) )[0]
        self.gThisLo = N.where(self.vThisPixel < self.ThreshBright)[0]

        #        print N.size(self.vThisPixel), N.size(self.gThisHi)
        print "DEBUG INFO:", N.shape(self.gThisHi), N.shape(self.gThisLo)


    def WriteSelectedPixel(self):

        """Writes the selected pixel to disk"""

        sName = 'SinglePix_c%i_r%i_Data.pickle' % (self.ThisCol, self.ThisRow)
        DSingle={}
        DSingle['vThisPixel'] = N.copy(self.vThisPixel)
        DSingle['vTimes'] = N.copy(self.vTimes)
        DSingle['gThisHi'] = N.copy(self.gThisHi)
        DSingle['gThisLo'] = N.copy(self.gThisLo)
        DSingle['vTheseTimes'] = N.copy(self.vTheseTimes)
    

        pickle.dump(DSingle, open(sName, 'w'))
        

    def FitGrowthSelected(self):

        """Fits linear growth-curve to the pixel-values below threshold"""

        # Now includes safety check in case unpopulated row (e.g. edge
        # of detector)

        if N.size(self.gThisLo) < 4:
            self.ThisLinr=N.zeros(2)
            return

        self.ThisLinr = N.polyfit(self.vTimes[self.gThisLo], self.vThisPixel[self.gThisLo], 1)


        if N.size(self.gThisHi) < 4:
            self.ThisQuad=N.zeros(3)
            return

        self.ThisQuad = N.polyfit(self.vTimes[self.gThisHi], self.vThisPixel[self.gThisHi], 2)

    def GetDeltaCounts(self):

        """Fit what we actually want - the difference between linear and
        actual as a function of actual"""

        # Ensure the linear fit has been done!
        self.FitGrowthSelected()
        
        # Predict the counts under the linear model, find the difference (Observed minus predicted)
        vCountsPredLin = N.polyval(self.ThisLinr, self.vTimes)
        self.vThisDeficit= vCountsPredLin - self.vThisPixel

        self.GetRoughRangesLinear()

    def GetRoughRangesLinear(self):

        """Gets count ranges roughly corresponding to the linear regime"""

        self.ThisStdLinear = N.std(self.vThisDeficit[self.gThisLo])
        self.bLinear = (N.abs(self.vThisDeficit) < self.ThisStdLinear * self.nSigLinear)

    def FitDeltaCounts(self):

        """Fits the growth-curve of delta-counts vs counts"""

        self.b2Fit = (~self.bLinear) & (self.vThisDeficit > 0)

        # can't do anything if bad
        if N.size(N.where(self.b2Fit)[0]) < 2:
            self.vParsDeficit=N.zeros(2)
            return

        # Be a bit smarter about this... use all points brighter than
        # the dimmest of these
        CountsMin = N.min(self.vThisPixel[N.where(self.b2Fit)[0]])
#        self.g2Fit = N.where( (self.vThisPixel > CountsMin) )[0]
        self.b2Fit = (self.vThisPixel > CountsMin) # & (self.vThisPixel < 59500.0)

        # The following should work:
        x2Fit = self.vThisPixel[self.b2Fit]
        y2Fit = N.log(N.abs(self.vThisDeficit[self.b2Fit]))
        if N.size(y2Fit) < 3:
            self.vParsDeficit=N.zeros(2)
            return
        self.vParsDeficit = N.polyfit(x2Fit, y2Fit, 1)
        
    def FitDeltasAgainstLinear(self):

        """Fits the deltas against the linear prediction"""

        xFit = N.polyval(self.ThisLinr, self.vTimes[self.b2Fit])
        yFit = self.vThisDeficit[self.b2Fit]
        self.vParsAgainstLinear = N.zeros(2)
        # WIC 2015-04-06 - return if badpix
        if N.size(yFit) < 3:
            return
        self.vParsAgainstLinear = N.polyfit(xFit, yFit,1)

        # Do one pass of clipping
        yPredPass1 = N.polyval(self.vParsAgainstLinear, xFit)
        fMed = N.median(yPredPass1 - yFit)
        fStd = N.std(yPredPass1 - yFit)
        gLo = N.where(N.abs(yPredPass1 - yFit - fMed)/fStd < 2.0)[0]
        if N.size(gLo) > 5:
            self.vParsAgainstLinear = N.polyfit(xFit[gLo], yFit[gLo],1)

    def InitFineCountsVector(self):

        """Sets up fine-grained vector of counts, for assessing the curve later on."""

        self.vFineCounts = N.linspace(100,N.max(self.aStack),5000)

        # Take the natural log of the fine counts and the threshold so
        # we don't need to re-compute it each time:
        self.vFineLnCounts = N.log(self.vFineCounts)

    def GetDeltaLimits(self):

        """Get the limits on the counts-deficit"""

        lnYPred = N.polyval(self.vParsDeficit, self.vFineCounts)
        YPred = N.exp(lnYPred)

        # we look for two points in particular: First, where the curve
        # goes above the range of points we're not fitting, and,
        # second, where deficit / Y rises above a threshold fraction
        # for the first time.

        # Since both are single-valued, we can just find the closest
        # point.
        iPokesOut = N.argmin(N.abs(YPred - self.ThisStdLinear * self.nSigLinear))
        iAtThresh = N.argmin(N.abs(YPred/self.vFineCounts - self.MinFracNonlin))

        # Now return those values
        self.ThisLowPoke = self.vFineCounts[iPokesOut]
        self.ThisLowThresh = self.vFineCounts[iAtThresh]

    def LoadExpon(self):

        """Loads the set of predictions for the exponential model"""

        if not os.access(self.FileExpon, os.R_OK):
            print "LoadExpon FATAL - cannot load exponential predictions"
            return

        self.aExpon = pickle.load(open(self.FileExpon,'r'))
        
    def LoadVsLinear(self):
        
        """Loads the set of predictions against linear input"""

        self.aLinrVsLinr = pickle.load(open(self.FileLL,'r'))

    def PredictDeficitExpon(self, CountsIn=50000.0):

        """Predicts the counts-deficit due to nonlinearity under the
exponential model."""

        self.aImgPred = self.aExpon[:,:,0]*CountsIn + self.aExpon[:,:,1]
        self.aImgPred = N.exp(self.aImgPred)/CountsIn

        # Do in percentage
        self.aImgPred *= 100.0

    def PredictDeficitVsLinear(self, CountsLinear=50000.0):

        self.aImgPred = self.aLinrVsLinr[:,:,1] \
                        + self.aLinrVsLinr[:,:,0] * CountsLinear

        # Find the measured counts, and scale by this
        aMeas = N.zeros(N.shape(self.aImgPred))+CountsLinear - self.aImgPred

        self.aImgPred /= aMeas
        self.aImgPred *= 100.0

        gLow = N.where(self.aImgPred < 0)
        self.aImgPred[gLow] = 0.0

    def SmoothDeficitPred(self):

        """Smooths the predicted deficit"""

        self.aImgSmooth = \
                          ndimage.filters.gaussian_filter(self.aImgPred, \
                                                          5.0)


    def ShowSelectedPixel(self):

        """Plot the curve for a single pixel"""

        vPixvals = self.vThisPixel
        gLo = self.gThisLo
        gHi = self.gThisHi

#        print N.size(vPixvals), N.size(self.vTimes), N.size(gLo), N.size(gHi)
        P.figure(1)
        P.clf()
        P.plot(self.vTimes[gLo], vPixvals[gLo], 'ko')
        P.plot(self.vTimes[gHi], vPixvals[gHi], 'ks')
        P.grid()

        # overplot the fits if we have them
        print self.ThisLinr
        if N.size(self.ThisLinr) > 0:
            P.plot(self.tFine, N.polyval(self.ThisLinr, self.tFine), 'g--')

        if N.size(self.ThisQuad) > 0:
            yQuad = N.polyval(self.ThisQuad, self.tFine)
            gBri = N.where(yQuad >= self.ThreshBright)[0]
            if N.size(gBri) > 1:                              
                P.plot(self.tFine[gBri], yQuad[gBri], 'r')
        P.ylim(0,70000.0)

    def ShowDeltaVsLinear(self):

        """Shows the delta as a function of the linear fit"""

        # This is useful when predicting what the nonlinearity would
        # be as a function of the *incident* light (so that we can
        # check for spatial patterns decoupled from patterns in
        # illumination)

        vX = N.polyval(self.ThisLinr, self.vTimes)
        gLo = ~self.b2Fit
        gHi = self.b2Fit # was b2fit
        vY = self.vThisDeficit
        P.figure(2)
        P.clf()
        P.subplot(221)
        P.plot(vX[gLo], vY[gLo], 'bo')
        P.plot(vX[gHi], vY[gHi], 'rs')

        P.plot(vX[gHi], N.polyval(self.vParsAgainstLinear, vX[gHi]), 'r-')

        P.subplot(223)
        P.semilogy(vX[gLo], vY[gLo], 'bo')
        P.semilogy(vX[gHi], vY[gHi], 'rs')

    def ShowCorrectionCurve(self, DoLog=True):
        
        """Shows the (observed - predicted) vs (observed) curve for a single pixel"""

        vX = self.vThisPixel
        gLo = ~self.b2Fit
        gHi = self.b2Fit # was b2fit
        vY = self.vThisDeficit

        P.figure(1)
        P.clf()

        P.subplots_adjust(hspace=0.25, wspace=0.35)

        # yes I know, mission bloat... make a nice figure here anyway!
        P.subplot(222)
        if not DoLog:
            P.plot(self.vTimes[gLo], self.vThisPixel[gLo], 'bo')
        else:
            P.semilogx(self.vTimes[gLo], self.vThisPixel[gLo], 'bo')
        P.plot(self.vTimes[gHi], self.vThisPixel[gHi], 'rs')
        P.plot(self.tFine, N.polyval(self.ThisLinr, self.tFine), 'b')
        P.grid(which='both')
        P.xlabel('Exposure Time t_exp (sec)')
        P.ylabel('Recorded counts')
        P.title('Observations: Recorded counts vs Exposure time')

        # show the equation
        sInc = 'Incident counts = %.3f + %.3f t_exp' \
               % (self.ThisLinr[1], self.ThisLinr[0])
        P.annotate(sInc,(0.05,0.90), xycoords='axes fraction', color='b')

        # BULB ANNOTATION COMES HERE 2015-04-06
#        sBulb = 'Flatfield bulb setting: 60V'
#        P.annotate(sBulb, (0.95, 0.07), color='k', xycoords='axes fraction', \
#                   horizontalalignment='right')

        yCopy = N.copy(P.ylim())

        if yCopy[-1] > self.ThreshBright:
            P.ylim(0, yCopy[-1])

        # Do this both ways
        P.subplot(224)
        P.semilogy(vX[gLo], N.abs(vY[gLo]), 'ko')
        P.semilogy(vX[gHi], N.abs(vY[gHi]), 'rs')
        P.grid()

        # Get the ranges of points that fall within the linear regime
        # (for our fitting)
        xLims = N.copy(P.xlim())
        Band = self.ThisStdLinear * self.nSigLinear 
        P.plot( [xLims[0], xLims[1]], [Band, Band], 'k--')
#        P.plot( [xLims[0], xLims[1]], [-Band, -Band], 'k--')

        # Overplot the growth-curve
#        xFine = N.linspace(N.min(vX[gHi]), N.max(vX)*1.1, 100)
        xFine = N.linspace(20000.0, N.max(vX)*1.02, 100, endpoint=True)
        yPred = N.polyval(self.vParsDeficit, xFine)
        P.semilogy(xFine, N.exp(yPred), 'r-')

        # Show the equation
        sAnno = 'y = %.2e exp(%.2e x)' % (N.exp(self.vParsDeficit[1]), self.vParsDeficit[0])
        P.annotate(sAnno,(0.10,0.90), xycoords='axes fraction', color='r')

        # Show the fitted limits
        yVec = N.copy(P.ylim())
        yLo = yVec[0]
        P.plot([self.ThisLowPoke, self.ThisLowPoke], [yLo, Band], \
               color='0.5', alpha=0.5)
        P.plot([self.ThisLowThresh, self.ThisLowThresh], \
               [yLo, self.ThisLowThresh * self.MinFracNonlin], \
               color='0.5', alpha=0.5)
        P.grid(which='both', axis='both')

        P.xlabel('Recorded counts')
        P.ylabel('| Incident minus Recorded |')
        P.title('Deficit vs recorded counts')
        P.subplot(223)

        # If the "vs-linear" has been defined, plot that instead.
        if N.size(self.vParsAgainstLinear) < 2:

            P.plot(vX[gLo], N.abs(vY[gLo]), 'bo')
            P.plot(vX[gHi], N.abs(vY[gHi]), 'rs')
            Band = self.ThisStdLinear * self.nSigLinear 
            P.plot( [xLims[0], xLims[1]], [Band, Band], 'k--')
            P.plot( [xLims[0], xLims[1]], [-Band, -Band], 'k--')
            yPred = N.polyval(self.vParsDeficit, xFine)
            P.plot(xFine, N.exp(yPred), 'r-')

            P.xlabel('Recorded counts')
            P.ylabel('Incident minus Recorded')
            P.title('Nonlinearity vs signal')
        
            
        else:

            # Otherwise, we plot the predicted curve.
            xShow = N.polyval(self.ThisLinr, self.vTimes)
            yShow = self.vThisDeficit
            yPred = N.polyval(self.vParsAgainstLinear, xShow)
            P.plot(xShow, yShow, 'ko')
            yLims = N.copy(P.ylim())
            P.plot(xShow, yPred, 'g-')
            P.ylim(yLims)

            # show the equation
            sEqn = 'y = %.3f + %.2f x' % (self.vParsAgainstLinear[1], \
                                          self.vParsAgainstLinear[0])

            sEqn = 'y = %.3f( x - %.2f )' % (self.vParsAgainstLinear[0], \
                                             N.abs(self.vParsAgainstLinear[1] \
                                                   / self.vParsAgainstLinear[0]))
            P.annotate(sEqn,(0.10,0.90), xycoords='axes fraction', color='g')


            P.title('Deficit vs incident counts')
            P.xlabel('Incident counts')
            P.ylabel('Incident minus Recorded, counts')
            
        P.grid()

        # Get the ranges of points that fall within the linear regime
        # (for our fitting)
        xLims = N.copy(P.xlim())

        
        P.suptitle('Linear regression for pixel (%i, %i)' % (self.ThisCol, self.ThisRow), fontsize=16)

        P.subplot(221)
        lSor = N.argsort(self.vTimes)
        iShow = lSor[-10]
        LExtent = [0, N.shape(self.aStack)[1], 0, N.shape(self.aStack)[0]]        

        zVals = N.copy(self.aStack[:,:,iShow])
        nStd = 4.0
        vMin = N.median(zVals) - nStd * N.std(zVals)
        vMax = N.median(zVals) + nStd * N.std(zVals)
        vMax = N.min([vMax, N.max(zVals) ])
        vMin = N.max([vMin, 0.])

        P.imshow(self.aStack[:,:,iShow], \
                 origin='lower', interpolation='nearest', extent=LExtent, aspect=None, cmap=CM.bone, \
                 vmin=vMin, vmax=vMax)  # WARN - HARDCODE!!
        P.plot([self.ThisCol], [self.ThisRow], 'rx', markersize=9)

#        P.plot([self.ThisCol, self.ThisCol], [LExtent[2], LExtent[3]], 'k-', lw=2,alpha=0.35)
#        P.plot([LExtent[0], LExtent[1]], [self.ThisRow, self.ThisRow], 'k-', lw=2,alpha=0.35)
        P.xlim(LExtent[0], LExtent[1])
        P.ylim(LExtent[-2], LExtent[-1])
        P.xlabel('X')
        P.ylabel('Y')
#        P.title('Dome-flat Image, pixel location')  # 2015-04-06
        P.title('Files "%s*.%s", pixel location' % (self.sStart, self.sExt))

#        P.grid()

        P.colorbar(fraction=0.05)
        
        # Save the image to disk
        P.savefig('UMDObs_PixelLinearity.png')



def TestAssembleStack(sStart='Illum__Delay', sExt='fits'):

    FS = FlatStack(sStart=sStart, sExt=sExt)
    FS.AssembleImageSet()
    print "Writing image set to stack..."
    FS.WriteImageSet()
    print "... done."
    FS.InitFineCountsVector()

    FS.RunStatsOnAllPixels()

    return

    # Show one pixel's curve
    FS.SelectPixelVector()
    FS.FitGrowthSelected()
    FS.ShowSelectedPixel()

def TestFittingDeficit(jCol=250, iRow=300, \
                       sStart='Illum__Delay', sExt='fits', \
                       DoAllStats=False, DoLog=False, tShort=0.00, \
                       DoWriteStack=False, \
                       UseSimple=True, DoingST8=True, \
                       LenRow=1, LenCol=1):

    """Test fitting the nonlinearity counts-deficit for a single pixel"""

    # Updated AGAIN to test clumping pixels by region. Will need to
    # update the loop through pixels when doing all stats, otherwise
    # will show very messy overlap.

    FF = FlatStack(sStart=sStart, sExt=sExt, MinTExp=tShort)
    FF.ThisCol = jCol
    FF.ThisRow = iRow
    FF.RegLenRow = LenRow
    FF.RegLenCol = LenCol
    FF.AssembleImageSet()

    if DoingST8:
        FF.nDiscrep=120

    if DoWriteStack:

        print "TestFittingDeficit WARN - writing the stack takes 20 minutes"
        print "adjust program if yuo really want this."
        
        ReallyWrite=False
        if ReallyWrite:

            print "TestFittingDeficit INFO - about to write stack to disk..."
            tZero = time.time()
            FF.WriteImageSet()
            print "TestFittingDeficit INFO - written stack: %.2f sec" \
            % (time.time() - tZero)

    FF.InitFineCountsVector()

    FF.InitSinglePixel()
    FF.SelectPixelVector()
    FF.WriteSelectedPixel()

    FF.GetDeltaCounts()
    FF.FitDeltaCounts()
    FF.FitDeltasAgainstLinear()
    FF.GetDeltaLimits()

    # Write to disk
    FF.WriteSelectedPixel()
    
    try:
        FF.ShowCorrectionCurve(DoLog=DoLog)
        FF.ShowDeltaVsLinear()
    except:
        NotOrigKind=True

    # Now that's done, run all the stats on all the pixels (may take a
    # while)
    ####

    print FF.ThisLowPoke
    print FF.ThisLowThresh

    if not DoAllStats:
        return
    print "Running statistics on all pixels..."
    FF.InitFineCountsVector()
    FF.RunStatsOnAllPixels(UseSimple=UseSimple)
    print "... done."


def TestPredict(CountsIn=55000.0, ColCut=631, RowCut=529):

    """Load the exponential prediction for the nonlinearity deficit, and
predict the counts-deficit map as a function of the signal counts"""

    DD = FlatStack()
    DD.LoadExpon()
    DD.LoadVsLinear()

    DD.PredictDeficitVsLinear(CountsIn)
#    DD.PredictDeficitExpon(CountsIn)
    DD.SmoothDeficitPred()
    
    print N.shape(DD.aImgPred)

#    return

    # Show this
    P.figure(5)
    P.clf()
    P.subplots_adjust(hspace=0.4, wspace=0.3)
    P.clf()
    P.subplot(222)
    vMed = N.median(DD.aImgPred)
    vStd = N.std(DD.aImgPred)
    nStd = 3.0
    P.imshow(DD.aImgPred, interpolation='nearest', origin='lower', \
             vmin=vMed - nStd*vStd, vmax=vMed + nStd*vStd, \
             cmap=CM.bone)
    P.title('Predicted nonlinearity, percent')
    P.xlabel('Column Number (X)')
    P.ylabel('Row Number (Y)')
    P.grid(color='y', alpha=0.75, lw=2)
    P.colorbar()

    aShape = N.shape(DD.aImgPred)

    # Do a cut through both
    P.subplot(221)
#    P.plot(DD.aImgPred[RowCut,:],color='0.5', marker='.', alpha=0.5, ls='none')
    P.plot(DD.aImgPred[RowCut,:],'k.', alpha=0.5)
    P.plot(DD.aImgPred[RowCut,:],color='0.15', alpha=0.25) 

    # Is this periodic?
#    lx = N.arange(N.shape(DD.aImgPred)[1])
#    bSho = (lx % 6 == 1)
#    lSho = lx[bSho]
#    ySho = DD.aImgPred[RowCut,:][bSho]
#    P.plot(lSho, ySho, 'g.')

    P.plot(DD.aImgSmooth[RowCut,:], 'r-', lw=2)
    P.xlabel('Column Number (X)')
    P.ylabel('Percent nonlinearity')
    P.title('Horizontal cut, row %i' % (RowCut))
    P.xlim(0,N.shape(DD.aImgPred)[1])

    P.subplot(224)
    
    l=N.arange(N.shape(DD.aImgPred)[0])
#    P.plot(DD.aImgPred[:,ColCut],l,color='0.5', marker='.', alpha=0.5, ls='none')
    P.plot(DD.aImgPred[:,ColCut],l,'k.', alpha=0.5) #color='0.5', alpha=0.5, marker='.')
    P.plot(DD.aImgPred[:,ColCut],l,color='0.15', alpha=0.25) 
    P.plot(DD.aImgSmooth[:,ColCut],l,'b-', lw=2)
    P.xlabel('Percent nonlinearity')
    P.ylabel('Row number (Y)')
    P.title('Vertical cut, column %i' % (ColCut))
    P.ylim(0,N.shape(DD.aImgPred)[0])

    P.subplot(223)
    P.imshow(DD.aImgSmooth, interpolation='nearest', origin='lower', \
             cmap=CM.bone)
    P.plot([ColCut, ColCut], [0, aShape[0]], 'b', lw=2, alpha=0.75)
    P.plot([0,aShape[1]], [RowCut, RowCut], 'r', lw=2, alpha=0.75)
    P.xlim(0,aShape[1])
    P.ylim(0,aShape[0])
    P.grid(color='y', alpha=1.00, lw=2)
    P.colorbar()
    P.title('Smoothed prediction (percent)')
    P.xlabel('Column Number (X)')
    P.ylabel('Row Number (Y)')

    P.suptitle('Prediction { (Incident - Recorded)/Recorded }, under incident flux %.1f counts' % (CountsIn), fontsize=14)

    # save the figure
    sFile = 'NonlinMap_Incident_%i.png' % (int(CountsIn))
    P.savefig(sFile)


#    # vector to fourier-transform
#    vRow = N.asarray(l,'float')
#    vDelta = (DD.aImgPred[:,800])

#    P.figure(6)
#    P.clf()
#    P.plot(vRow, vDelta)

#    pers = N.arange(1,200,1)
#    freq = 1.0/pers
#    print N.shape(freq), N.shape(vRow), N.shape(vDelta)
#    pgram = signal.lombscargle(vRow, vDelta, freq)
#    print N.shape(pgram)

##    spec = N.fft.fft(vDelta)
##    freq = N.fft.fftfreq(N.size(vDelta))

#    P.clf()
#    P.plot(freq, pgram)
