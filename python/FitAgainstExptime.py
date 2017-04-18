#
# FitAgainstExptime.py
#

#
# Given a set of exptimes, counts for a single pixel, fit a straight
# line and reject outliers.
#

import sys, os, shutil, time
import pickle
import numpy as np
import pylab as P
P.ion()

class PixTrend(object):

    """Class to hold the pixel data, good/bad objects, indices, and fit
parameters."""

    def __init__(self, vExptimes = np.array([]), vCounts=np.array([]) ):

        # Variable for exptime and counts
        self.vExptimes = np.copy(vExptimes)
        self.vCounts  = np.copy(vCounts)

        # which elements are being "included?"
        self.vInclude = np.copy(vCounts)

        # ID - which "clump" in exposure times each point belongs to
        self.vIDs = np.zeros(np.size(vCounts))  # initialise to zeros
        self.vUniqIDs = np.array([])

        # We might be reading in the curve of growth for a given pixel
        self.PathCanned='SinglePix_c250_r300_Data.pickle'

        # Minimum difference in exptime when distributing parts
        self.MinTDiff = 1.0e-3

        # Linear parameters for the central clump (so that we can
        # reject points)
        self.parsRough = np.zeros(2)

        # Control variable
        self.Verbose=True

    def ImportCannedPixel(self):

        """Import a pixel we made earlier"""

        if not os.access(self.PathCanned, os.R_OK):
            if self.Verbose:
                print "ImportCannedPixel FATAL - cannot read path %s" \
                    % (self.PathCanned)
            return

        # Get the dictionary
        DRead = pickle.load(open(self.PathCanned,'r'))
        self.vExptimes = np.copy(DRead['vTimes'])
        self.vCounts = np.copy(DRead['vThisPixel'])

        # Initialise the IDs
        self.vIDs = np.zeros(np.size(self.vCounts))

    def ClumpByExptime(self, ShowClumps=False):

        """Breaks the exptimes into groups by exptime"""

        # NOTE this has a very specific use-case: commanded exptimes
        # at least 1ms apart from each other (want to account for
        # rounding)
        nTimes = np.size(self.vExptimes)
        lTimes = np.argsort(self.vExptimes)
        ll = np.arange(nTimes)

        # New vector - times ordered in increasing exptime
        tOrdered = self.vExptimes[lTimes]
        tDiff = tOrdered - np.roll(tOrdered, 1)
        tDiff[0] = np.max(tDiff)+0.1  # set out the first one

        # Wherever tDiff > 0, then the corresponding exptime is the
        # FIRST in the ordered list at that exptime. Use this!
        gUniq = np.where(tDiff > self.MinTDiff)[0]
        if np.size(gUniq) < 1:
            return

        # Now generate array of IDs for each point
        self.vIDs = np.zeros(np.size(tOrdered))
        for iStep in range(np.size(gUniq)):
            ThisTime = tOrdered[gUniq[iStep]]
            gIsThisTime = np.where( np.abs(self.vExptimes - ThisTime) \
                                    < self.MinTDiff)[0]

            # Do nothing if a blank step (neat trick!!)
            if np.size(gIsThisTime) < 1:
                continue
 
            # Now fill in the ID for these datapoints
            self.vIDs[gIsThisTime] = iStep

        # Establish now the unique ID values
        self.FindUniqueIDs()

        # If not showing the debug plot, we're done here. Exit.
        if not ShowClumps:
            return

        P.figure(4)
        P.clf()
#        P.plot(ll, tOrdered)
        P.plot(ll, tDiff, 'k')
        P.scatter(ll, tOrdered, c=self.vIDs[lTimes], \
                  edgecolor='none')
        P.xlabel('lTimes')
        P.ylabel('vExptimes[lTimes]')
        P.colorbar()
        P.grid()

    def FindUniqueIDs(self):

        """Convenience-method - update the variable that gives unique ID
        values."""

        if np.size(self.vIDs) < 1:
            return

        self.vUniqIDs = np.unique(np.sort(self.vIDs))

    def InitInclude(self):

        """Initialises the "include" vectors"""

        self.vInclude = np.ones(np.size(self.vExptimes))

    def StatsPerExptime(self, nClip=2, nStd=3, nMin=3):

        """Converts the raw exptime measurements into mean, stddev,
        optionally clipping"""

        # Nothing to do if no unique points
        if np.size(self.vUniqIDs) < 1:
            return

        # initialise the "use" vectors (this might not be the best place
        # to do it?)
        self.InitInclude()

        # Set up vectors for exptime, mean, stddev
        self.vTimes = np.zeros(np.size(self.vUniqIDs))
        self.vMeans = np.copy(self.vTimes)
        self.vStdds = np.copy(self.vTimes)
        self.vNPts  = np.copy(self.vTimes)

        # When clipping, keep track of arrays in the ORIGINAL vectors
        # of which measurements are eliminated (will be useful to
        # check later on).
        for iID in self.vUniqIDs:
            
            # Which measurements are in this step?
            gThis = np.where(np.abs(self.vIDs - self.vUniqIDs[iID]) < 0.4)[0]
            
            # Need at least 2 measurements to get the stddev
            if np.size(gThis) < 2:
                continue

#            print self.vUniqIDs[iID]

            # Go through nClip + 1 times (first calculate with all,
            # then clip)
            for jClip in range(nClip + 1):
                
                # vector of points not yet rejected
                gUse = np.where(self.vInclude[gThis] > 0)[0]
                gKeep = gThis[gUse]
                if np.size(gKeep) < nMin:
                    continue

                # Perform the calculation
                ThisNPts = np.size(gKeep)
                ThisStd  = np.std(self.vCounts[gKeep])
                ThisTExp = np.mean(self.vExptimes[gKeep])

                # Use median if lots of points, mean otherwise.
                if ThisNPts < 6:
                    ThisMean = np.median(self.vCounts[gKeep])
                else:
                    ThisMean = np.median(self.vCounts[gKeep])

                # No point in clipping if we're in a place where
                # stddev = 0... This can happen if we're heavily
                # saturated, where all measurements read exactly the
                # same value.
                    
                # WATCHOUT - if converting to count rates, "1" would
                # no longer be a small stddev.
                if ThisStd < 5.0:   # use "small" value.
                    continue

                # Now evaluate which points are "close" to the "mean."
                # Update the "vInclude" values so that we have access
                # to the "include/don't include" data everywhere.
                if jClip < nClip:
                    vDist = np.abs(self.vCounts[gKeep] - ThisMean)/ThisStd
                    gFar = np.where(vDist > nStd)[0]
                    gOut = gKeep[gFar]
                    if np.size(gOut) > 0:
                        self.vInclude[gOut] = 0
                   
            # Now that the clipping passes are done, pass the estimate
            # up.
            self.vTimes[iID] = np.copy(ThisTExp)
            self.vMeans[iID] = np.copy(ThisMean)
            self.vStdds[iID] = np.copy(ThisStd)
            self.vNPts[iID] = np.copy(ThisNPts)

    def FitCentralLinear(self):

        """Having found the means and errors on each exptime point,
        estimate a straight line fit using the inner +/- 20% of points
        from the midpoint"""

        # (Could make that 20% tunable)

        # This won't work if fewer than six datapoints (can make
        # tunable)
        if np.size(self.vTimes) < 6:
            return
        
        # Argsort by exptime
        lTimes = np.argsort(self.vTimes)
        
        # Midpoint - keep as integers
        iMid = np.size(self.vTimes)/2
        iHalf = np.max([2, np.size(self.vTimes)/4])

        # Convenience:
        iLo = iMid - iHalf
        iHi = iMid + iHalf

        # Try this: (USE SOMETHING BETTER THAN POLYFIT!!! The whole
        # point is to maximize the use of a small number of points!!
        self.ParsRough = np.polyfit(\
            self.vTimes[iLo:iHi], self.vMeans[iLo:iHi], 1)


    def ShowPixel(self, ShowStats=False):

        """Convenience-method to plot the pixel curve read in"""

        if np.size(self.vExptimes) < 2:
            return

        # Show the IDs if we have assigned them already
        if np.size(self.vUniqIDs) > 1:
            colo = self.vIDs
        else:
            colo = 'g'

        P.figure(4)
        P.clf()

        if ShowStats:
            P.plot(self.vExptimes, self.vCounts, 'ks', alpha=0.25, \
                       zorder=5, ms=6)
            P.errorbar(self.vTimes, self.vMeans, self.vStdds, fmt='ro', \
                           zorder=10, linestyle='none', elinewidth=2, \
                           ms=4, ecolor='r')
            P.xlabel('Exptimes')
            P.ylabel('Counts')
            P.grid(which='both')
            P.title('Statistics')

            # If we have parameters, overplot the linear fit
            xFine = np.linspace(\
                np.min(self.vExptimes), \
                    np.max(self.vExptimes), 100)
            if np.sum(np.abs(self.ParsRough)) > 1e-2:
                P.plot(xFine, np.polyval(self.ParsRough, xFine), 'g--')

                print self.ParsRough
                P.plot(self.vTimes, \
                           np.polyval(self.ParsRough, self.vTimes) - self.vMeans, \
                           'gs')


            # exit
            return

        # Show ALL the datapoints
        P.scatter(self.vExptimes, self.vCounts, c=colo, edgecolor='0.95', \
                      zorder=5)

        # Show the datatpoints passing the "use" steps
        g = np.where(self.vInclude > 0)[0]
        if np.size(g) > 3:
            P.plot(self.vExptimes[g], self.vCounts[g], 'k.', zorder=10)

        P.xlabel('Exptimes')
        P.ylabel('Counts')
        P.grid(which='both')

def TestLoadCanned():

    """Tests loading a canned file"""

    PT = PixTrend()
    PT.ImportCannedPixel()

    tZero = time.time()
    PT.ClumpByExptime()
    PT.StatsPerExptime()
    
    PT.FitCentralLinear()

    print time.time() - tZero
#    print PT.vUniqIDs
    PT.ShowPixel(ShowStats=True)
