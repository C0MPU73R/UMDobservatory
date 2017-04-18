#
# ShowLinearity.py - shows the linearity predictions for the Dob
#

import pylab as P
P.ion()

import numpy as np
import pickle

import os, shutil, time

# Useful image-smoothing
from scipy import ndimage

# define some useful functions

def f_linear(x, m, c):

    """returns y = m.x + c"""

    return m * x + c

def f_quad(x, a, b, c):
    
    """returns y = a.x^2 + b.x + c"""

    return a* x**2 + b * x + c

def go(tObs=0.5, InFlux=50000.0,FiltSz=8., FiltTyp='gaussian', \
       DoSmooth=True, nStd=3.0, AbsFlux=False):

    """Runs the comparison"""

    DirPickles = os.getcwd()

    FilLinr='EstLinear.pickle'
    FilQuad='EstQuad.pickle'
    FilThre='LinearityThresholds.pickle'
    
    aLinr = pickle.load(open('%s/%s' % (DirPickles, FilLinr), 'r'))
    aQuad = pickle.load(open('%s/%s' % (DirPickles, FilQuad), 'r'))
    aThre = pickle.load(open('%s/%s' % (DirPickles, FilThre), 'r'))

    # "center of mass"
    tCM = aThre[:,:,1]

    # Apply the center of mass correction right away
    tCor = tObs - tCM

    # Work out the ObsTime from the input flux

    print "Predicting effective exposure time for input flux %.1f" % (InFlux)
    tCor = (InFlux - aLinr[:,:,1])/aLinr[:,:,0]

    print np.median(tCor), np.std(tCor), np.min(tCor), np.max(tCor)

    # Now predict the linear and quadratic at this tObs
    PredLinr = f_linear(tCor, aLinr[:,:,0], aLinr[:,:,1])
    PredQuad = f_quad(tCor, aQuad[:,:,0], aQuad[:,:,1], \
                        aQuad[:,:,2])

    PredDiff = 100.0 * ((PredQuad-PredLinr)/PredLinr)

    # if AbsFlux, this time report the flux deficit in counts.
    if AbsFlux:
        PredDiff = PredQuad - PredLinr

    pMed = np.median(PredDiff)
    pStd = np.std(PredDiff)
    print np.median(PredDiff), np.std(PredDiff)

    # try smoothing pred
#    PredSmooth = np.copy(PredDiff)
#    PredSmooth = ndimage.filters.gaussian_filter(PredDiff, 3.0)

# Dictionary of filter methods
    DF = {'uniform':ndimage.uniform_filter, \
           'gaussian':ndimage.gaussian_filter, \
           'median':ndimage.median_filter}


    PredSmooth = DF[FiltTyp](PredDiff, FiltSz)

    P.figure(1)
    P.clf()

    print np.median(PredLinr), np.std(PredLinr)

    if DoSmooth:
        P.imshow(PredSmooth, interpolation='nearest', origin='lower')
        P.title('Incident flux %.2f counts: Smoothed with %s, width: %.1f pix' \
                % (InFlux, FiltTyp, FiltSz))
        OutImg='Dob_NL_Incident_%i_Smooth.png' % (InFlux)
    else:
        P.imshow(PredDiff, \
                 interpolation='nearest', origin='lower', \
                 vmin=pMed - nStd*pStd, \
                 vmax=pMed + nStd*pStd)
#        print pMed - nStd*pStd
#        print pMed + nStd*pStd
#        print nStd
#        print pStd

        OutImg='Dob_NL_Incident_%i_Clipped.png' % (InFlux)

        P.title('Incident flux %.2f counts: clipped at %.2f sigma' \
                % (InFlux, nStd))

    P.colorbar()   
    P.grid(which='both')

    P.xlabel('X (pix)')
    P.ylabel('Y (pix)')

    P.annotate('Flux deficit (%)', (1.18, 0.5), xycoords='axes fraction',\
               rotation=270, ha='left', va='center', fontsize=14)

#    P.title('Incident flux %.2f counts: Smoothing %s, width: %.1f pix' \
#            % (InFlux, FiltTyp, FiltSz))

    P.savefig(OutImg)

    return

    nStd=0.5
    P.imshow(PredLinr - PredQuad, \
             interpolation='nearest', origin='lower', \
             vmin=pMed - nStd*pStd, \
             vmax=pMed + nStd*pStd)
    P.colorbar()
             
