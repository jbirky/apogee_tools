###############################################################################
# apogee.spec.lsf: Utilities to work with APOGEE LSFs
###############################################################################
import os, os.path
from functools import wraps
import warnings
import math
import numpy
from scipy import special, interpolate, sparse, ndimage
import scipy.sparse.linalg
import fitsio
import apogee_tools.apogee_hack.tools.read as apread
import apogee_tools.apogee_hack.tools.path as appath
from apogee_tools.apogee_hack.tools.download import _download_file
from apogee_tools.apogee_hack.spec.plot import apStarWavegrid

_SQRTTWO= numpy.sqrt(2.)
# Load wavelength solutions
_WAVEPIX_A= apread.apWave('a',ext=2)
_WAVEPIX_B= apread.apWave('b',ext=2)
_WAVEPIX_C= apread.apWave('c',ext=2)

def convolve(wav,spec,
             lsf=None,xlsf=None,dxlsf=None,fiber='combo',
             vmacro=6.):
    """
    NAME:
       convolve
    PURPOSE:
       convolve with the APOGEE LSF and resample to APOGEE's apStar wavelength grid
    INPUT:
       wav - wavelength array (linear in wavelength in \AA)
       spec - spectrum on wav wavelength grid [nspec,nwave]
       lsf= (None) pre-calculated LSF array from apogee.spec.lsf.eval
       Either:
          xlsf= (None) 1/integer equally-spaced pixel offsets at which the lsf=lsf input is calculated
          dxlsf= (None) spacing of pixel offsets
       fiber= if lsf is None, the LSF is calculated for this fiber
       vmacro= (6.) Gaussian macroturbulence smoothing to apply as well (FWHM or a [sparse] matrix like lsf on the same x grid; can be computed with apogee.modelspec.vmacro)
    OUTPUT:
       spectrum on apStar wavelength grid
    HISTORY:
       2015-03-14 - Written - Bovy (IAS)
    """
    # Parse LSF input
    if lsf is None:
        xlsf= numpy.linspace(-7.,7.,43)
        lsf= eval(xlsf,fiber=fiber)
    if not isinstance(lsf, sparse.dia_matrix):
        lsf= sparsify(lsf)
    if dxlsf is None:
        dx= xlsf[1]-xlsf[0]
    else:
        dx= dxlsf

    hires  = int(round(1./dx))
    l10wav = numpy.log10(apStarWavegrid())
    dowav  = l10wav[1]-l10wav[0]
    tmpwav = 10.**numpy.arange(l10wav[0],l10wav[-1]+dowav/hires,dowav/hires)
    tmp    = numpy.empty(len(l10wav)*hires)   

    # Interpolate the input spectrum, starting from a polynomial baseline
    if len(spec.shape) == 1: 
        spec = numpy.reshape(spec,(1,len(spec)))
    nspec = spec.shape[0]
    tmp   = numpy.empty((nspec,len(tmpwav)))
    for ii in range(nspec):
        baseline = numpy.polynomial.Polynomial.fit(wav,spec[ii], 4)
        ip       = interpolate.InterpolatedUnivariateSpline(wav,
                                                            spec[ii]/baseline(wav),
                                                            k=3)
        tmp[ii]  = baseline(tmpwav)*ip(tmpwav)

    # Use sparse representations to quickly calculate the convolution
    tmp= sparse.csr_matrix(tmp)
    return lsf.dot(tmp.T).T.toarray()[:,::hires]

def sparsify(lsf):
    """
    NAME:
       sparsify
    PURPOSE:
       convert an LSF matrix calculated with eval [ncen,npixoff] to a sparse [ncen,ncen] matrix with the LSF on the diagonals (for quick convolution with the LSF)
    INPUT:
       lsf - lsf matrix [ncen,npixoff] calculated by eval
    OUTPUT:
       sparse matrix with the lsf on the diagonals
    HISTORY:
       2015-03-14 - Written - Bovy (IAS)
    """
    nx= lsf.shape[1]
    diagonals= []
    offsets= []
    for ii in range(nx):
        offset= nx//2-ii
        offsets.append(offset)
        if offset < 0:
            diagonals.append(lsf[:offset,ii])
        else:
            diagonals.append(lsf[offset:,ii])
    return sparse.diags(diagonals,offsets)

def dummy(dx=1./3.,sparse=False):
    """
    NAME:
       dummy
    PURPOSE:
       return a 'dummy' LSF that is a delta function
    INPUT:
       dx= (1/3) spacing between LSF centers in the apStar grid
       sparse= (False) if True, return a sparse representation that can be passed to apogee.spec.lsf.convolve for easy convolution
    OUTPUT:
       LSF(x|pixel center);
       pixel centers are apStarWavegrid if dx=1, and denser 1/integer versions if dx=1/integer
    HISTORY:
       2015-03-23 - Written - Bovy (IAS)
    """
    # Are the x unit pixels or a fraction 1/hires thereof?
    hires= int(round(1./dx))
    # Setup output
    wav= apStarWavegrid()
    l10wav= numpy.log10(wav)
    dowav= l10wav[1]-l10wav[0]
    # Hi-res wavelength for output
    hireswav= 10.**numpy.arange(l10wav[0],l10wav[-1]+dowav/hires,dowav/hires)
    out= numpy.ones((len(hireswav),1))
    if sparse: out= sparsify(out)
    return out

def eval(x,fiber='combo',sparse=False):
    """
    NAME:
       eval
    PURPOSE:
       evaluate the LSF for a given fiber
    INPUT:
       x - Array of X values for which to compute the LSF, in pixel offset relative to pixel centers; the LSF is calculated at the x offsets for each pixel center; x need to be 1/integer equally-spaced pixel offsets
       fiber= ('combo') fiber number or 'combo' for an average LSF (uses the same one-based indexing as the APOGEE fibers [i.e., fibers range from 1 to 300])
       sparse= (False) if True, return a sparse representation that can be passed to apogee.spec.lsf.convolve for easy convolution
    OUTPUT:
       LSF(x|pixel center);
       pixel centers are apStarWavegrid if dx=1, and denser 1/integer versions if dx=1/integer
    HISTORY:
       2015-03-12 - Written based on Jon H's code (based on David N's code) - Bovy (IAS)
    """
    # Parse fiber input
    if isinstance(fiber,str) and fiber.lower() == 'combo':
        fiber= [50,100,150,200,250,300]
    elif isinstance(fiber,int):
        fiber= [fiber]
    elif not isinstance(fiber,list) and isinstance(fiber[0],int):
        raise ValueError('fiber input to apogee.spec.lsf.eval not understood ...')
    # Are the x unit pixels or a fraction 1/hires thereof?
    hires= int(round(1./(x[1]-x[0])))
    # Setup output
    wav= apStarWavegrid()
    l10wav= numpy.log10(wav)
    dowav= l10wav[1]-l10wav[0]
    # Hi-res wavelength for output
    hireswav= 10.**numpy.arange(l10wav[0],l10wav[-1]+dowav/hires,dowav/hires)
    out= numpy.zeros((len(hireswav),len(x)))
    for chip in ['a','b','c']:
        # Get pixel array for this chip, use fiber[0] for consistency if >1 fib
        pix= wave2pix(hireswav,chip,fiber[0])
        dx= numpy.roll(pix,-hires,)-pix
        dx[-1]= dx[-1-hires]
        dx[-2]= dx[-2-hires]
        dx[-3]= dx[-3-hires]
        xs= numpy.tile(x,(len(hireswav),1))\
            *numpy.tile(dx,(len(x),1)).T # nwav,nx     
        gd= numpy.bitwise_xor(True, numpy.isnan(pix))
        # Read LSF file for this chip
        lsfpars= apread.apLSF(chip,ext=0)
        # Loop through the fibers
        for fib in fiber:
            out[gd]+= raw(xs[gd],pix[gd],lsfpars[:,300-fib])
    out[out<0.]= 0.
    out/= numpy.tile(numpy.sum(out,axis=1),(len(x),1)).T
    if sparse: out= sparsify(out)
    return out
    
def raw(x,xcenter,params):
    """
    NAME:
       raw
    PURPOSE:
       Evaluate the raw APOGEE LSF (on the native pixel scale)
    INPUT:
       x - Array of X values for which to compute the LSF (in pixel offset relative to xcenter; the LSF is calculated at the x offsets for each xcenter if x is 1D, otherwise x has to be [nxcenter,nx]))
       xcenter - Position of the LSF center (in pixel units)
       lsfarr - the parameter array (from the LSF HDUs in the APOGEE data products)
    OUTPUT:
       LSF(x|xcenter))
    HISTORY:
       2015-02-26 - Written based on Nidever's code in apogeereduce - Bovy (IAS)
    """
    # Parse x
    if len(x.shape) == 1:
        x= numpy.tile(x,(len(xcenter),1))
    # Unpack the LSF parameters
    params= unpack_lsf_params(params)
    # Get the wing parameters at each x
    wingparams= numpy.empty((params['nWpar'],len(xcenter)))
    for ii in range(params['nWpar']):
        poly= numpy.polynomial.Polynomial(params['Wcoefs'][ii])       
        wingparams[ii]= poly(xcenter+params['Xoffset'])
    # Get the GH parameters at each x
    ghparams= numpy.empty((params['Horder']+2,len(xcenter)))
    for ii in range(params['Horder']+2):
        if ii == 1:
            # Fixed, correct for wing
            ghparams[ii]= 1.-wingparams[0]
        else:
            poly= numpy.polynomial.Polynomial(params['GHcoefs'][ii-(ii > 1)])
            ghparams[ii]= poly(xcenter+params['Xoffset'])
        # normalization
        if ii > 0: ghparams[ii]/= numpy.sqrt(2.*numpy.pi*math.factorial(ii-1))
    # Calculate the GH part of the LSF
    out= _gausshermitebin(x,ghparams,params['binsize'])
    # Calculate the Wing part of the LSF
    out+= _wingsbin(x,wingparams,params['binsize'],params['Wproftype'])
    return out

def _gausshermitebin(x,params,binsize):
    """Evaluate the integrated Gauss-Hermite function"""
    ncenter= params.shape[1]
    out= numpy.empty((ncenter,x.shape[1]))
    integ= numpy.empty((params.shape[0]-1,x.shape[1]))
    for ii in range(ncenter):
        poly= numpy.polynomial.HermiteE(params[1:,ii])
        # Convert to regular polynomial basis for easy integration
        poly= poly.convert(kind=numpy.polynomial.Polynomial)
        # Integrate and add up
        w1= (x[ii]-0.5*binsize)/params[0,ii]
        w2= (x[ii]+0.5*binsize)/params[0,ii]
        eexp1= numpy.exp(-0.5*w1**2.)
        eexp2= numpy.exp(-0.5*w2**2.)
        integ[0]= numpy.sqrt(numpy.pi/2.)\
            *(special.erf(w2/_SQRTTWO)-special.erf(w1/_SQRTTWO))
        out[ii]= poly.coef[0]*integ[0]
        if params.shape[0] > 1:
            integ[1]= -eexp2+eexp1
            out[ii]+= poly.coef[1]*integ[1]
        for jj in range(2,params.shape[0]-1):
            integ[jj]= (-w2**(jj-1)*eexp2+w1**(jj-1)*eexp1)\
                +(jj-1)*integ[jj-2]
            out[ii]+= poly.coef[jj]*integ[jj]
    return out

def _wingsbin(x,params,binsize,Wproftype):
    """Evaluate the wings of the LSF"""
    ncenter= params.shape[1]
    out= numpy.empty((ncenter,x.shape[1]))
    for ii in range(ncenter):
        if Wproftype == 1: # Gaussian
            w1=(x[ii]-0.5*binsize)/params[1,ii]
            w2=(x[ii]+0.5*binsize)/params[1,ii]
            out[ii]= params[0,ii]/2.*(special.erf(w2/_SQRTTWO)\
                                          -special.erf(w1/_SQRTTWO))
    return out

def unpack_lsf_params(lsfarr):
    """
    NAME:
       unpack_lsf_params
    PURPOSE:
       Unpack the LSF parameter array into its constituents
    INPUT:
       lsfarr - the parameter array
    OUTPUT:
       dictionary with unpacked parameters and parameter values:
          binsize: The width of a pixel in X-units
          Xoffset: An additive x-offset; used for GH parameters that vary globally
          Horder: The highest Hermite order
          Porder: Polynomial order array for global variation of each LSF parameter
          GHcoefs: Polynomial coefficients for sigma and the Horder Hermite parameters
          Wproftype: Wing profile type
          nWpar: Number of wing parameters
          WPorder: Polynomial order for the global variation of each wing parameter          
          Wcoefs: Polynomial coefficients for the wings parameters
    HISTORY:
       2015-02-15 - Written based on Nidever's code in apogeereduce - Bovy (IAS@KITP)
    """
    out= {}
    # binsize: The width of a pixel in X-units
    out['binsize']= lsfarr[0]
    # X0: An additive x-offset; used for GH parameters that vary globally
    out['Xoffset']= lsfarr[1]
    # Horder: The highest Hermite order
    out['Horder']= int(lsfarr[2])
    # Porder: Polynomial order array for global variation of each LSF parameter
    out['Porder']= lsfarr[3:out['Horder']+4]
    out['Porder']= out['Porder'].astype('int')
    nGHcoefs= numpy.sum(out['Porder']+1)
    # GHcoefs: Polynomial coefficients for sigma and the Horder Hermite parameters
    maxPorder= numpy.amax(out['Porder'])
    GHcoefs= numpy.zeros((out['Horder']+1,maxPorder+1))
    GHpar= lsfarr[out['Horder']+4:out['Horder']+4+nGHcoefs] #all coeffs
    CoeffStart= numpy.hstack((0,numpy.cumsum(out['Porder']+1)))
    for ii in range(out['Horder']+1):
        GHcoefs[ii,:out['Porder'][ii]+1]= GHpar[CoeffStart[ii]:CoeffStart[ii]+out['Porder'][ii]+1]
    out['GHcoefs']= GHcoefs
    # Wproftype: Wing profile type
    wingarr= lsfarr[3+out['Horder']+1+nGHcoefs:]
    out['Wproftype']= int(wingarr[0])
    # nWpar: Number of wing parameters
    out['nWpar']= int(wingarr[1])
    # WPorder: Polynomial order for the global variation of each wing parameter
    out['WPorder']= wingarr[2:2+out['nWpar']]
    out['WPorder']= out['WPorder'].astype('int')
    # Wcoefs: Polynomial coefficients for the wings parameters
    maxWPorder= numpy.amax(out['WPorder'])
    Wcoefs= numpy.zeros((out['nWpar'],maxWPorder+1))
    Wpar= wingarr[out['nWpar']+2:]
    WingcoeffStart= numpy.hstack((0,numpy.cumsum(out['WPorder']+1)))
    for ii in range(out['nWpar']):
        Wcoefs[ii,:out['WPorder'][ii]+1]= Wpar[WingcoeffStart[ii]:WingcoeffStart[ii]+out['WPorder'][ii]+1]
    out['Wcoefs']= Wcoefs
    return out

def scalarDecorator(func):
    """Decorator to return scalar outputs for wave2pix and pix2wave"""
    @wraps(func)
    def scalar_wrapper(*args,**kwargs):
        if numpy.array(args[0]).shape == ():
            scalarOut= True
            newargs= (numpy.array([args[0]]),)
            for ii in range(1,len(args)):
                newargs= newargs+(args[ii],)
            args= newargs
        else:
            scalarOut= False
        result= func(*args,**kwargs)
        if scalarOut:
            return result[0]
        else:
            return result
    return scalar_wrapper

@scalarDecorator
def wave2pix(wave,chip,fiber=300):
    """
    NAME:
       wave2pix
    PURPOSE:
       convert wavelength to pixel
    INPUT:
       wavelength - wavelength (\AA)
       chip - chip to use ('a', 'b', or 'c')
       fiber= (300) fiber to use the wavelength solution of
    OUTPUT:
       pixel in the chip
    HISTORY:
        2015-02-27 - Written - Bovy (IAS)
    """
    if chip == 'a':
        wave0= _WAVEPIX_A[300-fiber]
    if chip == 'b':
        wave0= _WAVEPIX_B[300-fiber]
    if chip == 'c':
        wave0= _WAVEPIX_C[300-fiber]
    pix0= numpy.arange(len(wave0))
    # Need to sort into ascending order
    sindx= numpy.argsort(wave0)
    wave0= wave0[sindx]
    pix0= pix0[sindx]
    # Start from a linear baseline
    baseline= numpy.polynomial.Polynomial.fit(wave0,pix0,1)
    ip= interpolate.InterpolatedUnivariateSpline(wave0,pix0/baseline(wave0),
                                                 k=3)
    out= baseline(wave)*ip(wave)
    # NaN for out of bounds
    out[wave > wave0[-1]]= numpy.nan
    out[wave < wave0[0]]= numpy.nan
    return out

@scalarDecorator
def pix2wave(pix,chip,fiber=300):
    """
    NAME:
       pix2wave
    PURPOSE:
       convert pixel to wavelength
    INPUT:
       pix - pixel
       chip - chip to use ('a', 'b', or 'c')
       fiber= (300) fiber to use the wavelength solution of
    OUTPUT:
       wavelength in \AA
    HISTORY:
        2015-02-27 - Written - Bovy (IAS)
    """
    if chip == 'a':
        wave0= _WAVEPIX_A[300-fiber]
    if chip == 'b':
        wave0= _WAVEPIX_B[300-fiber]
    if chip == 'c':
        wave0= _WAVEPIX_C[300-fiber]
    pix0= numpy.arange(len(wave0))
    # Need to sort into ascending order
    sindx= numpy.argsort(pix0)
    wave0= wave0[sindx]
    pix0= pix0[sindx]
    # Start from a linear baseline
    baseline= numpy.polynomial.Polynomial.fit(pix0,wave0,1)
    ip= interpolate.InterpolatedUnivariateSpline(pix0,wave0/baseline(pix0),
                                                 k=3)
    out= baseline(pix)*ip(pix)
    # NaN for out of bounds
    out[pix < 0]= numpy.nan
    out[pix > 2047]= numpy.nan
    return out

def lsf_gh(x,xcenter,par,dlsfgh=None,globalderiv=True,double=True,stp=False,nogauss=False,nowings=False):
    """
    NAME:
        lsf_gh
    PURPOSE:
        This returns the LSF using a superposition of 
        Gauss-Hermite functions.  The LSF is *always* normalized.
        Therefore, the coefficient for H0 is fixed to 1.

        NOTE: There are three ways to input the X/Xcenter values.
            See below.

    INPUTS:
        x       The array of X-values for which to compute the LSF.
                This must be a 1 or 2-dimensional array (see below).
        xcenter The position of the LSF center (in X-units).  This
                must be a scalar (and the X-values 1-dim) or a
                1-dimensional array (see below).

        NOTE: There are three ways to input X/Xcenter:
            1) X as a 1-dimensional array, and Xcenter as a scalar
               so the same center for all X.
            2) X as a 2-dimensional array [Ncenter,Nlsf] and X-center
               as a 1-dimensional array [Ncenter].  
            3) X as a 1-dimensional array and X-center as a 1-dimensional
               array with the same number of elements as X (i.e. there
               is an Xcenter for every X.  
     
        par     The input parameters:

        binsize  The width of a pixel in X-units.  If this is non-zero
                    then a "binned" Gauss-Hermite function is used.  If
                    binsize=0 then a "normal, unbinned" Gauss-Hermite
                   function is used.
        X0       An additive x-offset.  This is only used to
                   evaluate the GH parameters that vary globally
                   with X.
        Horder   The highest Hermite order, Horder=0 means
                   only a constant term (i.e. only Gaussian).
                   There are Horder Hermite coefficients (since we fix H0=1).
        Porder   This array gives the polynomial order for the
                   global variation (in X) of each LSF parameter.
                   That includes sigma and the Horder Hermite
                   coefficients (starting with H1 because we fix H0=1)
                   There will be Porder[i]+1 coefficients for
                   parameter i.
        GHcoefs  The polynomial coefficients for sigma and the
                   Horder Hermite parameters.  There are Porder[i]+1
                   coefficients for parameter i.  The Hermite parameters
                   start with H1 since we fix H0=1.
        Wproftype  The Wing profile type
        nWpar      The number of wing parameters.
        WPorder    An array similar to Porder to give the polynomial
                     order for the global varition (in X) of each wing parameter.
        Wcoefs     The polynomial coefficients for the wings parameters. 
                     There are WPorder[i]+1 coefficients for wing
                     parameters i.

      /nogauss  This sets H0=0 which removes the constant Gaussian
                   term and makes the sum=0.
      /stp      Stop at the end of the program.

    OUTPUTS:
        lsf       The LSF function
        dlsf      The partial derivative of lsf [Nx,Npar].

    USAGE:
        IDL>lsf=lsf_gh(x,xcenter,par)

    BY D. Nidever  March/April 2010
        June 2012, merged lsf_gh, lsf_gh1d, and lsf_gh2d into one
    PORTED TO PYTHON BY C. Theissen  February 2023
    """

    nx       = len(x)
    nxcenter = len(xcenter)
    npar     = len(par)

    #Not enough input parameters
    if nx == 0 or nxcenter == 0 or npar == 0:
      raise Exception('Syntax - lsf = lsf_gh(x,xcenter,par,dp)')
      return -1

    if len(dbl) == 0: dbl=1  # double-precision by default

    # Determine the type of X/Xcenter input
    inptype    = 0
    x_sz       = size(x)
    xcenter_sz = size(xcenter)
    if x_sz[0] == 1 and xcenter_sz[0] == 0: inptype=1
    if x_sz[0] == 2 and xcenter_sz[0] == 1: inptype=2
    if x_sz[0] == 1 and xcenter_sz[0] == 1: inptype=3
    if inptype == 0:
      print('X/XCENTER dimensions not compatible')
      print('1) X [Nlsf] and Xcenter (scalar)')
      print('2) X [Ncenter,Nlsf] and Xcenter [Ncenter]')
      print('3) X [N] and Xcenter [N]')
      return -1

    # Make internal X/Xcenter arrays depending on the
    # input type
    
    # normal
    if inptype == 1:
      XX   = x
      Xcen = [Xcenter + i for i in range(x_sz[1])]
    # 2D
    #  X [Ncenter,Nlsf]
    #  Xcenter [Ncenter]
    elif inptype == 2:
      XX   = np.multiply(x, np.ones_like(x))
      Xcen = [Xcenter] * x_sz[2]#( Xcenter; replicate(1.0d0,x_sz[2]) )(*)
    # 1D
    elif inptype == 3:
      XX   = x
      Xcen = Xcenter

    # Breaking up the parameters
    binsize  = par[0]
    Xoffset  = par[1]   # Additive Xoffset
    Horder   = par[2]
    Porder   = par[3:Horder+3]   # Horder+1 array
    nGHcoefs = np.total(Porder+1)

    # Getting the GH parameters that vary globally
    cpar = par[Horder+4:Horder+4+nGHcoefs-1]
    if double: 
        coefarr = dblarr(Horder+1,max(Porder)+1) 
    else:
        coefarr = fltarr(Horder+1,max(Porder)+1)
    cstart = [0] + [sum(Porder[:i+1]) for i in range(len(Porder))] # extra one at the end
    # Coefarr might have extra zeros at the end, but it shouldn't
    #  make a difference.
    for i in range(Horder):
        coefarr[i,0:Porder[i]] = cpar[cstart[i]:cstart[i]+Porder[i]]

    # XX and Xcen are now a 1D array
    # There is a separate center for each X value.
    # Now we need parameters for EACH X-element
    sz1 = size(XX)
    if double: 
        GHpar = np.zeros((sz1[1],Horder+1))
    else:
        GHpar = np.zeros((sz1[1],Horder+1), dtype=np.float32)

    nGHpar = Horder+1

    # To evaluate the GH parmeters that vary with X we evaluate
    # them at Center+Xoffset
    Xcenter1 = Xcen+Xoffset
    for i in range(Horder):
        GHpar[:,i] = np.array([np.poly1d(x, coefarr[i,:]) for x in Xcenter1])

    # Wing parameters
    if npar > (3+Horder+1+nGHcoefs):
        wcpar = par[3+Horder+1+nGHcoefs:] # the wing part of "par"

        # Nwpar     number of W parameters
        # WPorder   the polynomial order for each
        # Wing coefficients
        wproftype = wcpar[0]
        nWpar     = wcpar[1]
        wPorder   = wcpar[2:2+nWpar-1]
        nWcoefs   = total(wPorder+1)

        # Getting the Wing parameters that vary globally
        wcoef     = wcpar[nWpar+2:]
        wcoefarr  = np.zeros((nWpar, max(wPorder)+1))
        wcstart   = [0] + [sum(wPorder[:i+1]) for i in range(len(wPorder))] # extra one at the end
        # wcoefarr might have extra zeros at the end, but it shouldn't
        #  make a difference.
        for i in range(nWpar):    
            wcoefarr[i,0:wPorder[i]] = wcoef[wcstart[i]:wcstart[i]+wPorder[i]]

        # To evaluate the Wing parmeters that vary with X we evalute
        # them at the Center+Xoffset
        Wparam = np.zeros((sz1[1], nWpar))
        for i in range(nWpar):
            Wparam[:,i] = np.array([np.poly1d(x, wcoefarr[i,:]) for x in Xcenter1])


    if nowings==True: wparam = []

    # "Binned" Gauss-Hermite function
    if binsize > 0:

      # Derivatives requested
      if dlsfgh is not None:

        # Wings
        if len(wparam) > 0:
            inpar = np.zeros((sz1[1], nGHpar+nWpar+3))
            inpar[:,0] = binsize
            inpar[:,1] = Horder
            inpar[:,2:nGHpar+1] = GHpar
            inpar[:,nGHpar+2] = wproftype
            inpar[:,nGHpar+3:nGHpar+nWpar+2] = Wparam
    
            lsf = apgprofile(XX, Xcen, inpar, dp=dlsf)

        # NO wings
        else:

            inpar = np.zeros((sz1[1], nGHpar+2))
            inpar[:,0] = binsize
            inpar[:,1] = Horder
            inpar[:,2:nGHpar+1] = GHpar
    
            lsf = apgprofile(XX, Xcen, inpar, dp=dlsf)

      # No derivatives
      else:

        # Wings
        if len(wparam) > 0:

            inpar                            = np.zeros(sz1[1], nGHpar+nWpar+3)
            inpar[:,0]                       = binsize
            inpar[:,1]                       = Horder
            inpar[:,2:nGHpar+1]              = GHpar
            inpar[:,nGHpar+2]                = wproftype
            inpar[:,nGHpar+3:nGHpar+nWpar+2] = Wparam

            lsf = apgprofile(XX,Xcen,inpar)

        # NO wings
        else:

            inpar               = np.zeros(sz1[1], nGHpar+2)
            inpar[:,0]          = binsize
            inpar[:,1]          = Horder
            inpar[:,2:nGHpar+1] = GHpar

            lsf = apgprofile(XX,Xcen,inpar)

    # No binning
    else:

      raise Exception('No-binning not supported yet')
      return -1

      # Derivatives requested
      #if arg_present(dlsfgh) then begin
      #  lsf = GAUSSHERMITE(x,inpar,dp=dlsf)
      #endif else begin
      #  lsf = GAUSSHERMITE(x,inpar)
      #endelse

    #if keyword_set(nogauss) then inpar[*,3]=0   ; No constant Gaussian term, H0

    #stop

    # Derivatives requested
    if dlsfgh is not None:

      # Convert to global input parameter derivatives
      if globalderiv == True:

        # Breaking up the parameters
        # binsize - 0  FIXED
        # Xoffset - 1  FIXED
        # Horder  - 2  FIXED
        # Porder  - 3:Horder+3
        # GHcoefs - Horder+4:npar-1


        # We will return the derivatives of the height, center
        #  and the GHcoefs.
        nGHcoefs = total(Porder+1)
        nWcoefs  = total(WPorder+1)

        # Initialize the output array
        dlsfgh   = make_array(sz1[1],nGHcoefs+2+nWcoefs,value=lsf[0]*0)

        # The first one is the height
        dlsfgh[:,0] = dlsf[:,0]
        # The first one is the center
        dlsfgh[:,1] = dlsf[:,1]

        derivparcntr = 2   # counter, where the next parameters will start in the array

        # The parameters in dlsf are: [height, center, sigma, H0, H1, H2, H3, H4]
        # H0 is always fixed at 1. we want H1-H4
        # so this is now the derivative in [sigma, H1, H2, H3, H4]
        szlsf = size(dlsf)
        #dlsf_ghpar = make_array(szlsf[1],szlsf[2]-3,value=lsf[0]*0)
        dlsf_ghpar = make_array(szlsf[1],Horder+1,value=lsf[0]*0)
        dlsf_ghpar[:,0] = dlsf[:,2]
        if Horder > 0:
            dlsf_ghpar[:,1:Horder] = dlsf[:,4:Horder+3]

        # Loop through the GH functions
        for i in range(Horder):

          GHparind  = indgen(Porder[i]+1) + derivparcntr
          numGHpar  = long(Porder[i]+1)
          GHparind0 = derivparcntr
          GHparind1 = GHparind0 + numGHpar-1

          # We have the derivative wrt the GH coefficient
          # but we want it wrt the polynomial coefficient
          # GHpar = c_0 + c_1*x + c_2*x^2
          # So we need to multiply by the derivative of GHpar wrt
          # that parameter, e.g. d GHpar/ d c_i =  x^i
          # All that matters it the center, Xcenter1
         
          # dlsf_ghpar is [Nlines,Nlsfpix,5]
          # Xcenter1 is [Nlines,Nlsfpix]
          #num = [szlsf[1], szlsf[2], numGHpar]
          #dlsfgh[*,*,GHparind0:GHparind1] = REBIN(dlsf_ghpar[*,*,i],num) * $
          #                                  REBIN(Xcenter1,num)^REBIN(indgen(1,1,numGHpar),num)

          # This is 2x faster
          for j in range(numGHpar-1):
             dlsfgh[:,GHparind0+j] = dlsf_ghpar[:,i] * Xcenter1**j


          # increment the counter
          derivparcntr += numGHpar


        # The parameters in dlsf are: [height, center, sigma, H0, H1, H2, H3, H4,
        #                                W1, W2, ...]
        dlsf_wpar = dlsf[:,Horder+4::]

        # Loop through the Wing parameters
        for i in range(nWpar-1):

          Wparind  = indgen(WPorder[i]+1) + derivparcntr
          numWpar  = long(WPorder[i]+1)
          Wparind0 = derivparcntr
          Wparind1 = Wparind0 + numWpar-1

          # We have the derivative wrt the GH coefficient
          # but we want it wrt the polynomial coefficient
          # GHpar = c_0 + c_1*x + c_2*x^2
          # So we need to multiply by the derivative of GHpar wrt
          # that parameter, e.g. d GHpar/ d c_i =  x^i
          # All that matters it the center, Xcenter1

          #dlsfgh[*,Wparind0:Wparind1] = dlsf_wpar[*,i] # (Xcenter1^indgen(WPorder[i]+1))

          # This is 2x faster
          for j in range(numWpar-1):
             dlsfgh[:,Wparind0+j] = dlsf_wpar[:,i] * Xcenter1**j

          # increment the counter
          derivparcntr += long( WPorder[i]+1 )


      # Return normal LSF derivatives
      #  This is the derivative of the LSF+wing parameters evaluated at
      #  each Xcenter.  This is *not* of the global polynomial coefficients
      else:
        dlsfgh = dlsf

      # Format differently for INPTYPE=2
      if inptype == 2:
        dlsfgh_temp = dlsfgh
        dlsfgh_sz   = size(dlsfgh)
        dlsfgh      = make_array(x_sz[1],x_sz[2],dlsfgh_sz[2],value=dlsf[0]*0)
        for i in range(dlsfgh_sz[2]-1):
            dlsfgh[:,:,i]=dlsfgh_temp[:,i]

    # Change format to 2D for INPTYPE=2
    if inptype == 2:
      lsf_temp = lsf
      lsf      = dblarr(x_sz[1],x_sz[2])
      lsf[:,:] = lsf_temp

    if stp == True: return 0

    return lsf


