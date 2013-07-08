'''
TurbulenceTools.py

Collection of tools/functions to do spectral analysis from time series extracted from
a turbulent flow field

Functions included:
    nextpow2(i)
    dofft(sig,samplefrq)
    movingave(x, window_len)
    movingavesmart(x,window_len,window='flat')  #disabled!!
    xcorr(x, y=None, maxlags=None, norm='ceoff')
'''

#===========================================================================#
# load modules
#===========================================================================#
#standard modules
import math
#scientific modules
import numpy as np
import scipy.fftpack as spfft
import scipy.signal as spsig
import scipy as sp
import pylab as pl

#===========================================================================#
# functions
#===========================================================================#

def nextpow2(i):
    """
    Find the next power of two. Python shell output:
    
    Examples:
        >>>nextpow2(250)
        256
        >>>rtn = nextpow2(i)
        
    Arguments:
        i = [float or int] any number...
        
    Returns:
        rtn = [int] the next power of two
    """
    # do not use numpy here, math is much faster for single values
    buf = math.ceil(math.log(i)/math.log(2))
    return int(math.pow(2,buf))


def dofft(sig,samplefrq):
    """
    Do fft. Outputs are the frq and amp, used by 99% of people doing fft.
    
    Example:
        >>> frq,amp = dofft(sig=mySignal, samplefrq=mySampleFrq)
        
    Arguments:
        sig: [numpy.array]signal
        samplefrq: [int or float]sample frequency
        
    Returns:
        frq: [numpy.array]frequencies
        amp: [numpy.array]amplitudes
    """
    siglength = sig.shape[0]
    NFFT = nextpow2(siglength)
    #    NFFT = nextpow2(sig.shape(-1,)[0])
    amp = spfft.fft(sig.reshape(-1,),NFFT)/siglength
    frq = samplefrq/2*np.linspace(0,1,NFFT/2+1)
    #    frq = spfft.fftfreq(siglength, (samplefrq)**-1)
    #    frq = frq[:NFFT/2+1]
    amp = 2*abs(amp[:NFFT/2+1])   
    return frq,amp


def movingave(x, window_len):
    '''
    moving average of the signal x. This function has important side effects 
    (from http://stackoverflow.com/a/11352216/2219303)
    
    Example:
        >>> ave = movingave(x,window_len)
    
    Arguments:
        x: [numpy.array] signal
        window_len: [integer] windows length
    
    Returns:
        ave: [numpy.array]
    '''
    w = np.ones(int(window_len))/float(window_len)
    return np.convolve(x,w,'same')


    #def movingavesmart(x,window_len,window='flat'):
    #"""smooth the data using a window with requested size.
    #(from http://www.scipy.org/Cookbook/SignalSmooth)
       
       #This method is based on the convolution of a scaled window with the signal.
       #The signal is prepared by introducing reflected copies of the signal 
       #(with the window size) in both ends so that transient parts are minimized
       #in the begining and end part of the output signal.
       
       #input:
           #x: the input signal 
           #window_len: the dimension of the smoothing window; should be an odd integer
           #window: the type of window from 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'
               #flat window will produce a moving average smoothing.

       #output:
           #the smoothed signal
           
       #see also: 
       #numpy.hanning, numpy.hamming, numpy.bartlett, numpy.blackman, numpy.convolve
       #scipy.signal.lfilter
       #"""
       ##test input
       #if x.ndim != 1:
           #raise ValueError, "smooth only accepts 1 dimension arrays."
       #if x.size < window_len:
           #raise ValueError, "Input vector needs to be bigger than window size."
       #if window_len<3:
           #return x
       #if not window in ['flat', 'hanning', 'hamming', 'bartlett', 'blackman']:
           #raise ValueError, "Window is on of 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'"
       
       ##the interressant stuff
       #s=np.r_[x[window_len-1:0:-1],x,x[-1:-window_len:-1]]
       #print(s.shape)
       #if window == 'flat': #moving average
           #w=np.ones(window_len,'d')
       #else:
           #w=eval('numpy.'+window+'(window_len)')
       #y=np.convolve(w/w.sum(),s,mode='valid')
       #return y


def xcorr(x, y=None, maxlags=None, norm='ceoff',doDetrend=False):
    '''
    Cross-correlation using scipy.correlate
    
    modified 
    copy from
    http://subversion.assembla.com/svn/PySpectrum/trunk/src/spectrum/correlation.py
    
    Estimates the cross-correlation (and autocorrelation) sequence of a random
    process of length N. By default, there is no normalisation and the output
    sequence of the cross-correlation has a length 2*N+1. 
    
    :param array x: first data array of length N
    :param array y: second data array of length N. If not specified, computes the 
        autocorrelation. 
    :param int maxlags: compute cross correlation between [-maxlags:maxlags]
        when maxlags is not specified, the range of lags is [-N+1:N-1].
    :param str option: normalisation in ['biased', 'unbiased', None, 'coeff']
     
    The true cross-correlation sequence is
    
    .. math:: r_{xy}[m] = E(x[n+m].y^*[n]) = E(x[n].y^*[n-m])

    However, in practice, only a finite segment of one realization of the 
    infinite-length random process is available.
    
    The correlation is estimated using numpy.correlate(x,y,'full'). 
    Normalisation is handled by this function using the following cases:

        * 'biased': Biased estimate of the cross-correlation function
        * 'unbiased': Unbiased estimate of the cross-correlation function
        * 'coeff': Normalizes the sequence so the autocorrelations at zero 
           lag is 1.0.

    :return:
        * a numpy.array containing the cross-correlation sequence (length 2*N-1)
        * lags vector
        
    .. note:: If x and y are not the same length, the shorter vector is 
        zero-padded to the length of the longer vector.
    '''
    N = len(x)
    if y == None:
        y = x
        
    if doDetrend:
        x=spsig.detrend(x)
        y=spsig.detrend(y)
        
    assert len(x) == len(y), 'x and y must have the same length. Add zeros if needed'
    assert maxlags <= N, 'maxlags must be less than data length'
    
    if maxlags == None:
        maxlags = N-1
        lags = np.arange(0, 2*N-1)
    else:
        assert maxlags < N
        lags = np.arange(N-maxlags-1, N+maxlags)
    
    res = sp.correlate(x, y, mode='full')
    
    if norm == 'biased':
        Nf = float(N)
        res = res[lags] / float(N)    # do not use /= !! 
    elif norm == 'unbiased':
        res = res[lags] / (float(N)-abs(np.arange(-N+1, N)))[lags]
    elif norm == 'coeff':        
        Nf = float(N)
        rms = pl.rms_flat(x) * pl.rms_flat(y)
        #rms = (np.mean(x**2)*np.mean(y**2))**(0.5)
        res = res[lags] / rms / Nf
    else:
        res = res[lags]

    lags = np.arange(-maxlags, maxlags+1)        
    return res, lags
    

def xcorr_fft(x, y=None, maxlags=None, norm='ceoff',doDetrend=False):
    '''
    test to use fft for cross correleation
    '''

    N = len(x)
    if y == None:
        y = x
        
    if doDetrend:
        x=spsig.detrend(x)
        y=spsig.detrend(y)
        
    assert len(x) == len(y), 'x and y must have the same length. Add zeros if needed'
    assert maxlags <= N, 'maxlags must be less than data length'
    
    if maxlags == None:
        maxlags = N-1
        lags = np.arange(0, 2*N-1)
    else:
        assert maxlags < N
        lags = np.arange(N-maxlags-1, N+maxlags)
    
    res = spsig.fftconvolve(x, y[::-1], mode="full")
    
    if norm == 'biased':
        Nf = float(N)
        res = res[lags] / float(N)    # do not use /= !! 
    elif norm == 'unbiased':
        res = res[lags] / (float(N)-abs(np.arange(-N+1, N)))[lags]
    elif norm == 'coeff':        
        Nf = float(N)
        rms = pl.rms_flat(x) * pl.rms_flat(y)
        if rms==0:
            rms=1
        #rms = (np.mean(x**2)*np.mean(y**2))**(0.5)
        res = res[lags] / rms / Nf
    else:
        res = res[lags]
    
    res=res[(len(res)-1)/2:-1]
    lags = np.arange(0, maxlags)        
    return res, lags
