'''
TurbulenceTools.py

Collection of tools/functions to do spectral analysis from time series
extracted from a turbulent flow field

Functions included:
    *nextpow2*
     calculate the next power of 2. 
      
    *dofft*
     do an fft.
      
    *movingave*
     Moving average on a signal x with a user defined window size.
      
    *xcorr*
     Cross-corelation of two signal x and y. If x=y, then the auto correlation
     is computed.
    
    *xcorr_fft*
     Same as xcorr but much faster by using a fft.
    
    *twoPointCorr*
    
    *func_exp_correlation*
    
    *func_gauss_correlation*
    
    *fit_exp_correlation*
    
    *calcInegarlScale_expFit*
    
    *calcInegarlScale_trapz*
    
    *calcInegarlScale_simps*
    
    *calcIntegarlScale*
    
    *bandpass*
     Butterworth-Bandpass Filter.
     
    
    *bandstop*
     Butterworth-Bandstop Filter.
    
    *lowpass*
     Butterworth-Lowpass Filter.
    
    *highpass*
    Butterworth-Lowpass Filter.
    

Notes:
    *The functions "bandpass", "bandstop", "lowpass" and "highpass" are copied
     from the python library ObsPy". ObsPy is an open-source project dedicated
     to provide a Python framework for processing seismological data. See the
     following links:
     opspy website: https://github.com/obspy/obspy/wiki
     opspy.filter: http://docs.obspy.org/_modules/obspy/signal/filter.html
'''

#===========================================================================#
# load modules
#===========================================================================#
#standard modules
import warnings
import math

#scientific modules
import numpy as np
import scipy.fftpack as spfft
import scipy.signal as spsig
import scipy as sp
#import pylab as pl
from scipy.optimize import curve_fit
from scipy.integrate import simps

from pyFlowStat import Statistics as stat

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
        * i: [float or int] any number...

    Returns:
        * rtn: [int] the next power of two
    """
    # do not use numpy here, math is much faster for single values
    buf = math.ceil(math.log(i)/math.log(2))
    return int(math.pow(2,buf))


def dofft(sig,samplefrq,nperseg=512,detrend='constant'):
    """
    Estimate power spectral density using Welch's method. For more informations
    on the Welch's method implemented in python, visit:

    http://docs.scipy.org/doc/scipy-dev/reference/generated/scipy.signal.welch.html

    Example:
        >>> frq,psd = dofft(sig=mySignal, samplefrq=mySampleFrq)

    Arguments:
        * sig: [numpy.array] signal
        * samplefrq: [int or float] sample frequency

    Returns:
        * frq: [numpy.array]frequencies
        * psd: [numpy.array]amplitudes
    """
    # old crap...
    #siglength = sig.shape[0]
    #NFFT = nextpow2(siglength)
    ##    NFFT = nextpow2(sig.shape(-1,)[0])
    #amp = spfft.fft(sig.reshape(-1,),NFFT)/siglength
    #frq = samplefrq/2*np.linspace(0,1,NFFT/2+1)
    ##    frq = spfft.fftfreq(siglength, (samplefrq)**-1)
    ##    frq = frq[:NFFT/2+1]
    #amp = 2*abs(amp[:NFFT/2+1])
    #return frq,amp

    frq, psd = spsig.welch(sig,fs=samplefrq, window='hanning', noverlap=None, nperseg=nperseg, return_onesided=True, scaling='density', axis=-1,detrend=detrend)
    return frq, psd

def movingave(x, window_len):
    '''
    moving average of the signal x. This function has important side effects
    (from http://stackoverflow.com/a/11352216/2219303)

    Example:
        >>> ave = movingave(x,window_len)

    Arguments:
        * x: [numpy.array] signal
        * window_len: [integer] windows length

    Returns:
        * ave: [numpy.array]
    '''
    w = np.ones(int(window_len))/float(window_len)
    return np.convolve(x,w,'same')

def smooth(x,window_len=5,window='hanning'):
    if x.ndim != 1:
            raise ValueError, "smooth only accepts 1 dimension arrays."
    if x.size < window_len:
            raise ValueError, "Input vector needs to be bigger than window size."
    if window_len<3:
            return x
    if not window in ['flat', 'hanning', 'hamming', 'bartlett', 'blackman']:
            raise ValueError, "Window is on of 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'"
    s=np.r_[2*x[0]-x[window_len-1::-1],x,2*x[-1]-x[-1:-window_len:-1]]
    if window == 'flat': #moving average
            w=np.ones(window_len,'d')
    else:
            w=eval('np.'+window+'(window_len)')
    y=np.convolve(w/w.sum(),s,mode='same')
    return y[window_len:-window_len+1]
#def movingavesmart(x,window_len,window='flat'):
#    '''
#    smooth the data using a window with requested size.
#    (from http://www.scipy.org/Cookbook/SignalSmooth)
#
#    This method is based on the convolution of a scaled window with the signal.
#    The signal is prepared by introducing reflected copies of the signal
#    (with the window size) in both ends so that transient parts are minimized
#    in the begining and end part of the output signal.
#
#    input:
#    x: the input signal
#    window_len: the dimension of the smoothing window; should be an odd integer
#    window: the type of window from 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'
#    flat window will produce a moving average smoothing.
#
#    output:
#    the smoothed signal
#
#    see also:
#    numpy.hanning, numpy.hamming, numpy.bartlett, numpy.blackman, numpy.convolve
#    scipy.signal.lfilter
#    '''
#    #test input
#    if x.ndim != 1:
#        raise ValueError, "smooth only accepts 1 dimension arrays."
#    if x.size < window_len:
#        raise ValueError, "Input vector needs to be bigger than window size."
#    if window_len<3:
#        return x
#    if not window in ['flat', 'hanning', 'hamming', 'bartlett', 'blackman']:
#        raise ValueError, "Window is on of 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'"
#
#    #the interressant stuff
#    s=np.r_[x[window_len-1:0:-1],x,x[-1:-window_len:-1]]
#    print(s.shape)
#    if window == 'flat': #moving average
#        w=np.ones(window_len,'d')
#    else:
#        w=eval('numpy.'+window+'(window_len)')
#    y=np.convolve(w/w.sum(),s,mode='valid')
#    return y


def xcorr(x, y=None, maxlags=None, norm='ceoff',doDetrend=False):
    '''
    Cross-correlation using scipy.correlate

    copy from http://subversion.assembla.com/svn/PySpectrum/trunk/src/spectrum/correlation.py
    and futher modified for flow analysis.

    Estimates the cross-correlation (and autocorrelation) sequence of a random
    process of length N. By default, there is no normalisation and the output
    sequence of the cross-correlation has a length 2*N+1.

    Arguments:
        * x: first data array of length N
        * y: second data array of length N. If not specified, computes the
        autocorrelation.
        * maxlags: compute cross correlation between [-maxlags:maxlags]
        when maxlags is not specified, the range of lags is [-N+1:N-1].
        * norm: ['biased', 'unbiased', None, 'coeff'] normalisation
        * doDetrend: [bool()] do a data detrend. Useful from data from mesurment

    The true cross-correlation sequence is:
        * r_{xy}[m] = E(x[n+m].y^*[n]) = E(x[n].y^*[n-m])

    However, in practice, only a finite segment of one realization of the
    infinite-length random process is available.

    The correlation is estimated using numpy.correlate(x,y,'full').
    Normalisation is handled by this function using the following cases:

        * 'biased': Biased estimate of the cross-correlation function
        * 'unbiased': Unbiased estimate of the cross-correlation function
        * 'coeff': Normalizes the sequence so the autocorrelations at zero lag is 1.0.

    returns:
        * xcorr: [np.array, shape=(N-1,1)]a numpy.array containing the cross-correlation sequence
        * lags: [np.array, shape=(N-1,1)] lag vector

    notes:
        * If x and y are not the same length, the shorter vector is
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
        rms = stat.rms(x) * stat.rms(y)
        #rms = (np.mean(x**2)*np.mean(y**2))**(0.5)
        res = res[lags] / rms / Nf
    else:
        res = res[lags]

    lags = np.arange(-maxlags, maxlags+1)
    return res, lags


def xcorr_fft(x, y=None, maxlags=None, norm='coeff',doDetrend=False,oneSided=True):
    '''
    Cross-correlation using scipy.fftconvolve. Similar returns as TurbulenceTools.xcorr()
    but faster with fftconvolve.

    copy from http://subversion.assembla.com/svn/PySpectrum/trunk/src/spectrum/correlation.py
    and futher modified for flow analysis.

    Estimates the cross-correlation (and autocorrelation) sequence of a random
    process of length N. By default, there is no normalisation and the output
    sequence of the cross-correlation has a length 2*N+1.

    Arguments:
        * x: first data array of length N
        * y: second data array of length N. If not specified, computes the
        autocorrelation.
        * maxlags: compute cross correlation between [-maxlags:maxlags]
        when maxlags is not specified, the range of lags is [-N+1:N-1].
        * norm: ['biased', 'unbiased', None, 'coeff'] normalisation
        * doDetrend: [bool()] do a data detrend. Useful from data from mesurment

    The true cross-correlation sequence is:
        * r_{xy}[m] = E(x[n+m].y^*[n]) = E(x[n].y^*[n-m])

    However, in practice, only a finite segment of one realization of the
    infinite-length random process is available.

    The correlation is estimated using numpy.correlate(x,y,'full').
    Normalisation is handled by this function using the following cases:

        * 'biased': Biased estimate of the cross-correlation function
        * 'unbiased': Unbiased estimate of the cross-correlation function
        * 'coeff': Normalizes the sequence so the autocorrelations at zero lag is 1.0.

    returns:
        * xcorr: [np.array, shape=(N-1,1)]a numpy.array containing the cross-correlation sequence
        * lags: [np.array, shape=(N-1,1)] lag vector

    notes:
        * If x and y are not the same length, the shorter vector is
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

    res = spsig.fftconvolve(x, y[::-1], mode="full")

    if norm == 'biased':
        Nf = float(N)
        res = res[lags] / float(N)    # do not use /= !!
    elif norm == 'unbiased':
        res = res[lags] / (float(N)-abs(np.arange(-N+1, N)))[lags]
    elif norm == 'coeff':
        Nf = float(N)
        rms = stat.rms(x) * stat.rms(y)
        if rms==0:
            rms=1
        #rms = (np.mean(x**2)*np.mean(y**2))**(0.5)
        res = res[lags] / rms / Nf
    else:
        res = res[lags]
    if oneSided:
        res=res[(len(res)-1)/2:-1]
        lags = np.arange(0, maxlags)
    return res, lags
    
def twoPointCorr(x,y,subtractMean=True,norm=False):
    '''
    dot product of two vectors to claculate two point correlation
    '''
    #return np.dot(x,y)/(rms(x)*rms(y))/len(x)
    if subtractMean==True:
        x_prime=x-np.mean(x)
        y_prime=y-np.mean(y)
        #x_prime=scipy.signal.detrend(np.nan_to_num(x))
        #y_prime=scipy.signal.detrend(np.nan_to_num(y))
    else:
        x_prime=x
        y_prime=y
    
    #cc=np.correlate(x_prime,y_prime)
    #print x_prime.shape
    
    #cc=np.dot(smooth(x_prime,window_len=11),smooth(y_prime,window_len=11))
    cc=np.dot(x_prime,y_prime)
    #print cc.shape
    if norm==True:
        return cc/np.std(x_prime)/np.std(y_prime)/len(x)
    else:
        return cc
    #return np.correlate(x_prime,y_prime)
    
def twoPointCorr_fast(x_prime,y_prime,cc_max):
    '''
    dot product of two vectors to claculate normalized two point correlation
    '''
    cc=np.dot(x_prime,y_prime)
    return cc/cc_max
    
def func_exp_correlation(x, a):
    np.seterr('ignore')
    res = np.exp(-x/a)

    #print res
    return res
    
def func_gauss_correlation(x, a):
    np.seterr('ignore')
    res = np.exp(-np.pi*x**2/(a**2*4))

    #print res
    return res
    
def func_dblExp_correlation(x,a):
    if len(np.unique(np.diff(x)))>1:
        warnings.warn('x has to be evenly spaced')
    def _func_dblExp_correlation(x,a):
                if len(x)>1:
                    x=np.r_[x,x[-1]+np.diff(x)[-1]]
                    k_exp=func_exp_correlation(x,a)
                    kk_exp=spsig.fftconvolve(k_exp,k_exp, mode="full")
                    res,_=xcorr_fft(kk_exp)
                    return res[::2]
                else:
                    return []
                
    if np.any(x<0):
        x_plus=x[x>=0]
        x_neg=np.abs(x[x<=0])[::-1]

        res_plus=_func_dblExp_correlation(x_plus,a)
        res_neg=_func_dblExp_correlation(x_neg,a)
        res=np.r_[res_neg[::-1],res_plus[1:]]
    else:
        res=_func_dblExp_correlation(x,a)
    return res
    
def fit_exp_correlation(xdata,ydata):
    '''
    Fits an exponential function of shape exp(-x/a) to the data and returns a
    
    Arguments:
        * xdata: x-values (e.g lags)
        * ydata: y-values (e.g auto correlation coefficient)
        
    returns:
        * a: fitter parameter a
        * pcov: The estimated covariance of a.
    
    '''
    
    popt, pcov = curve_fit(func_exp_correlation,xdata,ydata)
    a=popt[0]
    return a,pcov
    
def fit_gauss_correlation(xdata,ydata):
    '''
    Fits an exponential function of shape exp(-x/a) to the data and returns a
    
    Arguments:
        * xdata: x-values (e.g lags)
        * ydata: y-values (e.g auto correlation coefficient)
        
    returns:
        * a: fitter parameter a
        * pcov: The estimated covariance of a.
    
    '''
    
    popt, pcov = curve_fit(func_gauss_correlation,xdata,ydata)
    a=popt[0]
    return a,pcov
    
def calcInegralScale_expFit(rho_i,dx):
    lags=np.arange(len(rho_i))
    try:
        popt, pcov = fit_exp_correlation(lags,rho_i)
        #self.data['Txx']=abs(popt[0])*np.sqrt(np.pi)*0.5*self.data['dt']
        L=popt*dx
    except RuntimeError:
        print("Error - curve_fit failed")
        L=np.nan
    return L,(pcov)
    
def calcInegralScale_trapz(rho_i,dx):
    L=np.trapz(y=rho_i, x=None, dx=dx, axis=-1)
    return L,()
    
def calcInegralScale_simps(rho_i,dx):
    L=simps(rho_i,dx=dx)
    return L,()
    
def calcIntegralScale(rho_i,dx=1.0,method=None):
    '''
    calculates the integral scale
    
    Arguments:
        * rho_i: [array] Normalized autocorrelation coefficients, starting with lag=0, and rho_i=1.0
        * method: [string] Speciefies the method to be used.
        
    returns:
        * L: integral scale
        * params: [tuple] additional output parameters
    '''
    rho_i=np.array(rho_i)
    #lags = np.linspace(0,(len(rho_i)-1.0)*dx,len(rho_i))
    L,params=method(rho_i,dx)
    return L,params


#-----------------------------------------------------------------------------#
def bandpass(data, freqmin, freqmax, df, corners=4, zerophase=False,axis=-1):
    """
    Butterworth-Bandpass Filter.

    Filter data from ``freqmin`` to ``freqmax`` using ``corners`` corners.

    :param data: Data to filter, type numpy.ndarray.
    :param freqmin: Pass band low corner frequency.
    :param freqmax: Pass band high corner frequency.
    :param df: Sampling rate in Hz.
    :param corners: Filter corners / orders.
    :param zerophase: If True, apply filter once forwards and once backwards.
        This results in twice the number of corners but zero phase shift in
        the resulting filtered trace.
    :return: Filtered data.
    """
    fe = 0.5 * df
    low = freqmin / fe
    high = freqmax / fe
    # raise for some bad scenarios
    if high > 1:
        high = 1.0
        msg = "Selected high corner frequency is above Nyquist. " + \
              "Setting Nyquist as high corner."
        warnings.warn(msg)
    if low > 1:
        msg = "Selected low corner frequency is above Nyquist."
        raise ValueError(msg)
    [b, a] = spsig.iirfilter(corners, [low, high], btype='band',
                       ftype='butter', output='ba')
    if zerophase:
        firstpass = spsig.lfilter(b, a, data,axis=axis)
        return spsig.lfilter(b, a, firstpass[::-1],axis=axis)[::-1]
    else:
        return spsig.lfilter(b, a, data,axis=axis)


def bandstop(data, freqmin, freqmax, df, corners=4, zerophase=False,axis=-1):
    """
    Butterworth-Bandstop Filter.

    Filter data removing data between frequencies ``freqmin`` and ``freqmax``
    using ``corners`` corners.

    :param data: Data to filter, type numpy.ndarray.
    :param freqmin: Stop band low corner frequency.
    :param freqmax: Stop band high corner frequency.
    :param df: Sampling rate in Hz.
    :param corners: Filter corners / orders.
    :param zerophase: If True, apply filter once forwards and once backwards.
        This results in twice the number of corners but zero phase shift in
        the resulting filtered trace.
    :return: Filtered data.
    """
    fe = 0.5 * df
    low = freqmin / fe
    high = freqmax / fe
    # raise for some bad scenarios
    if high > 1:
        high = 1.0
        msg = "Selected high corner frequency is above Nyquist. " + \
              "Setting Nyquist as high corner."
        warnings.warn(msg)
    if low > 1:
        msg = "Selected low corner frequency is above Nyquist."
        raise ValueError(msg)
    [b, a] = spsig.iirfilter(corners, [low, high],
                       btype='bandstop', ftype='butter', output='ba')
    if zerophase:
        firstpass = spsig.lfilter(b, a, data,axis=axis)
        return spsig.lfilter(b, a, firstpass[::-1],axis=axis)[::-1]
    else:
        return spsig.lfilter(b, a, data,axis=axis)


def lowpass(data, freq, df, corners=4, zerophase=False,axis=-1):
    """
    Butterworth-Lowpass Filter.

    Filter data removing data over certain frequency ``freq`` using ``corners``
    corners.

    :param data: Data to filter, type numpy.ndarray.
    :param freq: Filter corner frequency.
    :param df: Sampling rate in Hz.
    :param corners: Filter corners / orders.
    :param zerophase: If True, apply filter once forwards and once backwards.
        This results in twice the number of corners but zero phase shift in
        the resulting filtered trace.
    :return: Filtered data.
    """
    fe = 0.5 * df
    f = freq / fe
    # raise for some bad scenarios
    if f > 1:
        f = 1.0
        msg = "Selected corner frequency is above Nyquist. " + \
              "Setting Nyquist as high corner."
        warnings.warn(msg)
    [b, a] = spsig.iirfilter(corners, f, btype='lowpass', ftype='butter',
                       output='ba')
    if zerophase:
        firstpass = spsig.lfilter(b, a, data,axis=axis)
        return spsig.lfilter(b, a, firstpass[::-1],axis=axis)[::-1]
    else:
        return spsig.lfilter(b, a, data,axis=axis)


def highpass(data, freq, df, corners=4, zerophase=False,axis=-1):
    """
    Butterworth-Highpass Filter.

    Filter data removing data below certain frequency ``freq`` using
    ``corners`` corners.

    :param data: Data to filter, type numpy.ndarray.
    :param freq: Filter corner frequency.
    :param df: Sampling rate in Hz.
    :param corners: Filter corners / orders.
    :param zerophase: If True, apply filter once forwards and once backwards.
        This results in twice the number of corners but zero phase shift in
        the resulting filtered trace.
    :return: Filtered data.
    """
    fe = 0.5 * df
    f = freq / fe
    # raise for some bad scenarios
    if f > 1:
        msg = "Selected corner frequency is above Nyquist."
        raise ValueError(msg)
    [b, a] = spsig.iirfilter(corners, f, btype='highpass', ftype='butter',
                       output='ba')
    if zerophase:
        firstpass = spsig.lfilter(b, a, data,axis=axis)
        return spsig.lfilter(b, a, firstpass[::-1],axis=axis)[::-1]
    else:
        return spsig.lfilter(b, a, data,axis=axis)
#-----------------------------------------------------------------------------#