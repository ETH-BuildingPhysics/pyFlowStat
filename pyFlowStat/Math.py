#===========================================================================#
# load modules
#===========================================================================#
#standard modules
import sys

#scientific modules
import numpy as np

def interpy_lin_1d(x1,x2,y1,y2,yi):
    xi=(float(x2)-float(x1))/(float(y2)-float(y1))*(float(yi)-float(y1))+float(x1)
    return xi
    
def interpx_lin(yval,xi):
    xval=range(len(yval))
    yi=interpx_lin_1d(xval[int(np.floor(xi))],xval[int(np.ceil(xi))],yval[int(np.floor(xi))],yval[int(np.ceil(xi))],xi)
    return yi
    
def interpx_lin_1d(x1,x2,y1,y2,xi):
    yi=(float(y2)-float(y1))/(float(x2)-float(x1))*(float(xi)-float(x1))+float(y1)
    return yi
    
def interp1exp(x,y,xi):
    '''
    linear interpolation with linear extrapolation
    '''
    val=np.interp(xi,x,y,left=np.inf*-1.0,right=np.inf)
    if np.isneginf(val):
        return np.polyval(np.polyfit(x[:2],y[:2],1),xi)
    elif np.isposinf(val):
        return np.polyval(np.polyfit(x[-2:],y[-2:],1),xi)
    return val