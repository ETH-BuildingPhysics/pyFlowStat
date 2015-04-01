#===========================================================================#
# load modules
#===========================================================================#
#standard modules
import sys

#scientific modules
import numpy as np

def sample_autocovariance(x_in,maxlag):
    #to test the FFT implementation in TurbulenceTools
    x=scipy.signal.detrend(x_in)
    N=len(x)
    rho=[]
    for k in range(maxlag):
        tmpsum=[]
        for i in range(N-k):
            tmpsum.append(x[i+k]*x[i])
        rho.append(sum(tmpsum)/N)
    return rho/rho[0]
    
def NeffFactor(r1):
    return (1.0-r1)/(1.0+r1)
    
def VarRk(N,r11,z=1.96):
    rsq=np.array(r11)
    rsq=rsq**2
    varR=[]
    varR.append(np.nan)
    
    tmpsum=0
    for i in range(1,len(r11)):
        tmpsum=tmpsum+rsq[i]
        #tmpsum=np.sum(rsq[:i])
        varR.append(z*np.sqrt(1.0/N*(1.0+2.0*tmpsum)))
    return np.array(varR)
    
def VarR_k(N,r11,k,z=1.96):
    if k==0:
        return np.nan
    rsq=np.array(r11)
    rsq=rsq**2
    tmpsum=np.sum(rsq[1:(k+1)])
    #tmpsum=sum(r**2 for r in r11[:k])
    return z*np.sqrt(1.0/N*(1.0+2.0*tmpsum))
    
def r_conv(N_eff,z=1.96):
    return z/np.sqrt(N_eff)
    
    
def SE_r(N,r11):
    rsq=np.array(r11)
    rsq=rsq**2
    varR=[]
    varR.append(np.nan)
    varR.append(1.0/np.sqrt(N))
    tmpsum=0
    for i in range(2,len(r11)):
        tmpsum=tmpsum+rsq[i-1]
        #tmpsum=np.sum(rsq[:i])
        varR.append(np.sqrt(1.0/N*(1.0+2.0*tmpsum)))
    return np.array(varR)
    
def SE_r_k(N,r11,k):
    #Bartlett's formula for MA(l) processes, http://en.wikipedia.org/wiki/Correlogram 
    if k==0:
        return np.nan
    elif k==1:
        return 1.0/np.sqrt(N)
    rsq=np.array(r11)
    rsq=rsq**2
    tmpsum=np.sum(rsq[1:(k-1)+1])
    #tmpsum=sum(r**2 for r in r11[:k])
    return np.sqrt(1.0/N*(1.0+2.0*tmpsum))
   
def rms(x):
    return np.sqrt(np.mean(x**2))