import numpy as np
from pyFlowStat import Statistics
from pyFlowStat import Math
from pyFlowStat import TurbulenceTools as tt


'''
Methods to calculate the intergal scale from correlation coefficients
'''    
def calcTii_fitExp(dt,rii):
    t=fitExp(rii)
    return t*dt
    
def calcTii_exp(dt,rii):
    t=xcut_value(rii,1.0/np.exp(1.0))
    return t*dt
    
def calcTii_intMin(dt,rii):
    x_cut=xcut_min(rii)
    t=integrateUpTo(rii,x_cut)
    return t*dt
    
def calcTii_intZero(dt,rii):
    x_cut=xcut_value(rii,value=0.0)
    t=integrateUpTo(rii,x_cut)
    return t*dt
    
def calcTii_intFull(dt,rii):
    t=integrateUpTo(rii)
    return t*dt
    
def calcTii_intConvTail(dt,rii,z=1.96):
    x_cut=xcut_conv(len(rii),rii,z=z)
    t_tail=intExpTail(rii,x_cut)
    t_int=integrateUpTo(rii,x_cut)
    return (t_tail+t_int)*dt
    
def calcTii_intTail(dt,rii,value=0.0):
    x_cut=xcut_value(rii,value=0.0)
    t_tail=intExpTail(rii,x_cut)
    t_int=integrateUpTo(rii,x_cut)
    return (t_tail+t_int)*dt
'''
Helper functions
'''
def fitExp(r11,maxlag=None):
    if maxlag==None:
        maxlag=len(r11)
    xdata=np.arange(len(r11[:maxlag]))
    ydata=np.array(r11[:maxlag])
    try:
        popt, pcov = tt.fit_exp_correlation(xdata,ydata)
        T=popt
    except RuntimeError:
        print("Error - curve_fit failed")
        T=0
    return T
    
def remove(N,r11,z=1.96):
    r_out=[]
    r_out.append(r11[0])
    for i in range(1,len(r11)):
        varr=Statistics.VarR_k(N,r11,i,z=z)
        if r11[i]<varr:
            break
        else:
            r_out.append(r11[i])
    return r_out
    
def removeZero(r11):
    r_out=[]
    r_out.append(r11[0])
    for i in range(1,len(r11)):
        if r11[i]<0:
            break
        else:
            r_out.append(r11[i])
    return r_out
    
def removeMin(r11):
    r_out=[]
    r_out.append(r11[0])
    for i in range(1,len(r11)):
        if r11[i-1]<r11[i]:
            break
        else:
            r_out.append(r11[i])
    return r_out
    
def xcut_min(r11):
    for i in range(1,len(r11)):
        if r11[i-1]<r11[i]:
            return i-1
    return len(r11)-1
    
def xcut_value(r11,value=0.0):
    #lag=range(len(r11))
    for i in range(1,len(r11)):
        if r11[i]<value:
            xi=Math.interpy_lin_1d(i-1,i,r11[i-1],r11[i],value)
            return xi
    return len(r11)-1
    
def xcut_conv(N,r11,z=1.96):
    #lag=range(len(r11))
    for i in range(1,len(r11)):
        varr=Statistics.VarR_k(N,r11,i,z=z)
        if r11[i]<varr:
            xi=Math.interpy_lin_1d(i-1,i,r11[i-1],r11[i],varr)
            return xi
    return np.nan
    
def removeAdd(N,r11,z=1.96):
    r_out=remove(N,r11,z=z)
    x_cut=len(r_out)
    a=-len(r_out)/np.log(r_out[-1])
    print a
    for i in range(len(r_out)+1,N):
        r_out.append(np.exp(-float(i)/a))
    
    L_int=np.trapz(r_out)
    L_exp=a*np.exp(-float(x_cut)/a)
    return r_out,L_int,L_exp,x_cut
    
def intExpTail(rii,x_cut):
    if x_cut >= (len(rii)-1):
        return 0.0
    y_cut=Math.interpx_lin(rii,x_cut)
    a=-x_cut/np.log(y_cut)
    L_exp=a*np.exp(-float(x_cut)/a)
    return L_exp
    
def integrateUpTo(rii,x_max=None):
    if x_max==None:
        return np.trapz(rii)
    else:
        x,y=cutUpTo(rii,x_max)
        return np.trapz(y=y,x=x)

def cutUpTo(rii,x_max):
    if x_max==None:
        x_max=len(rii)-1
        
    nrpoints=np.floor(x_max)+1
    endtime=np.floor(x_max)
    xi_end=int(np.floor(x_max))
    xtmp=np.linspace(0,endtime,nrpoints)
    r_temp=rii[:xi_end+1]
    if x_max>np.floor(x_max):
        
        r_last=Math.interpx_lin_1d(xi_end,xi_end+1,rii[xi_end],rii[xi_end+1],x_max)
        xtmp=np.append(xtmp,x_max)
        r_temp=np.append(r_temp,r_last)
        
    return xtmp,r_temp