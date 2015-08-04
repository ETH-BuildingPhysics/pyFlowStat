# -*- coding: utf-8 -*-
"""
Created on Mon Dec 01 10:13:26 2014

@author: mimmer
"""

#from . import Decomposition
import matplotlib as mpl
import matplotlib.pyplot as plt
#import matplotlib.gridspec as gridspec
import numpy as np

def plotEigenvals(myDMD,ax,k=1,idx=None,mark=False,cmap=plt.cm.jet,min_size=10.0,max_size=1000.0):
    ritz = myDMD.result['ritz_vals']**k
    norms = myDMD.result['mode_norms']
    norms_val=np.sqrt(norms)
    
    if idx==None:
        idx_sorted=np.argsort(norms)
    else:
        idx_sorted=idx
    
    circle2=plt.Circle((0,0),1.0,color='k',fill=False,zorder=1)
    ax.add_artist(circle2)

    A=np.min(norms_val[idx_sorted])
    B=np.max(norms_val[idx_sorted])
    s=min_size + (norms_val[idx_sorted]-A)*(max_size-min_size)/(B-A)
    print np.max(s),np.min(s)
    
    if mark:
        A=np.min(norms_val)
        B=np.max(norms_val)
        s=min_size + (norms_val-A)*(max_size-min_size)/(B-A)
        ax.scatter(np.real(ritz)[idx_sorted],np.imag(ritz)[idx_sorted],facecolors='none',edgecolors='orange',s=s[idx_sorted]+10,lw = 1.0,zorder=3)
        idx_sorted2=np.argsort(norms)
        ax.scatter(np.real(ritz)[idx_sorted2],np.imag(ritz)[idx_sorted2],c=norms_val[idx_sorted2],edgecolors='none',s=s[idx_sorted2],cmap=cmap,zorder=2,lw=0.5)
    else:
        ax.scatter(np.real(ritz)[idx_sorted],np.imag(ritz)[idx_sorted],c=norms_val[idx_sorted],edgecolors='none',s=s,zorder=2,cmap=cmap)
    
    ax.set_aspect('equal')
    ax.set_xlabel(r'$Re\{\lambda_i\}$')
    ax.set_ylabel(r'$Im\{\lambda_i\}$')
    ax.set_xlim([-1.2,1.2])
    ax.set_ylim([-1.2,1.2])


    #plt.tight_layout()
    
def plotSpectrum(myDMD,ax,k=1,idx=None,width=1.0,alpha=0.5,cmap=plt.cm.jet):
    norms = np.sqrt(myDMD.result['mode_norms'])

    stab = np.log(myDMD.getEigAbs())
    f=myDMD.getFrqList()
    
    if idx==None:
        idx_sorted=np.argsort(norms)
    else:
        idx_sorted=idx
        
    #norm=mpl.colors.Normalize(vmin=np.min(stab),vmax=np.max(stab))
    color_val=stab
    max_size=1.0
    min_size=0.0
    A=np.min(color_val[idx_sorted])
    B=np.max(color_val[idx_sorted])
    s=min_size + (color_val[idx_sorted]-A)*(max_size-min_size)/(B-A)
    mycolors=cmap(s)

    for i in range(len(idx_sorted)):
        ax.bar(f[idx_sorted[i]],norms[idx_sorted[i]],width=width,facecolor=mycolors[i],alpha=alpha)
        
    ax.set_xlabel(r'$f [Hz]$')
    ax.set_ylabel(r'$||\Phi||$')

def plotGrowth(myDMD,ax,k=1,dt=1.0,idx=None,reverseAxis=False):
    '''
    plots the real and imaginary part of lambda=log(ritz_val)/dt
    '''
    log_ritz = np.log(myDMD.result['ritz_vals'])/dt#/(2.0*np.pi)
    norms = np.sqrt(myDMD.result['mode_norms'])

    norms_val=np.sqrt(norms)

    
    if idx==None:
        idx_sorted=np.argsort(norms)
    else:
        idx_sorted=idx
    #cmax=np.sort(norms_val[idx_sorted])[::-1][3]
    cmin=0
    cmax=np.max(norms_val[idx_sorted])

    #ax.scatter(np.imag(log_ritz)[idx_sorted],np.real(log_ritz)[idx_sorted],c=norms_val[idx_sorted],edgecolors='none',s=np.sqrt(norms_val[idx_sorted])*10)
    if not reverseAxis:
        ax.scatter(np.imag(log_ritz)[idx_sorted]/(2.0*np.pi),np.real(log_ritz)[idx_sorted],c=norms_val[idx_sorted],norm=mpl.colors.Normalize(vmin=cmin,vmax=cmax),edgecolors='none',s=np.sqrt(norms_val[idx_sorted])*10)
    else:
        ax.scatter(np.real(log_ritz)[idx_sorted],np.imag(log_ritz)[idx_sorted]/(2.0*np.pi),c=norms_val[idx_sorted],norm=mpl.colors.Normalize(vmin=cmin,vmax=cmax),edgecolors='none',s=np.sqrt(norms_val[idx_sorted])*10)
    
    #ax.set_aspect('equal')
    if not reverseAxis: 
        ax.set_xlabel(r'$ln(\lambda)_i$')
        ax.set_ylabel(r'$ln(\lambda)_r$')
    else:
        ax.set_xlabel(r'$\gamma$ $[1/s]$')
        #ax.set_ylabel(r'$Im(\lambda)$ $[\omega/2\pi]$')
        ax.set_ylabel(r'$f$ $[Hz]$')
    #ax.set_xlim([-1.2,1.2])
    #ax.set_ylim([-1.2,1.2])

    plt.tight_layout()
    
def plotModes(myDMD,idxList,nRows=8,nCols=6,iStart=0,vmin=-0.5,vmax=0.5,extent=None):
    testMode=np.real(myDMD.getMode(0))
    sizef=testMode.shape[0]/float(testMode.shape[1])
    width=20.0
    plt.figure(figsize=(width,width/nCols*nRows*sizef))
    
    for i in range(min((nCols*nRows)-1,len(idxList))):
        j=idxList[i+iStart]
        ax=plt.subplot(nRows,nCols,i+1)
        ax.imshow(np.real(myDMD.getMode(j)),vmin=vmin,vmax=vmax,interpolation='nearest',extent=extent)
        #ax.imshow(np.real(FullFieldDMD.getMode(j)*np.real(FullFieldDMD.result['ritz_vals'][j])),vmin=vmin_mode,vmax=vmax_mode,interpolation='nearest')
        myDMD.getFrqStr(j,2)
        ax.set_title(str(j)+' '+(myDMD.getFrqStr(j,2) + r' |$\lambda$|=' + str(np.round(myDMD.getEigAbs(idx=j),6))+r' |$\phi$| '+ str(np.round(np.sqrt(myDMD.getNorm(j)),2))))
        #ax.set_title(str(j)+' '+(myDMD.getFrqStr(j,2) + r' |$\lambda$|=' + str(np.round(myDMD.getEigAbs(idx=j),6))))
                
        #ax.set_title(myDMD.getFrqStr(j,2))
    plt.tight_layout()
    
def plotModesPiv(myDMD,idxList,nRows=8,nCols=6,iStart=0,vmin=-0.5,vmax=0.5,component=0,extent=None):
    testMode=np.real(myDMD.getMode(0,component=0))
    sizef=testMode.shape[0]/float(testMode.shape[1])
    width=20.0
    plt.figure(figsize=(width,width/nCols*nRows*sizef))
    for i in range(min((nCols*nRows)-1,len(idxList))):
        j=idxList[i+iStart]
        ax=plt.subplot(nRows,nCols,i+1)
        ax.imshow(np.real(myDMD.getMode(j,component=component)),vmin=vmin,vmax=vmax,interpolation='nearest',extent=extent)
        #ax.imshow(np.real(FullFieldDMD.getMode(j)*np.real(FullFieldDMD.result['ritz_vals'][j])),vmin=vmin_mode,vmax=vmax_mode,interpolation='nearest')
        myDMD.getFrqStr(j,2)
        ax.set_title(str(j)+' '+(myDMD.getFrqStr(j,2) + r' |$\lambda$|=' + str(np.round(myDMD.getEigAbs(idx=j),6))+r' |$\phi$| '+ str(np.round(myDMD.getNorm(j),2))))
        #ax.set_title(myDMD.getFrqStr(j,2))
    plt.tight_layout()
