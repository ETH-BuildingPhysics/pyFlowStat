# -*- coding: utf-8 -*-
"""
Created on Mon Dec 01 08:56:19 2014

@author: mimmer
"""
from . import Decomposition
from . import DecompositionPlotting
from IPython.html.widgets import interact
from IPython.html import widgets
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import numpy as np
import scipy.signal as sp

def modeAnalyzer(myDMD_Uy,k=(0,100),m=(0,500),vmin=-0.1,vmax=0.1,xscale='log'):
    def plotTime(k,m,vmin=vmin,vmax=vmax,plotSpectrum=True,magSort=True):
        idx_pos,idx_neg,f=myDMD_Uy.getFrqSortedIdx()
        normList=myDMD_Uy.result['mode_norms']
        stabList=np.absolute(myDMD_Uy.result['ritz_vals'])
    
        fig=plt.figure(figsize=(15,5))
        gs = gridspec.GridSpec(1, 5)
        ax1 = plt.subplot(gs[0,0])
        ax2 = plt.subplot(gs[0,1])
        ax3 = plt.subplot(gs[0,2])
        ax4 = plt.subplot(gs[0,3:5])
        
        #idx=[myDMD_Uy.getIdxforFrq(65.0)]
        if magSort:
            #idx=[myDMD_Uy.getGrowthSortedIdx()[::-1][m]]
            idx_stab=np.argsort(stabList**k)[::-1]
            idx_norm=np.argsort(normList)[::-1]
            

            #idx=[idx_stab[m]]
            idx=[idx_norm[m]]
        else:
            #idx_pos,idx_neg,fc=myDMD_Uy.getFrqSortedIdx()
            idx=[idx_pos[m]]
        
        #idx=[idx_pos[m]]
            
        
        out=myDMD_Uy.reconstructDMD(k=k,idx=idx_stab[1:])

        ax1.imshow(np.real(out),vmin=vmin,vmax=vmax)
        ax1.set_title(str(idx[0]))
        md=myDMD_Uy.getMode(idx[0])

        ax2.imshow(np.real(md),vmin=vmin,vmax=vmax)

        ax3.imshow(np.imag(md),vmin=vmin,vmax=vmax)

        if(plotSpectrum):
            f1=np.abs(myDMD_Uy.getFrq(idx[0]))
            #ax4.plot(f,normList)
            flist=myDMD_Uy.getFrqList()
            DecompositionPlotting.plotEigenvals(myDMD_Uy,ax4,k)
            ritz_sel=myDMD_Uy.result['ritz_vals'][idx[0]]**k
            ax4.plot(np.real(ritz_sel),np.imag(ritz_sel),'rx')
            #ax4.plot(flist[idx_pos],10*np.log10(stabList[idx_pos]))
            #ax4.plot(flist[idx_pos[m]],10*np.log10(stabList[idx_pos[m]]),'ro')
            
            #ax4.plot(f1,myDMD_Uy.result['mode_norms'][idx[0]],'ro')
            #print f1,myDMD_Uy.result['mode_norms'][idx[0]]
            #ax4.set_xscale(xscale)
            #ax4.set_yscale('log')
            ax4.set_title(str(f1)+ r' Hz |$\lambda$|=' + str(myDMD_Uy.getEigAbs()[idx[0]]))
        plt.show()
    kw=widgets.IntSliderWidget()
    kw.min=k[0]
    kw.max=k[1]
    kw.value=1
    
    mw=widgets.IntSliderWidget()
    mw.min=m[0]
    mw.max=m[1]
    mw.value=0
    
    ws=widgets.interaction.interactive(plotTime,k=kw,m=mw,vmin=widgets.FloatTextWidget(value=vmin),vmax=widgets.FloatTextWidget(value=vmax),plotSpectrum=True,magSort=True)
       
    return ws

def modeAnalyzerGrowth(myDMD_Uy,k=(0,100),m=(0,500),vmin=-0.1,vmax=0.1,xscale='log',min_norm=3.0,min_abs=0.0003,component=0,extent=None):
    def plotTime(k,m,sorting,vmin=vmin,vmax=vmax,plotSpectrum=True):
        idx_pos,idx_neg,f=myDMD_Uy.getFrqSortedIdx()
        
        normList=myDMD_Uy.result['mode_norms']
        stabList=np.absolute(myDMD_Uy.result['ritz_vals'])
    
        fig=plt.figure(figsize=(15,5))
        gs = gridspec.GridSpec(1, 5)
        ax1 = plt.subplot(gs[0,0])
        ax2 = plt.subplot(gs[0,1])
        ax3 = plt.subplot(gs[0,2])
        ax4 = plt.subplot(gs[0,3:5])
        
        #idx=[myDMD_Uy.getIdxforFrq(65.0)]
        
        if sorting==0:
            idx_sel=myDMD_Uy.getFilteredIndex(min_norm=min_norm,min_abs=min_abs)[::-1]
            idx_sel=[idx for idx in idx_pos if idx in idx_sel]
            idx=[idx_sel[m]]
            
        elif sorting==1:
            #idx=[myDMD_Uy.getGrowthSortedIdx()[::-1][m]]
            idx_norm=np.argsort(normList)[::-1]
            idx_sel=idx_norm

            #idx=[idx_stab[m]]
            idx=[idx_norm[m]]
        elif sorting==2:
            idx_stab=np.argsort(np.abs(np.log(stabList)))
            idx_sel=idx_stab
            idx=[idx_stab[m]]
        elif sorting==3:
            idx_none=np.linspace(0,len(normList),len(normList)+1)
            idx_sel=idx_none
            print idx_sel
            idx=[idx_none[m]]
        elif sorting==4:
            idx_sel=myDMD_Uy.getFilteredIndex(min_norm=min_norm,min_abs=min_abs)[::-1]
            idx=[idx_sel[m]]
        
        #idx=[idx_pos[m]]
        out=myDMD_Uy.reconstructDMD(k=k,idx=idx)
        if len(out.shape)>2:
            out=out[component]
        ax1.imshow(np.real(out),vmin=vmin,vmax=vmax,extent=extent)
        ax1.set_title(str(idx[0]))
        #md=myDMD_Uy.getMode(idx[0],component=component)
        md=myDMD_Uy.getMode(idx[0],component=component)
        print md.shape
        if len(md.shape)>2:
            md=md[component]
            
        ax2.imshow(np.real(md)/np.linalg.norm(md),vmin=vmin,vmax=vmax,extent=extent)
        ax2.grid()
        ax3.imshow(np.imag(md)/np.linalg.norm(md),vmin=vmin,vmax=vmax,extent=extent)
        ax3.grid()

        if(plotSpectrum):
            f1=np.abs(myDMD_Uy.getFrq(idx[0]))
            #ax4.plot(f,normList)
            flist=myDMD_Uy.getFrqList()
            if sorting==4 or sorting==0:
                DecompositionPlotting.plotGrowth(myDMD_Uy,ax4,k,idx=idx_sel[::-1],dt=myDMD_Uy.result['dt'])
            else:
                #DecompositionPlotting.plotGrowth(myDMD_Uy,ax4,k,idx=idx_sel[::-1])
                DecompositionPlotting.plotGrowth(myDMD_Uy,ax4,k,dt=myDMD_Uy.result['dt'])
            ritz_sel=np.log(myDMD_Uy.result['ritz_vals'][idx[0]]**k)/myDMD_Uy.result['dt']
            ax4.plot(np.imag(ritz_sel)/(2.0*np.pi),np.real(ritz_sel),'rx')
            #ax4.plot(flist[idx_pos],10*np.log10(stabList[idx_pos]))
            #ax4.plot(flist[idx_pos[m]],10*np.log10(stabList[idx_pos[m]]),'ro')
            
            #ax4.plot(f1,myDMD_Uy.result['mode_norms'][idx[0]],'ro')
            #print f1,myDMD_Uy.result['mode_norms'][idx[0]]
            #ax4.set_xscale(xscale)
            #ax4.set_yscale('log')
            ax4.grid()
            ax4.set_title(str(np.round(f1,5))+ r' Hz |$\lambda$|=' + str(np.round(myDMD_Uy.getEigAbs()[idx[0]],6))+r' |$\phi$| '+ str(np.round(normList[idx[0]],2)))
        plt.show()
    kw=widgets.IntSliderWidget()
    kw.min=k[0]
    kw.max=k[1]
    kw.value=1
    
    mw=widgets.IntSliderWidget()
    mw.min=m[0]
    mw.max=m[1]
    mw.value=0
    sw=widgets.Dropdown(options={'Frequency':0,'Norm':1,'Growth':2,'None':3,'Importance':4})
    #sw.values={'Frequency':0,'Norm':1,'Growth':2,'None':3,'Importance':4}
    ws=widgets.interaction.interactive(plotTime,k=kw,m=mw,vmin=widgets.FloatTextWidget(value=vmin),vmax=widgets.FloatTextWidget(value=vmax),plotSpectrum=True,sorting=sw)
       
    return ws
    
def modeAnalyzer2(myDMD_Uy,k=(0,100),m=(0,500),vmin=-0.1,vmax=0.1,xscale='linear'):
    def plotTime(k,m,plotSpectrum=True,magSort=True):
        #idx_pos,idx_neg,f=myDMD_Uy.getFrqSortedIdx()
        normList=myDMD_Uy.result['mode_norms']
        growthList=np.absolute(myDMD_Uy.result['ritz_vals'])
        f=myDMD_Uy.getFrqList()
    
        fig=plt.figure(figsize=(15,5))
        gs = gridspec.GridSpec(1, 5)
        ax1 = plt.subplot(gs[0,0])
        ax2 = plt.subplot(gs[0,1])
        ax3 = plt.subplot(gs[0,2])
        ax4 = plt.subplot(gs[0,3:5])

        #idx=[myDMD_Uy.getIdxforFrq(65.0)]
        if magSort:
            idx=[myDMD_Uy.getGrowthSortedIdx()[::-1][m]]
            #idx=[myDMD_Uy.getNormSortedIdx()[::-1][m]]
        else:
            #idx_pos,idx_neg,fc=myDMD_Uy.getFrqSortedIdx()
            #idx=[idx_pos[m]]
            idx=[range(len(myDMD_Uy.result['ritz_vals']))[m]]
            
        
        out=myDMD_Uy.reconstructDMD(k=k,idx=idx)

        ax1.imshow(np.real(out),vmin=vmin,vmax=vmax)
        ax1.set_title(str(idx[0]))
        md=myDMD_Uy.getMode(idx[0])

        ax2.imshow(np.real(md),vmin=vmin,vmax=vmax)

        ax3.imshow(np.imag(md),vmin=vmin,vmax=vmax)

        print idx,myDMD_Uy.getFrq(idx[0]),myDMD_Uy.result['ritz_vals'][idx[0]],myDMD_Uy.result['mode_norms'][idx[0]]
        if(plotSpectrum):
            f1=np.abs(myDMD_Uy.getFrq(idx[0]))
            ax4.plot(f1,np.absolute(myDMD_Uy.result['ritz_vals'][idx[0]]),'x')
            ax4.scatter(f,growthList,marker='d',c=normList,edgecolor='none')
            
            ax4.set_xscale(xscale)
            ax4.set_yscale('log')
            ax4.set_title(str(f1)+ ' Hz')
        plt.show()
    kw=widgets.IntSliderWidget()
    kw.min=k[0]
    kw.max=k[1]
    kw.value=0
    
    mw=widgets.IntSliderWidget()
    mw.min=m[0]
    mw.max=m[1]
    mw.value=0
    
    ws=widgets.interaction.interactive(plotTime,k=kw,m=mw,plotSpectrum=True,magSort=True)
    return ws
    
def wavelengthAnalyzer(myDMD_Uy,dx,L,iRow,vmin=-0.1,vmax=0.1,maxOnly=False):    
    def on_change2(pt):
        #maxi=sp.argrelmax(normList)[0]
        #pt=maxi[pt]
        
        fig=plt.figure(figsize=(20,10))
        
        gs = gridspec.GridSpec(2, 2)
        ax1 = plt.subplot(gs[:, 0])
        ax2 = plt.subplot(gs[0,1])
        ax3 = plt.subplot(gs[1,1])
        
        #ax=plt.subplot(1,2,1)
        ax1.plot(f,normList)
        
        
        
        ax1.plot(f[pt],normList[pt],'ko')
        #ax1.text(f[pt],normList[pt],str(f[pt])+ 'Hz')
        string='f={:.3f} Hz\nMode={:.0f}'.format(f[pt],pt)
        ax1.text(0.05, 0.95, string, transform=ax1.transAxes, fontsize=14,
            verticalalignment='top')
        
        ax1.set_xscale('log')
        ax1.set_yscale('log')
        
        #ax=plt.subplot(1,2,2)
        idxMode=myDMD_Uy.getIdxforFrq(f[pt])
        mode=myDMD_Uy.getMode(idxMode)
        ax2.imshow(np.real(mode),vmin=vmin,vmax=vmax,interpolation='nearest')
        
        uy=np.array(np.real(mode)[iRow,:])
        uy_imag=np.array(np.imag(mode)[iRow,:])
        ax3.plot(uy)
        ax3.plot(uy_imag,'r')
        maxi=sp.argrelmax(uy)[0]
        mini=sp.argrelmin(uy)[0]
        exti=np.sort(np.r_[maxi,mini])
        
        maxi_imag=sp.argrelmax(uy_imag)[0]
        mini_imag=sp.argrelmin(uy_imag)[0]
        exti_imag=np.sort(np.r_[maxi_imag,mini_imag])        
        
        print np.diff(exti)
        ax3.scatter(maxi,uy[maxi],marker=2)
        ax3.scatter(mini,uy[mini],marker=3)
        ax3.scatter(maxi_imag,uy_imag[maxi_imag],marker=2)
        ax3.scatter(mini_imag,uy_imag[mini_imag],marker=3)       

        ax3.set_xlim([0,np.real(mode).shape[1]])
        gamma=0
        print 'n=',L/(np.diff(maxi)*dx)+gamma
        print 'n=',L/(np.diff(mini)*dx)+gamma
        print 'n=',L/(np.diff(exti)*dx*2.0)+gamma
        
        print 'n=',L/(np.diff(maxi_imag)*dx)+gamma
        print 'n=',L/(np.diff(mini_imag)*dx)+gamma
        print 'n=',L/(np.diff(exti_imag)*dx*2.0)+gamma
    
    idx_pos,idx_neg,f=myDMD_Uy.getFrqSortedIdx()
    normList=myDMD_Uy.result['mode_norms'][idx_pos]
    ws=widgets.IntSliderWidget()
    ws.min=0
    if maxOnly==True:
        ws.max=len(sp.argrelmax(normList)[0])
    else:
        ws.max=len(normList)
    ws.value=1
    
    ws=widgets.interaction.interactive(on_change2,pt=ws)
    return ws