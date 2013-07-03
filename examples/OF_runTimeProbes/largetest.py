#! /usr/bin/env python

#standard modules
import sys

# scientific modules
import numpy as np
import scipy.fftpack as spfft
import scipy.stats as spstats
import scipy.signal as spsig

# plot modules
import matplotlib as mpl
#mpl.use('QT4Agg')         #load qt4 backend (optional)
import matplotlib.pyplot as plt
import pylab as pl

# special import to save as pdf (backend modification)
#from matplotlib.backends.backend_pdf import PdfPages

#special modules
from pyFlowStat.PointProbe import PointProbe as pp
#from pyFlowStat.TurbulenceTools import TurbulenceTools as tt
import pyFlowStat.TurbulenceTools as tt

# parameters for plots
params = {'font.family': 'serif',
          'axes.labelsize': 14,
          'text.fontsize': 14,
          'legend.fontsize': 14,
          'xtick.labelsize': 14,
          'ytick.labelsize': 14,}
pl.rcParams.update(params)

#file to load and probe location
probeFile1 = '10/U'
probeFile2 = '19.82/U'
probeLoc0 = [-0.1,-0.3,0.1]
probeLoc1 = [-0.1,-0.3,0.3]
probeLoc2 = [1.1,-0.3,0.1]
probeLoc3 = [1.1,-0.3,0.3]
probeLoc4 = [1.9,-0.3,0.1]
probeLoc5 = [1.9,-0.3,0.3]


# local function to load a point from various files 
def loadAndGenStat(probeName,probeLoc,ofFiles):
    tseries = []
    varseries = []
    for ofFile in ofFiles:
        t,var = pp(probeLoc=probeLoc,ofFile=ofFile)
        tseries.append(t)
        varseries.append(var)
    pt_t = tseries[0]
    pt_var = varseries[0]
    if len(tseries)>1:
        for i in range(len(tseries)-1):
            pt_t = np.hstack((pt_t,tseries[i+1]))
            pt_var = np.vstack((pt_var,varseries[i+1]))
    return pp.genPtStat(probeName=probeName,probeLoc=probeLoc,tserie=pt_t,Userie=pt_var)
            
pt0 = loadAndGenStat('pt0',probeLoc0,ofFiles=[probeFile1,probeFile2])
pt1 = loadAndGenStat('pt1',probeLoc1,ofFiles=[probeFile1,probeFile2])
pt2 = loadAndGenStat('pt2',probeLoc2,ofFiles=[probeFile1,probeFile2])
pt3 = loadAndGenStat('pt3',probeLoc3,ofFiles=[probeFile1,probeFile2])
pt4 = loadAndGenStat('pt4',probeLoc4,ofFiles=[probeFile1,probeFile2])
pt5 = loadAndGenStat('pt5',probeLoc5,ofFiles=[probeFile1,probeFile2])

# <markdowncell>

# <font size="5">Plots</font>

# <markdowncell>

# Generate plots

# <codecell>

#plot pt2: U and Uoo
fig0 = plt.figure(figsize=[10,8])
plt.subplot(311)
plt.subplots_adjust(hspace=0.5)
plt.grid(True)
plt.title(r'pt2: Velocity $U_{1}$ and mean velocity $\langle U_{1} \rangle$')
plt.xlabel('t [s]')
plt.ylabel(r'$U_{1}$ , $\langle U_{1} \rangle$ [m/s]')
plt.plot(pt2['t'],pt2['U'][:,0],label=r'$U_{1}$',color='blue')
plt.plot(pt2['t'],pt2['Uoo'][:,0],label=r'$\langle U_{1} \rangle$',color='darkblue',lw=2)
plt.legend()

plt.subplot(312)
plt.grid(True)
plt.title(r'pt2: Velocity $U_{2}$ and mean velocity $\langle U_{2} \rangle$')
plt.xlabel('t [s]')
plt.ylabel(r'$U_{2}$ , $\langle U_{2} \rangle$ [m/s]')
plt.plot(pt2['t'],pt2['U'][:,1],label=r'$U_{2}$',color='green')
plt.plot(pt2['t'],pt2['Uoo'][:,1],label=r'$\langle U_{2} \rangle$',color='darkgreen',lw=2)
plt.legend()

plt.subplot(313)
plt.grid(True)
plt.title(r'pt2: Velocity $U_{3}$ and mean velocity $\langle U_{3} \rangle$')
plt.xlabel('t [s]')
plt.ylabel(r'$U_{3}$ , $\langle U_{3} \rangle$ [m/s]')
plt.plot(pt2['t'],pt2['U'][:,2],label=r'$U_{3}$',color='red')
plt.plot(pt2['t'],pt2['Uoo'][:,2],label=r'$\langle U_{3} \rangle$',color='darkred',lw=2)
plt.legend()

# <codecell>

#plot pt4: U and Uoo
fig1 = plt.figure(figsize=[10,8])
plt.subplot(311)
plt.subplots_adjust(hspace=0.5)
plt.grid(True)
plt.title(r'pt4: Velocity $U_{1}$ and mean velocity $\langle U_{1} \rangle$')
plt.xlabel('t [s]')
plt.ylabel(r'$U_{1}$ , $\langle U_{1} \rangle$ [m/s]')
plt.plot(pt4['t'],pt4['U'][:,0],label=r'$U_{1}$',color='blue')
plt.plot(pt4['t'],pt4['Uoo'][:,0],label=r'$\langle U_{1} \rangle$',color='darkblue',lw=2)
plt.legend()

plt.subplot(312)
plt.grid(True)
plt.title(r'pt4: Velocity $U_{2}$ and mean velocity $\langle U_{2} \rangle$')
plt.xlabel('t [s]')
plt.ylabel(r'$U_{2}$ , $\langle U_{2} \rangle$ [m/s]')
plt.plot(pt4['t'],pt4['U'][:,1],label=r'$U_{2}$',color='green')
plt.plot(pt4['t'],pt4['Uoo'][:,1],label=r'$\langle U_{2} \rangle$',color='darkgreen',lw=2)
plt.legend()

plt.subplot(313)
plt.grid(True)
plt.title(r'pt4: Velocity $U_{3}$ and mean velocity $\langle U_{3} \rangle$')
plt.xlabel('t [s]')
plt.ylabel(r'$U_{1}$ , $\langle U_{3} \rangle$ [m/s]')
plt.plot(pt4['t'],pt4['U'][:,2],label=r'$U_{3}$',color='red')
plt.plot(pt4['t'],pt4['Uoo'][:,2],label=r'$\langle U_{3} \rangle$',color='darkred',lw=2)
plt.legend()

# <codecell>

#plot pt5: U and Uoo
fig2 = plt.figure(figsize=[10,8])
plt.subplot(311)
plt.subplots_adjust(hspace=0.5)
plt.grid(True)
plt.title(r'pt5: Velocity $U_{1}$ and mean velocity $\langle U_{1} \rangle$')
plt.xlabel('t [s]')
plt.ylabel(r'$U_{1}$ , $\langle U_{1} \rangle$ [m/s]')
plt.plot(pt5['t'],pt5['U'][:,0],label=r'$U_{1}$',color='blue')
plt.plot(pt5['t'],pt5['Uoo'][:,0],label=r'$\langle U_{1} \rangle$',color='darkblue',lw=2)
plt.legend()

plt.subplot(312)
plt.grid(True)
plt.title(r'pt5: Velocity $U_{2}$ and mean velocity $\langle U_{2} \rangle$')
plt.xlabel('t [s]')
plt.ylabel(r'$U_{2}$ , $\langle U_{2} \rangle$ [m/s]')
plt.plot(pt5['t'],pt5['U'][:,1],label=r'$U_{2}$',color='green')
plt.plot(pt5['t'],pt5['Uoo'][:,1],label=r'$\langle U_{2} \rangle$',color='darkgreen',lw=2)
plt.legend()

plt.subplot(313)
plt.grid(True)
plt.title(r'pt5: Velocity $U_{3}$ and mean velocity $\langle U_{3} \rangle$')
plt.xlabel('t [s]')
plt.ylabel(r'$U_{1}$ , $\langle U_{3} \rangle$ [m/s]')
plt.plot(pt5['t'],pt5['U'][:,2],label=r'$U_{3}$',color='red')
plt.plot(pt5['t'],pt5['Uoo'][:,2],label=r'$\langle U_{3} \rangle$',color='darkred',lw=2)
plt.legend()

# <codecell>

#plot: auto-correlation
fig3 = plt.figure(figsize=[16,8])
#pt0
plt.subplot(234)
plt.subplots_adjust(hspace=0.5)
plt.grid(True)
plt.title('pt0: auto-correlation coefficients')
plt.xlabel(r'$\tau$')
plt.ylabel(r'$r_{ij}(\tau)$')
plt.plot(pt0['taur11'],pt0['r11'],label=r'$r_{11}(\tau)$')
plt.plot(pt0['taur22'],pt0['r22'],label=r'$r_{22}(\tau)$')
plt.plot(pt0['taur33'],pt0['r33'],label=r'$r_{33}(\tau)$')
plt.xlim(0,500)
plt.legend()
#pt1
plt.subplot(231)
plt.grid(True)
plt.title('pt1: auto-correlation coefficients')
plt.xlabel(r'$\tau$')
plt.ylabel(r'$r_{ij}(\tau)$')
plt.plot(pt1['taur11'],pt1['r11'],label=r'$r_{11}(\tau)$')
plt.plot(pt1['taur22'],pt1['r22'],label=r'$r_{22}(\tau)$')
plt.plot(pt1['taur33'],pt1['r33'],label=r'$r_{33}(\tau)$')
plt.xlim(0,500)
plt.legend()
#pt2
plt.subplot(235)
plt.grid(True)
plt.title('pt2: auto-correlation coefficients')
plt.xlabel(r'$\tau$')
plt.ylabel(r'$r_{ij}(\tau)$')
plt.plot(pt2['taur11'],pt2['r11'],label=r'$r_{11}(\tau)$')
plt.plot(pt2['taur22'],pt2['r22'],label=r'$r_{22}(\tau)$')
plt.plot(pt2['taur33'],pt2['r33'],label=r'$r_{33}(\tau)$')
plt.xlim(0,500)
plt.legend()
#pt3
plt.subplot(232)
plt.grid(True)
plt.title('pt3: auto-correlation coefficients')
plt.xlabel(r'$\tau$')
plt.ylabel(r'$r_{ij}(\tau)$')
plt.plot(pt3['taur11'],pt3['r11'],label=r'$r_{11}(\tau)$')
plt.plot(pt3['taur22'],pt3['r22'],label=r'$r_{22}(\tau)$')
plt.plot(pt3['taur33'],pt3['r33'],label=r'$r_{33}(\tau)$')
plt.xlim(0,500)
plt.legend()
#pt4
plt.subplot(236)
plt.grid(True)
plt.title('pt4: auto-correlation coefficients')
plt.xlabel(r'$\tau$')
plt.ylabel(r'$r_{ij}(\tau)$')
plt.plot(pt4['taur11'],pt4['r11'],label=r'$r_{11}(\tau)$')
plt.plot(pt4['taur22'],pt4['r22'],label=r'$r_{22}(\tau)$')
plt.plot(pt4['taur33'],pt4['r33'],label=r'$r_{33}(\tau)$')
plt.xlim(0,500)
plt.legend()
#pt5
plt.subplot(233)
plt.grid(True)
plt.title('pt5: auto-correlation coefficients')
plt.xlabel(r'$\tau$')
plt.ylabel(r'$r_{ij}(\tau)$')
plt.plot(pt5['taur11'],pt5['r11'],label=r'$r_{11}(\tau)$')
plt.plot(pt5['taur22'],pt5['r22'],label=r'$r_{22}(\tau)$')
plt.plot(pt5['taur33'],pt5['r33'],label=r'$r_{33}(\tau)$')
plt.xlim(0,500)
plt.legend()

# <codecell>

#plot pt4: scatter plot uxuy
fig4 = plt.figure(figsize=[15,5])
gs = mpl.gridspec.GridSpec(1, 3, width_ratios=[1,1,1])
ptc = 'blue'
#uxuy
plt.subplot(gs[0])
plt.subplots_adjust(left=0.05,bottom=0.1,right=0.95,top=0.9,wspace=0.3,hspace=None)
plt.grid(True)
plt.xlabel(r'$u_{1}$')
plt.ylabel(r'$u_{2}$')
#plt.xlim(-0.2,0.2)
#plt.ylim(-0.2,0.2)
plt.scatter(pt4['U'][:,0],pt4['U'][:,1],label=r'$u_{1}u_{2}$',s=5,c=ptc,marker='o')
plt.axhline(0, color='black', lw=2)
plt.axvline(0, color='black', lw=2)
plt.legend()
#uxuz
plt.subplot(gs[1])
plt.grid(True)
plt.xlabel(r'$u_{1}$')
plt.ylabel(r'$u_{3}$')
#plt.xlim(-0.2,0.2)
#plt.ylim(-0.2,0.2)
plt.scatter(pt4['U'][:,0],pt4['U'][:,2],label=r'$u_{1}u_{3}$',s=5,c=ptc,marker='o')
plt.axhline(0, color='black', lw=2)
plt.axvline(0, color='black', lw=2)
plt.legend()
#uyuz
plt.subplot(gs[2])
plt.grid(True)
plt.xlabel(r'$u_{2}$')
plt.ylabel(r'$u_{3}$')
#plt.xlim(-0.2,0.2)
#plt.ylim(-0.2,0.2)
plt.scatter(pt4['U'][:,1],pt4['U'][:,2],label=r'$u_{2}u_{3}$',s=5,c=ptc,marker='o')
plt.axhline(0, color='black', lw=2)
plt.axvline(0, color='black', lw=2)
plt.legend()

# <codecell>

#plot pt4-5: u in Fourier space
fig5 = plt.figure(figsize=[12,5])
plt.subplot(121)
plt.subplots_adjust(wspace=0.3)
plt.title(r'pt4: $u_{i}$ in frequency domain')
plt.xlabel(r'$\omega$')
plt.ylabel(r'$u_{i}(\omega)$')
plt.grid(True)
plt.loglog()
plt.plot(pt4['u1frq'],pt4['u1amp'])
plt.plot(pt4['u2frq'],pt4['u2amp'])
plt.plot(pt4['u3frq'],pt4['u3amp'])
plt.xlim(1e-1,1e2)
plt.ylim(1e-6,1e-1)

plt.subplot(122)
plt.title(r'pt5: $u_{i}$ in frequency domain')
plt.xlabel(r'$\omega$')
plt.ylabel(r'$u_{i}(\omega)$')
plt.grid(True)
plt.loglog()
plt.plot(pt5['u1frq'],pt5['u1amp'])
plt.plot(pt5['u2frq'],pt5['u2amp'])
plt.plot(pt5['u3frq'],pt5['u3amp'])
plt.xlim(1e-1,1e2)
plt.ylim(1e-6,1e-1)

# <markdowncell>

# <font size="5">Energy Sectrum and Kolmogorov scales</font> ()

# <codecell>

help(tt.dofft)

# <codecell>

fig6 = plt.figure(figsize=[12,5])
# pt4
plt.subplot(121)
plt.subplots_adjust(wspace=0.3)
plt.grid(True)
plt.loglog()
plt.title(r'pt4: energy spectrum $Se_{ii}(\omega$)')
plt.plot(pt4['Se11frq'],pt4['Se11'],label=r'$Se_{11}(\omega)$')
plt.plot(pt4['Se22frq'],pt4['Se22'],label=r'$Se_{22}(\omega)$')
plt.plot(pt4['Se33frq'],pt4['Se33'],label=r'$Se_{33}(\omega)$')
plt.xlabel(r'$\omega$ $[Hz]$')
plt.ylabel(r'$Se_{ii}(\omega)$')
plt.xlim(1e-1,3e2)
plt.ylim(1e-7,1e1)
plt.legend()
# pt5
plt.subplot(122)
plt.grid(True)
plt.loglog()
plt.title(r'pt5: energy spectrum $Se_{ii}(\omega$)')
plt.plot(pt5['Se11frq'],pt5['Se11'],label=r'$Se_{11}(\omega)$')
plt.plot(pt5['Se22frq'],pt5['Se22'],label=r'$Se_{22}(\omega)$')
plt.plot(pt5['Se33frq'],pt5['Se33'],label=r'$Se_{33}(\omega)$')
plt.xlabel(r'$\omega$ $[Hz]$')
plt.ylabel(r'$Se_{ii}(\omega)$')
plt.xlim(1e-1,3e2)
plt.ylim(1e-7,1e1)
plt.legend()

# <markdowncell>

# Save plots

# <codecell>

fig0.savefig('pt2 - U-Uoo.pdf')
fig1.savefig('pt4 - U-Uoo.pdf')
fig2.savefig('pt5 - U-Uoo.pdf')
fig3.savefig('pt0-1-2-3-4-5 - auto-corelation.pdf')
fig4.savefig('pt4 - scatter.pdf')
fig5.savefig('pt4-5 - uii spectrum.pdf')
fig6.savefig('pt4-5 - energy spectrum Seii.pdf')

# <codecell>


