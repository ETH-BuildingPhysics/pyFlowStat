'''
PointProbe.py

Collection of tools/functions to load and/or get spectral analysis of time series extracted from
a turbulent flow field.

functions included:
    readOfRuntime(probeLoc,ofFile,rtnType='numpy')
    genPtStat(probeName,probeLoc,tserie,Userie)
'''


#===========================================================================#
# load modules
#===========================================================================#
#standard modules
import sys
import re

#scientific modules
import numpy as np
import scipy as sp
from scipy import signal
from scipy.optimize import curve_fit
# special modules

#from pyFlowStat.TurbulenceTools import TurbulenceTools as tt
import pyFlowStat.TurbulenceTools as tt

class PointProbe(object):
    """PointProbe Class"""
    
    def __init__(self,probeLoc=[],ofFile=''):
        if len(probeLoc)>0 and ofFile != '':
            self.probeLoc=probeLoc
            self.probeTimes,self.probeVar = self.readFromOpenFoam(probeLoc,ofFile)
            self.createDataDict()
        else:
            self.probeLoc=[]
            self.probeTimes=[]
            self.probeVar=[]
        return

    #===========================================================================#
    # functions
    #===========================================================================#
    def t(self):
        return self.data['t']
    def Ux(self):
        return self.data['U'][:,0]
    def Uy(self):
        return self.data['U'][:,1]
    def Uz(self):
        return self.data['U'][:,2]
    def Umag(self):
        return self.data['Umag']
    def ux(self):
        return self.data['u'][:,0]
    def uy(self):
        return self.data['u'][:,1]
    def uz(self):
        return self.data['u'][:,2]
    def Umean(self):
        return np.mean(self.data['Umag'])
     
    def uu_bar(self):
        return np.mean(pow(signal.detrend(self.ux()),2))
    def vv_bar(self):
        return np.mean(pow(signal.detrend(self.uy()),2))
    def ww_bar(self):
        return np.mean(pow(signal.detrend(self.uz()),2))
    def uv_bar(self):
        return np.mean(signal.detrend(self.ux())*signal.detrend(self.uy()))
    def uw_bar(self):
        return np.mean(signal.detrend(self.ux())*signal.detrend(self.uz()))
    def vw_bar(self):
        return np.mean(signal.detrend(self.uy())*signal.detrend(self.uz()))
    def TKE_bar(self):
        return 0.5*(self.uu_bar()+self.vv_bar()+self.ww_bar())
        
    def readFromOpenFoam(self,probeLoc,filepath):
        '''
        Read runtime probe generate by OpenFOAM. return time t and variable var.
        
        Arguments:
            probeLoc: [numpy.array or list with shape=(3)] Coordinate of probe (must be included in ofFile)
            ofFile:   [path] Path to OpenFOAM probe file
            rtnType:  ['numpy', 'array'. default='numpy'] Type of returned array
        
        Returns:
            t:   [numpy.array or list with shape=(N)] The time with N, the lenght of the serie.
            var: [numpy.array or list with shape=(a,N)] The variable with a the size of the variable (scalar, vector or tensor)
        '''
        probeTimes = []
        probeVar = []
        # read file 
        crs = open(filepath, 'r')
        lineno = 0
        for line in crs:
            # This regex finds all numbers in a given string.
            # It can find floats and integers writen in normal mode (10000) or with power of 10 (10e3).
            match = re.findall('[-+]?\d*\.?\d+e*[-+]?\d*', line)
            #print(match)
            if lineno==0:
                allXs = match
            if lineno==1:
                allYs = match
            if lineno==2:
                allZs = match
                ptFound = False
                for i in range(len(allXs)):
                    if (float(allXs[i])==probeLoc[0] and float(allYs[i])==probeLoc[1] and float(allZs[i])==probeLoc[2]):
                        ptPos = i
                        ptFound = True
                if ptFound==True:
                    #print('Probe found!')
                    pass
                else:
                    print('Probe not found!')
                    break          
            if lineno>3 and len(match)>0:
                if lineno==4:
                    varSize = int((len(match)-1)/(len(allXs)))  #check if probe var is scalar, vec or tensor
                    srtindex = 1+ptPos*varSize
                    endindex = srtindex+varSize
                probeTimes.append(float(match[0]))
                #probeTimes.append(float(match[0]))
                probeVar.append([float(var) for var in match[srtindex:endindex]])
            else:
                pass
            lineno = lineno+1
        crs.close()
        
        return (np.array(probeTimes),np.array(probeVar))
    def readFromLDA(self,probeLoc,filepath):
        probeVar = []
        probeTimes = []
        
        crs = open(filepath, 'r')
        lineno = 0
        for line in crs:
            if lineno>=5:
                data=line.split()
                probeTimes.append(float(data[0])/1000.0)
                probeVar.append([float(data[2]),0,0])
            lineno = lineno+1
        crs.close()
        self.probeLoc=probeLoc
        self.probeVar=np.array(probeVar)
        self.probeTimes=np.array(probeTimes)
        self.createDataDict()
        
    def cutData(self,indices):
        self.probeVar=self.probeVar[np.array(indices),:]
        self.probeTimes=self.probeTimes[np.array(indices),:]
        self.createDataDict()
        
    def createDataDict(self):
        '''
        Creates the "data" dcit from the following class variables:
            probeLoc:  [numpy.array with shape=(3)] Coordinate (x,y,z) of probe
            probeTimes:    [numpy.array with shape=(N)] List of time t.
            probeVar:    [numpy.array with shape=(N,3)] List of velocity vector U.
            
        Populates the "data" python dict with with the following keys:
            A python dict with the following keys:

                pos:  [numpy.array of shape=(3)] Probe location
                frq:  [float] Sample frequence
                
                U:    [numpy.array of shape=(N,3)] Velocity U
                t:    [numpy.array of shape=(N)] Time t
                u:    [numpy.array of shape=(N,3)] Velocity fluctuation u
                Umag: [numpy.array of shape=(N)] Velocity magnitute Umag
                umag: [numpy.array of shape=(N)] Fluctuating velocity magnitute umag   
                Uoo:  [numpy.array of shape=(N,3)] Mean velocity with infinit window size
        '''
        
        self.data = dict()
        self.data['pos'] = self.probeLoc
        # velocity and time
        self.data['U'] = self.probeVar
        self.data['t'] = self.probeTimes
        
        # timestep, sample frequence
        self.data['dt']=self.data['t'][1]-self.data['t'][0]
        self.data['frq'] = 1/self.data['dt']
        
        #Umag
        Umag = np.zeros(self.data['U'].shape[0])
        for i in range(self.data['U'].shape[0]):
            a=self.data['U'][i,:]
            nrm2, = sp.linalg.get_blas_funcs(('nrm2',), (a,))
            Umag[i] = nrm2(a)
            #Umag[i] = np.linalg.norm(self.data['U'][i,:])
        self.data['Umag']=Umag
        #mean
        Uoo = np.zeros((self.data['U'].shape))
        Uoo[:,0] = np.mean(self.data['U'][:,0])
        Uoo[:,1] = np.mean(self.data['U'][:,1])
        Uoo[:,2] = np.mean(self.data['U'][:,2])
        self.data['Uoo'] = Uoo
        # fluctuation
        self.data['u'] = self.data['U']-self.data['Uoo']
        #umag
        umag = np.zeros(self.data['u'].shape[0])
        for i in range(self.data['u'].shape[0]):
            a=self.data['u'][i,:]
            nrm2, = sp.linalg.get_blas_funcs(('nrm2',), (a,))
            umag[i] = nrm2(a)
            #umag[i] = np.linalg.norm(self.data['u'][i,:])

    def generateStatistics(self,doDetrend=True):
        '''
        Generate statistic from velocity vector U and time t for a given point P.
        The lenght of series U and t is N. This function is a convinent to generate some useful
        analysis in one line. This function may grow in the future.
        
        Arguments:
            doDetrend: detrend data bevor sigbal processing
            
        Populates the "data" python dict with with the following keys:

                rii:    [numpy.array of shape=(?)] Auto-correlation coefficent rii. For i=1,2,3 
                taurii: [numpy.array of shape=(?)] Time lags for rii. For i=1,2,3
                Rii:    [numpy.array of shape=(?)] Auto-correlation Rii. For i=1,2,3 
                tauRii: [numpy.array of shape=(?)] Time lags for Rii. For i=1,2,3
                
                uifrq:  [numpy.array of shape=(?)] u1 in frequency domain. For i=1,2,3
                uiamp:  [numpy.array of shape=(?)] amplitude of u1 in frequency domain. For i=1,2,3
                Seiifrq:[numpy.array of shape=(?)] Frequencies for energy spectrum Seii. For i=1,2,3
                Seii:   [numpy.array of shape=(?)] Energy spectrum Seii derived from Rii. For i=1,2,3
        '''   
        # auto correlation corefficient of u
        if doDetrend:
            ux=signal.detrend(self.ux());
            uy=signal.detrend(self.uy());
            uz=signal.detrend(self.uz());
        else:
            ux=self.ux();
            uy=self.uy();
            uz=self.uz();
        #ux=ux[-samples:-1]
        #uy=uy[-samples:-1]
        #uz=uz[-samples:-1]
        self.data['r11'],self.data['taur11'] = tt.xcorr_fft(ux, maxlags=None, norm='coeff')
        self.data['r22'],self.data['taur22'] = tt.xcorr_fft(uy, maxlags=None, norm='coeff')
        self.data['r33'],self.data['taur33'] = tt.xcorr_fft(uz, maxlags=None, norm='coeff')
        # auto correlation of u
        self.data['R11'],self.data['tauR11'] = tt.xcorr_fft(ux, maxlags=None, norm='none')
        self.data['R22'],self.data['tauR22'] = tt.xcorr_fft(uy, maxlags=None, norm='none')
        self.data['R33'],self.data['tauR33'] = tt.xcorr_fft(uz, maxlags=None, norm='none')
        
        
        #u in frequency domain
        self.data['u1frq'],self.data['u1amp'] = tt.dofft(sig=ux,samplefrq=self.data['frq'])
        self.data['u2frq'],self.data['u2amp'] = tt.dofft(sig=uy,samplefrq=self.data['frq'])
        self.data['u3frq'],self.data['u3amp'] = tt.dofft(sig=uz,samplefrq=self.data['frq'])
        #Time energy sectrum Se11 (mean: Rii in frequency domain...)
        self.data['Se11frq'],self.data['Se11'] = tt.dofft(sig=self.data['R11'],samplefrq=self.data['frq'])
        self.data['Se22frq'],self.data['Se22'] = tt.dofft(sig=self.data['R22'],samplefrq=self.data['frq'])
        self.data['Se33frq'],self.data['Se33'] = tt.dofft(sig=self.data['R33'],samplefrq=self.data['frq'])
        
    def lengthScale(self):
        def func_exp(x, a):
            np.seterr('ignore')
            res = np.exp(-x/a)
            #print res
            return res
            
        corr_keys=['taur11','taur22','taur33','r11','r22','r33']
        if len(set(corr_keys) & set(self.data.keys()))!=len(corr_keys):
            self.generateStatistics()
            print "Generating Missing Statistics"
            
        xdata=self.data['taur11']
        ydata=self.data['r11']
        popt, pcov = curve_fit(func_exp,xdata,ydata)
        self.data['Txx']=popt[0]*self.data['dt']
        self.data['Lxx']=self.data['Txx']*self.Umean()
        
        xdata=self.data['taur22']
        ydata=self.data['r22']
        popt, pcov = curve_fit(func_exp,xdata,ydata)
        self.data['Tyy']=popt[0]*self.data['dt']
        self.data['Lyy']=self.data['Tyy']*self.Umean()
        
        xdata=self.data['taur33']
        ydata=self.data['r33']
        popt, pcov = curve_fit(func_exp,xdata,ydata)
        self.data['Tzz']=popt[0]*self.data['dt']
        self.data['Lzz']=self.data['Tzz']*self.Umean()
        
    def detrend_periodic(self):
        def func_sin_u(Umean):
            f=Umean/2.0;
            def func_sin(x, a, b):
                return a*np.sin(2*np.pi*f*x+b)
            return func_sin
    
        def plot_sin(x, a, b,Umean):
            f=Umean/2.0;
            return a*np.sin(2*np.pi*f*x+b)
            
        xdata=self.t()
        Umean=self.Umean()
        self.data['u'][:,0]=signal.detrend(self.data['u'][:,0])
        ydata=self.ux()
        popt, pcov = curve_fit(func_sin_u(Umean),xdata,ydata)
        self.data['ux_old']=self.data['u'][:,0]
        self.data['ux_sin']=plot_sin(xdata,popt[0],popt[1],Umean)
        self.data['u'][:,0]=self.data['u'][:,0]-plot_sin(xdata,popt[0],popt[1],Umean)
        
    def detrend_sin(self):

        def func_sin(x, a, b, f):
                    return a*np.sin(2*np.pi*f*x+b)
            
        xdata=self.t()
        Umean=self.Umean()
        self.data['u'][:,0]=signal.detrend(self.data['u'][:,0])
        ydata=self.ux()
        popt, pcov = curve_fit(func_sin,xdata,ydata)
        self.data['ux_old']=np.copy(self.data['u'][:,0])
        self.data['ux_sin']=func_sin(xdata,popt[0],popt[1],popt[2])
        self.data['u'][:,0]=self.data['u'][:,0]-func_sin(xdata,popt[0],popt[1],popt[2])
        
def getDataPoints(filepath):
    f = open(filepath, 'r')
    xlist=[]
    ylist=[]
    zlist=[]
    lineno=0
    for line in f:
        # This regex finds all numbers in a given string.
        # It can find floats and integers writen in normal mode (10000) or with power of 10 (10e3).
        match = re.findall('[-+]?\d*\.?\d+e*[-+]?\d*', line)
        #print(match)
        if lineno==0:
            xlist = match
        if lineno==1:
            ylist = match
        if lineno==2:
           zlist = match
           break
        lineno+=1
    f.close()
    pointlist = np.zeros(shape=(len(xlist),3))
    pointlist[:,0]=xlist
    pointlist[:,1]=ylist
    pointlist[:,2]=zlist

    return pointlist        

def getVectorPointProbeList(filename):
    pointlist=getDataPoints(filename)
    pts=[PointProbe()]*len(pointlist)
    
    for i in range(0,len(pointlist)):
        pts[i]=PointProbe()
        pts[i].probeLoc=pointlist[i]
        
    # read file 
    crs = open(filename, 'r')
    lineno = 0
    for line in crs:
        # This regex finds all numbers in a given string.
        # It can find floats and integers writen in normal mode (10000) or with power of 10 (10e3).
        match = np.array(re.findall('[-+]?\d*\.?\d+e*[-+]?\d*', line))
        match=match.astype(np.float)
        #print(match)

        if lineno>3 and len(match)>0:
            for i in range(0,len(pointlist)):
                pts[i].probeTimes.append(match[0])
                pts[i].probeVar.append(match[1+i*3:4+i*3])
        else:
            pass
        lineno = lineno+1
    crs.close()
    
    for i in range(0,len(pointlist)):
        pts[i].probeTimes = np.array(pts[i].probeTimes)
        pts[i].probeVar = np.array(pts[i].probeVar)
        pts[i].createDataDict()
        
    return pts
        