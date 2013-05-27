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

# special modules

from pyFlowStat.TurbulenceTools import TurbulenceTools as tt

class PointProbe(object):
    """PointProbe Class"""
    
    def __init__(self,probeLoc,ofFile):
        self.probeLoc=probeLoc
        self.probeTimes,self.probeVar = self.readFromOpenFoam(probeLoc,ofFile)
        self.createDataDict()
        
        return

    #===========================================================================#
    # functions
    #===========================================================================#

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

    def createDataDict(self):
        self.data = dict()
        self.data['pos'] = self.probeLoc
        # velocity and time
        self.data['U'] = self.probeVar
        self.data['t'] = self.probeTimes
        
        # sample frequence
        self.data['frq'] = 1/(self.data['t'][1]-self.data['t'][0])
        #Umag
        Umag = np.zeros(self.data['U'].shape[0])
        for i in range(self.data['U'].shape[0]):
            Umag[i] = np.linalg.norm(self.data['U'][i,:])
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
            umag[i] = np.linalg.norm(self.data['u'][i,:])

    def generateStatistics(self,doDetrend=True):
        '''
        Generate statistic from velocity vector U and time t for a given point P.
        The lenght of series U and t is N. This function is a convinent to generate some useful
        analysis in one line. This function may grow in the future.
        
        Arguments:
            probeName: [string] Name of probe
            probeLoc:  [numpy.array with shape=(3)] Coordinate (x,y,z) of probe
            tserie:    [numpy.array with shape=(N)] List of time t.
            Userie:    [numpy.array with shape=(N,3)] List of velocity vector U.
            
        Returns:
            A python dict with the following keys:
                name: [string] Probe name
                loc:  [numpy.array of shape=(3)] Probe location
                frq:  [float] Sample frequence
                
                U:    [numpy.array of shape=(N,3)] Velocity U
                t:    [numpy.array of shape=(N)] Time t
                u:    [numpy.array of shape=(N,3)] Velocity fluctuation u
                Umag: [numpy.array of shape=(N)] Velocity magnitute Umag
                umag: [numpy.array of shape=(N)] Fluctuating velocity magnitute umag   
                Uoo:  [numpy.array of shape=(N,3)] Mean velocity with infinit window size
                
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
        
        self.data['r11'],self.data['taur11'] = tt.xcorr_fft(self.data['u'][:,0], maxlags=None, norm='coeff',doDetrend=doDetrend)
        self.data['r22'],self.data['taur22'] = tt.xcorr_fft(self.data['u'][:,1], maxlags=None, norm='coeff',doDetrend=doDetrend)
        self.data['r33'],self.data['taur33'] = tt.xcorr_fft(self.data['u'][:,2], maxlags=None, norm='coeff',doDetrend=doDetrend)
        # auto correlation of u
        self.data['R11'],self.data['tauR11'] = tt.xcorr_fft(self.data['u'][:,0], maxlags=None, norm='none',doDetrend=doDetrend)
        self.data['R22'],self.data['tauR22'] = tt.xcorr_fft(self.data['u'][:,1], maxlags=None, norm='none',doDetrend=doDetrend)
        self.data['R33'],self.data['tauR33'] = tt.xcorr_fft(self.data['u'][:,2], maxlags=None, norm='none',doDetrend=doDetrend)
        
        
        #u in frequency domain
        self.data['u1frq'],self.data['u1amp'] = tt.dofft(sig=self.data['u'][:,0],samplefrq=self.data['frq'])
        self.data['u2frq'],self.data['u2amp'] = tt.dofft(sig=self.data['u'][:,1],samplefrq=self.data['frq'])
        self.data['u3frq'],self.data['u3amp'] = tt.dofft(sig=self.data['u'][:,2],samplefrq=self.data['frq'])
        #Time energy sectrum Se11 (mean: Rii in frequency domain...)
        self.data['Se11frq'],self.data['Se11'] = tt.dofft(sig=self.data['R11'],samplefrq=self.data['frq'])
        self.data['Se22frq'],self.data['Se22'] = tt.dofft(sig=self.data['R22'],samplefrq=self.data['frq'])
        self.data['Se33frq'],self.data['Se33'] = tt.dofft(sig=self.data['R33'],samplefrq=self.data['frq'])
        
