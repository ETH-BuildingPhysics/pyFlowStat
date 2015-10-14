'''
PointProbeFunctions.py

Collection of functions for the PointProbe class.

Functions included:
    *
'''


#=============================================================================#
# load modules
#=============================================================================#
import sys
#import re
#import os
#import csv
#import collections
import h5py

#scientific modules
import numpy as np
import scipy as sp
#from scipy import signal
#from scipy.optimize import curve_fit
# special modules

#from pyFlowStat.TurbulenceTools import TurbulenceTools as tt
import pyFlowStat.PointProbe as PointProbe
import pyFlowStat.TurbulenceTools as tt
import pyFlowStat.Surface as Surface


#=============================================================================#
# functions
#=============================================================================#

def savePPlist_hdf5(ppList,hdf5file,keyrange='raw'):
    '''
    Save a point probe list, generate py getVectorPointProbeList for example,
    in a hdf5 data file. The hdf5 file will have the following minimal
    structure:

    myData.hdf5:
        * pointProbe1  (GROUP)
            * 'probeVar'   (DATASET)
            * 'probeTimes' (DATASET)
            * 'probeLoc'   (DATASET)
        * pointProbei  (GROUP)
            * 'probeVar'   (DATASET)
            * 'probeTimes' (DATASET)
            * 'probeLoc'   (DATASET)


    Arguments:
        * ppList: [python List] List of PointPorbe object
        * hdf5file: [str] path to target file.
        * keyrange: [str] keys included in the pointProbe which will be
          saved in the hdf5 file.
              * 'raw' = only probeVar, probeTimes and probeLoc (default)
              * 'full' = U, t and pos, plus all the other keys included
              in ppList[i].data

    Returns:
        * None
    '''
    fwm = h5py.File(hdf5file, 'w-')
    try:
        for i in range(len(ppList)):
            # group name
            gName = 'pointProbe'+str(i)
            #print('save '+str(gName))
            gppi = fwm.create_group(gName)
            # iter dict keys
            gppi.create_dataset('probeVar',data=ppList[i].probeVar)
            gppi.create_dataset('probeTimes',data=ppList[i].probeTimes)
            gppi.create_dataset('probeLoc',data=ppList[i].probeLoc)
            if keyrange=='raw':
                pass
            elif keyrange=='full':
                for key in ppList[i].data.keys():
                    gppi.create_dataset(key,data=ppList[i][key])
    finally:
        fwm.close()


def loadPPlist_hdf5(hdf5file,keyrange='raw',createDict=False,step=1,indices=None,verbose=False):
    '''
    Load and return a point probe list from a hdf5 data file. eager evaluation
    only. The hdf5 file must have the following minimal structure:

    myData.hdf5:
        * pointProbe1  (GROUP)
            * 'probeVar'   (DATASET)
            * 'probeTimes' (DATASET)
            * 'probeLoc'   (DATASET)
        * pointProbei  (GROUP)
            * 'probeVar'   (DATASET)
            * 'probeTimes' (DATASET)
            * 'probeLoc'   (DATASET)

    Arguments:
        * hdf5file: [str] path to source file.
        * keyrange: [str] keys included in the pointProbe which will be
          saved in the hdf5 file.
              * 'raw' = only probeVar, probeTimes and probeLoc (default)
              * 'full' = 'raw', plus all the other keys included in ppList[i].data
        * createDict: [bool] create data dict. Usefull if the hdf5 contains
          only the raw data or if you load only the raw data from a full
          hdf5.

    Returns:
        * ppList: [python list] list of PointProbe object.
    '''
    ppList = []
    fr = h5py.File(hdf5file, 'r')
    try:
        if not indices:
            indices=range(len(fr.keys()))
        for i in indices:
            gName = 'pointProbe'+str(i)
            if verbose:
                print('load '+str(gName))
            pp=PointProbe.PointProbe()
            pp.probeVar = fr[gName]['probeVar'][::step,...]
            pp.probeTimes = fr[gName]['probeTimes'][::step]
            pp.probeLoc = fr[gName]['probeLoc'].value

            if keyrange=='raw':
                pass
            elif keyrange=='full':
                for key in fr[gName].keys():
                    if (key=='probeVar' or key=='probeTimes' or key=='probeLoc'):
                        pass
                    else:
                        pp.data[str(key)] = fr[gName][key].value

            if createDict==False:
                pass
            else:
                pp.createDataDict()
            ppList.append(pp)
    finally:
        fr.close()
    return ppList

def actionPPlist_hdf5(hdf5file,actionFunction):
    '''
    Load and return a point probe list from a hdf5 data file. eager evaluation
    only. The hdf5 file must have the following minimal structure:

    myData.hdf5:
        * pointProbe1  (GROUP)
            * 'probeVar'   (DATASET)
            * 'probeTimes' (DATASET)
            * 'probeLoc'   (DATASET)
        * pointProbei  (GROUP)
            * 'probeVar'   (DATASET)
            * 'probeTimes' (DATASET)
            * 'probeLoc'   (DATASET)

    Arguments:
        * hdf5file: [str] path to source file.
        * keyrange: [str] keys included in the pointProbe which will be
          saved in the hdf5 file.
              * 'raw' = only probeVar, probeTimes and probeLoc (default)
              * 'full' = 'raw', plus all the other keys included in ppList[i].data
        * createDict: [bool] create data dict. Usefull if the hdf5 contains
          only the raw data or if you load only the raw data from a full
          hdf5.

    Returns:
        * ppList: [python list] list of PointProbe object.
    '''
    ppList = []
    resultList=[]
    
    fr = h5py.File(hdf5file, 'r')
    try:
        for i in range(len(fr.keys())):
            gName = 'pointProbe'+str(i)
            #print('load '+str(gName))
            probe=PointProbe.PointProbe()
            probe.probeVar = fr[gName]['probeVar'].value
            probe.probeTimes = fr[gName]['probeTimes'].value
            probe.probeLoc = fr[gName]['probeLoc'].value
            resultList.append(actionFunction(probe))
            probe=None
    finally:
        fr.close()
    return resultList
  
def createPointProbeFromSurfaceTimeSeries(surfaceTimeSeries,frq,i,j,doDetrend=True,createDict=True,genStat=True):
    '''
    Create a PointProbe from time resolved field data for a selected location i,j.
    adds the velocity vector time series, the times, probeLoc and calls 
    createDataDict and generateStatistics by default.

    Arguments:
        * surfaceTimeSeries: [SurfaceTimeSeries] pyFlowStat.Surface.SurfaceTimeSeries object
        * frq: [float]  Sample frequency
        * i: [int] Row of the point (in pixel).
        * j: [int] Column of the point (in pixel).
        * doDetrend: [bool] apply detrending in generateStatistics(). Default=True
        * createDict: [bool] call createDataDict. Default=True
        * genStat: [bool] call generateStatistics. Default=True

    Returns:
        * pt: [PointProbe] PointProbe object.
    '''
    
    vel=np.column_stack((surfaceTimeSeries.vx[:,i,j],surfaceTimeSeries.vy[:,i,j],surfaceTimeSeries.vz[:,i,j]))
    #if not np.isnan(vel).any():
    try:
        pt=PointProbe.PointProbe()
        pt.probeVar=vel
        pt.probeTimes=surfaceTimeSeries.t
        pt.probeLoc=[i,j]
        pt.createDataDict(action=createDict)
        if genStat==True:
            pt.addVectorMagnitude()
            pt.generateStatistics(doDetrend=doDetrend)

        return pt
    except:
        print "Unexpected error:", sys.exc_info()[0]
        return None
 
def createFromArray(dt,u,v,w,doDetrend=True):
    '''
    Create a PointProbe from three scalar arrays (velocity components)
    adds the velocity vector time series, the times, and calls createDataDict and generateStatistics

    Arguments:
        * u,v,w:   [numpy.array shape=(N)] Coordinate of probe (must be included in ofFile)
        * dt: [float]  Time step
        * doDetrend: [bool] apply detrending in generateStatistics()

    Returns:
        * pt: [PointProbe] PointProbe object.
    '''
    pt=PointProbe.PointProbe()
    pt.probeVar=np.transpose(np.vstack((u,v,w)))
    nrpoints=len(u)
    endtime=(len(u)-1)*dt
    pt.probeTimes=np.linspace(0,endtime,nrpoints)
    pt.createDataDict()
    pt.addVectorMagnitude() 
    pt.generateStatistics(doDetrend=doDetrend)
    return pt