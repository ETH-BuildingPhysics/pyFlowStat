'''
PointProbeFunctions.py

Collection of functions for the PointProbe class.

Functions included:
    *
'''


#=============================================================================#
# load modules
#=============================================================================#
#import sys
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
import pyFlowStat.PointProbe as pp
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
              * 'raw' = only U, t and pos (default)
              * 'full' = U, t and pos, plus all the other keys included
              in ppList[i].data

    Returns:
        * None
    '''
    fwm = h5py.File(hdf5file, 'w-')
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
    fwm.close()


def loadPPlist_hdf5(hdf5file,keyrange='raw',createDict=False):
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
              * 'raw' = only U, t and pos (default)
              * 'full' = 'raw', plus all the other keys included in ppList[i].data
        * createDict: [bool] create data dict. Usefull if the hdf5 contains
          only the raw data or if you load only the raw data from a full
          hdf5.

    Returns:
        * ppList: [python list] list of PointProbe object.
    '''
    ppList = []
    fr = h5py.File(hdf5file, 'r')
    for i in range(len(fr.keys())):
        gName = 'pointProbe'+str(i)
        #print('load '+str(gName))
        ppList.append(pp.PointProbe())
        ppList[i].probeVar = fr[gName]['probeVar'].value
        ppList[i].probeTimes = fr[gName]['probeTimes'].value
        ppList[i].probeLoc = fr[gName]['probeLoc'].value

        if keyrange=='raw':
            pass
        elif keyrange=='full':
            for key in fr[gName].keys():
                if (key=='probeVar' or key=='probeTimes' or key=='probeLoc'):
                    pass
                else:
                    ppList[i].data[str(key)] = fr[gName][key].value

        if createDict==False:
            pass
        else:
            ppList[i].createDataDict()
    fr.close()
    return ppList