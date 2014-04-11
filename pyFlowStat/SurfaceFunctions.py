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
import pyFlowStat.Surface as sr


#=============================================================================#
# functions
#=============================================================================#

def saveSurfaceList_hdf5(surfaceList,hdf5file,keyrange='raw'):
    '''
    Save a surface list in a hdf5 data file. The hdf5 file will have the
    following minimal structure:

    myData.hdf5:
        * Surface1  (GROUP)
            * 'vx'   (DATASET)
            * 'vy' (DATASET)
            * 'vz'   (DATASET)
            * 'dim'  (DATASET)
            * 'dimExtent' (DATASET)
        * Surfacei  (GROUP)
            * 'vx'   (DATASET)
            * 'vy' (DATASET)
            * 'vz'   (DATASET)
            * 'dim'  (DATASET)
            * 'dimExtent' (DATASET)

    Arguments:
        * surfaceList: [python list] a list of surfaces
        * hdf5file: [str] path to the target file.
        * keyrange: [str] defines which keys will be be saved the hdf5 file:
              * 'raw' = only vx,vy and vz (default)
              * 'full' =  'raw' plus every keys in the surface.

    Returns:
        * surfaceList: [python list] list of Surface object.
    '''
    fwm = h5py.File(hdf5file, 'w-')
    for i in range(len(surfaceList)):
        # group name
        gName = 'Surface'+str(i)
        gsurfi = fwm.create_group(gName)

        # save minimal data
        gsurfi.create_dataset('vx',data=surfaceList[i].vx)
        gsurfi.create_dataset('vy',data=surfaceList[i].vy)
        gsurfi.create_dataset('vz',data=surfaceList[i].vz)
        
        gsurfi.create_dataset('dx',data=surfaceList[i].dx)
        gsurfi.create_dataset('dy',data=surfaceList[i].dy)

        dim=[surfaceList[i].minX, surfaceList[i].minY, surfaceList[i].maxX, surfaceList[i].maxY]
        gsurfi.create_dataset('dim',data=dim)
        gsurfi.create_dataset('dimExtent',data=surfaceList[i].extent)

        # save extra data if specified:
        #save nothing more
        if keyrange=='raw':
            pass
        #add all data from the dictionnary
        elif keyrange=='full':
            for key in surfaceList[i].data.keys():
                if (key=='dx' or key=='dy'):
                    pass
                else:
                    gsurfi.create_dataset(key,data=surfaceList[i].data[key])
    fwm.close()


def loadSurfaceList_hdf5(hdf5file,keyrange='raw',createDict=False):
    '''
    Load and return a surface list from a hdf5 data file. eager evaluation
    only. The hdf5 file must have the following minimal structure:

    myData.hdf5:
        * Surface1  (GROUP)
            * 'vx'   (DATASET)
            * 'vy' (DATASET)
            * 'vz'   (DATASET)
            * 'dim'  (DATASET)
            * 'dimExtent' (DATASET)
        * Surfacei  (GROUP)
            * 'vx'   (DATASET)
            * 'vy' (DATASET)
            * 'vz'   (DATASET)
            * 'dim'  (DATASET)
            * 'dimExtent' (DATASET)

    Arguments:
        * hdf5file: [str] path to source file.
        * keyrange: [str] defines which keys will be loaded from the hdf5 file:
              * 'raw' = only vx,vy, vz, (default)
              * 'full' = 'raw' plus every keys from the hdf5 file.
        * createDict: [bool] create data dict. Usefull if the hdf5 contains
          only the raw data (vx, vy and vz).

    Returns:
        * surfaceList: [python list] list of Surface object.
    '''
    surfaceList = []
    fr = h5py.File(hdf5file, 'r')
    for i in range(len(fr.keys())):
        gName = 'Surface'+str(i)
        surfaceList.append(sr.Surface())

        # load minimum data
        surfaceList[i].vx = fr[gName]['vx'].value
        surfaceList[i].vy = fr[gName]['vy'].value
        surfaceList[i].vz = fr[gName]['vz'].value
        
        surfaceList[i].dx = fr[gName]['dx'].value
        surfaceList[i].dy = fr[gName]['dy'].value

        dim = fr[gName]['dim'].value
        surfaceList[i].minX = dim[0]
        surfaceList[i].minY = dim[1]
        surfaceList[i].maxX = dim[2]
        surfaceList[i].maxY = dim[3]
        surfaceList[i].extent = fr[gName]['dimExtent'].value

        # load extra data if specified:
        #load nothing more
        if keyrange=='raw':
            pass
        elif keyrange=='full':
            for key in fr[gName].keys():
                if (key=='vx' or key=='vy' or key=='vz' or key=='dim' or key=='dimExtent'):
                    pass
                else:
                    surfaceList[i].data[str(key)] = fr[gName][key].value

        if createDict==False:
            pass
        else:
            surfaceList[i].createDataDict()
    fr.close()
    return surfaceList


def loadSurface_hdf5(hdf5fileObj,surfaceNo,keyrange='raw'):
    '''
    load a single surface from a hdf5 file object (http://docs.h5py.org).
    The hdf5 file will have the following minimal structure:

    myData.hdf5:
        * Surface1  (GROUP)
            * 'vx'   (DATASET)
            * 'vy' (DATASET)
            * 'vz'   (DATASET)
            * 'dim'  (DATASET)
            * 'dimExtent' (DATASET)
        * Surfacei  (GROUP)
            * 'vx'   (DATASET)
            * 'vy' (DATASET)
            * 'vz'   (DATASET)
            * 'dim'  (DATASET)
            * 'dimExtent' (DATASET)

    Arguments:
        * hdf5fileObj: an h5py file object
        * surfaceNo: [int] surface number
        * keyrange: [string] 'raw' or 'full' (default='raw')

    Returns:
        * surf: a surface object (see pyFlowStat.Surface.Surface)


    Examples:
    >>> import h5py
    >>> from pyFlowStat.Surface import surface
    >>> from pyFlowStat.SurfaceFunctions import SurfaceFunctions
    >>> h5obj = h5pyFile('mydata.hdf5','r')
    >>> surf = SurfaceFunctions.loadSurface_hdf5(h5obj,2,keyrange='raw')
    >>> h5obj.close()
    '''
    s = sr.Surface()
    gName = 'Surface'+str(surfaceNo)

    #load minimal data
    s.vx = hdf5fileObj[gName]['vx'].value
    s.vy = hdf5fileObj[gName]['vy'].value
    s.vz = hdf5fileObj[gName]['vz'].value
    
    s.dx = hdf5fileObj[gName]['dx'].value
    s.dy = hdf5fileObj[gName]['dy'].value

    dim = hdf5fileObj[gName]['dim'].value
    s.minX = dim[0]
    s.minY = dim[1]
    s.maxX = dim[2]
    s.maxY = dim[3]
    s.extent = hdf5fileObj[gName]['dimExtent'].value

    # load extra data if specified:
    #load nothing more
    if keyrange=='raw':
        pass
    elif keyrange=='full':
        for key in hdf5fileObj[gName].keys():
            if (key=='vx' or key=='vy' or key=='vz' or key=='dim' or key=='dimExtent'):
                pass
            else:
                s.data[str(key)] = hdf5fileObj[gName][key].value

    return s



