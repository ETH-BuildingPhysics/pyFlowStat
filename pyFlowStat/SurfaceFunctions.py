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
        *surfaceList*: python list.
         A list of surfaces.
        *hdf5file*: python string.
         Path to the target file.
        *keyrange*: python string.
         defines which keys will be be saved the hdf5 file. if keyrange='raw',
         only vx,vy and vz are saved, if keyrange='full', 'raw' plus every keys
         in the surface are saved. Default='raw'.
         
    Returns:
        None
    '''
    fwm = h5py.File(hdf5file, 'w-')
    try:
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
                    if (key=='dx' or key=='dy' or key=='Ux' or key=='Uy' or key=='Uz'):
                        pass
                    else:
                        gsurfi.create_dataset(key,data=surfaceList[i].data[key])
    finally:
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
        *hdf5file*: string 
          path to the source file.
        *keyrange*: string
          Defines which keys will be loaded from the hdf5 file:
              * 'raw' = only vx,vy, vz, (default)
              * 'full' = 'raw' plus every keys from the hdf5 file.
        *createDict*: bool
         Create data dict. Usefull if the hdf5 contains
         only the raw data (vx, vy and vz). If keyrange='full', createDict
         is set to False, whatever the user selects.

    Returns:
        *surfaceList*: python list
         List of Surface object.
    '''
    surfaceList = []
    fr = h5py.File(hdf5file, 'r')
    try:
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
                    if (key=='dim' or key=='dimExtent'):
                        pass
                    elif (key=='vx'):
                        surfaceList[i].data[str('Ux')] = fr[gName][key].value
                    elif (key=='vy'):
                        surfaceList[i].data[str('Uy')] = fr[gName][key].value
                    elif (key=='vz'):
                        surfaceList[i].data[str('Uz')] = fr[gName][key].value
                    else:
                        surfaceList[i].data[str(key)] = fr[gName][key].value

            if createDict==True and keyrange=='raw':
                surfaceList[i].createDataDict()
            else:
                pass
    finally:
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
            if ( key=='dim' or key=='dimExtent'):
                pass
            elif (key=='vx'):
                s.data[str('Ux')] = hdf5fileObj[gName][key].value
            elif (key=='vy'):
                s.data[str('Uy')] = hdf5fileObj[gName][key].value
            elif (key=='vz'):
                s.data[str('Uz')] = hdf5fileObj[gName][key].value
            else:
                s.data[str(key)] = hdf5fileObj[gName][key].value

    return s


def loadPPfromSurf_hdf5(hdf5fileObj,pixloc,keyrange='raw',createDict=False,dt=1.0,ptloc=[0.0,0.0,0.0]):
    '''
    Generate a PointPorbe object from a list of surfaces saved in a hdf5.
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
        * hdf5fileObj: an h5py file object.
        * pixloc: [python list, shape=[2]] pixel location .
        * keyrange: [string] 'raw' or 'full' (default='raw'. 'full' not implemented).
        * createDict: [bool] create PointPorbe dictonnary (default=False).
        * dt: [float] time step between the surfaces.
        * ptloc: [python list, shape=[3]] point location (default=(0,0,0)).

    Returns:
        * pt: a PointProbe object (see pyFlowStat.PointProbe)


    Examples:
    >>> import h5py
    >>> from pyFlowStat.Surface import surface
    >>> from pyFlowStat.SurfaceFunctions import SurfaceFunctions
    >>> h5obj = h5pyFile('mydata.hdf5','r')
    >>> pt = SurfaceFunctions.loadPPfromSurf_hdf5(h5obj,(140,53),keyrange='raw',dt=0.04)
    >>> h5obj.close()
    '''
    pt = pp.PointProbe()

    # get number of surfaces
    nbsurf = len(hdf5fileObj.keys())
    print nbsurf
    # generate probeVar and probeTimes
    probeVar = np.zeros((nbsurf,3))
    probeTimes = np.zeros(nbsurf)


    # append velocity v and time t to PointProbe. Add extra variable if requiered
    for i in range(len(hdf5fileObj.keys())):
        gName = 'Surface'+str(i)

        vx = hdf5fileObj[gName]['vx'][pixloc[0],pixloc[1]]
        vy = hdf5fileObj[gName]['vy'][pixloc[0],pixloc[1]]
        vz = hdf5fileObj[gName]['vz'][pixloc[0],pixloc[1]]
        probeVar[i,0] = vx
        probeVar[i,1] = vy
        probeVar[i,2] = vz

        probeTimes[i] = i*dt

#       crappy... needs improvement
#        if keyrange=='raw':
#            pass
#        elif keyrange=='full':
#            for key in hdf5fileObj[gName].keys():
#                if ( key=='dim' or key=='dimExtent' or key=='vx' or key=='vy' or key=='vz'):
#                    pass
#                else:
#                    pt.data[key].append( hdf5fileObj[gName][key][pixloc[0],pixloc[1]] )

    # add location to PointPorbe
    pt.probeLoc.append(ptloc[0])
    pt.probeLoc.append(ptloc[1])
    pt.probeLoc.append(ptloc[2])

    # add probeVar and probeTimes to pt
    pt.probeVar = probeVar
    pt.probeTimes = probeTimes


    if createDict==True:
        pt.createDataDict()

    return pt

    
def corrField(f1,f2=None,i_ref=0,j_ref=0,norm=True):
    '''
    Two point correlation of an entire field with the point pt(i_ref,j_ref) as
    reference. The field has the dimension of NxM with T time realization.
    
    Arguments:
        *f1*: np.array of shape (T,N,M).
         In general, one of the fields (vx,vy or vz) of a 
         pyFlowStat.SurfaceTimeSeries object.
    
        *f2*: np.array of shape (T,N,M).
         Same as f1. If None, then f1 is used has second field. Default=None
         
        *i_ref*: python int.
         Horizontal location of the reference point. Default=0
         
        *j_ref*: python int.
         Vertical location of the reference point. Default=0
         
        *Norm*: python bool.
         Normalization of the correlation. Default=True.
         
         
    Returns:
        *res*: np.array of shape (N,M).
         the two point horizontal correlation. 
    '''
    f1sum=np.sum(f1,axis=0)
    if f2==None:
        f2=f1
        f2sum=f1sum
    else:
        f2sum=np.sum(f2,axis=0)
        
    res=np.empty(f1[0].shape, dtype=float)
    res[:] = np.nan
    k,dx,dy=f1.shape
    if np.isnan(f1sum[i_ref,j_ref]):
        res[i_ref,j_ref] = np.nan
    else:
        for i in range(dx):
            for j in range(dy):
                if np.isnan(f2sum[i,j]):
                    res[i,j]=np.nan
                else:
                    res[i,j]=tt.twoPointCorr(x=f1[:,i_ref,j_ref],y=f2[:,i,j],norm=True)
    if norm==True:
        return res/res[i_ref,j_ref]
    else:
        return res
        
def corrFieldHorz(f1,f2=None,j_ref=0,norm=True):
    '''
    Do the two point horizontal correlation of an entire field along a vertical
    line. Each points of the line is used as a reference for the two point 
    correlation. the field has dimension of NxM with T time realization.
    
    Arguments:
        *f1*: np.array of shape (T,N,M).
         In general, one of the fields (vx,vy or vz) of a 
         pyFlowStat.SurfaceTimeSeries object.
    
        *f2*: np.array of shape (T,N,M).
         Same as f1. If None, then f1 is used has second field. Default=None
         
        *j_ref*: python int.
         Horizontal location of the vertical line. Default=0
         
        *Norm*: python bool.
         Normalization of the correlation. Default=True.
         
         
    Returns:
        *res*: np.array of shape (N,M).
         the two point horizontal correlation. 
    '''
    f1sum=np.sum(f1,axis=0)
    if f2==None:
        f2=f1
        f2sum=f1sum
    else:
        f2sum=np.sum(f2,axis=0)
        
    res=np.empty(f1[0].shape, dtype=float)
    res[:] = np.nan
    k,dx,dy=f1.shape
    for i in range(dx):
        if np.isnan(f1sum[i,j_ref]):
            res[i,:]=np.nan
        else:
            for j in range(dy):
                if np.isnan(f1sum[i,j]) or np.isnan(f2sum[i,j]):
                    res[i,j]=np.nan
                else:
                    res[i,j]=tt.twoPointCorr(x=f1[:,i,j_ref],y=f2[:,i,j],norm=True)
            if norm==True:
                res[i,:]=res[i,:]/res[i,j_ref]
      
    return res
    
def corrFieldVert(f1,f2=None,i_ref=0,norm=True):
    '''
    Do the two point vertical correlation of an entire field along an
    horizontal line. Each points of the line is used as a reference for the two
    point correlation. the field has dimension of NxM with T time realization.
    
    Arguments:
        *f1*: np.array of shape (T,N,M).
         In general, one of the fields (vx,vy or vz) of a 
         pyFlowStat.SurfaceTimeSeries object.
    
        *f2*: np.array of shape (T,N,M).
         Same as f1. If None, then f1 is used has second field. Default=None
         
        *i_ref*: python int.
         Vertical location of the horizontal line. Default=0
         
        *Norm*: python bool.
         Normalization of the correlation. Default=True.
         
         
    Returns:
        *res*: np.array of shape (N,M).
         the two point vertical correlation. 
    '''
    f1sum=np.sum(f1,axis=0)
    if f2==None:
        f2=f1
        f2sum=f1sum
    else:
        f2sum=np.sum(f2,axis=0)
       
    res=np.empty(f1[0].shape, dtype=float)
    res[:] = np.nan
    k,dx,dy=f1.shape
    for j in range(dy):
        if np.isnan(f1sum[i_ref,j]):
            res[:,j]=np.nan
        else:
            for i in range(dx):
                if np.isnan(f1sum[i,j]) or np.isnan(f2sum[i,j]):
                    res[i,j]=np.nan
                else:
                    res[i,j]=tt.twoPointCorr(x=f1[:,i_ref,j],y=f2[:,i,j],norm=True)
            if norm==True:
                res[:,j]=res[:,j]/res[i_ref,j]
            
    return res
    
def corrVert(f1,f2=None,i_ref=0,j_ref=0,norm=True):
    '''
    Do the two point vertical correlation from a reference point
    pt(i_ref,j_ref). The field has dimension of NxM with T time realization.
    
    Arguments:
        *f1*: np.array of shape (T,N,M).
         In general, one of the fields (vx,vy or vz) of a 
         pyFlowStat.SurfaceTimeSeries object.
    
        *f2*: np.array of shape (T,N,M).
         Same as f1. If None, then f1 is used has second field. Default=None
         
        *i_ref*: python int.
         Vertical location of the horizontal line. Default=0
         
        *j_ref*: python int.
         Vertical location of the horizontal line. Default=0
         
        *Norm*: python bool.
         Normalization of the correlation. Default=True.
         
         
    Returns:
        *lags_l*: numpy array.
         Array of left/negative lags.
        
        *res_l*: numpy array.
         Array of left/negative results.
        
        *lags_r*: numpy array.
         Array of right/positive lags.
        
        *res_r*: numpy array.
         Array of right/positive results.       
    '''
    f1sum=np.sum(f1,axis=0)
    if f2==None:
        f2=f1
        f2sum=f1sum
    else:
        f2sum=np.sum(f2,axis=0)
    
    res=np.empty(f1[0].shape[0], dtype=float)
    res[:] = np.nan
    k,dx,dy=f1.shape
    for i in range(dx):
        if np.isnan(f1sum[i_ref,j_ref]) or np.isnan(f2sum[i_ref,j_ref]):
            res[i]=np.nan
        else:
            res[i]=tt.twoPointCorr(x=f1[:,i_ref,j_ref],y=f2[:,i,j_ref],norm=True)
    if norm==True:
        res=res/res[i_ref]
    res_l=res[:i_ref+1][::-1]
    res_r=res[i_ref:]
    lags_l=np.arange(len(res_l))
    lags_r=np.arange(len(res_r))
    return lags_l,res_l,lags_r,res_r