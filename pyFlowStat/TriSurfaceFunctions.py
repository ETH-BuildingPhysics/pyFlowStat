'''
TriSurfaceFunctions.py

Collection of functions for the following classes:
    * TriSurfaceMesh
    * TriSurfaceScalar
    * TriSurfaceVector
    * TriSurfaceSymmTensor
    * TriSurfaceTensor (to be implemented)
    
Functions included:
    * getSubTriSurfaceMesh
    * compressArray
    * getSubTriSurfaceVector
    * getSubTriSurfaceVectorList
    * saveTriSurfaceList_hdf5
    * loadTriSurfaceMesh_hdf5Parser
    * loadTriSurfaceVectorList_hdf5Parser
    * loadTriSurfaceVector_hdf5Parser
    * parseFoamFile_sampledSurface
    * parseVTK_ugly_sampledSurface
    
For documentation, see the docstring included in each function.
'''
#=============================================================================#
# load modules
#=============================================================================#
import h5py

import numpy as np

import matplotlib.tri as tri
import matplotlib.path as mplPath
import matplotlib.pyplot as plt

import pyFlowStat

#import pyFlowStat.CoordinateTransformation as coorTrans
import pyFlowStat.TriSurface as TriSurface
import pyFlowStat.TriSurfaceMesh as TriSurfaceMesh
import pyFlowStat.TriSurfaceScalar as TriSurfaceScalar
import pyFlowStat.TriSurfaceVector as TriSurfaceVector
import pyFlowStat.TriSurfaceSymmTensor as TriSurfaceSymmTensor
import pyFlowStat.Functions as func
import pyFlowStat.TriSurfaceContainer as TriSurfaceContainer

from pyFlowStat import TurbulenceTools as tt




#=============================================================================#
# functions
#=============================================================================#
def getSubTriSurfaceMesh(tsmSource, poly, op='in', mode='mid'):
    '''
    Return a sub TriSurfaceMesh defined inside or outside a polygon. This is
    a helper function. Use better getSubTriSurfaceVector.    
    
    Arguments:
        *tsmSource*: TriSurfaceMesh Object
        
        *poly*: numpy array of shape (N,2)
         list of N points of the polygon. Example for a square made by the
         points pt1, pt2, pt3 and pt4:
         >>> poly = np.array([[x1,y1],[x2,y2],[x3,y3],[x4,y4]])
         If poly ist set to None, it uses an already existing masked trinagulation.
         
        *op*: python string ('in' or 'out'). Default='in'
         Keep the triangles inside or outside the polygon

    Returns:
        *tsm*: TriSurfaceMesh object.

        *node_renum*: 
         Node renumbering. Useful for the compression of the data
    '''
    if poly:
        # create a boundBox from the polygon
        maskedTri = np.zeros(tsmSource.triangles.shape[0])
        bb = mplPath.Path(poly)
        
        # define the operation
        inPoly = None
        outPoly = None
        if op=='in':
            inPoly = 0
            outPoly = 1
        elif op=='out':
            inPoly = 1
            outPoly = 0
        else:
            raise ValueError('Argument "op" must be "in" or "out".')

        # define the mask and apply it.
        if mode=='each':
            xtri = tsmSource.x[tsmSource.triangles]
            ytri = tsmSource.y[tsmSource.triangles]
            for i in range(maskedTri.shape[0]):
                if np.sum(bb.contains_points(np.vstack((xtri[i],ytri[i])).T))>0:
                    maskedTri[i] = inPoly
                else:
                    maskedTri[i] = outPoly
        elif mode=='mid':
            xmid = tsmSource.x[tsmSource.triangles].mean(axis=1)
            ymid = tsmSource.y[tsmSource.triangles].mean(axis=1)
            for i in range(maskedTri.shape[0]):
                if bb.contains_point([xmid[i],ymid[i]])==1:
                    maskedTri[i] = inPoly
                else:
                    maskedTri[i] = outPoly
        else:
            raise ValueError('Argument "mode" must be "each" or "mid".')
        
        tsmSource.triangulation.set_mask(maskedTri)
    else:
        print 'no poly specified. Using existing triangulation'
    
    #compress the triangles and create a new TriSurface mesh
    trianalyzer = tri.TriAnalyzer(tsmSource.triangulation)
    (comp_triangles,comp_x,comp_y,node_renum) = trianalyzer._get_compressed_triangulation(return_node_renum=True)

    comp_z = compressArray(tsmSource._TriSurfaceMesh__z,node_renum)

    subTsm = TriSurfaceMesh.TriSurfaceMesh(x=comp_x,
                                           y=comp_y,
                                           z=comp_z,
                                           triangles=comp_triangles,
                                           mask=None,
                                           affTrans=tsmSource.affTrans,
                                           linTrans=tsmSource.linTrans)

    # remove the mash from the source
    tsmSource.triangulation.set_mask(None)   
    
    return subTsm,node_renum


def compressArray(z,node_renum):
    '''
    Compress an array z according the node renumbering list produced by the
    function getSubTriSurfaceMesh. See "getSubTriSurfaceMesh" and
    "getSubTriSurfaceVector" to learn about node_renum.
    
    Arguments:
        *z*: numpy array of shape (N,)
         The array to compress.
         
        *node_renum*: numpy array
         Node renumbering list generated by getSubTriSurfaceMesh or
         getSubTriSurfaceVector

    Returns:
        *comp_z*: numpy array of shape (M,)
         Compressed array    
    '''
    node_mask = (node_renum == -1)
    z[node_renum[~node_mask]] = z
    z = z[~node_mask]
    return z

   
def getSubTriSurfaceVector(tsvSource, poly, op='in', mode='mid',return_node_renum=False):
    '''
    Return a sub TriSurfaceVector (subTsv) from a TriSurfaceVector source.
    The sub part is cut out of a polygon. The part inside or ouside the 
    polygon can be kept. getSubTriSurfaceVector does not extract a sub-data
    stored in tsvSource.data. To do so, switch return_node_renum to True, and
    use the node_renum output to feed the function compressArray.
    
    Arguments:
        *tsvSource*: TriSurfaceVector Object
        
        *poly*: numpy array of shape (N,2)
         list of N points of the polygon. Example for a square made by the
         points pt1, pt2, pt3 and pt4:
         >>> poly = np.array([[x1,y1],[x2,y2],[x3,y3],[x4,y4]])
         
        *op*: python string ('in' or 'out'). Default='in'
         Keep the triangles inside or outside the polygon
         
        *return_node_renum*: python bool. Default=False
         If True, returns the node renumbering. Useful to compress array with
         TriSurfaceFunctions.compressArray()
         
    Returns:
        *subTsm*: TriSurfaceMesh object
        
        *subTsv*: TriSurfaceVector object
        
        *node_renum*: numpy array. Returned only if return_node_renum=True        
    '''
    subTsm,node_renum = getSubTriSurfaceMesh(tsmSource=tsvSource.triSurfaceMesh,
                                             poly=poly,
                                             op=op,
                                             mode=mode)
    
    comp_vx = compressArray(tsvSource.vx,node_renum)
    comp_vy = compressArray(tsvSource.vy,node_renum)
    comp_vz = compressArray(tsvSource.vz,node_renum)

    subProjectedField = tsvSource.projectedField

    subTsv = TriSurfaceVector.TriSurfaceVector(vx=comp_vx,
                                               vy=comp_vy,
                                               vz=comp_vz,
                                               time=tsvSource.time,
                                               triSurfaceMesh=subTsm,
                                               projectedField=subProjectedField,
                                               interpolation=None,
                                               kind=None)
 
    if return_node_renum==False:
        return subTsm, subTsv
    elif return_node_renum==True:
        return subTsm, subTsv, node_renum

def _getSubTriSurfaceContainer(tscSource,subTsm,node_renum):

    subTsc = TriSurfaceContainer.TriSurfaceContainer(subTsm)
    for fname,fdata in tscSource.fields.iteritems():
        if isinstance(fdata,TriSurfaceVector.TriSurfaceVector):
            comp_vx = compressArray(fdata.vx,node_renum)
            comp_vy = compressArray(fdata.vy,node_renum)
            comp_vz = compressArray(fdata.vz,node_renum)
            subProjectedField = fdata.projectedField
            subTs = TriSurfaceVector.TriSurfaceVector(vx=comp_vx,
                                                      vy=comp_vy,
                                                      vz=comp_vz,
                                                      time=fdata.time,
                                                      triSurfaceMesh=subTsm,
                                                      projectedField=subProjectedField,
                                                      interpolation=None,
                                                      kind=None)                                              
            subTsc.addTriSurface(subTs,fname)
           
        if isinstance(fdata,TriSurfaceScalar.TriSurfaceScalar):
            comp_s = compressArray(fdata.s,node_renum)
            subProjectedField = fdata.projectedField
            subTs = TriSurfaceScalar.TriSurfaceScalar(s=comp_s,
                                                      time=fdata.time,
                                                      triSurfaceMesh=subTsm,
                                                      projectedField=subProjectedField,
                                                      interpolation=None,
                                                      kind=None)                                              
            subTsc.addTriSurface(subTs,fname)
            
        if isinstance(fdata,TriSurfaceSymmTensor.TriSurfaceSymmTensor):
            comp_txx = compressArray(fdata.txx,node_renum)
            comp_txy = compressArray(fdata.txy,node_renum)
            comp_txz = compressArray(fdata.txz,node_renum)
            comp_tyy = compressArray(fdata.tyy,node_renum)
            comp_tyz = compressArray(fdata.tyz,node_renum)
            comp_tzz = compressArray(fdata.tzz,node_renum)
            subProjectedField = fdata.projectedField
            subTs = TriSurfaceSymmTensor.TriSurfaceSymmTensor(txx=comp_txx,
                                                              txy=comp_txy,
                                                              txz=comp_txz,
                                                              tyy=comp_tyy,
                                                              tyz=comp_tyz,
                                                              tzz=comp_tzz,
                                                              time=fdata.time,
                                                              triSurfaceMesh=subTsm,
                                                              projectedField=subProjectedField,
                                                              interpolation=None,
                                                              kind=None)                                              
            subTsc.addTriSurface(subTs,fname)
    return subTsc

def getSubTriSurfaceContainer(tscSource, poly, op='in', mode='mid',return_node_renum=False):
    '''
    Return a sub TriSurfaceContainer (subTsc) from a source TriSurfaceContainer.
    The sub part is cut out of a polygon. The part inside or ouside the 
    polygon can be kept. IF you want to compressed other fields based on the
    mesh stored in tscSource, switch "return_node_renum" to True, and
    use the output to feed the function compressArray.
    
    Arguments:
        *tscSource*: TriSurfaceContainer Object
        
        *poly*: numpy array of shape (N,2)
         list of N points of the polygon. Example for a square made by the
         points pt1, pt2, pt3 and pt4:
         >>> poly = np.array([[x1,y1],[x2,y2],[x3,y3],[x4,y4]])
         
        *op*: python string ('in' or 'out'). Default='in'
         Keep the triangles inside or outside the polygon
         
        *return_node_renum*: python bool. Default=False
         If True, returns the node renumbering. Useful to compress array with
         TriSurfaceFunctions.compressArray()
         
    Returns:
        *subTsc*: TriSurfacContainer object
        
        *node_renum*: numpy array. Returned only if return_node_renum=True        
    '''
    subTsm,node_renum = getSubTriSurfaceMesh(tsmSource=tscSource.triSurfaceMesh,
                                             poly=poly,
                                             op=op,
                                             mode=mode)
    subTsc = _getSubTriSurfaceContainer(tscSource,subTsm,node_renum)

    if return_node_renum==False:
        return subTsc
    elif return_node_renum==True:
        return subTsc, node_renum

def getSubTriSurfaceContainerList(tscSourceList, poly, op='in', mode='mid',return_node_renum=False):
    '''
    Return a sub TriSurfaceContainerList (subTscList) from a source TriSurfaceContainerList.
    The sub part is cut out of a polygon. The part inside or ouside the 
    polygon can be kept. IF you want to compressed other fields based on the
    mesh stored in tscSource, switch "return_node_renum" to True, and
    use the output to feed the function compressArray.
    
    Arguments:
        *tscSource*: TriSurfaceContainer Object
        
        *poly*: numpy array of shape (N,2)
         list of N points of the polygon. Example for a square made by the
         points pt1, pt2, pt3 and pt4:
         >>> poly = np.array([[x1,y1],[x2,y2],[x3,y3],[x4,y4]])
         
        *op*: python string ('in' or 'out'). Default='in'
         Keep the triangles inside or outside the polygon
         
        *return_node_renum*: python bool. Default=False
         If True, returns the node renumbering. Useful to compress array with
         TriSurfaceFunctions.compressArray()
         
    Returns:
        *subTsc*: TriSurfacContainer object
        
        *node_renum*: numpy array. Returned only if return_node_renum=True        
    '''
    tscSource=tscSourceList[0]
    subTsm,node_renum = getSubTriSurfaceMesh(tsmSource=tscSource.triSurfaceMesh,
                                             poly=poly,
                                             op=op,
                                             mode=mode)
    
    
    subTscList=[]
    for tscSource in tscSourceList:
        subTsc = _getSubTriSurfaceContainer(tscSource,subTsm,node_renum)
        subTscList.append(subTsc)

    if return_node_renum==False:
        return subTscList
    elif return_node_renum==True:
        return subTscList, node_renum

def getSubTriSurfaceVectorList(tsvListSource,
                               poly,
                               op='in',
                               mode='mid',
                               return_node_renum=False):
    '''
    Same as getSubTriSurfaceVector, but with a list of TriSurfaceVector list.
    
    Arguments:
        *tsvSource*: python list of lense N.
         List of source TriSurfaceVector Object.
        
        *poly*: numpy array of shape (N,2)
         list of N points of the polygon. Example for a square made by the
         points pt1, pt2, pt3 and pt4:
         >>> poly = np.array([[x1,y1],[x2,y2],[x3,y3],[x4,y4]])
         
        *op*: python string ('in' or 'out'). Default='in'
         Keep the triangles inside or outside the polygon
         
        *return_node_renum*: python bool. Default=False
         If True, returns the node renumbering. Useful to compress array with
         TriSurfaceFunctions.compressArray()
         
    Returns:
        *subTsm*: TriSurfaceMesh object
        
        *subTsvList*: python list of lense N.
         List of sub TriSurfaceVector object
        
        *node_renum*: numpy array. Returned only if return_node_renum=True      
    '''    
    subTsm,node_renum = getSubTriSurfaceMesh(tsmSource=tsvListSource[0].triSurfaceMesh,
                                             poly=poly,
                                             op=op,
                                             mode=mode)

    subTsvList = []    
    
    for tsvSource in tsvListSource:
        comp_vx = compressArray(tsvSource.vx,node_renum)
        comp_vy = compressArray(tsvSource.vy,node_renum)
        comp_vz = compressArray(tsvSource.vz,node_renum)
    
        subProjectedField = tsvSource.projectedField
    
        subTsv = TriSurfaceVector.TriSurfaceVector(vx=comp_vx,
                                                   vy=comp_vy,
                                                   vz=comp_vz,
                                                   time=tsvSource.time,
                                                   triSurfaceMesh=subTsm,
                                                   projectedField=subProjectedField,
                                                   interpolation=None,
                                                   kind=None)
        subTsvList.append(subTsv)
        
    if return_node_renum==False:
        return subTsm, subTsvList
    elif return_node_renum==True:
        return subTsm, subTsvList, node_renum

        
def saveTriSurfaceContainerList_hdf5(triSurfaceContainerList,hdf5fileName,names=[],indexingMode='time'):
    '''
    Save a list of TriSurface<type> in a hdf5 file.
    
    Arguments:
        *triSurfaceList*: python list
         Python list of TriSurface<type> objects.
         
        *varName* python string.
         Name of the variable <type>.
        
        *hdf5fileName*: python string
         Name of the target hdf5 file. It can also include the path.
         
        *extraVar*: python list of string.
         TriSurface<type> can holds extra data in the dict TriSurface<type>.data.
         Use extraVar to include them. Default=[] (empty list).
         
        *indexingMode*: python string
         Defines the key of each surfaces saved in the hdf5. Can be 'time' or
         'index'. Delault='time'.
         
    Returns:
        None
    '''
    tCont0=triSurfaceContainerList[0]
    
    fwm = h5py.File(hdf5fileName, 'w')
    try:
        #save the surface mesh
        gName = 'mesh'
        gMesh = fwm.create_group(gName)

        gMesh.create_dataset('points',data=tCont0.triSurfaceMesh.rawPoints())
        gMesh.create_dataset('faces',data=tCont0.triSurfaceMesh.triangles)
        
        #TODO: make nicer, get time from triSurfaceContainer, or check surfaces
        timeList=[[tCont[ts].time for ts in tCont.fields][0] for tCont in triSurfaceContainerList]
        
        for i,triSurfaceCont in enumerate(triSurfaceContainerList):
            # group name
            gName = str()
            if indexingMode=='time':
                gName = str(timeList[i])
            elif indexingMode=='index':
                gName = str(i)
            else:
                raise ValueError('Argument "indexingMode" must be "time" or "index".')
            
            gsurfi = fwm.create_group(gName)

            # save  data
            gsurfi.create_dataset('time',data=timeList[i])

            if len(names)!=0:
                for name in names:
                    gsurfi.create_dataset(name,data=triSurfaceCont[name].rawVars())
            else:
                for name in triSurfaceCont.fields.keys():
                    gsurfi.create_dataset(name,data=triSurfaceCont[name].rawVars())
        
    finally:
        fwm.close()
    


def saveTriSurfaceList_hdf5(triSurfaceList,varName,hdf5file,extraVar=[],indexingMode='time'):
    '''
    Save a list of TriSurface<type> in a hdf5 file.
    
    Arguments:
        *triSurfaceList*: python list
         Python list of TriSurface<type> objects.
         
        *varName* python string.
         Name of the variable <type>.
        
        *hdf5file*: python string
         Name of the target hdf5 file. It can also include the path.
         
        *extraVar*: python list of string.
         TriSurface<type> can holds extra data in the dict TriSurface<type>.data.
         Use extraVar to include them. Default=[] (empty list).
         
        *indexingMode*: python string
         Defines the key of each surfaces saved in the hdf5. Can be 'time' or
         'index'. Delault='time'.
         
    Returns:
        None
    '''
    fwm = h5py.File(hdf5file, 'w-')
    try:
        #save the surface mesh
        gName = 'mesh'
        gMesh = fwm.create_group(gName)

        gMesh.create_dataset('points',data=triSurfaceList[0].rawPoints())
        gMesh.create_dataset('faces',data=triSurfaceList[0].triangles)
        
        for i in range(len(triSurfaceList)):
            # group name
            gName = str()
            if indexingMode=='time':
                gName = str(triSurfaceList[i].time)
            elif indexingMode=='index':
                gName = str(i)
            else:
                raise ValueError('Argument "indexingMode" must be "time" or "index".')
            
            gsurfi = fwm.create_group(gName)

            # save  data
            gsurfi.create_dataset('time',data=triSurfaceList[i].time)
            gsurfi.create_dataset(varName,data=triSurfaceList[i].rawVars())
            
            if len(extraVar)!=0:
                for var in extraVar:
                    gsurfi.create_dataset(var,data=triSurfaceList[i][var])
            
    finally:
        fwm.close()
 
def loadTriSurfaceVectorList_hdf5Parser(hdf5Parser,
                                        varName,
                                        TriSurfaceMesh,
                                        projectedField=False):
    '''
    Helper function. See loadTriSurfaceVectorList_hdf5.
    '''
    # create the TriSurface list
    triSurfaceList = []
    keys = hdf5Parser.keys()
    keys.sort()
    try:
        keys.pop(keys.index('mesh'))
    except:
        pass
    
    for key in keys:
        tsv = TriSurfaceVector.TriSurfaceVector.readFromHdf5(hdf5Parser=hdf5Parser,
                                                             varName=varName,
                                                             triSurfaceMesh=TriSurfaceMesh,
                                                             time=key,
                                                             projectedField=projectedField)
        triSurfaceList.append(tsv)
        
    return triSurfaceList


def loadTriSurfaceVectorList_hdf5(hdf5file,
                                  varName,
                                  xViewBasis,
                                  yViewBasis=None,
                                  viewAnchor=(0,0,0),
                                  srcBasisSrc=[[1,0,0],[0,1,0],[0,0,1]],
                                  projectedField=False):
    '''
    Load all (N) TriSurfaceVectors stored in "hdf5file". the TriSurfaceMesh
    object associated to the surfaces is also returned.
    
    Arguments:
        *hdf5file*: python string.
         Name or path to the source hdf5 file.
         
        *varName*: python string.
         Name of the variable to load.
         
        *xViewBasis*: python array of shape=3.
         X direction of the surface, defined in the source base.
         
        *viewAnchor*: python array of shape=3.
         Origin of the surface coordinate system, defined in the source base.
         Default=(0,0,0)
         
        *srcBasisSrc*: python array of shape=3x3.
         Default=[[1,0,0],[0,1,0],[0,0,1]].
         
        *projectedField*: bool.
         Default=False.
    
    Returns:
        *tsvList*: python list with N entries
         List filled with TriSurfaceVector objects.
         
        *tsm* TriSurfaceMesh object
    '''
    # open the hdf5 parser
    fr = h5py.File(hdf5file, 'r')
    try:
        # load the mesh
        tsm = TriSurfaceMesh.TriSurfaceMesh.readFromHdf5(hdf5Parser=fr,
                                                         xViewBasis=xViewBasis,
                                                         yViewBasis=yViewBasis,
                                                         viewAnchor=viewAnchor,
                                                         srcBasisSrc=srcBasisSrc)
        
       
        tsvList = loadTriSurfaceVectorList_hdf5Parser(hdf5Parser=fr,
                                                      varName=varName,
                                                      TriSurfaceMesh=tsm,
                                                      projectedField=projectedField)  
       
    finally:
        fr.close()    
    return tsvList, tsm

def loadTriSurfaceContainerList_hdf5(hdf5FileName,
                                           varNames,
                                           xViewBasis,
                                           yViewBasis=None,
                                           viewAnchor=(0,0,0),
                                           srcBasisSrc=[[1,0,0],[0,1,0],[0,0,1]],
                                           projectedField=False,
                                           minVal=None,maxVal=None,step=1):
    '''
    Load all (N) TriSurfaceVectors stored in "hdf5file". the TriSurfaceMesh
    object associated to the surfaces is also returned.
    
    Arguments:
        *hdf5FileName*: python string.
         Name or path to the source hdf5 file.
         
        *varNames*: python list of strings.
         Name of the variable to load. set to [] to load all fields
         
        *xViewBasis*: python array of shape=3.
         X direction of the surface, defined in the source base.
         
        *viewAnchor*: python array of shape=3.
         Origin of the surface coordinate system, defined in the source base.
         Default=(0,0,0)
         
        *srcBasisSrc*: python array of shape=3x3.
         Default=[[1,0,0],[0,1,0],[0,0,1]].
         
        *projectedField*: bool.
         Default=False.
    
    Returns:
        *tsContainerList*: python list with N entries
         List filled with TriSurfaceContainer objects.
    '''
    
    tsContainerList=[]
    # open the hdf5 parser
    fr = h5py.File(hdf5FileName, 'r')
    try:
        tsContainerList=loadTriSurfaceContainerList_hdf5Parser(hdf5Parser=fr,
                                           varNames=varNames,
                                           xViewBasis=xViewBasis,
                                           yViewBasis=yViewBasis,
                                           viewAnchor=viewAnchor,
                                           srcBasisSrc=srcBasisSrc,
                                           projectedField=projectedField,
                                           minVal=minVal,maxVal=maxVal,step=step)

    finally:
        fr.close()    
    return tsContainerList
 
    
def loadTriSurfaceContainerList_hdf5Parser(hdf5Parser,
                                           varNames,
                                           xViewBasis,
                                           yViewBasis=None,
                                           viewAnchor=(0,0,0),
                                           srcBasisSrc=[[1,0,0],[0,1,0],[0,0,1]],
                                           projectedField=False,
                                           minVal=None,maxVal=None,step=1):
    '''
    Load all (N) TriSurfaceVectors stored in "hdf5file". the TriSurfaceMesh
    object associated to the surfaces is also returned.
    
    Arguments:
        *hdf5Parser*: h5File object.
         
        *varNames*: python list of strings.
         Name of the variable to load. set to [] to load all fields
         
        *xViewBasis*: python array of shape=3.
         X direction of the surface, defined in the source base.
         
        *viewAnchor*: python array of shape=3.
         Origin of the surface coordinate system, defined in the source base.
         Default=(0,0,0)
         
        *srcBasisSrc*: python array of shape=3x3.
         Default=[[1,0,0],[0,1,0],[0,0,1]].
         
        *projectedField*: bool.
         Default=False.
    
    Returns:
        *tsContainerList*: python list with N entries
         List filled with TriSurfaceContainer objects.
    '''
    # create the TriSurfaceContainer list
    tscList = []
    
    # get a valid list of all time step 
    allTs = []
    allTs = hdf5Parser.keys()
    try:
        allTs.pop(allTs.index('mesh'))
    except:
        pass
    allTsorted = func.sortNumStrList(allTs,minVal=minVal,maxVal=maxVal,step=step)
    
    # load the mesh
    tsm = TriSurfaceMesh.TriSurfaceMesh.readFromHdf5(hdf5Parser=hdf5Parser,
                                                     xViewBasis=xViewBasis,
                                                     yViewBasis=yViewBasis,
                                                     viewAnchor=viewAnchor,
                                                     srcBasisSrc=srcBasisSrc)
    
    # TriSurfaceContainer list
    for ts in allTsorted:
        tsc = TriSurfaceContainer.TriSurfaceContainer(tsm)                                  
        res=tsc.addFieldFromHdf5(hdf5Parser,names=varNames,key=str(ts),projectedField=projectedField)
        if not res:
            return tscList
        tscList.append(tsc)
        
    return tscList


def computeVectorStatistics(tscl,field,start=None,stop=None,step=None,fieldOnly=False):
    '''
    Compute the mean and the covariance of a time dependent vector field.
    The snapshot of the vector field are stored into a TriSrufaceContainerList.
    The range and the number of snapshot used for
    the statistics can be defined with start, stop and step.
    
    Arguments:
        *tsvl*: list of triSurfaceContainer
         Such list can be easily generate with the functions
         "loadTriSurfaceContainerList_hdf5".
        
        *field*: string
         Name of the field.
         
        *start*: integer
         Begin of the list used for the statistics. Default: start=None.
         
        *stop*: integer
         End of the list used for the statistics. Default: end=None.
         
        *step*: integer
         Subsampling of the list. Default: step=None.
         
        *fieldOnly*: bool
         If True, only the Mean field and the covariance field are returned as
         two numpy array. If False, the mean field and the covariance field are
         returned in a TriSurfaceVector and TriSurfaceSymmTensor object
         respectively. Default: fieldOnly=False
    '''
    # create a big multi dimentional array with the following structure
        # 1st dim = time
        # 2nd dim = data
        # 3rd dim = componant
    tsvTimeSeries = np.zeros([ len(tscl[start:stop:step]), tscl[0][field](0).shape[0] ,  3  ])
    for i in range(tsvTimeSeries.shape[0]):
        tsvTimeSeries[i,:,0] = tscl[i][field](0)
        tsvTimeSeries[i,:,1] = tscl[i][field](1)
        tsvTimeSeries[i,:,2] = tscl[i][field](2)
        
    # compute Mean and covariances.
    # np.cov has no axis parameter, therefore a for loop is needed...
    vMean = np.mean(tsvTimeSeries,axis=0)
    vCov = np.zeros([vMean.shape[0],6])
    for pt in range(tsvTimeSeries.shape[1]):
        Uxt = tsvTimeSeries[:,pt,0]
        Uyt = tsvTimeSeries[:,pt,1]
        Uzt = tsvTimeSeries[:,pt,2]
        data = np.vstack((Uxt,Uyt,Uzt))
        covdata = np.cov(data)
        vCov[pt,0] =  covdata[0,0]
        vCov[pt,1] =  covdata[0,1]
        vCov[pt,2] =  covdata[0,2]
        vCov[pt,3] =  covdata[1,1]
        vCov[pt,4] =  covdata[1,2]
        vCov[pt,5] =  covdata[2,2]
        
    #create a new triSurfaceContainer with the new data vMean and vCov
    if fieldOnly==False:
        tsvMean = TriSurfaceVector.TriSurfaceVector(
                                   vx=vMean[:,0],
                                   vy=vMean[:,1],
                                   vz=vMean[:,2],
                                   time=tscl[0][field].time,
                                   triSurfaceMesh=tscl[0][field].triSurfaceMesh)
        tsstCov = TriSurfaceSymmTensor.TriSurfaceSymmTensor(
                                        txx=vCov[:,0],
                                        txy=vCov[:,1],
                                        txz=vCov[:,2],
                                        tyy=vCov[:,3],
                                        tyz=vCov[:,4],
                                        tzz=vCov[:,5],
                                        time=tscl[0][field].time,
                                        triSurfaceMesh=tscl[0][field].triSurfaceMesh)
        return tsvMean, tsstCov
    else:
        return vMean,vCov


def getSortedTimes_hdf5(hdf5fileName,asFloat=True):
    '''
    returns times of the surfaces stored in h5 file, by reading the keys
    
    Arguments:
        *asFloat*: bool
        return original keys if false, else convert to float np.array
    '''
    
    hdf5Parser = h5py.File(hdf5fileName, 'r')
    allTsorted=[]
    try:
        allTs = []
        allTs = hdf5Parser.keys()
        try:
            allTs.pop(allTs.index('mesh'))
        except:
            pass
        allTsorted = func.sortNumStrList(allTs)
        
        if asFloat:
            allTsorted=np.array(allTsorted,dtype=float)
            print '#:',len(allTsorted),'min:',np.min(allTsorted),'max:',np.max(allTsorted)
            return allTsorted
        else:
            return allTsorted
    except Exception as e:
        print e
    finally:
        hdf5Parser.close()

def getIndex(tCont,x_ref,y_ref):
    r_list=np.abs(tCont.triSurfaceMesh.x-x_ref)+np.abs(tCont.triSurfaceMesh.y-y_ref)
    i_ref= np.argmin(r_list)
    x=tCont.triSurfaceMesh.x[i_ref]
    y=tCont.triSurfaceMesh.y[i_ref]
    return i_ref,(x,y)

def twoPointCorr(tContList,field,x_ref,y_ref,comp=0,idx_lst=[]):
    '''
    Given a list of TriSurfaceContainers, computes the two point correleation wrt to a reference position x_ref,y_ref.
    
    Arguments:
        *tContList*: list of TriSurfaceContainer
        
        *field*: string, key to field to compute cross correlation on
        
        *x_ref*: float
        
        *y_ref*: float
        
        *comp*: int
        component to use
        
        *idx_lst*: list of int
        indices of where the two-point correleation will be computed. Has to include i_ref
        
    Returns:
        *CCorr_u*: numpy array
        
        *i_ref*: int, index of reference position
        
        *x*:     float, actual x reference position
        
        *y*:     float, actual y reference position
    '''
    
    tCont=tContList[0]
    i_ref,(x,y)=getIndex(tCont,x_ref,y_ref)
    
    
    U_ref=[tc[field](comp)[i_ref] for tc in tContList]
    CCorr_u=[]
    
    if len(idx_lst)<1:
        idx_lst=range(len(tCont.triSurfaceMesh.x))
    if i_ref not in idx_lst:
        raise ValueError('i_ref not in list')
    
    for i in idx_lst:
        U_i=[tc[field](comp)[i] for tc in tContList]
        cc=pyFlowStat.TurbulenceTools.twoPointCorr(U_ref,U_i,subtractMean=True,norm=True)
        CCorr_u.append(cc)
    CCorr_u=np.array(CCorr_u)
    return CCorr_u,i_ref,x,y
    
    
def getCCorrHorizontal(tContList,comp,x_ref,y_ref,field='U'):
    tCont=tContList[0]

    i_ref,(x,y)=tCont.triSurfaceMesh.getIndex(x_ref,y_ref)
    x_pos,idx_x=tCont.triSurfaceMesh.getHorizontalLine(x_ref,y_ref)
    
    idx_x_r=idx_x[x_pos>=x]
    idx_x_l=idx_x[x_pos<=x]
    x_pos_r=tCont.triSurfaceMesh.x[idx_x_r]
    x_pos_l=tCont.triSurfaceMesh.x[idx_x_l]
    
    ccorr_x_r,_,_,_=twoPointCorr(tContList,field,x_ref,y_ref,comp=comp,idx_lst=idx_x_r)
    ccorr_x_l,_,_,_=twoPointCorr(tContList,field,x_ref,y_ref,comp=comp,idx_lst=idx_x_l)
    #ccorr_y,idx,x,y=twoPointCorr(tContList,'U',x_ref,y_ref,comp=comp,idx_lst=idx_y)

    x_r=np.abs(x_pos_r-x)
    x_l=np.abs(x_pos_l[::-1]-x)
    #x_r=x_pos[len(x_pos)//2:]-x_ref
    #x_l=x_pos[:len(x_pos)//2+1][::-1]-x_ref
    #ccorr_x_r=ccorr_x[len(ccorr_x)//2:]
    #ccorr_x_l=ccorr_x[:len(ccorr_x)//2+1][::-1]
    return x_l,ccorr_x_l[::-1],x_r,ccorr_x_r

def getMeanCCorrHorizontal(tContList,comp,x_ref_list,y_ref,doPlot=False,field='U'):
    x_l_lst=[]
    ccorr_x_l_lst=[]
    x_r_lst=[]
    ccorr_x_r_lst=[]
    for x_ref in x_ref_list:
        x_l,ccorr_x_l,x_r,ccorr_x_r=getCCorrHorizontal(tContList,comp,x_ref,y_ref,field=field)
        
        x_l_lst.append(x_l)
        ccorr_x_l_lst.append(ccorr_x_l)
        x_r_lst.append(x_r)
        ccorr_x_r_lst.append(ccorr_x_r)

    x_l_lst=np.array(x_l_lst)
    ccorr_x_l_lst=np.array(ccorr_x_l_lst)
    x_r_lst=np.array(x_r_lst)
    ccorr_x_r_lst=np.array(ccorr_x_r_lst)
    
    minidx_l=np.min([len(x) for x in x_l_lst])
    minidx_r=np.min([len(x) for x in x_r_lst])
    minidx=np.min([minidx_l,minidx_r])

    ccorr_m=np.mean(np.array([cc[:minidx] for cc in ccorr_x_l_lst]),axis=0)
    x_m=x_l_lst[0][:minidx]
    
    if doPlot:
        fig=plt.figure(figsize=(8.27,11.69/2))
        ax=plt.subplot(1,1,1)
        for i,x in enumerate(x_l_lst):
            ax.plot(x_l_lst[i][:minidx],ccorr_x_l_lst[i][:minidx],lw=0.5)
            ax.plot(x_r_lst[i][:minidx],ccorr_x_r_lst[i][:minidx],lw=0.5)
            ax.plot(x_m,ccorr_m,'r--',lw=4.0)
        plt.tight_layout()

    return x_m,ccorr_m#x_l_lst,ccorr_x_l_lst,x_r_lst,ccorr_x_r_lst

def getMeanHorizontalScale(tContList,comp,x_lst,ylist,field='U',scale=1000.0):
    '''
    Compute the vertical profile of the horizontal lengthscale of the velocity 
    componant "comp". A list of triSurfaceContainer is used as the source
    data. The vertical locations are defined with "ylist". To get a better 
    estimation of the lengthscale, sevral x position can be defined with
    "x_lst".
    
    The lengthscale is defined as the area under the two-point correlation
    function. As the intergral can be tricky to compute with real data, an
    gaussian function is fitted to the data and used to compute the intergral.
    
    A figure is generated with the two-point correlation and the fitted function
    to check the fitting quality.
    
    Arguments:
        *tContList*: a list TriSurfaceContainer object of length N.
         The list of TriSurfaceContainer used to compute the lengthscale
         
        *comp*: int.
         Componant of the velocity.
         
        *x_lst*: list of int.
         x location of the vertical profile. x_lst are the index, not the
         position in meter. If the list is longer than one, several vertical 
         profile are averaged togther to get a better/smoother estimation of 
         the lengthscale.
         
        *ylist*: list of int.
         vertical location used for the vertical profile.
         
        *field*: string
         The name of the field used to compute its lengthscale.
         
        *scale*: float.
         Scaling factor for the plots
         
    Returns:
        *Lz*: list of float
         The horizontal lengthscale. ylist are the vertical positions.        
    '''
    #cmap=Plotting.getColorMap(ylist[0],ylist[-1],'parula')
    Lz=[]

    fig=plt.figure(figsize=(8.27,11.69/4))
    ax=plt.subplot(1,1,1)
    for h in ylist:
        tCont=tContList[0]
        y_ref=h
        
        x_m,ccorr_m=getMeanCCorrHorizontal(tContList,comp,x_lst,y_ref,field=field)
        #x_l,ccorr_x_l,x_r,ccorr_x_r=getCCorrHorizontal(tContList,x_ref,y_ref)

        l=tt.fit_gauss_correlation(x_m/scale,ccorr_m)
        Lz.append(l[0]*scale)
        ax.plot(x_m,ccorr_m)
        ax.plot(x_m,tt.func_gauss_correlation(x_m,l[0]*scale),'k--')
        #ax.plot(x_pos-1500,ccorr_x,label='h='+str(h)+' m',c=cmap.to_rgba(h))

    ax.legend(loc=2)
    #ax.set_xlim([0,400])
    #ax.set_ylim([0,1.1])
    plt.tight_layout()
    return Lz
    
def getCCorrVertical(tContList,comp,x_ref,y_ref):
    tCont=tContList[0]
    
    i_ref,(x,y)=tCont.triSurfaceMesh.getIndex(x_ref,y_ref)
    y_pos,idx_y=tCont.triSurfaceMesh.getVerticalLine(x_ref,y_ref)
    
    idx_y_up=idx_y[y_pos>=y]
    idx_y_down=idx_y[(y_pos<=y) & (y_pos>0)]
    
    y_pos_up=y_pos[y_pos>=y]
    y_pos_down=y_pos[(y_pos<=y) & (y_pos>0)]
    
    ccorr_y_up,_,_,_=twoPointCorr(tContList,'U',x_ref,y_ref,comp=comp,idx_lst=idx_y_up)
    ccorr_y_down,_,_,_=twoPointCorr(tContList,'U',x_ref,y_ref,comp=comp,idx_lst=idx_y_down)
    #ccorr_y,idx,x,y=twoPointCorr(tContList,'U',x_ref,y_ref,comp=comp,idx_lst=idx_y)
    
    y_up=np.abs(y_pos_up-y)
    y_down=np.abs(y_pos_down[::-1]-y)
    
    return y_down,ccorr_y_down[::-1],y_up,ccorr_y_up

def getMeanCCorrVertical(tContList,comp,x_ref_list,y_ref,doPlot=False):
    x_l_lst=[]
    ccorr_x_l_lst=[]
    x_r_lst=[]
    ccorr_x_r_lst=[]
    for x_ref in x_ref_list:
        x_l,ccorr_x_l,x_r,ccorr_x_r=getCCorrVertical(tContList,comp,x_ref,y_ref)
        
        x_l_lst.append(x_l)
        ccorr_x_l_lst.append(ccorr_x_l)
        x_r_lst.append(x_r)
        ccorr_x_r_lst.append(ccorr_x_r)

    x_l_lst=np.array(x_l_lst)
    ccorr_x_l_lst=np.array(ccorr_x_l_lst)
    x_r_lst=np.array(x_r_lst)
    ccorr_x_r_lst=np.array(ccorr_x_r_lst)
    
    minidx_l=np.min([len(x) for x in x_l_lst])
    minidx_r=np.min([len(x) for x in x_r_lst])

    #l = down
    #r = up
    
    ccorr_m_down=np.mean(np.array([cc[:minidx_l] for cc in ccorr_x_l_lst]),axis=0)
    ccorr_m_up=np.mean(np.array([cc[:minidx_r] for cc in ccorr_x_r_lst]),axis=0)
    
    
    x_m_down=x_l_lst[0][:minidx_l]
    x_m_up=x_r_lst[0][:minidx_r]
    
    if doPlot:
        fig=plt.figure(figsize=(8.27,11.69/2))
        ax=plt.subplot(1,2,1)
        for i,x in enumerate(x_l_lst):
            ax.plot(x_l_lst[i][:minidx_l],ccorr_x_l_lst[i][:minidx_l],lw=0.5)
            ax.plot(x_m_down,ccorr_m_down,'r--',lw=4.0)
            
        ax=plt.subplot(1,2,2)
        for i,x in enumerate(x_r_lst):
            ax.plot(x_r_lst[i][:minidx_r],ccorr_x_r_lst[i][:minidx_r],lw=0.5)
            ax.plot(x_m_up,ccorr_m_up,'r--',lw=4.0)
        plt.tight_layout()

    return x_m_down,ccorr_m_down,x_m_up,ccorr_m_up

def getMeanVerticalScale(tContList,comp,x_lst,ylist):
    #cmap=Plotting.getColorMap(ylist[0],ylist[-1],'parula')
    Ly=[]

    fig=plt.figure(figsize=(8.27,11.69/4))
    ax=plt.subplot(1,1,1)
    for h in ylist:
        tCont=tContList[0]
        y_ref=h
        
        x_m_down,ccorr_m_down,x_m_up,ccorr_m_up=getMeanCCorrVertical(tContList,comp,x_lst,y_ref,doPlot=False)
        #x_m,ccorr_m=getMeanCCorrHorizontal(tContList,comp,x_lst,y_ref)
        #x_l,ccorr_x_l,x_r,ccorr_x_r=getCCorrHorizontal(tContList,x_ref,y_ref)

        l=tt.fit_gauss_correlation(x_m_up[:-1]/1000.0,ccorr_m_up[:-1])
        Ly.append(l[0]*1000)
        ax.plot(x_m_up,ccorr_m_up)
        ax.plot(x_m_up,tt.func_gauss_correlation(x_m_up,l[0]*1000.0),'k--')
        #ax.plot(x_pos-1500,ccorr_x,label='h='+str(h)+' m',c=cmap.to_rgba(h))

    ax.legend(loc=2)
    #ax.set_xlim([0,400])
    #ax.set_ylim([0,1.1])
    plt.tight_layout()
    return Ly
    
def getMeanVerticalLine(tContList):
    tCont=tContList[0]
    tsm=tCont.triSurfaceMesh
    x_pos_top,_=tsm.getHorizontalLine(0,np.max(tsm.y))

    U_lst=[]
    cov_list=[]
    y_pos,idx_y=tsm.getVerticalLine(0,0)
    for y in y_pos:
        U_lst_t=[]
        for x in x_pos_top:
            j,_=tsm.getIndex(x,y)
            Ux=np.array([tc['U'](0)[j] for tc in tContList])
            Uy=np.array([tc['U'](1)[j] for tc in tContList])
            Uz=np.array([tc['U'](2)[j] for tc in tContList])
            U=np.vstack([Ux,Uy,Uz]).T
            U_lst_t.append(U)
        U_lst_t=np.array(U_lst_t)
        U_lst_t=U_lst_t.reshape((U_lst_t.shape[0]*U_lst_t.shape[1],U_lst_t.shape[2]))

        c_tmp=np.cov(U_lst_t.T)
        cov_list.append([c_tmp[0,0],c_tmp[0,1],c_tmp[0,2],c_tmp[1,1],c_tmp[1,2],c_tmp[2,2]])
        #cov_list_t=np.array(cov_list_t)

        #cov_list.extend(cov_list_t)
        U_lst.append(U_lst_t)
    U_lst=np.array(U_lst)
    cov_list=np.array(cov_list)
    
    UMean=np.mean(U_lst,axis=1)
    
    return y_pos,UMean,cov_list