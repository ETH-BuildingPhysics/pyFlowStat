'''
TriSurfaceFunctions.py

Collection of functions for the following classes:
    * TriSurfaceVector
    * TriSurfaceTensor (to be implemented)
    * TriSurfaceScalar (to be implemented)

Functions included:
    *
'''
#=============================================================================#
# load modules
#=============================================================================#
import h5py
import re

import numpy as np

import matplotlib.tri as tri
import matplotlib.path as mplPath

import CoordinateTransformation as coorTrans
import TriSurfaceVector as TriSurfaceVector
import TriSurfaceMesh as TriSurfaceMesh



#=============================================================================#
# functions
#=============================================================================#
def getTransformation(viewAnchor,
                      xViewBasis,
                      yViewBasis,
                      srcBasisSrc=[[1,0,0],[0,1,0],[0,0,1]]):
    '''
    Return the affine and linear transfomation obecjt to move from the
    coordinate system of the 3D CFD simulation to the 2D surface.
    
    Arguments:
    
    Returns:
        *affTrans*: AffineTransformation object.
            
        *linTrans*: LinearTransformation object.
    '''
    # check and convert arguments
    srcBasisSrc = np.array(srcBasisSrc,dtype=float)
    if srcBasisSrc.shape!=(3,3):
        raise ValueError('srcBasis must be a 3x3 matrix')
        
    xViewBasis = np.array(xViewBasis,dtype=float)
    yViewBasis = np.array(yViewBasis,dtype=float)
    if xViewBasis.shape!=(3,) or yViewBasis.shape!=(3,):
        raise ValueError('xViewBasis.shape and yViewBasis. ',
                         'shape must be equal to (3,)')
        
    # get the basis and the transformation object
    tgtBasisSrc = np.zeros((3,3))
    tgtBasisSrc[:,0] = xViewBasis
    tgtBasisSrc[:,1] = yViewBasis
    tgtBasisSrc[:,2] = np.cross(xViewBasis,yViewBasis)
    afftrans = coorTrans.AffineTransfomation(srcBasisSrc,tgtBasisSrc,viewAnchor)
    lintrans = coorTrans.LinearTransformation(srcBasisSrc,tgtBasisSrc)
    
    return afftrans, lintrans


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
         
        *op*: python string ('in' or 'out'). Default='in'
         Keep the triangles inside or outside the polygon

    Returns:
        *tsm*: TriSurfaceMesh object.

        *node_renum*: 
         Node renumbering. Useful for the compression of the data
    '''
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

    return subTsm,node_renum


def compressArray(z,node_renum):
    '''
    Compress an array according the node renumbering list produced by the
    function getSubTriSurfaceMesh.
    
    Arguments:
    
    
    Returns:
    
    '''
    node_mask = (node_renum == -1)
    z[node_renum[~node_mask]] = z
    z = z[~node_mask]
    return z
    
def getSubTriSurfaceVector(tsvSource, poly, op='in', mode='mid',return_node_renum=False):
    '''
    Do stuff...    
    
    Arguments:
        *tsmSource*: TriSurfaceMesh Object
        
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
        
        *node_renum*: numpy array. Returned only if return_node_renum=Trus
         
    '''
    subTsm,node_renum = getSubTriSurfaceMesh(tsmSource=tsvSource.triSurfaceMesh,
                                             poly=poly,
                                             op=op,
                                             mode=mode)
    
    comp_vx = compressArray(tsvSource.vx,node_renum)
    comp_vy = compressArray(tsvSource.vy,node_renum)
    comp_vz = compressArray(tsvSource.vz,node_renum)

    
    subProjectedField = tsvSource._TriSurfaceVector__projectedField
    subInterpolation = tsvSource._TriSurfaceVector__interType
    subKind = tsvSource._TriSurfaceVector__interKind

    subTsv = TriSurfaceVector.TriSurfaceVector(vx=comp_vx,
                                               vy=comp_vy,
                                               vz=comp_vz,
                                               time=tsvSource.time,
                                               triSurfaceMesh=subTsm,
                                               projectedField=subProjectedField,
                                               interpolation=subInterpolation,
                                               kind=subKind)
                                               
                                              
    
    if return_node_renum==False:
        return subTsm, subTsv
    elif return_node_renum==True:
        return subTsm, subTsv, node_renum

   
def saveTriSurfaceList_hdf5(triSurfaceList,hdf5file):
    '''
    '''
    fwm = h5py.File(hdf5file, 'w-')
    
    #save the surface mesh
    gName = 'mesh'
    gMesh = fwm.create_group(gName)

    gMesh.create_dataset('points',data=triSurfaceList[0].rawPoints)
    gMesh.create_dataset('triangles',data=triSurfaceList[0].triangles)
    
    for i in range(len(triSurfaceList)):
        # group name
        gName = 'TriSurface'+str(i)
        gsurfi = fwm.create_group(gName)

        # save  data
        gsurfi.create_dataset('time',data=triSurfaceList[i].time)
        gsurfi.create_dataset('vars',data=triSurfaceList[i].rawVars)
        
    fwm.close()


def loadTriSurfaceList_hdf5(hdf5file,
                            viewAnchor,
                            xViewBasis,
                            yViewBasis,
                            srcBasisSrc=[[1,0,0],[0,1,0],[0,0,1]],
                            projectedField=True):
    '''
    '''
    # open the hdf5 parser
    fr = h5py.File(hdf5file, 'r')  
   
    # create the transformation objects
    afftrans, lintrans = getTransformation(viewAnchor=viewAnchor,
                                           xViewBasis=xViewBasis,
                                           yViewBasis=yViewBasis,
                                           srcBasisSrc=srcBasisSrc)
  
    # get mest data
    gName = 'mesh'
    points = fr[gName]['points'].value
    triangles = fr[gName]['triangles'].value
    
    ptsSrc = points
    ptsTgt = np.zeros((ptsSrc.shape[0],ptsSrc.shape[1]))
    for i in range(ptsSrc.shape[0]):
        ptsTgt[i,:] = afftrans.srcToTgt(ptsSrc[i,:])
    
    tsm = TriSurfaceMesh.TriSurfaceMesh(x=ptsTgt[:,0],
                                        y=ptsTgt[:,1],
                                        z=ptsTgt[:,2],
                                        triangles=triangles,
                                        mask=None,
                                        affTrans=afftrans,
                                        linTrans=lintrans)

    # create the TriSurface list
    triSurfaceList = []
    for i in range(len(fr.keys())-1):
        gName = 'TriSurface'+str(i)
        time = fr[gName]['time'].value
        data = fr[gName]['vars'].value 
        
        # get x and y vector (in ptsTgt)
        ptsSrc = points
        ptsTgt = np.zeros((ptsSrc.shape[0],ptsSrc.shape[1]))
        for i in range(ptsSrc.shape[0]):
            ptsTgt[i,:] = afftrans.srcToTgt(ptsSrc[i,:])

        #get vectors (in vecsTgt)
        vecsSrc = data
        vecsTgt = np.zeros((vecsSrc.shape[0],vecsSrc.shape[1]))
        if projectedField==True:
            for i in range(vecsSrc.shape[0]):
                vecsTgt[i,:] = lintrans.srcToTgt(vecsSrc[i,:])
        else:
            vecsTgt = vecsSrc
        
        tsv = TriSurfaceVector.TriSurfaceVector(vx=vecsTgt[:,0],
                                                vy=vecsTgt[:,1],
                                                vz=vecsTgt[:,2],
                                                time=time,
                                                triSurfaceMesh=tsm,
                                                projectedField=projectedField,
                                                interpolation=None,
                                                kind=None,)
        triSurfaceList.append(tsv)
    
    fr.close()    
    return triSurfaceList, tsm

    
def parseFoamFile_sampledSurface(foamFile):
    '''
    Parse a foamFile generated by the OpenFOAM sample tool or sampling library.
    
    Note:
        * It's a primitiv parser, do not add header in your foamFile!
        * Inline comment are allowed only from line start. c++ comment style.
        * It's REALLY a primitive parser!!!
        
    Arguments:
        * foamFile: [str()] Path of the foamFile

    Returns:
        * output: [numpy.array()] data store in foamFile
    '''
    output = []
    catchFirstNb = False
    istream = open(foamFile, 'r')
    for line in istream: 
        # This regex finds all numbers in a given string.
        # It can find floats and integers writen in normal mode (10000) or
        # with power of 10 (10e3).
        match = re.findall('[-+]?\d*\.?\d+e*[-+]?\d*', line)
        if (line.startswith('//')):
            pass
        if (catchFirstNb==False and len(match)==1):
            catchFirstNb = True
        elif (catchFirstNb==True and len(match)>0):
            matchfloat = list()
            for nb in match:                
                matchfloat.append(float(nb))
            if len(matchfloat)==1:
                output.append(matchfloat[0])
            else:
                output.append(matchfloat)
        else:
            pass
    istream.close()
    return np.array(output)

  
def parseVTK_ugly_sampledSurface(vtkfile):
    '''
    Parse a VTK file generate by the surface sampling tool of OpenFOAM. The
    surface has N grid points and M triangles. The data stored at each grid
    points has a dimension D.
    
    Warnings: This is a VERY primitive and ugly parser!! the python-to-vtk
    binding should be used instead of the following shitty code! 
    Nevertheless, this shit works :-)!!
    
    Arguments:
        *vtkfile*: python string
         Path to the vtk file
         
    Returns:
        *points*: numpy array of shape (N,3)
         List of points composing the grid.
         
        *polygon*: numpy array of shape (M,3)
         List of triangles. Technically, this parser can return a List of any
         type of ploygon, e.g: triangle, square, pentagon...
 
        *pointData* numpy array of shape (N,D)
         List of data associate with each point of the grid.
    '''
    pointsOut = []
    polyOut = []
    pointDataOut = []
    
    istream = open(vtkfile, 'r')
    line = istream.readline()

    # catch the begin of the list of points
    # -------------------------------------
    catchLine = False
    while catchLine==False:
        if (line.startswith('DATASET POLYDATA')):
            catchLine = True
        line = istream.readline()
        
    # catch the number of points
    nbpoints = int(re.findall('[-+]?\d*\.?\d+e*[-+]?\d*', line)[0])
    pti = 0
    line = istream.readline()
    
    #store the points in pointsOut
    while (pti<nbpoints):
        match = re.findall('[-+]?\d*\.?\d+e*[-+]?\d*', line)
        pointsOut.append(match)
        line = istream.readline()
        pti = pti+1
    pointsOut = np.asarray(pointsOut,dtype=float)
    
    # catch the begin of the list of polygons and the number of polygon
    # -----------------------------------------------------------------
    catchLine = False
    nbpoly = 0
    while catchLine==False:
        if (line.startswith('POLYGONS')):
            catchLine = True
            nbpoly = int(re.findall('[-+]?\d*\.?\d+e*[-+]?\d*', line)[0])
        line = istream.readline()
    polyi = 0
    
    #store the polygons in polyOut
    while polyi<nbpoly:
        match = re.findall('[-+]?\d*\.?\d+e*[-+]?\d*', line)
        polyOut.append(match[1:])
        line = istream.readline()
        polyi = polyi+1
    polyOut = np.asarray(polyOut,dtype=float)
    
    # catch the begin of the list of point data
    # -----------------------------------------
    catchLine = False
    nbptdata = 0
    while catchLine==False:
        if (line.startswith('POINT_DATA')):
            catchLine = True
            nbptdata = int(re.findall('[-+]?\d*\.?\d+e*[-+]?\d*', line)[0])
        line = istream.readline()
    ptdatai = 0
    
    # jump the line starting with "FIELD attributes"
    line = istream.readline()
    
    # catch the dimension of point data and the number of point data
    match = re.findall('[-+]?\d*\.?\d+e*[-+]?\d*', line)
    dimptdata = int(match[0])
    line = istream.readline()
    
    #store the point data in pointDataOut
    if dimptdata==1:
        pointDataOut = re.findall('[-+]?\d*\.?\d+e*[-+]?\d*', line)
    else: 
        while ptdatai<nbptdata:
            match = re.findall('[-+]?\d*\.?\d+e*[-+]?\d*', line)
            pointDataOut.append(match)
            line = istream.readline()
            ptdatai = ptdatai+1
    pointDataOut = np.asarray(pointDataOut,dtype=float)
    
    istream.close()
    
    return pointsOut, polyOut, pointDataOut
