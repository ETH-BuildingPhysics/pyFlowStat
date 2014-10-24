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

import numpy as np

import CoordinateTransformation as coorTrans
import TriSurfaceVector as TriSurfaceVector



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
        gsurfi.create_dataset('data',data=triSurfaceList[i].rawData)
        
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

    # create the TriSurface list
    triSurfaceList = []
    for i in range(len(fr.keys())):
        gName = 'TriSurface'+str(i)
        time = fr[gName]['time'].value
        data = fr[gName]['data'].value 
        
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
        
        tsv = TriSurfaceVector.TriSurfaceVector(x=ptsTgt[:,0],
                                                y=ptsTgt[:,1],
                                                z=ptsTgt[:,2],
                                                vx=vecsTgt[:,0],
                                                vy=vecsTgt[:,1],
                                                vz=vecsTgt[:,2],
                                                time=time,
                                                triangles=triangles,
                                                mask=None,
                                                projectedField=projectedField,
                                                interpolation=None,
                                                kind=None,
                                                affTrans=afftrans,
                                                linTrans=lintrans)
        triSurfaceList.append(tsv)
    
    fr.close()    
    return triSurfaceList