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
import numpy as np

import CoordinateTransformation as coorTrans



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