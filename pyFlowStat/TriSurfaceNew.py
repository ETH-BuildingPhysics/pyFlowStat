'''
TriSurfaceNew.py

New version of TriSurface. TriSurfaceNew will use matplotlib.tri classes
extensively.

!!! Still in alpha version !!!

See domcumention in class definition
'''

import re

import numpy as np
import matplotlib.tri as tri




class TriSurfaceNew(tri.Triangulation):
    '''
    class TriSurfaceNew
    '''
    
    def __init__(self, x, y, z, triangles=None, mask=None):
        '''
        base constructor from a list of x, y and z. list of triangles and mask 
        optional.
        '''
        tri.Triangulation.__init__(self, x, y, triangles=None, mask=None)
        
        self.z = z
        self.data = dict()
        
    @classmethod
    def readFromFoamFile(cls,
                         varsFile,
                         pointsFile,
                         facesFile,
                         viewAnchor,
                         xViewBasis,
                         yViewBasis,
                         srcBasisSrc=[[1,0,0],[0,1,0],[0,0,1]],
                         mask=None):
        '''
        Construct from a surface saved  by OpenFOAM in foamFile format.
        '''
        # check and convert arguments
        srcBasisSrc = np.array(srcBasisSrc,dtype=float)
        if srcBasisSrc.shape!=(3,3):
            raise ValueError('srcBasis must be a 3x3 matrix')
            
        xViewBasis = np.array(xViewBasis,dtype=float)
        yViewBasis = np.array(yViewBasis,dtype=float)
        if xViewBasis.shape!=(3,) or yViewBasis.shape!=(3,):
            raise ValueError('xViewBasis.shape and yViewBasis.shape must be equal to (3,)')
            
        # get the basis and the transformation object
        tgtBasisSrc = np.zeros((3,3))
        tgtBasisSrc[:,0] = xViewBasis        
        tgtBasisSrc[:,1] = yViewBasis
        tgtBasisSrc[:,2] = np.cross(xViewBasis,yViewBasis)
        afftrans = affineTransfomation(srcBasisSrc,tgtBasisSrc,viewAnchor)
        lintrans = linearTransformation(srcBasisSrc,tgtBasisSrc)

        # get x and y vector (in ptTrt)
        ptsSrc = parseFoamFile(pointsFile)
        ptsTgt = np.zeros((ptsSrc.shape[0],ptsSrc.shape[1]))
        for i in range(ptsSrc.shape[0]):
            ptsTgt[i,:] = afftrans.srcToTgt(ptsSrc[i,:])
            
        #get triangles
        triangles = parseFoamFile(facesFile)[:,1:4]

        # feed x, y and triangles to the base constructor
        return cls(ptsTgt[:,0],ptsTgt[:,1],parseFoamFile(varsFile),triangles,mask)
 
       
    @classmethod
    def readFromVTK(cls):
        '''
        Construct from a surface saved by OpenFOAM in VTK format.
        '''
        raise NotImplementedError('The method is still not implemented')
    
        
#    @classmethod
#    def test(cls, x, y, z, triangles):
#        ts = cls(x, y, z, triangles, mask=None)
#
#        return ts

        
    
    
    
    
class affineTransfomation(object):
    '''
    An affine transfomation is define as ya = A*xa, with A the affine matrix
    (shape=(4,4)), xa (shape=(4)) the augmented vector defined in the source
    basis and ya (shape=(4)) the augmented vector defined in the target basis.
    
    The transformation in the other direction is defined as xa = invA*ya, with
    invA (shape=(4,4)) the inverse of A.
    
    Arguments:
        * srcBasisSrc: [numpy.array, shape=(3,3)] The source coordinate system
          e1, e2 and e3 packed in matrix [e1,e2,e3].
        * tgtBasisSrc: [numpy.array, shape=(3,3)] The target coordinate system
          f1, f2 and f3, defined in the source basis and packed in matrix 
          [f1,f2,f3.]
        * 
    '''
    
    def __init__(self,srcBasisSrc,tgtBasisSrc,tgtTransSrc):  
        self.A = self.affineMat(tgtBasisSrc,tgtTransSrc)
        self.invA = np.linalg.inv(self.A)

    def affineVec(self,vec):
        '''
        Return affine vector (shape=4) from standard vector (shape=3).
        
        Arguments:
            * vec: [np.array, shape=)] a vector
        
        Return:
            * affineVec [np.array, shape=4] the affine vector. Same as vec, 
              but with a trailing 1.
        '''
        affineVec = np.zeros(len(vec)+1)
        affineVec[-1] = 1
        affineVec[0:len(vec)] = affineVec[0:len(vec)]+vec
        return affineVec
    
    def affineMat(self,mat,trans):
        '''
        Return affine matrix (shape=(4,4)) from  a standard matrix 
        (shape=(3,3)) + a translation vector (shape=(3)).
        
        Arguments:
            * mat: [np.array, shape=(3,3)] a matrix.
            * trans: [np.array, shape=(3)] translation vector from srcBasisTgt 
              to tgtBasisSrc.
        
        Return:
            * affineMat [np.array. affineMat.shape=(4,4)] the affine Matrix.
        '''
        nbl = mat.shape[0]  #number of line
        nbr = mat.shape[0]  #number of row
        affineMat = np.zeros((nbl+1,nbr+1))
        affineMat[0:nbl,0:nbl] = affineMat[0:nbl,0:nbl]+mat
        affineMat[0:nbr,-1] = affineMat[0:nbr,-1] + trans
        affineMat[-1,-1] = 1
        return affineMat
        
    def tgtToSrc(self,vec):
        '''
        From a vector defined in the target basis, return the same vector
        defined in the source basis.
        '''
        return np.dot(self.A,self.affineVec(vec))[0:3]
        
    
    def srcToTgt(self,vec):
        '''
        From a vector defined in the source basis, return the same vector
        defined in the target basis.
        '''
        return np.dot(self.invA,self.affineVec(vec))[0:3]


class linearTransformation(object):
    '''
    A class which defines a linear transformation
    '''    
    def __init__(self,srcBasisSrc,tgtBasisSrc):  
        self.A = tgtBasisSrc
        self.invA = np.linalg.inv(self.A)
        
    def tgtToSrc(self,vec):
        return np.dot(self.A,vec)
         
    def srcToTgt(self,vec):
        return np.dot(self.invA,vec)
        
            


def parseFoamFile(foamFile):
    '''
    Parse a foamFile generated by the OpenFOAM sample tool or sampling library.
    
    Note:
        * It's a primitiv parse, do not add header in your foamFile!
        * Inline comment are allowed only from line start. c++ comment style.
        
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
            output.append(matchfloat)
        else:
            pass
    istream.close()
    return np.array(output)  