'''
CoordinateTransfomation.py

Collection of tools/classes to do coordinate transformation.
'''

import numpy as np


class AffineTransformation(object):
    '''
    An affine transfomation is define as ya = A*xa, with A the affine matrix
    (shape=(4,4)), xa (shape=(4)) the augmented vector defined in the source
    basis and ya (shape=(4)) the augmented vector defined in the target basis.
    
    The transformation in the other direction is defined as xa = invA*ya, with
    invA (shape=(4,4)) the inverse of A.
    
    Arguments:
        *srcBasisSrc*: numpy array of shape=(3,3)
         The source coordinate system e1, e2 and e3 packed in matrix [e1,e2,e3].
        *tgtBasisSrc*: numpy array of shape=(3,3)
         The target coordinate system f1, f2 and f3, defined in the source
         basis and packed in matrix [f1,f2,f3]
        *tgtTransSrc*: numpy array of shape=(3,)
         Translation vector of the target basis regarding the source basi
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


class LinearTransformation(object):
    '''
    A class which defines a linear transformation y = A*x, with x a vector
    defined in the source basis, y the same vector x but defined in the target
    basis and A the linear tranformation matrix to go from the source to the
    taget
    '''    
    def __init__(self,srcBasisSrc,tgtBasisSrc):  
        self.A = tgtBasisSrc
        self.invA = np.linalg.inv(self.A)
        
    def tgtToSrc(self,vec):
        return np.dot(self.A,vec)
         
    def srcToTgt(self,vec):
        return np.dot(self.invA,vec)
        