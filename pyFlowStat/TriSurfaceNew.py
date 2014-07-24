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

import pyFlowStat.triZinterpolator as triz




class TriSurfaceNew(tri.Triangulation):
    '''
    class TriSurfaceNew. New implentation of the class trisurface, which
    derived from the matplotlib class matplotlib.tri.Triangulation.
    
    variables:
        *x*: array of shape (npoints).
         x-coordinates of grid points. From the parent class.

        *y*: array of shape (npoints).
         y-coordinates of grid points. From the parent class.

        *triangles*: integer array of shape (ntri,3).
         For each triangle, the indices of the three points that make
         up the triangle, ordered in an anticlockwise manner. From the parent 
         class.
         
        *mask*: optional boolean array of shape (ntri).
         Which triangles are masked out. From the parent class.

        *edges*: integer array of shape (?,2).
         All edges of non-masked triangles.  Each edge is the start
         point index and end point index.  Each edge (start,end and
         end,start) appears only once. From the parent class.

        *neighbors*: integer array of shape (ntri,3).
         For each triangle, the indices of the three triangles that
         share the same edges, or -1 if there is no such neighboring
         triangle.  neighbors[i,j] is the triangle that is the neighbor
         to the edge from point index triangles[i,j] to point index
         triangles[i,(j+1)%3]. From the parent class.
         
        *z*: array of shape (npoints,dim).
         value at the coordinate (x,y).
             * dim=1 for scalar
             * dim=3 for vector
             * etc...
         
        *interpolator*: list of length (dim).
         list of object from class matplotlib.tri.CubicTriInterpolator. One 
         obejct triCubicInterpolator per element of z[i]. Examples: 3 objects
         for a vector...
         
        *data*: dictionnary.
         Storage for all kind of stuff.
        
         
     Methods:
         *__init__*: base constructor
          Base constructor of the class. The constructor "readFrom*" should
          be prefered to create a TriSurfaceNew object
     
         *readFromFoamFile*: constructor.
          Construct from a surface saved  by OpenFOAM in foamFile format.
          
         *readFromVTK*: constructor.
          Construct from a surface saved by OpenFOAM in VTK format.
         
         *rawGrad*:
          Computer the gradient at the all the x,y location of the triangulation.
          In opposition to the methode gradient, no interpolation is done.
          
          
         *gradient*:
          Compute the gradient of z at location (x,y). x and y can be arrays.
          
         *interpolate*:
          Compute z at location (x,y). x and y can be arrays.

    '''
    
    def __init__(self, x, y, z, triangles=None, mask=None):
        '''
        base constructor from a list of x, y and z. list of triangles and mask 
        optional.
        
        Arguments:
            *x*: array of shape (npoints).
             x-coordinates of grid points.
    
            *y*: array of shape (npoints).
             y-coordinates of grid points.
    
            *triangles*: integer array of shape (ntri,3).
             For each triangle, the indices of the three points that make
             up the triangle, ordered in an anticlockwise manner. If no
             triangles is passed, a Delauney triangulation is computed. 
             Default=None.
             
            *mask*: optional boolean array of shape (ntri).
             Which triangles are masked out.
             
        '''
        if mask==None:
            tri.Triangulation.__init__(self, x, y, triangles=triangles, mask=None)
        else:
            triang = tri.Triangulation(x, y, triangles=triangles, mask=mask)
            trianalyzer = tri.TriAnalyzer(triang)
            (comp_triangles,
             comp_x,
             comp_y,
             tri_renum,
             node_renum) = trianalyzer._get_compressed_triangulation(True, True)
            node_mask = (node_renum == -1)
            z[node_renum[~node_mask]] = z
            z = z[~node_mask]
            tri.Triangulation.__init__(self, comp_x, comp_y, triangles=comp_triangles, mask=None)
            

        
        self.z = z
        self.data = dict()
        
        # interpolator. An oject of class ArrayTriInterpolator. 
        # Created if needed
        self.interpolator = None
        
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
            raise ValueError('xViewBasis.shape and yViewBasis. ',
                             'shape must be equal to (3,)')
            
        # get the basis and the transformation object
        tgtBasisSrc = np.zeros((3,3))
        tgtBasisSrc[:,0] = xViewBasis        
        tgtBasisSrc[:,1] = yViewBasis
        tgtBasisSrc[:,2] = np.cross(xViewBasis,yViewBasis)
        afftrans = AffineTransfomation(srcBasisSrc,tgtBasisSrc,viewAnchor)
        #lintrans = LinearTransformation(srcBasisSrc,tgtBasisSrc)

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
        raise NotImplementedError('The method is not implemented')
  
    def rawGrad(self):
        '''
        Calculate and save the gradient at all (self.x,self.y) locations.
        '''
        self.create_interpolator()
        self.data['grad'] = np.array([self.interpolator[i].gradient(self.x,self.y) for i in range(len(self.z[0,:]))]).T
        
        
    def gradient(self,x,y):
        '''
        Return gradient at location (x,y). x,y can be arrays
        '''
        self.create_interpolator()
        return np.array([self.interpolator[i].gradient(x,y) for i in range(len(self.z[0,:]))]).T
        
        
    def interpolate(self,x,y):
        '''
        Return interpolated value at location (x,y). x and y can be arrays. The
        member variable interpolation can also be used directly.
        '''
        self.create_interpolator()
        return np.array([self.interpolator[i](x,y) for i in range(len(self.z[0,:]))]).T

        

    def create_interpolator(self):
        '''
        Create the list of interpolator object.
        '''
        if self.interpolator==None:
            self.interpolator  = list()
            for i in range(len(self.z[0,:])):
                self.interpolator.append(triz.CubicTriZInterpolator(self, self.z[:,i]))
            
        
    def __iter__(self):
        '''
        Iterable on member "data".
        '''
        return self.data.itervalues()

    def __getitem__(self, key):
        '''
        Getter for key "key" on member dictionnary "data"
        '''
        return self.data[key]

    def __setitem__(self, key, item):
        '''
        Setter for key "key" on member dictionnary "data"
        '''
        self.data[key] = item
           
    # New implemtations of some parent methods
    #----------------------------------------#    
    def set_mask(self, mask):
        raise NotImplementedError('set_mask needs a new implementation!')

        
    
class AffineTransfomation(object):
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
        
            


def parseFoamFile(foamFile):
    '''
    Parse a foamFile generated by the OpenFOAM sample tool or sampling library.
    
    Note:
        * It's a primitiv parse, do not add header in your foamFile!
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
            output.append(matchfloat)
        else:
            pass
    istream.close()
    return np.array(output)  