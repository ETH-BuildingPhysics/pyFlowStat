'''
TriSurface.py
'''


import numpy as np
#import matplotlib.tri as tri

import pyFlowStat.CoordinateTransformation as coorTrans
import pyFlowStat.ParserFunctions as ParserFunctions


class TriSurface(object):
    '''
 
    '''
    
    # constructors #
    #--------------#
    
    def __init__(self,
                 time,
                 triSurfaceMesh,
                 projectedField=True,
                 interpolation=None,
                 kind=None):
        '''
        base constructor.
        '''
        self.triSurfaceMesh = triSurfaceMesh
        
        self.data = dict()
        self.data_i = dict()
        
        self.time = float(time)

        # "private" member variable. Don't play with them if you are not sure...        
        self.interType = interpolation
        self.interKind  = kind
        self.projectedField = projectedField
        
    @classmethod
    def readFromFoamFile(cls,
                         varsFile,
                         triSurfaceMesh,
                         time,
                         projectedField=True):
        '''    
        '''
        raise NotImplementedError('TriSurface subclasses should implement readFromFoamFile.')

 
    @classmethod   
    def readFromVTK(cls,
                    vtkFile,
                    triSurfaceMesh,
                    time,
                    projectedField=True):
        '''
        '''     
        raise NotImplementedError('TriSurface subclasses should implement readFromVTK.')

  
    # getters #
    #---------#
    @property
    def x(self):
        return self.triSurfaceMesh.x
        
    @property    
    def y(self):
        return self.triSurfaceMesh.y
        
    @property
    def triangulation(self):
        return self.triSurfaceMesh.triangulation
        
    @property
    def triangles(self):
        return self.triSurfaceMesh.triangles
        
    @property
    def affTrans(self):
        return self.triSurfaceMesh.affTrans
    
    @property
    def linTrans(self):
        return self.triSurfaceMesh.linTrans


    # setters #
    #---------#



    # class methods #
    #---------------#
    def rawPoints(self):
        '''
        Return the grid points in the source coordinate system.
        
        Returns:
            *rawPoints*: numpy array of shape (N,3)
        '''
        return self.triSurfaceMesh.rawPoints()
            
     
    def rawVars(self):
        '''
        '''
        raise NotImplementedError('TriSurface subclasses should implement rawVars.')
        

    def surfaceVars(self):
        '''
        '''
        raise NotImplementedError('TriSurface subclasses should implement surfaceVars.')




    def gradientxy(self,x,y):
        '''
        '''
        raise NotImplementedError('TriSurface subclasses should implement gradientxy.')


    def mat(self,varName):
        '''
        '''
        return mat(self.data[varName])
        
        
    def unmat(self,varName):
        '''
        '''
        return unmat(self.data[varName])
    
    
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
        
    # class methods - adders #
    #------------------------#
        
    def addInterpolator(self,interpolation='cubic', kind='geom'):
        '''
        '''
        raise NotImplementedError('TriSurface subclasses should implement addInterpolator.')
        
            
    
    def addField(self,field,fieldname):
        '''
#        Add a field F of dimension d (e.g: d=3 for a vector filed) to the
#        current TriSurfaceVector object TSV. The grid of F (N points) must be
#        identical, in term of number of points and their location, as the grid
#        of TSV. F will be stored in TSV.data['fieldName'] or TSV['fieldName'].
#        
#        Arguments:
#            *field*: numpy array of shape (N,d).
#            
#            *fieldName*: python string.
        '''
        fieldSrc = field
        fieldShape = fieldSrc.shape
        fieldTgt = np.zeros(fieldShape)

        if (self.projectedField==True and len(fieldShape)>1):            
            if fieldShape[1]==3:
                for i in range(fieldShape[0]):
                    fieldTgt[i,:] = self.linTrans.srcToTgt(fieldSrc[i,:])
            if fieldShape[1]==6:
                for i in range(fieldShape[0]):
                    a = fieldSrc[i,:]
                    A = np.array([[a[0],a[1],a[2]],
                                  [a[1],a[3],a[4]],
                                  [a[2],a[4],a[5]]])
                    B = self.linTrans.srcToTgt(A)
                    fieldTgt[i,:] = np.array([B[0,0],B[0,1],B[0,2],B[1,1],B[1,2],B[2,2]])
                    
                
        else:
            fieldTgt = fieldSrc
        self.data[fieldname] = fieldTgt
            
        
    def addFieldFromFoamFile(self,fieldFile,fieldname):
        '''
#        Add a field F (shape d) stored in a foamFile to the current
#        TriSurfaceVector object. See docstring from self.addField() for more
#        information.
        '''
        #get field
        fieldSrc = ParserFunctions.parseFoamFile_sampledSurface(fieldFile)
        self.addField(fieldSrc,fieldname)
        
        
    def addFieldFromVTK(self,fieldFile,fieldname):
        '''
#        Add a field F (shape d) stored in a VTK file to the current
#        TriSurfaceVector object. See docstring from self.addField() for more
#        information.
        '''
        #get field
        points, polygon, fieldSrc = ParserFunctions.parseVTK_ugly_sampledSurface(fieldFile)
        self.addField(fieldSrc,fieldname)
    
    
    def addGradient(self):
        '''
        Calculate and save the gradient at all point of the grid. As expected,
        the dvidz does not exist.
        '''   
        raise NotImplementedError('TriSurface subclasses should implement addGradient.')
        
        
def getTransformation(viewAnchor,
                      xViewBasis,
                      yViewBasis,
                      srcBasisSrc=[[1,0,0],[0,1,0],[0,0,1]]):
    '''
    Return the affine and linear transfomation object for coordinate
    transformation from a source coordinate system S to a target coordinate 
    system T.
    
    Arguments:
        *viewAnchor*: python list or numpy array of shape (3,).
         Location of the origin of T, defined in S.
         
        *xViewBasis*: python list or numpy array of shape (3,).
         x axis of T, defined in S.
        
        *yViewBasis*: python list or numpy array of shape (3,).
         y axis of T, defined in S.
        
        *srcBasisSrc*: python list or numpy array of shape (3,3).
         x,y,z axis of S, defined in S. For advenced user only.
         Default=[[1,0,0],[0,1,0],[0,0,1]]
          
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
    afftrans = coorTrans.AffineTransformation(srcBasisSrc,tgtBasisSrc,viewAnchor)
    lintrans = coorTrans.LinearTransformation(srcBasisSrc,tgtBasisSrc)
    
    return afftrans, lintrans

    
def mat(field):
        '''
        return a symmTensor field (shape=[N,6]) or Tensor field (shape=[N,9])
        as a "real" tensor field (shape=[N,3,3]). if "field" is a vector or a
        scalar, nothing is done and return as given.
        
        In the TriSurface ecosystem, a symmTenor field is stored as a line to
        save memory. For consistency, a Tensor field is also stored as a line.
        Nevertheless this memory efficient storing is incompatible with linear
        algebra calculus. Therefore the function mat() converts a tensor like
        field into "real" Tensor field.
        
        Arguments:
            *field*: numpy array. Shape=[N,d], with d an int.
             Field to convert.
             
        Returns:
            *realField* numpy array. shape=[N,3,3]
             Return the field as a "real" tensor. If "field" is not tensor
             like, the field is returned unmodified.
        '''
        tensorType = len(field[0])
        
        if tensorType==6:
            tgt = np.zeros([field.shape[0],3,3])   # target field
            # fill line by line
            tgt[:,0,:] = field[:,0:3]  # add index 11,12,13 to tgt
            tgt[:,1,:] = np.array([field[:,1],field[:,3],field[:,4]]).T  # add index 21,22,23 to tgt
            tgt[:,2,:] = np.array([field[:,2],field[:,4],field[:,5]]).T # add index 31,32,33 to tgt
            return tgt
        elif tensorType==9:
            tgt = np.zeros([field.shape[0],3,3])   # target field
            # fill line by line
            tgt[:,0,:] = field[:,0:3]  # add index 11,12,13 to tgt
            tgt[:,1,:] = field[:,3:6]  # add index 21,22,23 to tgt
            tgt[:,2,:] = field[:,6:9] # add index 31,32,33 to tgt
            return tgt
        else:
            return field
   
         
def unmat(field):
    '''
    Opposite of mat()....
    '''
    # get field type
    fieldType = None
    if (len(field.shape)==3 and field.shape[1]==3 and field.shape[2]==3):  #field is a 3x3 tensor
        if [field[0,0,1],field[0,0,2],field[0,1,2]]==[field[0,1,0],field[0,2,0],field[0,2,1]]: #field is symmTensor
            fieldType = 'symmTensor'
        else:
           fieldType = 'tensor'
    else:
       fieldType = 'notTensor'
           
    # convert field and return it
    if fieldType=='symmTensor':
        return np.hstack(( field[:,0,:] , field[:,1,1:3] , field[:,2,2].reshape(field.shape[0],1) ))
    elif fieldType=='tensor':
        return np.hstack(( field[:,0,:] , field[:,1,:] , field[:,2,:]))
    elif fieldType=='notTensor':
        return field            
    
    
 
