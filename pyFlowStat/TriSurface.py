'''
TriSurface.py
'''


import numpy as np
#import matplotlib.tri as tri

import TriSurfaceFunctions


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
        return TriSurfaceFunctions.mat(self.data[varName])
        
        
    def unmat(self,varName):
        '''
        '''
        return TriSurfaceFunctions.unmat(self.data[varName])
    
    
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
        fieldSrc = TriSurfaceFunctions.parseFoamFile_sampledSurface(fieldFile)
        self.addField(fieldSrc,fieldname)
        
        
    def addFieldFromVTK(self,fieldFile,fieldname):
        '''
#        Add a field F (shape d) stored in a VTK file to the current
#        TriSurfaceVector object. See docstring from self.addField() for more
#        information.
        '''
        #get field
        points, polygon, fieldSrc = TriSurfaceFunctions.parseVTK_ugly_sampledSurface(fieldFile)
        self.addField(fieldSrc,fieldname)
    
    
    def addGradient(self):
        '''
        Calculate and save the gradient at all point of the grid. As expected,
        the dvidz does not exist.
        '''   
        raise NotImplementedError('TriSurface subclasses should implement addGradient.')
        
        
    
            
    
    
 
