'''
TriSurfaceVector.py
'''

#import re

import numpy as np
import matplotlib.tri as tri

import pyFlowStat.TriSurface as TriSurface
import pyFlowStat.ParserFunctions as ParserFunctions


class TriSurfaceScalar(TriSurface.TriSurface):
    '''
    class TriSurfaceScalar.
    
    Class which describes a scalar field s on a 2D (flat)
    triangulated surface S of N points. The triangulation is a TriSurfaceMesh
    object of M triangles.
    
    The coordinate system of S is the following:
        * x = the horizontal coordinate, pointing East (in-plane)
        * y = the vertical coordinate, pointing North (in-plane)
        * z = the off-plane coordinate, pointing torward you     
    '''
    
    # constructors #
    #--------------#
    
    def __init__(self,
                 s,
                 time,
                 triSurfaceMesh,
                 projectedField=False,
                 interpolation=None,
                 kind=None):
        '''
        base constructor.
        
        Arguments:          
            *s*: numpy array of shape (npoints).
             scalar values.
         
            *time*: python float
             timestep of the surface. If this information does not matter,
             use 0.
             
            *triSurfaceMesh* :  TriSurfaceMesh object.
             TriSurfaceMesh object, which holds the mesh information.
             
            *projectedField* python bool (default=False)
             Unused for the scalar field, but might be used for the fields
             added with the methods "addField", "addFieldFromFoamFile" and
             "addFieldFromVtk"
             
            *interpoation*: python string. 
             type of interpolation used. Value: "cubic" or "linear".
             Default="cubic". 

            *kind*: python string.
             Definition of the cubic interpolation type. Value: "geom" or
             "min_E". Default="geom". "min_E" is supposed to be the more 
             accurate, but "geom" is way faster.
        '''
        super(TriSurfaceScalar,self).__init__(time=time,
                                        triSurfaceMesh=triSurfaceMesh,
                                        projectedField=projectedField,
                                        interpolation=interpolation,
                                        kind=kind)

        self.s = np.asarray(s)     
        self.s_i = None

        
    @classmethod
    def readFromFoamFile(cls,
                         varsFile,
                         triSurfaceMesh,
                         time,
                         projectedField=False):
        '''
        Construct from a surface saved  by OpenFOAM in foamFile format.
        
        Arguments:
            *varsFile*: python string.
             Path to the file holding the scalar field.
             
            *time*: python float
             timestep of the surface. If this information does not matter,
             use 0.
             
            *triSurfaceMesh* :  TriSurfaceMesh object.
             TriSurfaceMesh object, which holds the mesh information.
             
            *projectedField* python bool (default=False)
             Unused for the scalar field, but might be used for the fields
             added with the methods "addField", "addFieldFromFoamFile" and
             "addFieldFromVtk"
        '''

        #get scalars
        slrsTgt = ParserFunctions.parseFoamFile_sampledSurface(varsFile)

        # update class member variables
        return cls(s=slrsTgt,
                   time=time,
                   triSurfaceMesh=triSurfaceMesh,
                   projectedField=projectedField,
                   interpolation=None,
                   kind=None)
 
 
    @classmethod 
    def readFromHdf5(cls,
                     hdf5Parser,
                     varName,
                     triSurfaceMesh,
                     time,
                     projectedField=False):
        '''
        '''
        gTime = str(time)
        time = hdf5Parser[gTime]['time'].value
        slrsTgt = hdf5Parser[gTime][varName].value
        
        # update class member variables
        return cls(s=slrsTgt,
                   time=time,
                   triSurfaceMesh=triSurfaceMesh,
                   projectedField=projectedField,
                   interpolation=None,
                   kind=None)
    
    
    @classmethod   
    def readFromVTK(cls,
                    vtkFile,
                    triSurfaceMesh,
                    time,
                    projectedField=False):
        '''
        Construct from a surface saved by OpenFOAM in VTK format.
        
        Arguments:
            *varsFile*: python string.
             Path to the vtk-file.
             
            *time*: python float
             timestep of the surface. If this information does not matter,
             use 0.
             
            *triSurfaceMesh* :  TriSurfaceMesh object.
             TriSurfaceMesh object, which holds the mesh information.
             
            *projectedField* python bool (default=False)
             Defines if the data fields has to be projected in the basis of the
             surface. 
        '''     
        # read VTK file
        ptsSrc, triangles, slrsTgt = ParserFunctions.parseVTK_ugly_sampledSurface(vtkFile)


        # update class member variables
        return cls(s=slrsTgt,
                   time=time,
                   triSurfaceMesh=triSurfaceMesh,
                   projectedField=projectedField,
                   interpolation=None,
                   kind=None)
  
    # getters #
    #---------#      


    # setters #
    #---------#


    # class methods #
    #---------------#
    def __call__(self,dim=0):
        return self.component(dim=dim)
        
    def component(self,dim=0):
        return self.s
        
    def interpolate(self,x,y,dim=0):
        try:
            return self.s_i(x,y)
        except:
            raise ValueError('this method needs interpolators. Please run',
                                 'method "addInterpolator" first.')
        
    def rawVars(self):
        '''
        Return the scalar field defined in the source coordinate system.
        
        Returns:
            *rawData*: numpy array of shape (N,)
        '''
        return self.s
        

    def surfaceVars(self):
        '''
        Return the scalar field as saved in the TriSurfaceScalar object.
        
        Returns:
            *surfaceVars*: numpy array of shape (N,)
        '''
        return self.s


    def gradientxy(self,x,y):
        '''
        Calculate the gradient at the point pt(x,y) and return it. The method
        works also for a list of N points. In such case, x and y are arrays.
        
        Arguments:
            *x*: python float or numpy array (shape=N).
             x coordinate of the point pt.
            
            *y*: python float or numpy array (shape=N).
             y coordinate of the point pt.
             
        Returns:
            *dsdx, dsdy*: python tuple of two float or numpy array.
             The gradient at point pt(x,y).
        '''
        if self.s_i!=None: 
            dsdx, dsdy = self.s_i.gradient(x,y)
            return dsdx, dsdy
        else:
            raise ValueError('this method needs interpolators. Please run',
                             'method "addInterpolator" first.')

  
    # class methods - adders #
    #------------------------#
        
    def addInterpolator(self,interpolation='cubic', kind='geom'):
        '''
        Add interpolator Object.
        '''
        self.interType = interpolation
        self.interKind = kind
        if self.interType=='cubic':
            self.s_i = tri.CubicTriInterpolator(self.triangulation, self.s, kind=self.interKind)
        elif self.interType=='linear':
            self.s_i = tri.LinearTriInterpolator(self.triangulation, self.s)
        else:
            raise ValueError('Interpolation must be "cubic" or "linear".')
            

    def addGradient(self):
        '''
        Calculate and save the gradient at all point of the grid.
        '''   
        dsdx, dsdy = self.gradientxy(self.x,self.y)
        self.data['dsdx'] = dsdx
        self.data['dsdy'] = dsdy