'''
TriSurfaceVector.py
'''

import re

import numpy as np
import matplotlib.tri as tri

import CoordinateTransformation as coorTrans
import TriSurfaceFunctions


class TriSurfaceVector(object):
    '''
    class TriSurfaceVector.
    
    Class which describes a 3D vector field V=(vx,vy,vz) on a 2D (flat)
    triangulated surface S of N points. The triangulation is a grid of M 
    triangles. 
    The coordinate system of S is the following:
        * x = the horizontal coordinate, pointing East (in-plane)
        * y = the vertical coordinate, pointing North (in-plane)
        * z = the off-plane coordinate, pointing torward you   
    Therefore, vx and vy are the in-plane vector coordinate and vz the
    off-plane one.
    
    To construct a TriSurfaceVector object, use better the @classmethod
    "readFromFoamFile" or "readFromVTK".

    Member variables:
        *triangulation* : triangulation object. 
         Triangulation object, which holds the grid. Therefore, it also holds
         the grid points coordinates (x and y), the triangles and some other
         stuffs. See matplotlib.Tri.Trangulation.Triangulation for more
         details.
         
        *vx*: numpy array of shape (N,)
         
        *vy*: numpy array of shape (N,)
        
        *vz*: numpy array of shape (N,)
        
        *vx_i*: CubicTriInterpolator or LinearTriInterpolation object
        
        *vy_i*: CubicTriInterpolator or LinearTriInterpolation object
        
        *vz_i*: CubicTriInterpolator or LinearTriInterpolation object
        
        *data*: python dict
        
        *data_i*: python dict
        
        *__interType*: python string
        
        *__interKind*: python string
        
        *__affTrans*:
        
        *__linTrans*:
        
        *__projectedField*: python bool
        
    Member functions:
        *__init__*: base constructor.
        
        *readFromFoamFile*: constructor from a foamFile generate by OpenFOAM.
        
        *readFromVTK*: contructor from a VTK file generate by OpenFOAM.
        
        *x*: Returns the x coordinate of the grid points.
        
        *y*: Returns the y coordinate of the grid points.
        
        *trangles*
        
        *getInterpolator*
        
        *addFields*

        *addFieldFromFoamFile*        
        
        *gradient*
        
        *gradient_i*
        
        
    '''
    
    # constructors #
    #--------------#
    
    def __init__(self,
                 x,
                 y,
                 z,
                 vx,
                 vy,
                 vz,
                 time,
                 triangles=None,
                 mask=None,
                 projectedField=True,
                 interpolation=None,
                 kind=None,
                 affTrans=None,
                 linTrans=None):
        '''
        base constructor.
        
        Arguments:
            *x*: numpy array of shape (npoints).
             x-coordinates of grid points.
    
            *y*: numpy array of shape (npoints).
             y-coordinates of grid points.
             
            *vx*: numpy array of shape (npoints).
             "in-plan" vector componant. Labeled x TriSurfaceVector
             
            *vy*: numpy array of shape (npoints).
             "in-plan" vector componant. Labeled y in TriSurfaceVector
             
             *vz*: numpy array of shape (npoints).
             "off-plan" vector componant. Labeled z in TriSurfaceVector

            *triangles*: integer array of shape (ntri,3).
             For each triangle, the indices of the three points that make
             up the triangle, ordered in an anticlockwise manner. If no
             triangles is passed, a Delauney triangulation is computed. 
             Default=None.
             
            *mask*: optional boolean array of shape (ntri).
             Which triangles are masked out.
             
            *projectedField* python bool (default=True)
             Defines if the data fields has to be projected in the basis of the
             surface. Example where such option really matters:
             Let us consider a 3D laminar flow with the follow velocity field: 
             U=(Ux,0,0). Let us consider the surface S1 with the normal n1=(0,ny,0)
             and the surface S2 with the normal n2=(nx,ny,0). On S1, Ux will be
             in the surface plan, therefore S1.gradient() makes sense. On the 
             other hand, Ux in NOT in the plane of S2.
             
             
            *interpoation*: python string. 
             type of interpolation used. Value: "cubic" or "linear".
             Default="cubic". 

            *kind*: python string.
             Definition of the cubic interpolation type. Value: "geom" or
             "min_E". Default="geom". "min_E" is supposed to be the more 
             accurate, but "geom" is way faster.
             
            *affTrans*: AffineTransformation object.
            
            *linTrans*: LinearTransformation object.
        '''

        self.triangulation = tri.Triangulation(x, y, triangles=triangles, mask=mask)

        self.vx=np.asarray(vx)
        self.vy=np.asarray(vy)
        self.vz=np.asarray(vz)
        
        self.vx_i = None
        self.vy_i = None
        self.vz_i = None
        
        self.data = dict()
        self.data_i = dict()
        
        self.time = float(time)
        
        
        # "private" member variable. Don't play with them if you are not sure...
        self.__z = z
        
        self.__affTrans = affTrans
        self.__linTrans = linTrans
        
        self.__interType = interpolation
        self.__interKind  = kind
        
        self.__projectedField = projectedField
        
    @classmethod
    def readFromFoamFile(cls,
                         varsFile,
                         pointsFile,
                         facesFile,
                         time,
                         viewAnchor,
                         xViewBasis,
                         yViewBasis,
                         srcBasisSrc=[[1,0,0],[0,1,0],[0,0,1]],
                         projectedField=True):
        '''
        Construct from a surface saved  by OpenFOAM in foamFile format.
        '''
        afftrans, lintrans = TriSurfaceFunctions.getTransformation(viewAnchor,
                                                                   xViewBasis,
                                                                   yViewBasis,
                                                                   srcBasisSrc)

        # get x and y vector (in ptsTgt)
        ptsSrc = parseFoamFile_sampledSurface(pointsFile)
        ptsTgt = np.zeros((ptsSrc.shape[0],ptsSrc.shape[1]))
        for i in range(ptsSrc.shape[0]):
            ptsTgt[i,:] = afftrans.srcToTgt(ptsSrc[i,:])

        #get vectors (in vecsTgt)
        vecsSrc = parseFoamFile_sampledSurface(varsFile)
        vecsTgt = np.zeros((vecsSrc.shape[0],vecsSrc.shape[1]))
        if projectedField==True:
            for i in range(vecsSrc.shape[0]):
                vecsTgt[i,:] = lintrans.srcToTgt(vecsSrc[i,:])
        else:
            vecsTgt = vecsSrc
        
        #get triangles
        triangles = parseFoamFile_sampledSurface(facesFile)[:,1:4]
        

        # update class member variables
        return cls(x=ptsTgt[:,0],
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
 
    @classmethod   
    def readFromVTK(cls,
                    vtkFile,
                    time,
                    viewAnchor,
                    xViewBasis,
                    yViewBasis,
                    srcBasisSrc=[[1,0,0],[0,1,0],[0,0,1]],
                    projectedField=True):
        '''
        Construct from a surface saved by OpenFOAM in VTK format.
        '''
        afftrans, lintrans = TriSurfaceFunctions.getTransformation(viewAnchor,
                                                                   xViewBasis,
                                                                   yViewBasis,
                                                                   srcBasisSrc)
       
        # read VTK file
        ptsSrc, triangles, vecsSrc = parseVTK_ugly_sampledSurface(vtkFile)
        
        # Transform the points
        ptsTgt = np.zeros((ptsSrc.shape[0],ptsSrc.shape[1]))
        for i in range(ptsSrc.shape[0]):
            ptsTgt[i,:] = afftrans.srcToTgt(ptsSrc[i,:])
            
        # transform the data
        vecsTgt = np.zeros((vecsSrc.shape[0],vecsSrc.shape[1]))
        if projectedField==True:
            for i in range(vecsSrc.shape[0]):
                vecsTgt[i,:] = lintrans.srcToTgt(vecsSrc[i,:])
        else:
            vecsTgt = vecsSrc

        # update class member variables
        return cls(x=ptsTgt[:,0],
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
  
    # getters #
    #---------#

    @property
    def x(self):
        '''
        Get x coordinate of the grid points.
        '''
        return self.triangulation.x
        
    @property    
    def y(self):
        '''
        Get y coordinate of the grid points.
        '''
        return self.triangulation.y
        
    @property
    def triangles(self):
        '''
        Get triangles from the grid.
        '''
        return self.triangulation.triangles
        
    @property
    def affTrans(self):
        return self.__affTrans
    
    @property
    def linTrans(self):
        return self.__linTrans
        
        
    @property
    def rawPoints(self):
        '''
        '''
        surfacePoints = np.vstack((self.x,self.y,self.__z)).T
        rawPoints = np.zeros((surfacePoints.shape[0],surfacePoints.shape[1]))
        for i in range(surfacePoints.shape[0]):
            rawPoints[i,:] = self.affTrans.tgtToSrc(surfacePoints[i,:])
        return rawPoints
            
    @property        
    def rawData(self):
        '''
        '''
        surfaceData = np.vstack((self.vx,self.vy,self.vz)).T
        rawData = np.zeros((surfaceData.shape[0],surfaceData.shape[1]))
        if self.__projectedField==True:
            for i in range(surfaceData.shape[0]):
                rawData[i,:] = self.linTrans.tgtToSrc(surfaceData[i,:])
        else:
            rawData = surfaceData
        return rawData

        
        
        
    
        
    # setters #
    #---------#


    # class methods #
    #---------------#
    def getInterpolator(self,interpolation='cubic', kind='geom'):
        self.__interType = interpolation
        self.__interKind = kind
        if self.__interType=='cubic':
            self.vx_i = tri.CubicTriInterpolator(self.triangulation, self.vx, kind=self.__interKind)
            self.vy_i = tri.CubicTriInterpolator(self.triangulation, self.vy, kind=self.__interKind)
            self.vz_i = tri.CubicTriInterpolator(self.triangulation, self.vz, kind=self.__interKind)
        elif self.__interType=='linear':
            self.vx_i = tri.LinearTriInterpolator(self.triangulation, self.vx)
            self.vy_i = tri.LinearTriInterpolator(self.triangulation, self.vy)
            self.vz_i = tri.LinearTriInterpolator(self.triangulation, self.vz)
        else:
            raise ValueError('Interpolation must be "cubic" or "linear".')
    
    
    def addField(self,field,fieldname):
        '''
        Add a field F of dimension d (e.g: d=3 for a vector filed) to the
        current TriSurfaceVector object TSV. The grid of F (N points) must be
        identical, in term of number of points and their location, as the grid
        of TSV. F will be stored in TSV.data['fieldName'] or TSV['fieldName'].
        
        Arguments:
            *field*: numpy array of shape (N,d).
            
            *fieldName*: python string.
        '''
        fieldSrc = field
        fieldShape = fieldSrc.shape
        fieldTgt = np.zeros((fieldShape[0],fieldShape[1]))

        if (self.__projectedField==True and fieldShape[1]>1):
            for i in range(fieldShape[0]):
                fieldTgt[i,:] = self.__linTrans.srcToTgt(fieldSrc[i,:])
        else:
            fieldTgt = fieldSrc
        self.data[fieldname] = fieldTgt
            
        
    def addFieldFromFoamFile(self,fieldFile,fieldname):
        '''
        Add a field F of dimension d (e.g: d=3 for a vector filed) to the
        current TriSurfaceVector object TSV. The grid of F (N points) must be
        identical, in term of number of points and their location, as the grid
        of TSV. F will be stored in TSV.data['fieldName'] or TSV['fieldName'].
        
        Arguments:
            *field*: numpy array of shape (N,d).
            
            *fieldName*: python string.
        '''
        #get field
        fieldSrc = parseFoamFile_sampledSurface(fieldFile)
        self.addField(fieldSrc,fieldname)
        
    
    def gradient(self):
        '''
        Calculate and save the gradient at all point of the grid.
        '''
        self.data['dvxdx'],self.data['dvxdy'] = self.vx_i.gradient(self.x(),self.y())
        self.data['dvydx'],self.data['dvydy'] = self.vy_i.gradient(self.x(),self.y())
        self.data['dvzdx'],self.data['dvzdy'] = self.vz_i.gradient(self.x(),self.y())
 
       
    def gradient_i(self,x,y):
        '''
        Calculate the gradient at the point pt(x,y) and return it.
        
        Arguments:
            *x*: python float or numpy array.
             x coordinate of the point pt.
            
            *y*: python float or numpy array.
             y coordinate of the point pt.
             
        Returns:
            *dvxdx, dvxdy, dvydx, dvydy, dvzdx, dvzdy *: python tuple of four float.
             The gradient at point pt(x,y).
        '''
        dvxdx, dvxdy = self.vx_i.gradient(x,y)
        dvydx, dvydy = self.vy_i.gradient(x,y)
        dvzdx, dvzdy = self.vz_i.gradient(x,y)
        return dvxdx, dvxdy, dvydx, dvydy, dvzdx, dvzdy


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

    # catch the begin of the list ogf points of the grid
    # --------------------------------------------------
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

        
    