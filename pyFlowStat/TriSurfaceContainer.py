'''
TriSurfaceContainer.py
'''
import numpy as np
import os,sys
import glob
#import matplotlib.tri as tri
import matplotlib.tri as tri

from pyFlowStat.TriSurface import TriSurfaceDict
from pyFlowStat.TriSurfaceVector import TriSurfaceVector
from pyFlowStat.TriSurfaceScalar import TriSurfaceScalar
from pyFlowStat.TriSurfaceSymmTensor import TriSurfaceSymmTensor
from pyFlowStat.TriSurfaceMesh import TriSurfaceMesh


class TriSurfaceContainer(object):
    '''
 
    '''
    
    # constructors #
    #--------------#
    
    def __init__(self,triSurfaceMesh):
        '''
        base constructor.
        '''
        self.triSurfaceMesh = triSurfaceMesh
        self.fields=TriSurfaceDict()
        self.data=dict()
        self.time=0.0
        
    @classmethod
    def createFromTriSurface(cls,triSurface,name):
        '''
        '''
        c=cls(triSurfaceMesh=triSurface.triSurfaceMesh)
        c.fields[name]=triSurface
        return c
    
    @classmethod
    def createFromFoamFolder(cls,pathname,xViewBasis,yViewBasis=None,viewAnchor=(0,0,0),time=0.0,names=[]):
        '''
        This function creates a TriSurfaceContainer for an OpenFOAM sampled
        surface folder (e.g./postProcessing/surfaces/0/planeXY/) and loads all
        the fields.
        
        Arguments:
            *pathname*: string.
             Path ot the surface (e.g:/postProcessing/surfaces/0/planeXY/)
             
            *xViewBasis*: list or numpy.array of shape=3.
             x direction of the surface defined in the OpenFOAM coordinate
             system.
             
            *yViewBasis*: list or numpy.array of shape=3.
             y direction of the surface defined in the OpenFOAM coordinate
             system.Automatically cacluated with the normal of the surface if
             set to None. Default=None
             
            *viewAnchor*: list or numpy.array of shape=3.
             Anchor location of the surface defined in the OpenFOAM coordinate
             system. Default=(0,0,0)
             
            *time* float.
             timeStep of the surface. Default=0
             
            *names*: python list
             List of surface to load in the container.
             
        Returns:
            *tsc*: pyflowStat.TriSurfaceContainer object 
        '''
        pointsFile=os.path.join(pathname,'points')
        facesFile=os.path.join(pathname,'faces')

        tsm=TriSurfaceMesh.readFromFoamFile(pointsFile=pointsFile,
                                           facesFile=facesFile,
                                           viewAnchor=viewAnchor,
                                           xViewBasis=xViewBasis,
                                           yViewBasis=yViewBasis)
        
        c=cls(triSurfaceMesh=tsm)
        c.data['name']=os.path.basename(pathname)
        c.data['pathname']=pathname
        c.addFieldFromFoamFolder(names=names,time=time)
        return c

    @classmethod
    def createFromHdf5(cls,
                        hdf5Parser,
                        xViewBasis,
                        yViewBasis=None,
                        viewAnchor=(0,0,0),
                        srcBasisSrc=[[1,0,0],[0,1,0],[0,0,1]]):
        '''
        '''
        tsm = TriSurfaceMesh.readFromHdf5(hdf5Parser=hdf5Parser,
                                          xViewBasis=xViewBasis,
                                          yViewBasis=yViewBasis,
                                          viewAnchor=viewAnchor,
                                          srcBasisSrc=srcBasisSrc)
        return cls(triSurfaceMesh=tsm)


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

        
    # class methods #
    #---------------#
    def __getitem__(self, key):
        '''
        Getter for key "key" on member dictionary "fields"
        '''
        return self.fields[key]
        
    def __setitem__(self, key, item):
        '''
        Add an existing TriSurface<type> (for example: TriSurfaceScalar or 
        TriSurfaceVector) to TriSurfaceContainer.
        '''
        self.addTriSurface(item,key)

        
    def addTriSurface(self,triSurface,name):
        '''
        Add an existing TriSurface<type> (for example: TriSurfaceScalar or 
        TriSurfaceVector) to TriSurfaceContainer.
        
        Arguments:
            *triSurface*: TriSurface<type> object.
             The TrisSurface<type> to add.
            
            *name*: string.
             Name of the added field.
        '''
        if triSurface.triSurfaceMesh is not self.triSurfaceMesh:
            raise ValueError("triSurfaceMesh is not identical")
        self.fields[name]=triSurface
        

    def addFoamScalarField(self,fieldpath,name,time,projectedField=False):
        '''
        '''
        tSurf=TriSurfaceScalar.readFromFoamFile(fieldpath,self.triSurfaceMesh,time=time,projectedField=projectedField)
        self.fields[name] = tSurf
        
    
    def addFoamVectorField(self,fieldpath,name,time,projectedField=False):
        '''
        '''
        tSurf=TriSurfaceVector.readFromFoamFile(fieldpath,self.triSurfaceMesh,time=time,projectedField=projectedField)
        self.fields[name] = tSurf

                               
    def addFoamSymmTensorField(self,fieldpath,name,time,projectedField=False):
        '''
        '''
        tSurf=TriSurfaceSymmTensor.readFromFoamFile(fieldpath,self.triSurfaceMesh,time=time,projectedField=projectedField)           
        self.fields[name] = tSurf
        
    def addFieldFromFoamFolder(self,surfacepath='',names=[],time=0.0):
        '''
        '''
        if surfacepath=='':
            surfacepath=self.data['pathname']
        
        if os.path.exists(surfacepath):
            scalarFields=glob.glob(os.path.join(surfacepath,'scalarField','*'))
            vectorFields=glob.glob(os.path.join(surfacepath,'vectorField','*'))
            symmTensorFields=glob.glob(os.path.join(surfacepath,'symmTensorField','*'))
            namesToLoad = []
            if len(names)<1:
                for name in scalarFields:
                    namesToLoad.append(os.path.basename(name))
                for name in vectorFields:
                    namesToLoad.append(os.path.basename(name))
                for name in symmTensorFields:
                    namesToLoad.append(os.path.basename(name))
            else:
                namesToLoad=[n for n in names]           
            
            for field in scalarFields:
                key=os.path.basename(field)
                if key in namesToLoad:
                    self.addFoamScalarField(field,key,time)
            for field in vectorFields:
                key=os.path.basename(field)
                if key in namesToLoad:
                    self.addFoamVectorField(field,key,time)
            for field in symmTensorFields:
                key=os.path.basename(field)
                if key in namesToLoad:
                    self.addFoamSymmTensorField(field,key,time)
        else:
            raise IOError("Folder does not exist")
            
    def addFieldFromHdf5(self,hdf5Parser,key,names=[],projectedField=False):
        '''
        Add a list field stored in a hdf5 to the TriSurfaceContainer.
        
        Arguments:
            *hdf5Parser*: h5py parser object.
             Parser object of the source hdf5 file.
             
            *key*: python string.
             The time (as a key) to extract from the HDF5. If it does not exist, IOError
             is returned.
             
            *names*: python list of string.
             Name of the fields to extract from the HDF5. It can be a single
             field (names=['oneField']) or multiple. if names=[], all the
             fields are loaded. Default: names=[].
             
            *projectedField*: python bool
             Project the fields in the surface coordinate system.
             Default: projectedField=False.
        
        Usages:
            >>> import h5py
            >>> parser = h5py.File('myfile.h5','r')
            >>> tsc = TriSurfaceContainer.createFromHdf5(parser,xViewBasis=[1,0,0])
            >>> tsc.addFieldFromHdf5(parser,time=1.5,names=['U','S'])
        '''
        gTime = key
        time=float(key)
        if len(names)==0:
            names = hdf5Parser[gTime].keys()
            names.pop(names.index('time'))

        for name in names:
            try:
                
                dataShape = hdf5Parser[gTime][name].value.shape
                if len(dataShape)==1:  #data is a scalar
                    tss = TriSurfaceScalar.readFromHdf5(hdf5Parser=hdf5Parser,
                                                        varName=name,
                                                        triSurfaceMesh=self.triSurfaceMesh,
                                                        key=key,
                                                        projectedField=projectedField)
                    self.fields[name] = tss
                elif len(dataShape)==2 and dataShape[1]==3:  #data is a vector
                    tsv = TriSurfaceVector.readFromHdf5(hdf5Parser=hdf5Parser,
                                                        varName=name,
                                                        triSurfaceMesh=self.triSurfaceMesh,
                                                        key=key,
                                                        projectedField=projectedField)
                    self.fields[name] = tsv
                elif len(dataShape)==2 and dataShape[1]==6:  #data is a symmtensor
                    tsst = TriSurfaceSymmTensor.readFromHdf5(hdf5Parser=hdf5Parser,
                                                             varName=name,
                                                             triSurfaceMesh=self.triSurfaceMesh,
                                                             ley=key,
                                                             projectedField=projectedField)
                    self.fields[name] = tsst
                else:
                    raise IOError('variable of name "'+name+'" is not a'
                                  'scalar, not a vector, not a symmTensor.')
            except KeyError as e:
                print e
                print gTime,name
            except:
                print "Unexpected error:", sys.exc_info()[0]
                #raise IOError('time "'+gTime+'" and/or name "'+name+'" does not '
                              #'exist as key in the HDF5 parser.')