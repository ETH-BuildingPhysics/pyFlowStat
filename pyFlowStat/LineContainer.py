'''
LineContainer.py
'''
import numpy as np
import os
import glob

from pyFlowStat.Line import LineDict
from pyFlowStat.LineScalar import LineScalar
from pyFlowStat.LineVector import LineVector
from pyFlowStat.LineSymmTensor import LineSymmTensor


class LineContainer(object):
    '''
 
    '''
    
    # constructors #
    #--------------#
    
    def __init__(self):
        '''
        base constructor.
        '''
        self.lines=LineDict()
        self.data=dict()
        self.time=0.0
        
    @classmethod
    def createFromFoamFolder(cls,pathname,time=0.0,names=[]):
        c=cls()
        c.data['name']=os.path.basename(pathname)
        c.data['pathname']=pathname
        c.time=time
        #TODO: find names here?
        c.addLinesFromFoamFolder(names=names,time=time)
        return c
        
    # class methods #
    #---------------#
    def __getitem__(self, key):
        '''
        Getter for key "key" on member dictionary "fields"
        '''
        return self.lines[key]
        
    def __setitem__(self, key, item):
        '''
        Add an existing TriSurface<type> (for example: TriSurfaceScalar or 
        TriSurfaceVector) to TriSurfaceContainer.
        '''
        self.addLine(item,key)

        
    def addLine(self,line,name):
        '''
        Add an existing Line<type> (for example: LineScalar or 
        LineVector) to LineContainer.
        
        Arguments:
            *line*: Line<type> object.
             The Line<type> to add.
            
            *name*: string.
             Name of the added field.
        '''

        self.lines[name]=line
        

    def addFoamScalarLines(self,linePath,lnName='',underscoreHeaders=[]):
        '''
        '''
        lnName,headers = getxyfileInfo(linePath,lnName,underscoreHeaders)
        print(headers)
        fileData = np.loadtxt(linePath)
        xyz = fileData[:,0:3]
        for i,h in enumerate(headers):
            keyName = lnName+'_'+h
            s = fileData[:,3+i+0]
            ls = LineScalar(xyz,s)
            self.addLine(ls,keyName)
        
    
    def addFoamVectorLines(self,linePath,lnName='',underscoreHeaders=[]):
        '''
        '''
        lnName,headers = getxyfileInfo(linePath,lnName,underscoreHeaders)
        print(headers)
        fileData = np.loadtxt(linePath)
        xyz = fileData[:,0:3]
        for i,h in enumerate(headers):
            keyName = lnName+'_'+h
            vx = fileData[:,3+i*3+0]
            vy = fileData[:,3+i*3+1]
            vz = fileData[:,3+i*3+2]
            lv = LineVector(xyz,vx,vy,vz)
            self.addLine(lv,keyName)

                               
    def addFoamSymmTensorLines(self,linePath,lnName='',underscoreHeaders=[]):
        '''
        '''
        lnName,headers = getxyfileInfo(linePath,lnName,underscoreHeaders)
        print(headers)
        fileData = np.loadtxt(linePath)
        xyz = fileData[:,0:3]
        for i,h in enumerate(headers):
            keyName = lnName+'_'+h
            txx = fileData[:,3+i*6+0]
            txy = fileData[:,3+i*6+1]
            txz = fileData[:,3+i*6+2]
            tyy = fileData[:,3+i*6+3]
            tyz = fileData[:,3+i*6+4]
            tzz = fileData[:,3+i*6+5]
            lst = LineSymmTensor(xyz,txx,txy,txz,tyy,tyz,tzz)
            self.addLine(lst,keyName)
        
    def addLinesFromFoamFolder(self,lineFolder,lnName=[]):
        '''
        '''
        if os.path.exists(linepath):
            #def getxyfileType(filePath,headers):
#    '''
#    Get the type (scalar, vector, symmtensor, tensor) of xyfile gererated by
#    the sample tool of OpenFOAM
#    '''
#    fr = open(filePath,'r')
#    line0 = fr.readline()
#    fr.close()
#    
#    entries = line0.split('\t')
#    dataType = str()
#    if (len(entries)-3)==len(headers):      #filePath has scalars
#        dataType = 'scalar'
#    elif (len(entries)-3)==3*len(headers):  #filePath has vector
#        dataType = 'vector'
#    elif (len(entries)-3)==6*len(headers):  #filePath has symmtensor
#        dataType = 'symmtensor'
#    elif (len(entries)-3)==6*len(headers):  #filePath has tnesor
#        dataType = 'tensor'
#        
#    return dataType
            pass
        else:
            raise IOError("Folder does not exist")

         
def getxyfileInfo(filePath,lnName='',underscoreHeaders=[]):
    '''
    '''
    path,fileName = os.path.split(filePath)
    extType = '.xy'
    
    if len(lnName)==0:
        lnName = fileName.split('_')[0]

    print(lnName)
    headerOnlyRaw = fileName[(len(lnName)+1):-len(extType)]
    headerOnlyClean = headerOnlyRaw

    # check if some stupid named variables are in the header
    for h in underscoreHeaders:
        if headerOnlyRaw.find(h)!=0:
            headerOnlyClean = headerOnlyClean.replace(h,h.replace('_', ''))
            
    headers = headerOnlyClean.split('_')
    return lnName,headers