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
    def createFromFoamFolder(cls,pathname,time=0.0,names=[],underscoreHeaders=[],filterstr=''):
        '''
        Create a new LineContainer from the line stored in a foam folder. If
        the argument "name" is empty, all the lines are loaded. This 
        constructor only works with lines saved in the *.xy format, which
        coresponds to the "raw" OpenFOAM format. All the lines are saved in the
        member dict "lines" with the follow name:
            * lineName_variableName
        
        In general, the OpenFOAM lines are saved in:
            * "OpenFOAMcase/postProcessing/set/time"
        
        Arguments:
        *pathname*: string
         Path to the foam folder.
         
        *time*: float
         corresponding time of the lines.
         
        *names*: list string
         List of lines to load. If the list is empty, all the lines are loaded.
         If one of the line name has an underscore in it, lnNames must be used.
         Default=[]
         
        *underscoreHeaders*: list of string
         List of header which have an underscore. For example, if the file
         contains "p_rgh", use underscoreHeaders=["p_rgh"]. Default=[]
        '''
        c=cls()
        c.data['name']=os.path.basename(pathname)
        c.data['pathname']=pathname
        c.time=time
        c.addLinesFromFoamFolder(lineFolder=pathname,
                                 lnNames=names,
                                 underscoreHeaders=underscoreHeaders,filterstr=filterstr)
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
             Name of the added field. To follow the convention used in all
             the loader of this class, name should be defined as
             "lineName_variableName".
        '''

        self.lines[name]=line
        

    def addFoamScalarLines(self,linePath,lnName='',underscoreHeaders=[]):
        '''
        Add a line of scalars data to the current LineContainer.
        
        Arguments:
        *filePath*: string
         Path the *.xy file.
         
        *lnName*: string
         Name of the line. Mendatory only if the line name has an underscore.
         Default=''.
         
        *underscoreHeaders*: list of string
         List of header which have an underscore. For example, if the file
         contains "p_rgh", use underscoreHeaders=["p_rgh"]. Default=[]
        '''
        lnName,headers = getxyfileInfo(linePath,lnName,underscoreHeaders)
        fileData = np.loadtxt(linePath)
        xyz = fileData[:,0:3]
        for i,h in enumerate(headers):
            keyName = lnName+'_'+h
            s = fileData[:,3+i+0]
            ls = LineScalar(xyz,s)
            self.addLine(ls,keyName)
        
    
    def addFoamVectorLines(self,linePath,lnName='',underscoreHeaders=[]):
        '''
        Add a line of vectors data to the current LineContainer.
        
        Arguments:
        *filePath*: string
         Path the *.xy file.
         
        *lnName*: string
         Name of the line. Mendatory only if the line name has an underscore.
         Default=''.
         
        *underscoreHeaders*: list of string
         List of header which have an underscore. For example, if the file
         contains "p_rgh", use underscoreHeaders=["p_rgh"]. Default=[]
        '''
        lnName,headers = getxyfileInfo(linePath,lnName,underscoreHeaders)
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
        Add a line of symmTensors data to the current LineContainer.
        
        Arguments:
        *filePath*: string
         Path the *.xy file.
         
        *lnName*: string
         Name of the line. Mendatory only if the line name has an underscore.
         Default=''.
         
        *underscoreHeaders*: list of string
         List of header which have an underscore. For example, if the file
         contains "p_rgh", use underscoreHeaders=["p_rgh"]. Default=[]
        '''
        lnName,headers = getxyfileInfo(linePath,lnName,underscoreHeaders)
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
       
       
    def addFoamLines(self,linePath,lnName='',underscoreHeaders=[]):
        '''
        Add a line of data to the current LineContainer. Data can be scalar,
        vector or symmTensor.
        
        Arguments:
        *filePath*: string
         Path the *.xy file.
         
        *lnName*: string
         Name of the line. Mendatory only if the line name has an underscore.
         Default=''.
         
        *underscoreHeaders*: list of string
         List of header which have an underscore. For example, if the file
         contains "p_rgh", use underscoreHeaders=["p_rgh"]. Default=[]
        '''
        lnName,headers = getxyfileInfo(filePath=linePath,
                                       lnName=lnName,
                                       underscoreHeaders=underscoreHeaders)
        fr = open(linePath,'r')
        line0 = fr.readline()
        fr.close()
        entries = line0.split('\t')
        if (len(entries)-3)==len(headers):      #filePath has scalars
            self.addFoamScalarLines(linePath,lnName=lnName,underscoreHeaders=underscoreHeaders)
        elif (len(entries)-3)==3*len(headers):  #filePath has vector
            self.addFoamVectorLines(linePath,lnName=lnName,underscoreHeaders=underscoreHeaders)
        elif (len(entries)-3)==6*len(headers):  #filePath has symmtensor
            self.addFoamSymmTensorLines(linePath,lnName=lnName,underscoreHeaders=underscoreHeaders)
        


    def addLinesFromFoamFolder(self,lineFolder,lnNames=[],underscoreHeaders=[],filterstr=''):
        '''
        Add all the lines specified by lnNames, of format *.xy,  included in a
        foam folder. Example of foam folder:
            * "/path/to/my/OpenfoamCase/postProcessing/sets/120/".
        
        Arguments:
        *lineFolder*: string
         Path the *.xy file.
         
        *lnNames*: list string
         List of lines to load. If the list is empty, all the lines are loaded.
         If one of the line name has an underscore in it, lnNames must be used.
         Default=[]
         
        *underscoreHeaders*: list of string
         List of header which have an underscore. For example, if the file
         contains "p_rgh", use underscoreHeaders=["p_rgh"]. Default=[]
        '''
        # As it is coded now, this function works. But it is horrible code!!!
        if os.path.exists(lineFolder):
            allFilePath = glob.glob(lineFolder+'/*'+filterstr+'*.xy')
            
            for filePath in allFilePath:
                if len(lnNames)==0:
                    self.addFoamLines(linePath=filePath,
                                      lnName='',
                                      underscoreHeaders=underscoreHeaders)
                # load only the lines included in lnNames
                else:
                    path,fileName = os.path.split(filePath)
                    for lnToLoad in lnNames:
#                        print('working with line '+lnToLoad)
                        if fileName.startswith(lnToLoad)==True:
#                            print('file "'+fileName+'" starts with '+lnToLoad+' = True')
                            self.addFoamLines(linePath=filePath,
                                              lnName=lnToLoad,
                                              underscoreHeaders=underscoreHeaders)
                        
        else:
            raise IOError("Folder does not exist")
            
    def getSortedLineKeys(self,dir=0,filter=None):
        keys=sorted(self.lines, key=lambda key: self.lines[key].xyz[0,dir])
        if filter:
            keys=[k for k in keys if filter in k]
        return keys
        
    def getFilteredLineDict(self,filter):
        return {k:v for k,v in self.lines.iteritems() if filter in k}
    

         
def getxyfileInfo(filePath,lnName='',underscoreHeaders=[]):
    '''
    from an OpenFOAM *.xy, get the line name and headers (name of the data
    contained in the file). The *.xy file format is structured as follow: The 
    data are saved without header in the file, but the headers are part of the
    name of the file! Example for the "file myLine_U_UMean.xy":
        * "myLine"= name of the line
        * "U"= first header. Velocity
        * "UMean" = second header. Mean velocity
        
    This function processes and extracts the data of the file name.
    
    Arguments:
        *filePath*: string
         Path the *.xy file.
         
        *lnName*: string
         Name of the line. Mendatory only if the line name has an underscore.
         Default=''.
         
        *underscoreHeaders*: list of string
         List of header which have an underscore. For example, if the file
         contains "p_rgh", use underscoreHeaders=["p_rgh"]. Default=[]
         
    Returns:
        *lnName*: string
         Name of the line
         
        *headers*: string
         Header in the file name.
    '''
    path,fileName = os.path.split(filePath)
    extType = '.xy'
    
    if len(lnName)==0:
        lnName = fileName.split('_')[0]
        
    headerOnlyRaw = fileName[(len(lnName)+1):-len(extType)]
    headerOnlyClean = headerOnlyRaw

    # check if some stupid named variables are in the header
    for h in underscoreHeaders:
        if headerOnlyRaw.find(h)!=0:
            headerOnlyClean = headerOnlyClean.replace(h,h.replace('_', ''))
            
    headers = headerOnlyClean.split('_')
    return lnName,headers