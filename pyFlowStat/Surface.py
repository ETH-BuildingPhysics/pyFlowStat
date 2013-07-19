#===========================================================================#
# load modules
#===========================================================================#
#standard modules
import sys

#scientific modules
import numpy as np
import scipy as sp
import os

# special modules
from ctypes import *

TypeName = ["Image", "2D-PIV-Vector (header, 4x(Vx,Vy))",
            "2D-Vector (Vx,Vy)", "2D-PIV+p.ratio (header, 4x(Vx,Vy), peakratio)",
          "3D-Vector (Vx,Vy,Vz)", "3D-Vector+p.ratio (header, 4x(Vx,Vy), peakratio)"]

WORD=c_ushort

#typedef struct AttributeList 
#{
#   char*          name;
#   char*          value;
#   AttributeList* next;
#} AttributeList;
class AttributeList(Structure):
    pass
AttributeList._fields_=[("name",c_char_p),("value",c_char_p),("next",POINTER(AttributeList))]

#   union 
#	{
#      float*   floatArray;
#      Word*    wordArray;
#   };
class _bufarray(Union):
    _fields_=[("floatArray",POINTER(c_float)),("wordArray",POINTER(WORD))]
    


#typedef struct
#{
#   int         isFloat;
#   int         nx,ny,nz,nf;
#   int         totalLines;
#	int			vectorGrid;			// 0 for images
#	int			image_sub_type;	// BufferFormat_t
#   union 
#	{
#      float*   floatArray;
#      Word*    wordArray;
#   };
#	BufferScaleType	scaleX;		// x-scale
#	BufferScaleType	scaleY;		// y-scale
#	BufferScaleType	scaleI;		// intensity scale
#	bool*			bMaskArray;			// mask array, NULL if no mask exists
#} BufferType;
class BufferScaleType(Structure):
    _fields_=[("factor",c_float),("offset",c_float),("description",c_char*16),
              ("unit",c_char*16)]
    
#typedef struct
#{
#	float	factor;
#	float offset;
#	char	description[16];
#	char	unit[16];
#} BufferScaleType;
class BufferType(Structure):
    _anonymous_ = ("bufarray",)
    _fields_=[("isFloat",c_int),("ny",c_int),("nx",c_int),("nz",c_int),("nf",c_int),
             ("totalLines",c_int),("vectorGrid",c_int),("image_sub_type",c_int),
             ("bufarray",_bufarray),("scaleX",BufferScaleType),("scaleY",BufferScaleType),
             ("scaleI",BufferScaleType),("bMaskArray",POINTER(c_bool))]
             
def getMode(buf,theX_,theY_,width_,frameOffset):
    mode = int(buf.floatArray[theX_ + theY_*width_ + frameOffset])
    if mode<0:
        return -1
    elif mode>4:
        #// interpolated or filled vector
        mode = 4
    mode=mode-1
    return mode
             
#Read file of type IMG/IMX/VEC, returns error code ImReadError_t
#extern "C" int EXPORT ReadIMX ( const char* theFileName, BufferType* myBuffer, AttributeList** myList );
class Surface(object):
    def __init__(self):
        self.vx=[]
        self.vy=[]
        self.vz=[]
        return
        
    def createDataDict(self):
        '''
        Creates the "data" dictionnary from member variables vx, vy, vz, dx, dy
        
        Member variable data (python ditionary) is created.
        
        By default, the following keys are included in data:
            Ux:  [numpy.array.shape=(ny,nx)] Velocity Ux
            Uy:  [numpy.array.shape=(ny,nx)] Velocity Uy
            Uz:  [numpy.array.shape=(ny,nx)] Velocity Uz
            dx:    [float] spacing in x dirction
            dy:    [float] spacing in y dirction
        '''

        self.data = dict()
        self.data['Ux'] = self.vx
        self.data['Uy'] = self.vy
        self.data['Uz'] = self.vz
        self.data['dx'] = self.dx
        self.data['dy'] = self.dy
        
    def generateFields(self):
        '''
        Generates additional dictionary entries.
        '''
        Umag = np.zeros(self.data['Ux'].shape)
        Umag2D = np.zeros(self.data['Ux'].shape)
        Umag = np.sqrt(self.data['Ux']**2+self.data['Uy']**2+self.data['Uz']**2)
        Umag2D = np.sqrt(self.data['Ux']**2+self.data['Uy']**2)
        self.data['Umag']=Umag
        self.data['Umag2D']=Umag2D

        np_dudy,np_dudx=np.gradient(self.vx,-self.dy,self.dx)
        np_dvdy,np_dvdx=np.gradient(self.vy,-self.dy,self.dx)
        vort_z=np_dvdx-np_dudy
        self.data['VortZ']=vort_z
        self.data['TKE']=0.5*(self.vx**2+self.vy**2+self.vz**2)
        self.data['Div2D']=np_dudx+np_dvdy
        
    def generateStatistics(self,MeanFlowSurface):
        '''
        Generate statistics by loading a mean flow surface (of same size)
        '''
        self.data['ux']=self.data['Ux']-MeanFlowSurface.data['Ux']
        self.data['uy']=self.data['Uy']-MeanFlowSurface.data['Uy']
        self.data['uz']=self.data['Uz']-MeanFlowSurface.data['Uz']
        self.data['TKE_fluct']=0.5*(self.data['ux']**2+self.data['uy']**2+self.data['uz']**2)
        
    def readFromVC7(self,filename,v=False):
        '''
        reads PIV vector data in tha Davis format, using the 64bit windows DLL
        '''
        dllpath = os.path.dirname(os.path.realpath(__file__))
        self.ReadIMX64 = cdll.LoadLibrary(dllpath+"\ReadIMX64.dll")
        
        tmpBuffer = BufferType()
        self.attributeLst = AttributeList()
        self.vx=[]
        self.vy=[]
        self.vz=[]
        
        res = self.ReadIMX64.ReadIM7(filename, byref(tmpBuffer), byref(self.attributeLst))
        
        #print res
        if res>0:
            print "Error reading image"
            return
            
        if v:
            print "Size (ny, nx)"
            print tmpBuffer.ny
            print tmpBuffer.nx
            print TypeName[tmpBuffer.image_sub_type]
        
        theFrame=0
        frameOffset = theFrame * tmpBuffer.nx * tmpBuffer.ny * tmpBuffer.nz;
        width = tmpBuffer.nx;
        height = tmpBuffer.ny;
        componentOffset = width * height;
        
        
        self.vx=np.empty((width,height), dtype=float)
        self.vx[:] = np.NAN
        self.vy=np.empty((width,height), dtype=float)
        self.vy[:] = np.NAN
        self.vz=np.empty((width,height), dtype=float)
        self.vz[:] = np.NAN
        
        if tmpBuffer.image_sub_type == 3 or tmpBuffer.image_sub_type == 1:
            for theY in range(0,width):
                for theX in range(0,height):
                    mode = getMode(tmpBuffer,theX,theY,height,frameOffset)
                    if mode >= 0:
                        self.vx[theY,theX] = (tmpBuffer.floatArray[theX + theY*height + frameOffset + componentOffset*(mode*2+1)]*tmpBuffer.scaleI.factor+tmpBuffer.scaleI.offset)
                        self.vy[theY,theX] = (tmpBuffer.floatArray[theX + theY*height + frameOffset + componentOffset*(mode*2+2)]*tmpBuffer.scaleI.factor+tmpBuffer.scaleI.offset)
                    else:
                        pass
        if tmpBuffer.image_sub_type == 4:
            for theY in range(0,width):
                for theX in range(0,height):
                    self.vx[theY,theX] = (tmpBuffer.floatArray[theX + theY*height + frameOffset]*tmpBuffer.scaleI.factor+tmpBuffer.scaleI.offset)
                    self.vy[theY,theX] = -1*(tmpBuffer.floatArray[theX + theY*height + frameOffset + componentOffset]*tmpBuffer.scaleI.factor+tmpBuffer.scaleI.offset)
                    self.vz[theY,theX] = (tmpBuffer.floatArray[theX + theY*height + frameOffset + componentOffset*2]*tmpBuffer.scaleI.factor+tmpBuffer.scaleI.offset)                    
        if tmpBuffer.image_sub_type == 5:
            for theY in range(0,width):
                for theX in range(0,height):
                    mode = getMode(tmpBuffer,theX,theY,height,frameOffset)
                    if mode >= 0:
                        self.vx[theY,theX]=(tmpBuffer.floatArray[theX + theY*height + frameOffset + componentOffset*(mode*3+1)]*tmpBuffer.scaleI.factor+tmpBuffer.scaleI.offset)
                        self.vy[theY,theX]=-1*(tmpBuffer.floatArray[theX + theY*height + frameOffset + componentOffset*(mode*3+2)]*tmpBuffer.scaleI.factor+tmpBuffer.scaleI.offset)
                        self.vz[theY,theX]=(tmpBuffer.floatArray[theX + theY*height + frameOffset + componentOffset*(mode*3+3)]*tmpBuffer.scaleI.factor+tmpBuffer.scaleI.offset)
                    else:
                        pass
#        print tmpBuffer.scaleX.factor
#        print tmpBuffer.scaleX.offset
#        print tmpBuffer.scaleY.factor
#        print tmpBuffer.scaleY.offset
        #self.maxX=
        #self.maxY=
        if np.isnan(self.vz).all():
            self.vz.fill(0)
            
        self.dx=abs(tmpBuffer.scaleX.factor*tmpBuffer.vectorGrid)
        self.dy=abs(tmpBuffer.scaleY.factor*tmpBuffer.vectorGrid)
        self.minX=tmpBuffer.scaleX.factor*tmpBuffer.vectorGrid*(0.5)+tmpBuffer.scaleX.offset
        self.maxY=tmpBuffer.scaleY.factor*tmpBuffer.vectorGrid*(0.5)+tmpBuffer.scaleY.offset
        self.maxX=tmpBuffer.scaleX.factor*tmpBuffer.vectorGrid*(tmpBuffer.ny-0.5)+tmpBuffer.scaleX.offset
        self.minY=tmpBuffer.scaleY.factor*tmpBuffer.vectorGrid*(tmpBuffer.nx-0.5)+tmpBuffer.scaleY.offset
        self.extent=[self.minX-(self.dx/2),self.maxX+(self.dx/2),self.minY-(self.dy/2),self.maxY+(self.dy/2)]
        self.ReadIMX64.DestroyBuffer(tmpBuffer)
        self.createDataDict()
        #plot(vx)
        #plot(vy)
        #plot(vz)
        
#        			mode = (int) theBuffer.floatArray[ theX + theY*width + frameOffset ];
#			if (mode<=0)
#			{	// disabled vector
#				return true;
#			}
#			if (mode>4)
#			{	// interpolated or filled vector
#				mode = 4;
#			}
#			mode--;
#			vx = theBuffer.floatArray[ theX + theY*width + frameOffset + componentOffset*(mode*2+1) ];
#			vy = theBuffer.floatArray[ theX + theY*width + frameOffset + componentOffset*(mode*2+2) ];

def getVC7SurfaceList(directory,nr=0,step=1):
    '''
    Get a list of Surfaces read from PIV data
    '''
    filelist=getVC7filelist(directory,nr,step)
    
    surfaces=[]
    surfaces=[Surface()]*len(filelist)
    #os.chdir(directory)
    if nr==0:
        nr=len(filelist)
    for i in range(0,min(len(filelist),nr)):
        print("reading " + filelist[i])
        surfaces[i]=Surface()
        surfaces[i].readFromVC7(os.path.join(directory,filelist[i]))
    return surfaces
    
def getVC7filelist(directory,nr=0,step=1):
    '''
    Get a list of filenames of PIV vetor data files
    '''
    filelist=[]

    for files in os.listdir(directory):
        if files.endswith(".vc7"):
            filelist.append(files)
    filelist.sort()
    filelist=filelist[0::step]
    if nr==0:
        nr=len(filelist)
    filelist=filelist[0:min(len(filelist),nr)]
    
    return filelist