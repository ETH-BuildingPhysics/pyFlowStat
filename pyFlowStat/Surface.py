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
        pass
        
    def readFromVC7(self,filename):
        dllpath = os.path.dirname(os.path.realpath(__file__))
        self.ReadIMX64 = cdll.LoadLibrary(dllpath+"\ReadIMX64.dll")
        
        self.buffer = BufferType()
        self.attributeLst = AttributeList()
        self.vx=[]
        self.vy=[]
        self.vz=[]
        
        res = self.ReadIMX64.ReadIM7(filename, byref(self.buffer), byref(self.attributeLst))
        
        print res
        if res>0:
            print "Error reading image"
            return
            
        print "isFloat:" 
        print self.buffer.isFloat
        print "ny"
        print self.buffer.ny
        print self.buffer.nx
        print TypeName[self.buffer.image_sub_type]
        
        theFrame=0
        frameOffset = theFrame * self.buffer.nx * self.buffer.ny * self.buffer.nz;
        width = self.buffer.nx;
        height = self.buffer.ny;
        componentOffset = width * height;
        
        
        self.vx=np.empty((width,height), dtype=float)
        self.vx[:] = np.NAN
        self.vy=np.empty((width,height), dtype=float)
        self.vy[:] = np.NAN
        self.vz=np.empty((width,height), dtype=float)
        self.vz[:] = np.NAN
        if self.buffer.image_sub_type == 3:
            for theY in range(0,width):
                for theX in range(0,height):
                    mode = getMode(self.buffer,theX,theY,height,frameOffset)
                    if mode >= 0:
                        self.vx[theY,theX] = (self.buffer.floatArray[theX + theY*height + frameOffset + componentOffset*(mode*2+1)]*self.buffer.scaleI.factor+self.buffer.scaleI.offset)
                        self.vy[theY,theX] = (self.buffer.floatArray[theX + theY*height + frameOffset + componentOffset*(mode*2+2)]*self.buffer.scaleI.factor+self.buffer.scaleI.offset)
                    else:
                        pass
        if self.buffer.image_sub_type == 5:
            for theY in range(0,width):
                for theX in range(0,height):
                    mode = getMode(self.buffer,theX,theY,height,frameOffset)
                    if mode >= 0:
                        self.vx[theY,theX]=(self.buffer.floatArray[theX + theY*height + frameOffset + componentOffset*(mode*3+1)]*self.buffer.scaleI.factor+self.buffer.scaleI.offset)
                        self.vy[theY,theX]=(self.buffer.floatArray[theX + theY*height + frameOffset + componentOffset*(mode*3+2)]*self.buffer.scaleI.factor+self.buffer.scaleI.offset)
                        self.vz[theY,theX]=(self.buffer.floatArray[theX + theY*height + frameOffset + componentOffset*(mode*3+3)]*self.buffer.scaleI.factor+self.buffer.scaleI.offset)
                    else:
                        pass
        self.ReadIMX64.DestroyBuffer(self.buffer)
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

def getVC7SurfaceList(directory,nr):
    
    filelist=[]

    for files in os.listdir(directory):
        if files.endswith(".vc7"):
            filelist.append(files)
    filelist.sort()
    surfaces=[]
    surfaces=[Surface()]*min(len(filelist),nr)
    #os.chdir(directory)
    
    for i in range(0,min(len(filelist),nr)):
        print("reading " + filelist[i])
        surfaces[i]=Surface()
        surfaces[i].readFromVC7(os.path.join(directory,filelist[i]))
    return surfaces