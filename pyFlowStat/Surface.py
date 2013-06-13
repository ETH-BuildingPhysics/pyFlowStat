#===========================================================================#
# load modules
#===========================================================================#
#standard modules
import sys

#scientific modules
import numpy as np
import scipy as sp

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
        self.ReadIMX64 = cdll.LoadLibrary("ReadIMX64.dll")
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
        
        self.vx=[]
        self.vy=[]
        self.vz=[]
        for theY in range(0,width):
            theX=149
            mode = getMode(self.buffer,theX,theY,height,frameOffset)
            if mode >= 0:
                self.vx.append(self.buffer.floatArray[theX + theY*height + frameOffset + componentOffset*(mode*3+1)]*self.buffer.scaleI.factor+self.buffer.scaleI.offset)
                self.vy.append(self.buffer.floatArray[theX + theY*height + frameOffset + componentOffset*(mode*3+2)]*self.buffer.scaleI.factor+self.buffer.scaleI.offset)
                self.vz.append(self.buffer.floatArray[theX + theY*height + frameOffset + componentOffset*(mode*3+3)]*self.buffer.scaleI.factor+self.buffer.scaleI.offset)
            else:
                self.vx.append(0)
                self.vy.append(0)
                self.vz.append(0)
        
        #plot(vx)
        #plot(vy)
        #plot(vz)
    