'''
TriSurfaceNew.py

New version of TriSurface. TriSurfaceNew will use matplotlib.tri classes
extensively.

!!! Still in alpha version !!!

See domcumention in class definition
'''


class TriSurfaceNew(object):
    """
    class TriSurface
    
    A TriSurface object "triSurf" reads and holds a cutting plane "pl" of a 3D field. "pl" must be planar (...) but can be in
    any orientation. To plot "pl" as a contour plot "cplt", the x and y direction of "cplt" must be specified to "triSurf" with
    xViewBasis and yViewBasis. With xViewBasis and yViewBasis, you can entirely control how "triSurf" is displayed.
    
    Note:
        * Sbasis is the standard basis e1=(1,0,0), e1=(0,1,0), e1=(0,0,1). "triSurf" and the 3D fields are defined in this basis
        * Vbasis is the view basis. Vbasis definition in Sbasis is:v1=xViewBasis, v2=yViewBasis, v3=v1 x v2. v1 and v2 must be in the plane.
    
    See also:
        * matplotlib.pyplot.tricontourf: surface plot
        * matplotlib.pyplot.tricontour:  contour plot
        * matplotlib.pyplot.triplot:     plot mesh
        
    """
    
    def __init__(self):
        pass
    
    
    
    
class affineTransfomatrion(object):
    '''
    A affine transfomation is define as y = A*x + b.
    
    srcCoorSys = [[],]
    '''
    
    def __init__(self,srcBase_inSrc,tgtBase_inSrc): 
        pass
    