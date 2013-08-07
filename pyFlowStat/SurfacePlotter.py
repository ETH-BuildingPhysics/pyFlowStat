import matplotlib.pyplot as plt
import numpy as np
import pyFlowStat.Surface as Surface

def PlotField(ax,surface,field,vmin,vmax,offset=[0,0]):
    extent=[0,0,0,0]
    extent[0]=surface.extent[0]-offset[0]
    extent[1]=surface.extent[1]-offset[0]
    extent[2]=surface.extent[2]-offset[1]
    extent[3]=surface.extent[3]-offset[1]
    im=ax.imshow(surface.data[field],vmin=vmin,vmax=vmax,interpolation='nearest',extent=extent)
    return im
    
def PlotContour(ax,surface,field,vmin,vmax,offset=[0,0]):
    yrange = np.arange(surface.minY,surface.maxY+surface.dy,surface.dy)
    yrange=np.flipud(yrange)
    xrange = np.arange(surface.minX,surface.maxX+surface.dx,surface.dx)
#    self.OffsetXpos=xrange[self.xpos_left_wall]
#    self.OffsetYpos=yrange[self.ypos_rooftop]
    xrange=xrange-offset[0]
    yrange=yrange-offset[1]

    X,Y = np.meshgrid(xrange, yrange)
    contour_levels = np.linspace(vmin, vmax, 20)
    contour_levels_label = np.linspace(vmin, vmax, 10)
    print field
    cts=ax.contourf(X,Y,surface.data[field],contour_levels,alpha=.75)
#    ax.colorbar()

    C = ax.contour(X,Y,surface.data[field], contour_levels_label, colors='black', linewidth=.5)
    ax.clabel(C, inline=1, fontsize=10)
    return cts

#    xfill=[-100,0,0,100,100,200,200,-100]
#    yfill=[0,0,-100,-100,0,0,-120,-120]
#    fill(xfill,yfill,'k')
#    
#    xlim([-30,130])
#    ylim([-105,130])