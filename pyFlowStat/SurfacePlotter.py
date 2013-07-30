import matplotlib.pyplot as plt
import numpy as np
import pyFlowStat.Surface as Surface

def PlotField(ax,surface,field,vmin,vmax):
    im=ax.imshow(surface.data[field],vmin=vmin,vmax=vmax,interpolation='nearest',extent=surface.extent)
    return im
    
def PlotContour(ax,surface,field,vmin,vmax):
    yrange = np.arange(surface.minY,surface.maxY+surface.dy,surface.dy)
    yrange=np.flipud(yrange)
    xrange = np.arange(surface.minX,surface.maxX+surface.dx,surface.dx)
#    self.OffsetXpos=xrange[self.xpos_left_wall]
#    self.OffsetYpos=yrange[self.ypos_rooftop]
#    xrange=xrange-self.OffsetXpos
#    yrange=yrange-self.OffsetYpos

    X,Y = np.meshgrid(xrange, yrange)
    contour_levels = np.linspace(vmin, vmax, 20)
    contour_levels_label = np.linspace(vmin, vmax, 10)
    ax.contourf(X,Y,surface.data[field],contour_levels,alpha=.75)
#    ax.colorbar()
    C = ax.contour(X,Y,surface.data[field], contour_levels_label, colors='black', linewidth=.5)
    ax.clabel(C, inline=1, fontsize=10)


#    xfill=[-100,0,0,100,100,200,200,-100]
#    yfill=[0,0,-100,-100,0,0,-120,-120]
#    fill(xfill,yfill,'k')
#    
#    xlim([-30,130])
#    ylim([-105,130])