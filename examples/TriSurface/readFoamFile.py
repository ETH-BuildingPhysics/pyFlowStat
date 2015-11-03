#! /usr/bin/env python


# example.triSurface.readFoamFile

# Surface plot of an OpenFOAM sample plane saved as a foamFile

#See also:
#        * matplotlib.pyplot.tricontourf: surface plot
#        * matplotlib.pyplot.tricontour:  contour plot
#        * matplotlib.pyplot.triplot:     plot mesh


# plot modules
import matplotlib.pyplot as plt

# pyFlowStat modules
from pyFlowStat.TriSurfaceContainer import TriSurfaceContainer



# create an empty triSurface object


# fill object srf with a foamFile sample plane
#   *The openFOAM sample plane is extract from a 3D geometry and has an arbitrary orientation. Therefore, to plot
#    some orientation information are needed to plot the plane in 2D with matPlotLib.
#       *viewAnchor is a point included in the sample plane. This point will become the point (0,0) in the
#        matPlotLib 2D plot.
#       *xViewBasis is a vector coplanar with the sample plane. This axis will become the x-axis in the
#        matPlotLib 2D plot. The norm of xViewAxis is 1.
#       *yViewBasis is a vector coplanar with the sample plane AND perpendicular to xViewAxis. This axis will become the y-axis in the
#        matPlotLib 2D plot. The norm of yViewAxis is 1.
srf=TriSurfaceContainer.createFromFoamFolder('OFsamplePlane_foamFile/planeX/', viewAnchor=[0,0,-0.08], xViewBasis=[0,1,0], yViewBasis=[0,0,1])

print len(srf.triangulation.x)
print len(srf['U'](0))
# plot sample plane U(x) 
#======================
plt.figure('tricontourf')  # create a figure
plt.axis('equal')          # equal length axis
plt.tricontourf(srf.triangulation,srf['U'](0),100)  # fill the figure with the a contour plot with 100 contours
plt.colorbar() # draw colorbar
plt.xlabel('y [m]')
plt.ylabel('z [m]')
plt.title('surface plot with contourf')  # add figure title

plt.figure('shading = flat')   # create a figure
plt.axis('equal')          # equal length axis
plt.tripcolor(srf.triangulation,srf['U'](0),shading='flat')
plt.colorbar() # draw colorbar
plt.xlabel('y [m]')
plt.ylabel('z [m]')
plt.title('surface plot with tripcolor, flat shading')  # add figure title

plt.figure('shading = gouraud')   # create a figure
plt.axis('equal')          # equal length axis
plt.tripcolor(srf.triangulation,srf['U'](0),shading='gouraud')
plt.colorbar() # draw colorbar
plt.xlabel('y [m]')
plt.ylabel('z [m]')
plt.title('surface plot with tripcolor, gouraud shading')  # add figure title


plt.show()