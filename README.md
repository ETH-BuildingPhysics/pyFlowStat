pyFlowStat
==========

purpose:
    Python tools for statistical analyses of flow data.

version:
    3.0
    
Authors:
    Marc Immer   (aaa.aaa@aaa.com)
    Marcel Vonlanthen  (aaa.aaa@aaa.com)
    
python version:
    tested on python 2.7.x
    
python requirement:
    The following python modules must be installed:
        * numpy
        * scipy
        * matplotlib
        
    The following python modules should be installed for the advenced
    features:
        * h5py: for HDF5 save and load capabilites. See methods in
          "PointProbeFunctions.py", "SurfaceFunctions.py", and
          "TriSurfaceFunctions.py" for more information.
          (http://www.h5py.org/)
        * modred: for POD and DMD decomposition. See class POD and DMD. 
          (https://pypi.python.org/pypi/modred)
    
    
Installation:
    1. Download pyFlowStat on your computer
    2. go in the pyFlowStat folder
    3. in super user mode, type "python setup.py install"
    
