from distutils.core import setup

import glob, os

#scriptlist = glob.glob(os.path.join('bin', '*.py'))



setup(name='pyFlowStat',
      version="1.0",
      packages=['pyFlowStat'],
      description='Python tools for statistical analyses of flow data',
      url='http://www.carmeliet.arch.ethz.ch/',
      author='ETH/EMPA',
      author_email='',
	  package_data={'pyFlowStat': ['ReadIMX64.dll']},
      include_package_data=True
      #scripts=scriptlist,
      )

try:
    import numpy
except ImportError, e:
        print "\n\n"
        print "numpy python-package not installed. Please install numpy"
        print "\n\n"

try:
    import scipy
except ImportError, e:
        print "\n\n"
        print "scipy python-package not installed. Please install scipy"
        print "\n\n"

try:
    import matplotlib
except ImportError, e:
        print "\n\n"
        print "matplotlib python-package not installed. Please install matplotlib"
        print "\n\n"
