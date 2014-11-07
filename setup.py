from distutils.core import setup
from distutils.version import StrictVersion

import glob, os

#scriptlist = glob.glob(os.path.join('bin', '*.py'))



setup(name='pyFlowStat',
      version="3.0",
      packages=['pyFlowStat'],
      description='Python tools for statistical analyses of flow data',
      url='http://www.carmeliet.arch.ethz.ch/',
      author='ETH/EMPA',
      author_email='',
      package_data={'pyFlowStat': ['ReadIMX64.dll']}
      #scripts=scriptlist,
      )

errors=0
print "\n"
try:
    import numpy
    print "Found numpy",numpy.version.version
except ImportError, e:
	errors+=1
        print "\n"
        print "numpy python-package not installed. Please install numpy"
        print "\n"

try:
    import scipy
    print "Found scipy",scipy.version.version
    if not StrictVersion(scipy.version.version) >= StrictVersion('0.8.0'):
        errors+=1
	print "\n"
        print "scipy python-package too old. Please update scipy to at least version 0.8.x"
        print "\n"
except ImportError, e:
        errors+=1
        print "\n"
        print "scipy python-package not installed. Please install scipy"
        print "\n"

try:
    import matplotlib
    print "Found matplotlib",matplotlib.__version__
except ImportError, e:
        errors+=1
        print "\n"
        print "matplotlib python-package not installed. Please install matplotlib"
        print "\n"

print "\nInstallation finished.",errors,"Error(s)"
