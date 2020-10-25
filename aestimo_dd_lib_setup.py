from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext
import numpy

EXT_MODULES = [Extension("aestimo_dd_lib", ["aestimo_dd_lib.pyx"])]
setup(
    name = 'aestimo_dd_lib' ,
    cmdclass = {'build_ext': build_ext},
    ext_modules = EXT_MODULES,
    include_dirs=[numpy.get_include()]
)