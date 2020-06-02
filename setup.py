#!/bin/env python
# -*- coding: utf-8 -*-
"""
File Information:
-----------------
setuptools script for aestimo project
"""
from setuptools import setup
import os, sys

def read(fname):
    return open(os.path.join(os.path.dirname(__file__), fname)).read()


## cython stuff ################################################################
from setuptools.extension import Extension

cmdclass = { }
ext_modules = [ ]

# making sure that the compile Cython files in the distribution are up-to-date
from setuptools.command.sdist import sdist as _sdist

class sdist(_sdist):
    def run(self):
        # Make sure the compiled Cython files in the distribution are up-to-date
        from Cython.Build import cythonize
        cythonize(['psi_at_inf_cython.pyx'])
        cythonize(['aestimo_dd_lib.pyx'])
        _sdist.run(self)
        #print('compiling cython module into c source')

cmdclass['sdist'] = sdist

# using cython if it is installed on the users distribution.
try:
    from Cython.Distutils import build_ext as _build_ext
    #raise ImportError
    use_cython = True
    ext = '.pyx'
except ImportError:
    from setuptools.command.build_ext import build_ext as _build_ext
    use_cython = False
    ext = '.c' 
    
ext_modules += [
    Extension("aestimo.psi_at_inf_cython", [ "psi_at_inf_cython" +ext ]),
    Extension("aestimo.aestimo_dd_lib", [ "aestimo_dd_lib" +ext])
    ]

# if numpy needed to be installed to as a dependency of aestimo, this might enable
# the compilation of the cython extension to continue successfully

class build_ext(_build_ext):
    def finalize_options(self):
        _build_ext.finalize_options(self)
        # Prevent numpy from thinking it is still in its setup process:
        __builtins__.__NUMPY_SETUP__ = False
        try:
            import numpy
        except ImportError:
            #guess
            site_packages = [p for p in sys.path if sys.prefix in p and '-packages' in p][0]
            include_dir = os.path.join(site_packages,'numpy','core','include')
            #still see an error since bdist_wheel is attempted be built before numpy is installed.
            #this would all work except that 'setup_requires' parameter doesn't seem to work with pip.
        else:
            include_dir = numpy.get_include()
        self.include_dirs.append(include_dir)

cmdclass.update({ 'build_ext': build_ext })

################################################################################


setup(  name='aestimo',
        version='2.0.2',
        description='A bandstructure simulator of semiconductor nanostructures called quantum wells.',
        long_description= read('README.md'),
        classifiers=[
          "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
          "Programming Language :: Python :: 2",
          "Development Status :: 5 - Production/Stable",
          "Intended Audience :: Science/Research",
          "Natural Language :: English",
          "Operating System :: OS Independent",
          "Topic :: Scientific/Engineering :: Physics",
          "Topic :: Scientific/Engineering"
           ],
        author='robochat',
        author_email='rjsteed@talk21.com',
        url='http://www.aestimosolver.org',
        license='GPLv3',
        keywords='quantum well semiconductor nanostructure optical transitions',
        package_dir = {'aestimo': ''},
        packages=['aestimo'],
        package_data={'aestimo':['README.md','AUTHORS.md','COPYING.md','CHANGELOG.md',
                                 'psi_at_inf_cython.pyx','psi_at_inf_cython.c',
                                 'aestimo_dd_lib.pyx','aestimo_dd_lib.c','aestimo_dd_lib.pyd',
                                 'tutorials/*','examples/*.py']},
        scripts=['scripts/aestimo','scripts/aestimo_eh'],
        install_requires=['numpy>1.7.0','matplotlib','scipy'],
        zip_safe=False, #we want users to be able to easily see and edit the scripts
        #setup_requires=['numpy'], #causes problems with pip?
        cmdclass = cmdclass,
        ext_modules=ext_modules,
        )
