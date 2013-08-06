#!/bin/python
"""compile the cythonised version of the psi_at_inf function with the command
python setup_cython.py build_ext --inplace --compiler=mingw32
"""
from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext
import numpy

ext_modules = [Extension("psi_at_inf_cython", ["psi_at_inf_cython.pyx"])]

setup(
  name = 'psi_at_inf',
  cmdclass = {'build_ext': build_ext},
  ext_modules = ext_modules,
  include_dirs=[numpy.get_include()]
)