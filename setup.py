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
    filepath = os.path.join(os.path.dirname(__file__), fname)
    if os.path.exists(filepath):
        return open(filepath, encoding='utf-8').read()
    return ''


setup(
    name='aestimo',
    version='4.0.0',
    description='1D Schrödinger-Poisson, Drift-Diffusion and Quantum-Well Heterostructure Simulator with Modern GUI.',
    long_description=read('README.md'),
    long_description_content_type='text/markdown',
    classifiers=[
        "License :: OSI Approved :: GNU General Public License v3 or later (GPLv3+)",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3.11",
        "Programming Language :: Python :: 3.12",
        "Programming Language :: Python :: 3.13",
        "Development Status :: 5 - Production/Stable",
        "Intended Audience :: Science/Research",
        "Natural Language :: English",
        "Operating System :: OS Independent",
        "Topic :: Scientific/Engineering :: Physics",
        "Topic :: Scientific/Engineering",
    ],
    author='sblisesivdin and Aestimo Contributors',
    author_email='sblisesivdin@gmail.com',
    url='https://github.com/aestimosolver/aestimo',
    license='GPLv3',
    keywords='quantum well semiconductor nanostructure optical transitions drift-diffusion solar laser led',
    packages=['aeslibs'],
    py_modules=[
        'aestimo',
        'aestimo_gui',
        'database',
        'config',
        'characterize_solar',
    ],
    package_data={
        'aeslibs': ['*.py'],
    },
    include_package_data=True,
    install_requires=[
        'numpy>=1.20.0',
        'scipy>=1.7.0',
        'matplotlib>=3.4.0',
        'customtkinter>=5.0.0',
        'darkdetect>=0.8.0',
        'packaging>=20.0',
        'pillow>=8.0.0',
    ],
    entry_points={
        'console_scripts': [
            'aestimo = aestimo:main',
            'aestimo-gui = aestimo_gui:main',
        ],
        'gui_scripts': [
            'aestimo-gui-win = aestimo_gui:main',
        ],
    },
    python_requires='>=3.9',
    zip_safe=False,
)
