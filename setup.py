#!/usr/bin/env python
from setuptools import setup, Extension
from Cython.Distutils import build_ext
from Cython.Build import cythonize

setup(
    name='shand',
    version='0.1',
    description='A pipeline for investigating cospeciation in microbiomes',
    scripts=['scripts/shand'],
    url='http://github.com/ryneches/Shand',
    author='Russell Neches',
    author_email='ryneches@ucdavis.edu',
    license='BSD',
    packages=['shand'],
    install_requires=[
        'pandas',
        'screed',
        'hat_trie',
        'scikit-bio',
        'pyprind',
        'psutil',
        'cython'
    ],
    zip_safe=False,
    ext_modules = cythonize( 'shand/quicktree.pyx' ),
    test_suite = 'nose.collector'
)
