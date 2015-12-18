#!/usr/bin/env python

"""
setup.py file for SWIG example
"""

from distutils.core import setup, Extension


CentroidFold_module = Extension('_CentroidFold',
                                sources=['centroid_fold_wrap.cxx'],
                                include_dirs=['..', '../src', '../src/contrafold'],
                                define_macros=[('SWIG',1), ('HAVE_CONFIG_H',1)],
                                libraries=['centroid', 'contrafold', 'RNA', 'boost_system-mt'],
                                library_dirs=['../src', '../src/contrafold']
                                )

setup (name = 'CentroidFold',
       version = '0.0.11',
       author      = "CentroidFold",
       description = """CentroidFold wrapper""",
       ext_modules = [CentroidFold_module],
       py_modules = ["CentroidFold"],
       )
