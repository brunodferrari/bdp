# -*- coding: utf-8 -*-
"""
Created on Sat Dec 11 05:41:07 2021

@author: bferrari
"""


"""
setup.py file for SWIG example
"""

from distutils.core import setup, Extension


example_module = Extension('_bgraph',
                           sources=['bgraph_wrap.cxx', 'bgraph.cpp'],
                           )

setup (name = 'example',
       version = '0.1',
       author      = "SWIG Docs",
       description = """Simple swig example from docs""",
       ext_modules = [example_module],
       py_modules = ["bgraph"],
       )    