#!/usr/bin/env python

from distutils.core import setup


setup(name='nest_toolbox',
      version='1.0',
      description='A collection of tools to help building models and visualizing results using Nest simulator',
      author='Mikael Lindahl',
      author_email='lindahlm@csc.kth.se',
      url='http://www.csc.kth.se/~lindahlm',
      packages=['nest_toolbox'],
      package_dir = {'nest_toolbox': 'src'},
      package_data={'NeuroTools': ['README']},
     )
