#from distutils.core import setup
from setuptools import setup, find_packages

print('\nInstalling profitroll module...\n')

import profitroll

setup(name='profitroll',
      description='Python module for 2D evolution problems using NetCDF format applied to tropopause dynamics',
      version='0.1.0',
      author='Olivier Goux, Lucas Lange, Hugo Levy, Pablo Richard, Mathieu Roule, Maxence Seymat',
      packages=find_packages(),
      requires=['numpy', 'netCDF4', 'os', 'time', 'datetime', 'copy', 'matplotlib', 'IPython', 'base64', 'ipywidgets']
      )