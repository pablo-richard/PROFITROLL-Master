<div align="center"><img src="profitroll.png"></div>

# <div align="center">PROFITROLL README</div>

This is the GitHub repository of Profitroll: an Object Oriented Python Module for numerical simulations of the tropopause dynamics and satellite water-vapor images.

Quick-presentation:

Profitroll is a general-purpose numerical tool focused on meteorological applications. Its current version features the simulation of stratospheric intrusions into the troposphere under quasi-geostrophic assumptions. It is further capable of producing qualitative images of what could be seen from a meteorological satellite in the water-vapor spectral band (between 6 and 7 microns). The program's architecture is general enough to allow the implementation of new routines, as long as they are compliant with a given template (see the demonstration Jupyter notebook). Developers are free to extend the existing version as this is an opensource project.

docs directory:

Full documentation of the source code can be found in the 'docs' folder. Documentation has been compiled using Sphinx https://www.sphinx-doc.org/en/master/ for elegant and clear access to information regarding the module objects and methods.

profitroll directory:

This directory contains all the relevant python source files of the profitroll module. Its installation is straightforward:
- Download the profitroll directory from this git repository, or clone the full git repository.
- Open a terminal and navigate to the folder where the setup.py file is.
- From here, run the command "python setup.py install"
- You can then check the demonstration jupyter notebook which shows how to use profitroll in pratice.
- Enjoy ;)
